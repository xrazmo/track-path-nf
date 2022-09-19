import getopt, os, sys
from time import sleep
from xml.dom import INVALID_STATE_ERR
from Bio import Entrez
from tqdm.contrib.concurrent import process_map
from glob import glob
from xml.sax.handler import ContentHandler
from xml.sax import make_parser, handler
import pandas as pd
import sqlite3 as lit

class BioSampler(ContentHandler):

    def __init__(self):
        self.cur_tag = None
        self.cur_attr = None
        self.dataDic = {}
        self.dataCollect = []
        self.flag = False

    def startElement(self, tag, attrs):
        if tag == 'Attribute':
            self.cur_attr = attrs.get('display_name')
        elif tag == 'BioSample':
            self.dataDic["accession"]=attrs.get("accession")
            self.dataDic["publication_date"]=attrs.get("publication_date")
            self.dataDic["submission_date"]=attrs.get("submission_date")
        elif tag == 'Organism':
            self.dataDic["taxonomy_name"]=attrs.get("taxonomy_name")
            self.dataDic["taxonomy_id"]=attrs.get("taxonomy_id")

        self.cur_tag = tag

    def endElement(self, tag):
        if tag=="BioSample":
            self.dataCollect.append(self.dataDic)
            self.dataDic = {}

    def characters(self, content):
        if len(content.strip('\n ')) > 0:
            if self.cur_tag == 'Attribute':
                self.dataDic[self.cur_attr] = content
            if self.cur_tag == 'Title':
                self.dataDic["Title"] = content

    def parse(self, inDir,outDir,filename="BioSample_parseOut.csv"):
        for ff in glob(os.path.join(inDir,'*.xml')):
            with open(ff) as hdl:
                parser = make_parser()
                parser.setFeature(handler.feature_namespaces, 0)
                parser.setContentHandler(self)
                parser.parse(hdl)

        df = pd.DataFrame(self.dataCollect)
        df.to_csv(os.path.join(outDir,filename))

def fork_fetching(assembly_report_csv,params):

    df = pd.read_csv(assembly_report_csv,header=0,sep='\t')
    ids = df['biosample'].tolist()
    id_size = len(ids)
    ch_size = round(id_size/params["cpus"])
    chunks=[]
    for i in range(0, id_size, ch_size):
     chunks.append({**{"ids":ids[i:i+ch_size]},**params})

    process_map(fetchEtrez,chunks)


def fetchEtrez(params):
        
        Entrez.email = ""
        tmpgi = []
        k = 0
        tmpdir = os.path.join(params["out"],'.tmp')
        os.system(f'mkdir -p {tmpdir}')
        while k < len(params["ids"]):

            try:
                if len(tmpgi) > 10 or (k >= len(params["ids"])-1 and len(tmpgi)>0):
                    sleep(1)
                    outFile = os.path.join(tmpdir, f'{tmpgi[0]}.tmp')
                    print(f'Fetching in {outFile}', end='')
                    sleep(0.1)
                    handle = Entrez.efetch(db=params["db"], id=f"{','.join(tmpgi)}", 
                            rettype=params["rettype"], retmode=params["retmode"])

                    bin_mode = ['', 'b'][params["ext"] == 'xml']
                    with open(outFile, 'w' + bin_mode) as hdl:
                        hdl.write(handle.read())
                    os.rename(outFile, os.path.join(tmpdir, f'{tmpgi[0]}.{params["ext"]}'))
                    print('=> saved.')
                    tmpgi = []
                tmpgi.append(params["ids"][k])
                k += 1
            except Exception as e:
                print(e)
                sleep(1)
        print('Thread finished!')

def update_genome(db,in_tsv):

    def clean_col(content):
        content = str(content).strip().lower()
        rep_str = ["missing",'unknown','not applicable','not collected']
        eq_str = ['nan','none','na','n/a','null']
        for s in eq_str:
            if len(content)<1 or content==s:
                return ''

        for s in rep_str:
            content = content.replace(s,"")
        return content

    def fix_date(dt):
        dt =clean_col(dt)
        if len(dt)<4: # at least should contain the year
            return ''
        # remove time
        dt = dt.split('t')[0]
        y,m,d = '','',''
        for sep in ['-','/']:
            pp = dt.split(sep)
            if len(pp)== 3:
                y,m,d = [pp,[pp[2],pp[0],pp[1]]][sep=='/']
                break
            elif len(pp)== 2:
                y,m = pp
                break
            elif len(pp)== 1 and len(pp[0])==4:
                y = pp[0]
                break
        return '-'.join(list(filter(lambda x: len(x)>0,[y,m,d])))


    conn = lit.connect(db)
    cur = conn.cursor()
    

    df = pd.read_csv(in_tsv,header=0,index_col=0)
    df['src']= df.apply(lambda row: "; ".join(
                                              list(
                                                filter(
                                                    lambda it: len(it)>1,[clean_col(row["isolation source"]),clean_col(row["host"])]
                                                      )
                                                  )
                                            ),axis=1)
    df['date'] = df.apply(lambda row : fix_date(row["collection date"]),axis=1)
    df['loc'] = df.apply(lambda row : clean_col(row["geographic location"]),axis=1)

    for i,row in df.iterrows():
        print(f"\rprogress {i}/{len(df.index)}",end='')
        cur.execute("UPDATE genomes SET collection_date=?,source=?,geographic_location=? WHERE biosample = ?",
        (row["date"],row["src"],row["loc"],row["accession"]))

    conn.commit()
    conn.close()
    
def main(argv):
    inputfile = ''
    params = {}
    try:
        opts, args = getopt.getopt(argv, "hi:o:d:t:m:x:p:r:")
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            params["out"] = arg
        elif opt in ("-d"):
            params["db"] = arg
        elif opt in ("-t"):
            params["rettype"] = arg
        elif opt in ("-m"):
            params["retmode"] = arg
        elif opt in ("-x"):
            params["ext"] = arg
        elif opt in ("-p"):
            params["cpus"] = int(arg)
        elif opt in ("-r"):
            params["ref_db"] = arg

    # fork_fetching(inputfile,params)
    out_fname ="BioSample_parseOut.csv"
    # saxObj = BioSampler()
    # saxObj.parse(os.path.join(params["out"],".tmp"),params["out"],filename=out_fname)
    # os.system(f'rm {os.path.join(params["out"],".tmp/*.xml")}')
    update_genome(params["ref_db"],os.path.join(params["out"],out_fname))


if __name__ == '__main__':
    main(sys.argv[1:])