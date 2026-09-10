import getopt, os, sys
import pandas as pd

def parse_report(input_report):
    sample = os.path.basename(input_report).split('.')[0]
    df = pd.read_csv(input_report,header=None,
                    names=["perc_cov","num_cov","num_assign","rank","taxid","sci_name"]
                    ,index_col=None,sep='\t')
    df = df[df["perc_cov"]>10]
    rank_order = ['kingdom', 'superkingdom','phylum','class','order','family','genus','species']
    map_rank = {'kingdom':'K', 'superkingdom':'D','phylum':'P','class':'C','order':'O','family':'F','genus':'G','species':'S'}
    rank_order.reverse()
    out_dic = {"id":sample,"level":"na","name":'na',"taxid":"na","tag":'na'}
    for lvl in rank_order:
        tmp = df[df['rank']==map_rank[lvl]]
        if(tmp.size>0):
            idx = tmp["perc_cov"].idxmax()
            if idx>-1:
                name = str(df.loc[idx,"sci_name"]).strip()
                taxid = str(df.loc[idx,"taxid"]).strip()
                tag = name.lower().replace(' ','_')
                out_dic = {"id":sample,"level":lvl,"name":name,"tag":tag,"taxid":taxid}
                break
    pd.DataFrame([out_dic]).to_csv(f'{sample}.taxa.csv',index=None)        

def main(argv):
    input_report = ''
    
    try:
        opts, _ = getopt.getopt(argv, "i:",["input_report="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--input_report","-i"]:
            input_report = arg
        
    parse_report(input_report)

if __name__ == '__main__':
    main(sys.argv[1:])
