import getopt, os, sys
import pandas as pd
import sqlite3 as lit

def save_samples(db,incsv):
 
    conn = lit.connect(db)
    cur = conn.cursor()

    df = pd.read_csv(incsv,header=0)
    row = df.loc[0]
    cur.execute(f"INSERT INTO genomes (sample,tax_lvl,taxonomy,taxid) VALUES (?,?,?,?)"
    ,(str(row["id"]),row["level"],row["tag"],str(row["taxid"])))

    conn.commit()
    conn.close()


def save_seqtype(db,intsv,type):
    accession = os.path.basename(intsv).split('.')[0]
    accession = accession.replace('-v-','.')
    conn = lit.connect(db)
    cur = conn.cursor()

    df = pd.read_csv(intsv,header=0,index_col=0,sep='\t')
    data = df.iloc[0]
    if(type=='mlst'):
        cur.execute(f"UPDATE genomes set mlst_st=?,mlst_loci=?,"  
            f"mlst_alleles=?,clonal_complex=? where sample=?",(str(data.st),data.loci,data.alleles,data.clonal_complex,accession))
    else:
        cur.execute(f"UPDATE genomes set cgmlst_loci=?,cgmlst_coverage=?,"  
            f"cgmlst_alleles=? where sample=?",(data.loci,float(data.coverage),data.alleles,accession))

    conn.commit()
    conn.close()

def save_orfs(db,intsv):
    sample = os.path.basename(intsv).split('.')[0]
    sample = sample.replace('-v-','.')
  
    commands = []
    ge_acc = set()
    with open(intsv) as hdl:
        for line in hdl:
            if(line.startswith('#')):
                continue
            else:
                prt = line.split(';')[0].split('\t')
                genome_accession = prt[0]
                ge_acc.add(genome_accession)
                orf_cnt = prt[-1].split('_')[-1]   
                orf_accession = f"{genome_accession}_{orf_cnt}"       
                sidx,eidx,_,strand = prt[3:7]
                commands.append("INSERT INTO orfs (sample,accession,start_index,end_index,strand) "
                                f"VALUES ('{sample}','{orf_accession}',{sidx},{eidx},'{strand}')")

    # chech if it has been already stored
    conn = lit.connect(db)
    cur = conn.cursor()
    ge_acc = list(ge_acc)
    row = cur.execute("SELECT count(id) from CONTIGS WHERE sample=? and accession=?",(sample,ge_acc[0])).fetchone()

    if row[0]>0:
        return

    for acc in ge_acc:
        commands.append(f"INSERT INTO contigs (sample,accession) VALUES ('{sample}','{acc}')")


    for cmd in commands:
        cur.execute(cmd)

    conn.commit()
    conn.close()    

def save_annotation(db,intsv,idty_th=50,cov_th=50):
    filename = os.path.basename(intsv).split('.')[0]
    sample,ref_db = filename.split('__')
    sample = sample.replace('-v-','.')
    
    conn = lit.connect(db)
    cur = conn.cursor()
    
    row = cur.execute("SELECT count(id) from annotations WHERE sample=? and ref_db=?",(sample,ref_db)).fetchone()
    
    if row[0]>0:
        print(f'The genome {sample} has been already added to database.')
        return

    df = pd.read_csv(intsv,header=None,sep='\t')
    df.columns = ["qseqid","sseqid","pident","qcovhsp","scovhsp","mismatch","gaps","evalue",
    "bitscore","length","qlen","slen","qstart","qend","sstart","send","stitle"]
   
    df['qcov'] = df.apply(lambda row: round(min(100,100.0 * row.length/(row.qlen/3.)),2),axis=1 )
    df['scov'] = df.apply(lambda row: round(100.0 * row.length/row.slen,2),axis=1)
    df[['pident', 'qcov']]=df[['pident', 'qcov']].astype(float)
    filt_df = df.loc[(df["qcov"] > cov_th) & (df["pident"]> idty_th)]
   

    for _,row in filt_df.iterrows():
       
        cur.execute("INSERT INTO annotations (sample,orfs_accession,sseqid,pident,qcov,scov"
                    ",gaps,mismatch,stitle,ref_db,bitscore,length) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    (sample,row.qseqid,row.sseqid,row.pident,row.qcov,row.scov
                    ,row.gaps,row.mismatch,row.stitle,ref_db,row.bitscore,row.length))

    conn.commit()
    conn.close()

def main(argv):
    db = ''
    type = ''
    intsv = ''
    
    try:
        opts, _ = getopt.getopt(argv, "d:i:t:",["db=",'in=',"type="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--db","-d"]:
            db = arg
        if opt in ["--type","-t"]:
            type = arg
        if opt in ["--in","-b"]:
            intsv = arg
    
    if type == 'samples':
        save_samples(db,intsv)
    elif type in ['cgmlst','mlst']:
        save_seqtype(db,intsv,type)
    elif type=='orfs':
        save_orfs(db,intsv)
    elif type =='annotations':
        save_annotation(db,intsv)
    

if __name__ == '__main__':

    main(sys.argv[1:])