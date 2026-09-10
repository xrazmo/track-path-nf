import getopt, os, sys
from re import T
import pandas as pd
import sqlite3 as lit

def update_genome(db,in_tsv):

    conn = lit.connect(db)
    cur = conn.cursor()

    df = pd.read_csv(in_tsv,header=0,sep="\t")
   
    for i,row in df.iterrows():
        print(f"\rprogress {i}/{len(df.index)}",end='')
        cur.execute("UPDATE genomes SET bioproject=?,biosample=?,taxid=?,name=? WHERE assembly_accession = ?",
        (row["bioproject"],row["biosample"],row["species_taxid"],row["organism_name"],row["assembly_accession"]))

    conn.commit()
    conn.close()

def main(argv):
    db = ''
    in_tsv = ''
    
    try:
        opts, _ = getopt.getopt(argv, "d:i:t:",["db=",'in='])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--db","-d"]:
            db = arg
        if opt in ["--in","-b"]:
            in_tsv = arg
    
    
    update_genome(db,in_tsv)
    
if __name__ == '__main__':

    main(sys.argv[1:])