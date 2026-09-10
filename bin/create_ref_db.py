import getopt, os, sys
import pandas as pd
import sqlite3 as lit
from glob import glob

def creat_db(db_path,replace):
   
    if replace!="false":
        dbname,ext = os.path.basename(db_path).split('.')
        prv_dbs = glob(os.path.join(os.path.dirname(db_path),f"{dbname}*.{ext}"))
        prv_dbs.sort(reverse=True)
        ll = len(prv_dbs)
        for db in prv_dbs:
            nr = str(ll).zfill(3)
            os.rename(db,os.path.join(os.path.dirname(db_path),f"{dbname}_{nr}.{ext}"))
            ll-=1

    conn = lit.connect(db_path)
    tables_sql = []
    tables_sql.append("CREATE TABLE genomes (" \
        "id                  INTEGER PRIMARY KEY,"
        "sample  VARCHAR UNIQUE," \
        "taxonomy                VARCHAR," \
        "tax_lvl                VARCHAR," \
        "taxid               VARCHAR," \
        "source              VARCHAR," \
        "collection_date     VARCHAR," \
        "geographic_location VARCHAR," \
        "mlst_loci           VARCHAR," \
        "mlst_alleles        VARCHAR," \
        "mlst_st             INTEGER," \
        "clonal_complex      VARCHAR," \
        "cgmlst_loci         TEXT," \
        "cgmlst_alleles      TEXT," \
        "cgmlst_coverage     FLOAT);")

    tables_sql.append("CREATE TABLE orfs (" \
        "id   INTEGER PRIMARY KEY AUTOINCREMENT," \
        "sample VARCHAR," \
        "accession          VARCHAR," \
        "start_index        INTEGER," \
        "end_index          INTEGER," \
        "strand             CHAR)")

    tables_sql.append("CREATE TABLE contigs (" \
        "id       INTEGER PRIMARY KEY AUTOINCREMENT," \
        "sample VARCHAR," \
        "accession          VARCHAR);")    
    
    tables_sql.append("CREATE TABLE annotations (" \
        "id                 INTEGER PRIMARY KEY AUTOINCREMENT," \
        "sample VARCHAR," \
        "orfs_accession     VARCHAR," \
        "sseqid             VARCHAR," \
        "pident             DECIMAL," \
        "length             INTEGER," \
        "bitscore           REAL," \
        "qcov               DECIMAL," \
        "scov               DECIMAL," \
        "gaps               VARCHAR," \
        "mismatch           VARCHAR," \
        "stitle             VARCHAR," \
        "ref_db             VARCHAR);")

    for cmd in tables_sql:
        try:
            conn.execute(cmd)
        except Exception as e:
            print(e)
            return -1
    conn.commit()
    return 0

def main(argv):
    db_path = ''
    replace_prv = False
    try:
        opts, _ = getopt.getopt(argv, "d:r:",["db=","replace="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--db","-d"]:
            db_path = arg

        if opt in ["--replace","-r"]:
            replace_prv = arg

    creat_db(db_path,replace_prv)

if __name__ == '__main__':
    main(sys.argv[1:])