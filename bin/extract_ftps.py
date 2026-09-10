import getopt, os, sys
import pandas as pd
from random import sample

def extract_records(in_tsv,tax_ids,selected_columns,subsample=-1):
    
    def organism_name(names):
        name_set = list(map(lambda n: set(n.split()),names))
        common_name= list(name_set[0].intersection(*name_set))
        # find the right order
        
        tmp = names[0].split()
        ordered_name = [''  for _ in tmp]
        for cn in common_name:
           idx = tmp.index(cn)
           ordered_name[idx] = cn.replace('.','')
        ordered_name = [w for w in ordered_name if len(w)>0]
        return '_'.join(map(str.lower,ordered_name))

    print('Loading assembly report...')
    df = pd.read_csv(in_tsv,sep='\t',header=0,skiprows=1,low_memory=False)
    df.rename(columns={"# assembly_accession":"assembly_accession"},inplace=True)
    
    df_lst = []
    selected_columns.append("organism_name")
    selected_columns = list(set(selected_columns))
    print('Filtering selected taxa...')
    for tx in tax_ids:
        print(f'\t-Taxon:{tx}')
        tmp_df = df.loc[df["species_taxid"] == tx, selected_columns]
        sel_orgname = tmp_df["organism_name"].iloc[:20].to_list()
        orgname = organism_name(sel_orgname)
        
        print(f'\t\t-species_name:{orgname}')
        tmp_df["species_name"] = orgname
        
        if subsample>0 and len(tmp_df.index)>subsample:
            tmp_df = tmp_df.loc[sample(list(tmp_df.index),subsample)]
        df_lst.append(tmp_df)
    print('sp') 
    cdf = pd.concat(df_lst)
    cdf = cdf[selected_columns+["species_name"]]
    filename = f"selected_assembly_records_{['all','smpl-'+str(subsample)][subsample>0]}.tsv"
    cdf.to_csv(filename,sep='\t',index=False)
    print(f'Saved as {filename}')
   

def main(argv):
    assembly_report = ''
    tax_ids = []
    selected_col = ["assembly_accession","ftp_path"]  
    samples = -1
    try:
        opts, _ = getopt.getopt(argv, "r:x:c:s:",["report=","tax_id=","columns=","sample="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--report","-r"]:
            assembly_report = arg
        if opt in ["--tax_id","-x"]:
            tax_ids = list(map(int,arg.split(',')))
        if opt in ["--columns","-c"]:
            selected_col += arg.split(',')
        if opt in ["--sample","-s"]:
            samples = int(arg)
    
    extract_records(assembly_report,tax_ids,selected_col,samples)
if __name__ == '__main__':

   main(sys.argv[1:])
   