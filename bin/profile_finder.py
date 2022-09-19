import getopt, os, sys
import pandas as pd

def get_mlst_profile(in_tsv,blastout):
    # identifying the ST from blast results, considering analyzing one isolate at the time
    profiles = pd.read_csv(in_tsv,sep='\t',header=0,index_col=0)
    df = pd.read_csv(blastout,sep='\t',header=None)
    df.columns=["qseqid","sacc","pident","length","slen","qlen"]
    # report the perfect match
    
    df['cov'] = df.apply(lambda row: round(100.0 * row["length"]/row["slen"],2),axis=1)
    df['sacc'] = df.apply(lambda row: row.sacc+'*' if (row["pident"]>90 and row["pident"]<100 and row["cov"]>90 and row["cov"]<100) else row.sacc ,axis=1)

    perfect_match_cond = map(lambda locus: "({}=={})".format(*locus.split('_')),
                         df.loc[(df["cov"]==100) & (df["pident"]==100),"sacc"].to_list())
    perfect_match = profiles[profiles.eval('&'.join(perfect_match_cond))]
    st = ''
    cc = ''
    if len(perfect_match) == 1:
        st = perfect_match.index[0]     
        if "clonal_complex" in profiles.columns:
            cc = perfect_match['clonal_complex']
    
    pdict={}
    for sacc in df.loc[(df["cov"]==100) & (df["pident"]==100),"sacc"].to_list():
         locus,allele = sacc.split("_") 
         pdict[locus] = allele  

    result = {"st":[st],"loci":[],"alleles":[],'clonal_complex':cc}
    for col in profiles.columns:
        if col in ["clonal_complex"]:
            continue
        result["loci"].append(col)
        al = '?'
        if col in pdict:
           al=pdict[col]
        result["alleles"].append(al)
    
    result["loci"] = [','.join(result["loci"])]
    result["alleles"] = [','.join(result["alleles"])]

    df_o = pd.DataFrame.from_dict(result)
    fname = os.path.basename(blastout).split('.')[0]
    df_o.to_csv(f'{fname}.profile.txt',sep='\t')

def get_cgmlst_profile(in_tsv,blastout):
    
    loci = []
    with open(in_tsv) as hdl:
        for row in hdl:
            loci.append(row.strip('\t').split('\t')[0])
    #remove the header
    loci.pop(0)
    df = pd.read_csv(blastout,sep='\t',header=None)
    df.columns=["qseqid","sacc","pident","length","slen","qlen"]
   
    # report the perfect match
    df['cov'] = df.apply(lambda row: round(100.0 * row["length"]/row["slen"],2),axis=1)

    pdict={}
    for sacc in df.loc[(df["cov"]==100) & (df["pident"]==100),"sacc"].to_list():
         locus,allele = sacc.split("__") 
         pdict[locus] = allele  
    
    # having an order for profiles
    result = {"loci":[],"alleles":[]}
    cov = 0
    for locus in loci:
        result["loci"].append(locus)
        al = '?'
        if locus in pdict:
           al=pdict[locus]
           cov += 1.
        result["alleles"].append(al)
    
    result["coverage"] = [round(100. * cov / len(loci),2)]
    result["loci"] = [",".join(result["loci"])]
    result["alleles"] = [",".join(result["alleles"])]
 
    df_o = pd.DataFrame.from_dict(result)
    fname = os.path.basename(blastout).split('.')[0]
    df_o.to_csv(f'{fname}.profile.txt',sep='\t')

def main(argv):
    profile_tsv = ''
    type = ''
    blastout = ''
    
    try:
        opts, _ = getopt.getopt(argv, "p:b:t:",["profile=",'blastout=',"type="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--profile","-p"]:
            profile_tsv = arg
        if opt in ["--type","-i"]:
            type = arg
        if opt in ["--blastout","-b"]:
            blastout = arg
        
    if type=='mlst':
        get_mlst_profile(profile_tsv,blastout)
    elif type =='cgmlst':
        get_cgmlst_profile(profile_tsv,blastout)

if __name__ == '__main__':

    main(sys.argv[1:])
