


import getopt, os, sys
import json
from subprocess import Popen
from glob import glob

def mkdir(dir):
    try:
        os.mkdir(dir)
    except:
        pass

def dnld_card(dir,ignore):
   out_dir = os.path.join(dir,'dmnd')
   mkdir(out_dir)

   os.chdir(out_dir)
   prv_fasta = glob('card.fasta')  
   prv_dmnd = glob('card.dmnd')  
   
   if ignore or (len(prv_fasta)==0 and len(prv_dmnd)==0):
        tmpDir = "./_tmp"
        mkdir(tmpDir)
        os.chdir(tmpDir)
        cmd_lst = ['wget https://card.mcmaster.ca/latest/data', 
        'mv ./data ./card-data.tar.bz2', 
        'tar -xvf card-data.tar.bz2',
        'mv protein_fasta_protein_knockout_model.fasta ../card.fasta',
        'rm -rf ./*']
        for cmd in cmd_lst:
            Popen(cmd, shell=True).wait()
    
def dnld_vfdb(dir,ignore):

   out_dir = os.path.join(dir,'dmnd')
   mkdir(out_dir)
   
   os.chdir(out_dir)
   prv_fasta = glob('vfdb.fasta')  
   prv_dmnd = glob('vfdb.dmnd')  
   
   if ignore or (len(prv_fasta)==0 and len(prv_dmnd)==0):
        tmpDir = "./_tmp"
        mkdir(tmpDir)
        os.chdir(tmpDir)
        cmd_lst = ['wget --output-document=VFDB_setA_pro.fas.gz http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz', 
        'gzip -d VFDB_setA_pro.fas.gz', 
        'mv VFDB_setA_pro.fas ../vfdb.fasta',
        'rm -rf ./*']
        for cmd in cmd_lst:
            Popen(cmd, shell=True).wait()

def merge_fasta(fa_dir,name):
    ffs = glob(os.path.join(fa_dir,'*.fasta'))
    sep='__'
    content = []
    for i,f in enumerate(ffs):
        print(f"\r{i+1}/{len(ffs)}",end='') 
        filename = os.path.basename(f).split('.fasta')[0]
        with open(f) as hdl:
            for line in hdl:
                if line.startswith('>'):
                    line = line.replace('>',f">{filename}{sep}")
                content.append(line)
        
    with open(os.path.join(fa_dir,name),'w') as hdl:
        hdl.write(''.join(content))

    for f in ffs:
        os.remove(f)

def dnld_cgmlst(dir,conf_dic,ignore):
   
   out_dir = os.path.join(dir,'cgmlst')
   mkdir(out_dir)

   for dd in conf_dic["cgmlst"]:
        os.chdir(out_dir)
        db_name = dd["db_name"]
        fa_lnk = dd["url"]
        
        if not os.path.isdir(db_name) or ignore:
            mkdir(db_name)
            os.chdir(db_name)
            zipfilename = f"{db_name}-cgmlst.zip"
            tsvfilename = f"{db_name}.loci.tsv"

            cmd_lst = [f'wget --output-document={tsvfilename} {fa_lnk}/locus/?content-type=csv',
                    #    f'mv index.html?content-type=csv {tsvfilename}',
                       f'wget --output-document={zipfilename} {fa_lnk}/alleles', 
                    #    f'mv index.html {zipfilename}',
                       f'unzip {zipfilename}']

            for cmd in cmd_lst:
                print(cmd)
                Popen(cmd, shell=True).wait()

            allele_f = f'{db_name}.alleles.fna'
            merge_fasta(os.path.join(out_dir,db_name),allele_f)

def dnld_mlst(dir,conf_dic,ignore):        
    
    out_dir = os.path.join(dir,'mlst')
    mkdir(out_dir)
    for dd in conf_dic["mlst"]:
        os.chdir(out_dir)
        db_name = dd["db_name"]

        if not os.path.isdir(db_name) or ignore:
            mkdir(db_name)
            os.chdir(db_name)
            for locus in dd["loci"]:
                url = dd["url_template"].replace('***',locus)
                Popen(f"wget --output-document={locus}.fasta {url}", shell=True).wait()

            Popen(f"wget --output-document={db_name}.profiles.tsv {dd['profile']}", shell=True).wait()

            cmd = f"for f in *.fasta; do (cat \"${{f}}\"; echo) >> {db_name}.alleles.fna ; done"
            Popen(cmd, shell=True).wait()
            Popen("rm *.fasta", shell=True).wait()

def main(argv):
    out_dir = ''
    ignore = False
    dbs_config = {}
    try:
        opts, _ = getopt.getopt(argv, "d:i:r:",["out_dir=","ignore=","ref_json="])
       
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--out_dir","-d"]:
            out_dir = arg
        if opt in ["--ignore","-i"]:
            ignore = True
        if opt in ["--ref_json","-g"]:
            refdb_json = arg
            with open(refdb_json) as hdl:
                dbs_config = json.load(hdl)
               
    dnld_card(out_dir,ignore)
    dnld_vfdb(out_dir,ignore)
    dnld_cgmlst(out_dir,dbs_config,ignore)
    dnld_mlst(out_dir,dbs_config,ignore)

if __name__ == '__main__':

    main(sys.argv[1:])
