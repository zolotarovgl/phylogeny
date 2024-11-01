import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml



#### Tool functions #####
def align(fasta_file, output_file, ncpu, mafft_opt):

    logging.info(f"Aligning sequences with mafft")

    cmd = (f"mafft --reorder --thread {ncpu} {mafft_opt} {fasta_file} > {output_file}")
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f"Alignment done: {output_file}")

def trim(input_file,output_file,mode = "kpic-gappy", g = "0.7"):
    cmd = f"clipkit {input_file} -m {mode} -o {output_file} -g {g}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def align_and_trim(input_file,output_file,ncpu = 1,mafft_opt = "",clipkit_mode = "kpic-gappy",clipkit_g = 0.7, clean = True):
    if not os.path.exists(input_file):
        logging.error(f"{input_file} doesn't exist")
        sys.exit(1)
    tmpfile = input_file + '.tmp'
    
    align(input_file,tmpfile,ncpu = ncpu,mafft_opt = mafft_opt)
    trim(tmpfile,output_file,mode = clipkit_mode,g = clipkit_g)

    if clean:
        cmd = f"rm {tmpfile}"
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)

def phylogeny(fasta_file, output_prefix, cptime = 5000, nstop = 50, nm = 500, ntmax = 15, bb = 500, quiet = "",iqtree2 = "iqtree2"):
    # phylogeny(fasta_file, output_prefix, cptime = 1000, nstop = 100, nm = 10000, ntmax = 15, bb = 1000, quiet = "",iqtree2 = "iqtree2")
    logging.info(f"Phylogeny: {fasta_file} {output_prefix}")
    cmd = f"{iqtree2} -s {fasta_file} -m TEST -mset LG,WAG,JTT -nt AUTO -ntmax {ntmax} -bb {bb} -pre {output_prefix} -nm {nm} -nstop {nstop} -cptime {cptime} {quiet} --redo"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def possvm(treefile,output_prefix = None,reference_names = None, ogprefix = "OG", possvm = 'submodules/possvm-orthology/possvm.py'):
    logging.info(f"Possvm: {treefile}")
    # get the location of the possvm submodule 
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    possvm = scriptdir + '/../' + possvm
    
    if reference_names:
        reference_names = f"-r {reference_names}"
    else:
        reference_names = ""
    cmd = f"python {possvm} -ogprefix {ogprefix} -skipprint -method lpa -itermidroot 10 -i {treefile} {reference_names}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

# Phylo-search functions 

def blastp(query,target,db,outfile,ncpu=1,outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
    cmd = f"makeblastdb -in {target} -dbtype prot -out {db} > /dev/null 2>&1"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = f'blastp -query {query} -out {outfile} -db {db} -evalue 1e-5 -num_threads {ncpu} -outfmt "{outfmt}" > /dev/null 2>&1'
    #blastp -evalue 1e-5 -num_threads $NCPU -query $QUERY -db tmp/target -out search/${PREF}.blastp.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

def cluster(fasta_file,out_prefix,temp_dir):
    cmd = f"mmseqs easy-cluster {fasta_file} {out_prefix} {temp_dir} --cluster-reassign > /dev/null 2>&1"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

def get_fasta_names(fasta_file,out_file):
    cmd = f"grep '>' {fasta_file} | sed 's/>//g' | sort | uniq > {out_file}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def check_tempdir(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    else:
        print(f'Error: {dirpath} already exists. Delete it or set a different one!')
        quit()

def filter_clusters(cluster_file,fasta_file,query_ids_file,soi = None):
    # fasta_file should contain all the sequences 

    # Filtering :
    # The clusters should contain: the SOI and the query sequences
    logging.info('Filtering clusters:')
    print(cluster_file)
    print(fasta_file)
    print(query_ids_file)
    print(soi)

    print('filtering done')