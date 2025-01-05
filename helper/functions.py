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

    cmd = (f"mafft --reorder --quiet --thread {ncpu} {mafft_opt} {fasta_file} > {output_file}")
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f"Alignment done: {output_file}")

def trim(input_file,output_file,mode = "kpic-gappy", g = "0.7",logfile = '/dev/null 2>&1'):
    cmd = f"clipkit {input_file} -m {mode} -o {output_file} -g {g} > {logfile}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def align_and_trim(input_file,output_file,ncpu = 1,mafft_opt = "",clipkit_mode = "kpic-gappy",clipkit_g = 0.7, clean = True, logfile = '/dev/null'):
    if not os.path.exists(input_file):
        logging.error(f"{input_file} doesn't exist")
        sys.exit(1)
    tmpfile = input_file + '.tmp'
    
    align(input_file,tmpfile,ncpu = ncpu,mafft_opt = mafft_opt)
    trim(tmpfile,output_file,mode = clipkit_mode,g = clipkit_g,logfile = logfile)

    if clean:
        cmd = f"rm {tmpfile}"
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)

# Phylogeny wrappers
def phylogeny(fasta_file,output_prefix,ntmax = 1,method = 'iqtree2'):
    if method == 'iqtree2':
        phylogeny_iqtree(fasta_file,output_prefix,ntmax = ntmax)
    elif method == 'fasttree':
        outfile = output_prefix + ".treefile"
        phylogeny_fasttree(fasta_file,outfile)

def phylogeny_iqtree(fasta_file, output_prefix, cptime = 5000, nstop = 50, nm = 500, ntmax = 15, bb = 1000, quiet = "",iqtree2 = "iqtree2",logfile = '/dev/null'):
    # Main output: {output_prefix}.treeflie
    # phylogeny(fasta_file, output_prefix, cptime = 1000, nstop = 100, nm = 10000, ntmax = 15, bb = 1000, quiet = "",iqtree2 = "iqtree2")
    logging.info(f"Phylogeny: {fasta_file} {output_prefix}")
    cmd = f"{iqtree2} -s {fasta_file} -m TEST -mset LG,WAG,JTT -nt AUTO -ntmax {ntmax} -bb {bb} -pre {output_prefix} -nm {nm} -nstop {nstop} -cptime {cptime} {quiet} --redo > {logfile} 2>&1"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f'Created {output_prefix}.treefile')

def phylogeny_fasttree(fasta_file, output_file):
    logging.info(f"Phylogeny: {fasta_file} {output_file}")
    cmd = f"fasttree -gtr -quiet {fasta_file} > {output_file}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f'Created {output_file}')


def possvm(treefile,output_prefix = None,reference_names = None, ogprefix = "OG", possvm = 'submodules/possvm-orthology/possvm.py',logfile = '/dev/null'):
    logging.info(f"Possvm: {treefile}\nLog: {logfile}")
    
    # get the location of the possvm submodule 
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    possvm = scriptdir + '/../' + possvm
    
    if reference_names:
        reference_names = f"-r {reference_names}"
    else:
        reference_names = ""
    cmd = f"python {possvm} -ogprefix {ogprefix} -skipprint -method lpa -itermidroot 10 -i {treefile} {reference_names} >> {logfile} 2>&1"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

# Phylo-search functions 

def blastp(query,target,db,outfile,ncpu=1,evalue = "1e-5",outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",logfile = '/dev/null'):
    cmd = f"makeblastdb -in {target} -dbtype prot -out {db} > {logfile}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = f'blastp -query {query} -out {outfile} -db {db} -evalue {evalue} -num_threads {ncpu} -outfmt "{outfmt}" >> {logfile} 2>&1'
    #blastp -evalue 1e-5 -num_threads $NCPU -query $QUERY -db tmp/target -out search/${PREF}.blastp.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def cluster(fasta_file,out_prefix,temp_dir,logfile = '/dev/null',method = 'mmseqs2',ncpu = 1,mcl_inflation = "1.1",cluster_prefix = "HG",verbose = False):
    if method == 'mmseqs2':
        cmd = f"mmseqs easy-cluster -s 7.5 --cov-mode 0 --cluster-mode 2 {fasta_file} {out_prefix} {temp_dir} --cluster-reassign >> {logfile} 2>&1"
        if verbose:
             logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
    elif method == 'diamond_mcl':
        cmd = f"diamond makedb --in {fasta_file} -d {fasta_file} --quiet"
        if verbose: 
           logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)

        diamond_max_target_seqs = 30
        cmd = f"diamond blastp --more-sensitive --max-target-seqs {diamond_max_target_seqs} -d {fasta_file} -q  {fasta_file} -o {out_prefix}_diamond.csv --quiet --threads {ncpu}"
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"awk '{{ print $1,$2,$12 }}' {out_prefix}_diamond.csv > {out_prefix}_diamond.abc"
        if verbose:
             logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True) 
        cmd = f"mcl {out_prefix}_diamond.abc --abc -I {mcl_inflation} -o {out_prefix}_mcl.tsv 2> /dev/null"
        
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)  
        cmd = f"""
        cat {out_prefix}_mcl.tsv | awk '{{ for (i = 1; i <= NF; i++) print "{cluster_prefix}"NR"\\t"$i }}' > {out_prefix}_cluster.tsv
        """
        if verbose:
             logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)  

    else:
        logging.info(f'Unknown clustering method {method}!')
        quit()

def get_fasta_names(fasta_file,out_file):
    cmd = f"grep '>' {fasta_file} | sed 's/>//g' | sort | uniq > {out_file}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def check_dir(dirpath,force = False):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
        print(f'Directory created: {dirpath}')
    else:
        if not force:
            print(f'Error: {dirpath} already exists. Delete it or use --force!')
            quit()
        else:
            print(f'Directory {dirpath} exists, but --force is set! Continuing ...')


from Bio import SeqIO
def retrive_sequences(input_fasta,output_fasta,ids_to_keep):
    with open(output_fasta, "w") as outfile:
         for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in ids_to_keep:
                SeqIO.write(record, outfile, "fasta")
    print(f'Created {output_fasta}')


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
