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

def phylogeny(fasta_file, output_prefix, cptime = 1000, nstop = 100, nm = 10000, ntmax = 15, bb = 1000, quiet = "",iqtree2 = "iqtree2"):
    logging.info(f"Phylogeny: {fasta_file} {output_prefix}")
    cmd = f"{iqtree2} -s {fasta_file} -m TEST -mset LG,WAG,JTT -nt AUTO -ntmax {ntmax} -bb {bb} -pre {output_prefix} -nm {nm} -nstop {nstop} -cptime {cptime} {quiet} --redo"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def possvm(treefile,output_prefix = None,reference_names = None, possvm = 'submodules/possvm-orthology/possvm.py'):
    logging.info(f"Possvm: {treefile}")
    # get the location of the possvm submodule 
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    possvm = scriptdir + '/../' + possvm
    
    if reference_names:
        reference_names = f"-r {reference_names}"
    else:
        reference_names = ""
    cmd = f"python {possvm} -skipprint -method lpa -itermidroot 10 -i {treefile} {reference_names}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
