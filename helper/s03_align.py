import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml

# Input: fasta file 
# Output: trimmed alignment 


def align(fasta_file, output_file, ncpu, mafft_mode = ""):

    logging.info(f"Aligning sequences with mafft")

    cmd = (f"mafft {fasta_file} {mafft_mode} > {output_file}")
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f"Alignment done: {output_file}")

def trim(input_file,output_file,mode = "kpic-gappy", g = "0.7"):
    cmd = f"clipkit {input_file} -m {mode} -o {output_file} -g {g}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def align_and_trim(input_file,output_file,ncpu = 1,clipkit_mode = "kpic-gappy",clipkit_g = 0.7, clean = True):
    if not os.path.exists(input_file):
        logging.error(f"{input_file} doesn't exist")
        sys.exit(1)
    tmpfile = input_file + '.tmp'
    
    align(input_file,tmpfile,ncpu = ncpu)
    trim(tmpfile,output_file,mode = clipkit_mode,g = clipkit_g)

    if clean:
        cmd = f"rm {tmpfile}"
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
