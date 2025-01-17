import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
import ete3


#### Tool functions #####
def align(fasta_file, output_file, ncpu, mafft_opt,verbose):
    if verbose:
        n  = count_seqs(fasta_file,verbose)
        logging.info(f"Aligning {n} sequences with mafft")

    cmd = (f"mafft --reorder --quiet --thread {ncpu} {mafft_opt} {fasta_file} > {output_file}")
    if verbose:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    if verbose:
        logging.info(f"Alignment done: {output_file}")

def trim(input_file,output_file,mode = "kpic-gappy", g = "0.7",logfile = '/dev/null 2>&1',verbose = True):
    cmd = f"clipkit {input_file} -m {mode} -o {output_file} -g {g} > {logfile}"
    if verbose:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)


def align_and_trim(input_file,output_file,ncpu = 1,mafft_opt = "",clipkit_mode = "kpic-gappy",clipkit_g = 0.7, clean = True, logfile = '/dev/null',verbose = True):
    if not os.path.exists(input_file):
        logging.error(f"{input_file} doesn't exist")
        sys.exit(1)
    tmpfile = input_file + '.tmp'
    
    align(input_file,tmpfile,ncpu = ncpu,mafft_opt = mafft_opt,verbose = verbose)
    trim(tmpfile,output_file,mode = clipkit_mode,g = clipkit_g,logfile = logfile,verbose = verbose)

    if clean:
        cmd = f"rm {tmpfile}"
        #logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)

# Phylogeny wrappers
def get_node_support_range(treefile):
    # Load the tree
    #tree = Tree(treefile, format=1)  # format=1 assumes Newick with bootstrap values
    tree = ete3.PhyloTree("%s" % (treefile))
    # Collect support values from internal nodes

    support_values = [node.support for node in tree.traverse() if not node.is_leaf() and node.support is not None]
    min_support = min(support_values)
    max_support = max(support_values)
    return(min_support,max_support)

def phylogeny(fasta_file,output_file,ntmax = 1,method = 'iqtree2'):
    if method == 'iqtree2':
        phylogeny_iqtree(fasta_file,output_file,ntmax = ntmax)
    elif method == 'fasttree':
        #outfile = output_prefix + ".treefile"
        phylogeny_fasttree(fasta_file,output_file)
    if not os.path.isfile(output_file):
        logging.error('Phylogeny has failed. Can not find {output_file}! Aborting ...')
        quit()

def phylogeny_iqtree(fasta_file, output_file, cptime = 5000, nstop = 50, nm = 500, ntmax = 15, bb = 1000, quiet = "",iqtree2 = "iqtree2",logfile = '/dev/null',verbose = True):
    # Main output: {output_prefix}.treeflie
    # phylogeny(fasta_file, output_prefix, cptime = 1000, nstop = 100, nm = 10000, ntmax = 15, bb = 1000, quiet = "",iqtree2 = "iqtree2")

    # iqtree creates the files given a prefix {PREFIX}.treefile 
    # If the output file name provided and ends in .tree - use as a prefix 
    logging.info(f"Phylogeny: {fasta_file} {output_file}")
    print(output_file)
    if output_file.endswith('.treefile'):
        output_prefix = output_file.replace('.treefile','')
        if verbose:
            logging.info(f"IQTREE2: outputfile ends with .treefile => {output_prefix}")
    elif output_file.endswith('.tree'):
        output_prefix = output_file.replace('.tree','')
        if verbose:
            logging.info(f"IQTREE2: outputfile ends with .treefile => {output_prefix}")
    else:
        output_prefix = output_file
        if verbose:
            logging.info(f"IQTREE2: output prefix: {output_prefix}")
    cmd = f"{iqtree2} -s {fasta_file} -m TEST -mset LG,WAG,JTT -nt AUTO -ntmax {ntmax} -bb {bb} -pre {output_prefix} -nm {nm} -nstop {nstop} -cptime {cptime} {quiet} --redo > {logfile} 2>&1"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    #logging.info(f'IQTREE2: Created {output_prefix}.treefile')
    cmd = f"mv {output_prefix}.treefile {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def phylogeny_fasttree(fasta_file, output_file):
    logging.info(f"Phylogeny: {fasta_file} {output_file}")
    cmd = f"fasttree -gtr -quiet {fasta_file} > {output_file}"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f'Created {output_file}')


def possvm(treefile,output_prefix = None,reference_names = None, ogprefix = "OG", possvm = 'submodules/possvm-orthology/possvm.py',logfile = '/dev/null',refsps = None,min_support_transfer = 50):
    logging.info(f"Possvm: {treefile}\nLog: {logfile}")
    # Adjust min_support according to the value range in provided tree 
    print(min_support_transfer)
    nsr = get_node_support_range(treefile)
    if min_support_transfer:
        if min_support_transfer > nsr[1]:
            logging.info(f"Minimum node support ({min_support_transfer}) is bigger than the maximum observed support value ({nsr[1]}); Adjusting the threshold to {round(min_support_transfer/100,2)}") 
            min_support_transfer = min_support_transfer / 100
    # get the location of the possvm submodule 
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    possvm = scriptdir + '/../' + possvm
    
    if reference_names:
        reference_names = f"-r {reference_names}"
    else:
        reference_names = ""
    
    if refsps:
        reference_species = f"-refsps {refsps}" 
    else:
        reference_species = ""
    cmd = f"python {possvm} -ogprefix {ogprefix} -skipprint -method lpa -itermidroot 10 -min_support_transfer {min_support_transfer}  -i {treefile} {reference_names} {reference_species} >> {logfile} 2>&1"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

# Phylo-search functions 

def blastp(query,target,db,outfile,ncpu=1,evalue = "1e-5",min_perc = None,outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",logfile = '/dev/null',verbose = False):
    if not os.path.isfile(f'{db}.pdb'):
        cmd = f"makeblastdb -in {target} -dbtype prot -out {db} > {logfile}"
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
    else:
        logging.info(f'Found db files {db}. Skipping db building')
    cmd = f'blastp -query {query} -out {outfile} -db {db} -evalue {evalue} -num_threads {ncpu} -outfmt "{outfmt}" >> {logfile} 2>&1'
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    if min_perc:
        logging.info(f'Minimum BLASTP hit percetage is set to {min_perc}. Filtering')
        cmd = f"cat {outfile} | awk '$3>={min_perc}' > {outfile}.filtered; mv {outfile}.filtered {outfile}"
        if verbose:
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

def count_seqs(fasta_file, verbose=False):
    cmd = f"grep -c '>' {fasta_file}"
    if verbose > 1:
        logging.info(cmd)
    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    return int(result.stdout.strip())

def get_fasta_names(fasta_file,out_file,verbose = False):
    cmd = f"grep '>' {fasta_file} | sed 's/>//g' | sort | uniq > {out_file}"
    if verbose > 1:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

def check_tool(tool_name):
    try:
        # Try to run the tool with a harmless argument like --help
        cmd = [tool_name, '--help']
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # If return code is 0, it means the tool executed successfully
        if result.returncode == 0:
            print(f"{tool_name} is available and can be launched.")
        else:
            print(f"Error: {tool_name} is not functioning properly.")
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: {tool_name} is not available on your system.")
        sys.exit(1)

def check_dir(dirpath,force = False):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
        logging.info(f'Directory created: {dirpath}')
    else:
        if not force:
            print(f'Directory {dirpath} exists, but --force is not set! Continuing ...')
            


from Bio import SeqIO
def retrive_sequences(input_fasta,output_fasta,ids_to_keep):
    with open(output_fasta, "w") as outfile:
         for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in ids_to_keep:
                SeqIO.write(record, outfile, "fasta")
    logging.info(f'Created {output_fasta}')


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
