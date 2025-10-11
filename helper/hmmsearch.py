import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
from Bio import SeqIO

def activate_conda_environment(env_name):
    logging.info(f"Activating conda environment: {env_name}")
    subprocess.run(f"conda activate {env_name}", shell=True, check=True)

def extract_sequences(fasta_file, ids_file, output_file):
    # Read the list of IDs from the file
    with open(ids_file, 'r') as f:
        ids_to_extract = set(line.strip() for line in f)
    
    # Parse the FASTA file and extract sequences
    with open(output_file, 'w') as output_handle:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in ids_to_extract:
                SeqIO.write(record, output_handle, 'fasta')

def fetch_hmm(hmm,pfam_db,outfile):
    if not os.path.exists(pfam_db):
        logging.error(f"Specified PFAM database ({pfam_db}) doesn't exist. Specify using --pfam_db option!")
        sys.exit(1)
    cmd = f"hmmfetch {pfam_db} {hmm} > {outfile}"
    logging.info(cmd)
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        print(f"Return code: {e.returncode}")
        print(f"Command: {e.cmd}")

    if os.path.isfile(outfile) and os.path.getsize(outfile) == 0:
        logging.error(f"""
                {outfile} is empty. 
                This means PFAM db file doesn't contain {hmm}!
                Check the db file {pfam_db} or the .hmm model name spelling!""")
        sys.exit(1)
    else:
        logging.info(f"Fetched {hmm} from {pfam_db}: {outfile}")


def do_hmmsearch(hmm,hmm_dir, out, fasta, cpu, threshold, pfam_db = None,verbose = False):
    logging.info("Running HMMSEARCH")
    logging.info(f"hmm: {hmm}")
    logging.info(f"hmm_dir: {hmm_dir}")
    logging.info(f"out: {out}")
    logging.info(f"fasta: {fasta}")
    logging.info(f"cpu: {cpu}")
    logging.info(f"threshold: {threshold}")

    hmm_path = os.path.join(hmm_dir,hmm + '.hmm')

    if not os.path.exists(hmm_dir):
        os.makedirs(hmm_dir)
    if not os.path.exists(hmm_path) or os.path.getsize(hmm_path) == 0:
        logging.error(f"HMM model file {hmm_path} not found! Trying to fetch from {pfam_db} ...")
        if not pfam_db:
            logging.error(f"Please, specify Pfam-A.hmm location using --pfam_db option! I will try to fetch {hmm}.hmm from it. Thank you!")
            sys.exit(1)
        fetch_hmm(hmm,pfam_db,hmm_path)
    else:
        logging.info(f'Found {hmm_path}')
    
    if threshold == "GA":
        cmd = (
            f"hmmsearch --domtblout {out}.domtable --cut_ga "
            f"--cpu {cpu} {hmm_path} {fasta} 1> /dev/null"
        )
    else:
        cmd = (
            f"hmmsearch --domtblout {out}.domtable --domE {threshold} "
            f"--cpu {cpu} {hmm_path} {fasta} 1> /dev/null"
        )
    
    subprocess.run(cmd, shell=True, check=True)
    
    hmmsearch_outfile = f"{out}.domtable.csv.tmp"
    with open(f"{out}.domtable", 'r') as infile, open(hmmsearch_outfile, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                columns = line.split()
                outfile.write('\t'.join([columns[0], columns[17], columns[18], columns[3],columns[4],columns[11]]) + '\n')
    if verbose:
        print(f"Hmmsearch created: {hmmsearch_outfile}")

def merge_results(prefix,gene_family_name, searches_dir, fasta_file, domain_expand = 50,verbose = False, do_clean = False):
    #prefix = gene_family_name
    fam_id = gene_family_name
    logging.info(f'Merging results ...')
    # Check for any hits
    cmd = f"cat {searches_dir}/{prefix}.{fam_id}.hmmsearch.*.domtable.csv.tmp | grep -v '#' | cut -f 1 | sort -u > {searches_dir}/{prefix}.{fam_id}.genes.list"
    subprocess.run(cmd, shell=True, check=True)
    num_hit_genes_cmd = f"cat {searches_dir}/{prefix}.{fam_id}.genes.list | wc -l"
    num_hit_genes = int(subprocess.check_output(num_hit_genes_cmd, shell=True).strip().decode('utf-8'))
    logging.info(f"# {fasta_file}: {fam_id} | # genes found = {num_hit_genes}")

    if num_hit_genes == 0:
        logging.info(f"# {fasta_file}: {fam_id} | Omit downstream analyses")
    else:
        # Find most-inclusive region that includes all domain hits in the protein
        cmd = f"cat {searches_dir}/{prefix}.{fam_id}.hmmsearch.*.domtable.csv.tmp > {searches_dir}/{prefix}.{fam_id}.domtable.csv.tmp"
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
        cmd = (
            f"bedtools merge -i <(sort -k1,1 -k2,2n {searches_dir}/{prefix}.{fam_id}.domtable.csv.tmp) "
            f"-c 4 -o collapse -d 100000 > {searches_dir}/{prefix}.{fam_id}.domains.csv.tmp2 "
        )
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
        # 4.01.25 - issue with I/O on .fai
        # Extract complete sequences
        #cmd = f"samtools faidx {fasta_file} -r {searches_dir}/{prefix}.{fam_id}.genes.list > {searches_dir}/{prefix}.{fam_id}.seqs.fasta"
        #cmd = f'xargs  -P 1 -I {{}} samtools faidx {fasta_file} < {searches_dir}/{prefix}.{fam_id}.genes.list > {searches_dir}/{prefix}.{fam_id}.seqs.fasta'
        #cmd = f'python get_seq.py {fasta_file} {searches_dir}/{prefix}.{fam_id}.genes.list {searches_dir}/{prefix}.{fam_id}.seqs.fasta'
        
        #if verbose:
        #    logging.info(cmd)
        #subprocess.run(cmd, shell=True, check=True)
        
        # Solution that doesn't require samtools 
        extract_sequences(fasta_file,f"{searches_dir}/{prefix}.{fam_id}.genes.list",f"{searches_dir}/{prefix}.{fam_id}.seqs.fasta")
        # Dict sequence lengths
        cmd = f"samtools faidx {searches_dir}/{prefix}.{fam_id}.seqs.fasta"
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f"Expanding domain ranges: {domain_expand}")
        # Expand domain region by a fixed amount of aa
        cmd = (
            f"bedtools slop -i {searches_dir}/{prefix}.{fam_id}.domains.csv.tmp2 "
            f"-g {searches_dir}/{prefix}.{fam_id}.seqs.fasta.fai -b {domain_expand} "
            f"| awk 'BEGIN{{OFS=@\\t@}}{{ print $1, $2+1, $3, $4 }}' > {searches_dir}/{prefix}.{fam_id}.domains.csv "
        )
        cmd = cmd.replace("@",'"')
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
        # Extract domain region
        cmd = (
            f"bedtools getfasta -fi {searches_dir}/{prefix}.{fam_id}.seqs.fasta -bed "
            f"<(awk 'BEGIN{{OFS=@\\t@}}{{ print $1, $2, $3, $1 }}' {searches_dir}/{prefix}.{fam_id}.domains.csv) "
            f" | sed -E 's/:[0-9]+-[0-9]+$//g' > {searches_dir}/{prefix}.{fam_id}.domains.fasta"
        )
        cmd = cmd.replace("@",'"')
        if verbose:
            logging.info(cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
        # Strip domain coordinates from the domain file

        # Report
        num_unique_domains_cmd = f"cut -f1 {searches_dir}/{prefix}.{fam_id}.domains.csv | wc -l"
        num_unique_domains = subprocess.check_output(num_unique_domains_cmd, shell=True).strip().decode('utf-8')
        logging.info(f"{fasta_file}: {fam_id} | # unique domains = {num_unique_domains}")
        if do_clean:
            def list_temp_files(directory, prefix, endings):
                matching_files = []
                for file in os.listdir(directory):
                    if file.endswith(tuple(endings)) and file.startswith(prefix):
                        matching_files.append(directory + '/' + file)
                return matching_files
            

            endings = [".tmp","tmp2","dmnd","domtable","genes.list"]
            temp_files = list_temp_files(searches_dir,f'{prefix}.{fam_id}',endings)
            logging.info("Removing temporary files...")
            def remove_files(file_list):
                for file in file_list:
                    try:
                        os.remove(file)
                        #print(f"Removed: {file}")
                    except FileNotFoundError:
                        print(f"File not found: {file}")
                    except PermissionError:
                        print(f"Permission denied: {file}")
                    except Exception as e:
                        print(f"Error removing {file}: {e}")
            remove_files(temp_files)

def parse_gene_family_info(gene_family_info):
    gene_families = {}
    with open(gene_family_info, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t', fieldnames=["class_name", "hmms", "inflation", "min_seqs", "threshold", "pref1", "pref2"])
        for row in reader:
            class_name = row["class_name"]
            gene_families[class_name] = {
                "hmms": row["hmms"].split(','),
                "inflation": row["inflation"],
                "min_seqs": row["min_seqs"],
                "threshold": row["threshold"],
                "pref1": row["pref1"],
                "pref2": row["pref2"]
            }
    return gene_families

def hmmsearch(fasta_file, gene_family_info, gene_family_name, output_dir, pfam_db, hmm_dir,ncpu, domain_expand = 50,verbose = 1, do_clean = True):
    logging.info(f"# {fasta_file}: {gene_family_name} | HMM search")
    gene_families = parse_gene_family_info(gene_family_info)
   # if verbose: 
        #print(gene_families)
    if gene_family_name not in gene_families:
        logging.error(f"Gene family {gene_family_name} not found in gene family info file.")
        return

    gene_family = gene_families[gene_family_name]
    prefix = gene_family["pref2"]
    hmms = gene_family["hmms"]
    threshold = gene_family["threshold"]
    
    cpu = ncpu
    hmm_dir = hmm_dir
    domain_expand = int(domain_expand)

    searches_dir = output_dir
    os.makedirs(searches_dir, exist_ok=True)
    
    tmp_file = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.domains.csv.tmp")
    with open(tmp_file, 'w') as outfile:
        for hmm in hmms:
            logging.info(f"Running HMM search for {gene_family_name} with HMM: {hmm}")
            output_pref = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.hmmsearch.domain_{hmm}")
            do_hmmsearch(hmm,hmm_dir, output_pref, fasta_file, cpu, threshold,pfam_db,verbose = verbose)
            # Concatenate results
            output_file = f'{output_pref}.domtable'
            with open(output_file, 'r') as infile:
                for line in infile:
                    outfile.write(line)
    # Check for any hits
    genes_list_file = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.genes.list")
    cmd = f"grep -v '#' {tmp_file} | cut -f 1 -d ' ' | sort -u > {genes_list_file}"
    subprocess.run(cmd, shell=True, check=True)
    num_hit_genes_cmd = f"cat {genes_list_file} | wc -l"
    num_hit_genes = int(subprocess.check_output(num_hit_genes_cmd, shell=True).strip().decode('utf-8'))
    logging.info(f"# {fasta_file}: {gene_family_name} | # genes found = {num_hit_genes}")

    if num_hit_genes == 0:
        logging.info(f"# {fasta_file}: {gene_family_name} | Omit downstream analyses")
    else:
        merge_results(prefix,gene_family_name, searches_dir, fasta_file, domain_expand,verbose = verbose, do_clean = do_clean)
    logging.info(f'Hmmsearch done: {output_dir}/')

