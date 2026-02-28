import os
import logging
from Bio import SeqIO
from helper.functions import cluster


def top_n(file_path):
    counts = {}
    with open(file_path, 'r') as file:
        for line in file:
            first_col = line.split('\t')[0]  # Get the first column
            counts[first_col] = counts.get(first_col, 0) + 1
        # Sort by counts (descending) and get the most common value
    return (len(counts),max(counts.values()))

def cluster_counts(file_path):
    counts = {}
    with open(file_path, 'r') as file:
        for line in file:
            first_col = line.split('\t')[0]  # Get the first column
            counts[first_col] = counts.get(first_col, 0) + 1
        # Sort by counts (descending) and get the most common value
    return (counts)

def top_n(file_path):
    counts = cluster_counts(file_path)
    return (len(counts),max(counts.values()))

def extract_sequences(fasta_file, ids_file, output_file):
    # Read the list of IDs from the file
    with open(ids_file, 'r') as f:
        ids_to_extract = set(line.strip() for line in f)
    
    # Parse the FASTA file and extract sequences
    with open(output_file, 'w') as output_handle:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in ids_to_extract:
                SeqIO.write(record, output_handle, 'fasta')
def cluster_dict(file_path):
    clusters = {}
    with open(file_path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if not cols:
                continue
            cluster = cols[0]
            clusters.setdefault(cluster, []).append(cols[1])
    return clusters


def recluster_global(info_pref,fasta_file,out_prefix,temp_dir,logfile,ncpu,clustering_method,cluster_prefix,inflation, max_N,max_iterations = 30,inflation_step = 0.1, verbose = False,logging = logging,cluster_log = None,inflation_max = 100):
    inflation = float(inflation) + float(inflation_step)
    max_N = int(max_N)
    logging.info(f'{info_pref}: Trying to recluster with higher inflation...')
    cluster_file = out_prefix + '_cluster.tsv'
    iteration = 0
    cluster(fasta_file = fasta_file,out_prefix = out_prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = ncpu,method = clustering_method, cluster_prefix = cluster_prefix, mcl_inflation = inflation, verbose = verbose, logging = logging)
    N,max_N_obs = top_n(cluster_file)
    while max_N_obs > max_N and iteration < max_iterations and inflation <= inflation_max:
        inflation += inflation_step
        inflation = round(inflation,1)
        cluster(fasta_file = fasta_file,out_prefix = out_prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = ncpu,method = clustering_method, cluster_prefix = cluster_prefix, mcl_inflation = inflation, verbose = verbose,logging = logging)
        N,max_N_obs = top_n(cluster_file)
        logging.info(f'{info_pref}: Iteration: {iteration} / {max_iterations}; Inflation: {inflation}; N clusters: {N}; N max: {max_N_obs}')
        iteration += 1
    
    # Evaluate
    if max_N_obs > max_N:
        if iteration >= max_iterations:
            logging.warning(
                f'{info_pref}: Max iterations {max_iterations} reached '
                f'and max cluster size ({max_N_obs}) is still > allowed ({max_N}).'
            )
        elif inflation >= inflation_max:
            logging.warning(
                f'{info_pref}: Maximum inflation {inflation_max} reached '
                f'and max cluster size ({max_N_obs}) is still > allowed ({max_N}).'
            )
        else:
            logging.warning(
                f'{info_pref}: Reclustering stopped but max cluster size '
                f'({max_N_obs}) still exceeds allowed ({max_N}).'
            )

        status = False

    else:
        logging.info(
            f'{info_pref}: Iterative clustering finished successfully: '
            f'Iteration {iteration}; Inflation {inflation}; '
            f'Max cluster size {max_N_obs}'
        )
        status = True
    return(cluster_file,status)

def recluster_hg_local(args,hg_id,sequence_ids,fasta_file,temp_dir,ncpu,max_N,inflation = 1.1,inflation_step = 0.1, max_iterations = 30, verbose = False,clustering_method = 'diamond_mcl'):

    temp_hg = os.path.join(temp_dir, hg_id)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(temp_hg, exist_ok = True)
    ids_file = os.path.join(temp_hg,'ids.txt')
    hg_fasta = os.path.join(temp_hg,'ids.fasta')
    with open(ids_file,"w") as f:
        f.write("\n".join(sequence_ids))
    logging.info(f'ids written to: {ids_file}')
    extract_sequences(ids_file = ids_file,fasta_file = fasta_file,output_file = hg_fasta)    
    # recluster the hgs 
    out_pref = temp_hg + '/cluster'
    logging.info(f'Re-clustering: {temp_hg}: {out_pref}')
    
    subcluster_file,status = recluster_global(info_pref = hg_id,fasta_file = hg_fasta,out_prefix = out_pref,
                    temp_dir = temp_hg , logfile = None,ncpu = ncpu,
                    clustering_method = clustering_method, max_N = max_N, max_iterations = max_iterations, 
                    inflation = inflation,inflation_step = inflation_step,cluster_prefix = hg_id + ".",
                    verbose = verbose)
    #subcluster_file = out_pref + '_cluster.tsv'
    counts = cluster_counts(subcluster_file)
    logging.info(f'Iterative reclustering done: {hg_id}: {subcluster_file}')
    
    return(subcluster_file,status)