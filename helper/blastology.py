import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
from helper import functions 


from helper.functions import align_and_trim
from helper.functions import phylogeny
from helper.functions import possvm
from helper.functions import blastp
from helper.functions import cluster

def filter_clusters(query,temp_dir,cluster_file,soi,require_soi,min_n,refnames_file,cluster_prefix,output_directory,prefix,joint_fasta_fname,cluster_directory,verbose = False):
    logging.info('Cluster filtering')
    import csv
    
    
    query_ids_file = os.path.join(temp_dir,'query.ids') 
    functions.get_fasta_names(fasta_file = query,out_file = query_ids_file,verbose = verbose)
    #functions.filter_clusters(cluster_file = cluster_file )
    with open(query_ids_file, "r") as file:
        query_ids = [line.strip() for line in file]

    clusters = {}
    with open(cluster_file, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for cluster_name, sequence_name in reader:
            clusters.setdefault(cluster_name, []).append(sequence_name)

    logging.info(f'Cluster filtering: {cluster_file} \nN query sequences: {len(query_ids)}\nN clusters: {len(clusters)}')
    # report
    clusters_query = [k for k,v in clusters.items() if any(elem in query_ids for elem in v)]
    clusters_small = [k for k,v in clusters.items() if len(v) < min_n]
    clusters_soi = [k for k,v in clusters.items() if any(soi in elem for elem in v)]
    clusters_query_small = [x for x in clusters_query if len(clusters[x]) < min_n]
    clusters_soi_small = [x for x in clusters_soi if len(clusters[x]) < min_n]
    logging.info(f'\nClusters with query: {len(clusters_query)}\nClusters with soi: {len(clusters_soi)}\nClusters small: {len(clusters_small)}\nClusters with query & small: {len(clusters_query_small)}\nClusters with SOI & small: {len(clusters_soi_small)}')


    # Filter by query sequences
    clusters_filt = [k for k,v in clusters.items() if any(elem in query_ids for elem in v) and len(v) >= min_n]
    cluters_filt = [x for x in clusters_filt if any(soi in elem for elem in clusters[x])]
    clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}

    logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters passing filtering.')
    clusters = clusters_filt_d
    
    if require_soi:
        if not soi == "":
            logging.info('Filtering by species of interest {soi}')
            clusters_filt = [k for k,v in clusters.items() if any(soi in elem  for elem in v)]
            logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters with SOI ({soi}) sequences.')
            clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}
            clusters = clusters_filt_d

    # rename each cluster and store the sequences 
    clusters_renamed = {}
    for i,(k,v) in enumerate(clusters.items()):
        clusters_renamed.update({cluster_prefix+str(i):v})

    # if provided the refnames, read in and store for each cluster 
    if refnames_file: 
        refnames = {}
        with open(refnames_file, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            for key, value in reader:
                refnames[key] = value
        cl_to_refname = {}
        for k,v in clusters_renamed.items():
            cl_to_refname.update({k:",".join(sorted([refnames[x] for x in v if x in refnames.keys()]))})
    
    # write down the filtered clustering file 
    cluster_tabfile = os.path.join(output_directory,prefix + '.clusters.tsv')
    with open(cluster_tabfile, "w") as file:
        for k,v in clusters_renamed.items():
            if refnames_file:
                file.write(k + "\t" + ",".join(v) + "\t" + cl_to_refname[k] + "\n")
            else:
                file.write(k + "\t" + ",".join(v) +  "\n")
    logging.info(f'Created: {cluster_tabfile}')


    # For each cluster, create a separate file:
    for cl_id in clusters_renamed.keys():
        fasta_file = joint_fasta_fname
        ids_to_keep = clusters_renamed[cl_id]  # Replace with your list of IDs
        cluster_fasta = cluster_directory + "/" + cl_id +  ".fasta"
        functions.retrive_sequences(joint_fasta_fname, cluster_fasta, ids_to_keep)

def parse_args(args):
    min_n = 3 # minimal number of sequences in the cluster
    cluster_prefix = args.prefix + '.' + args.cluster_prefix # a prefix to add to the cluster 
    output_directory = args.output_dir
    cluster_directory = os.path.join(output_directory,'clusters')
    
    query = args.query
    target = args.target
    temp_dir = args.temp_dir
    prefix = args.prefix
    soi = args.soi
    require_soi = bool(soi)
    
    force = args.force
    refnames_file = args.refnames

    ncpu = args.ncpu
    mafft = args.mafft
    phy_method = args.phymethod
    #globals().update(locals())
    #print(query,target,temp_dir,prefix,soi,force,refnames_file,ncpu,min_n,cluster_prefix,output_directory,cluster_directory,require_soi)
    return(query,target,temp_dir,prefix,soi,force,refnames_file,ncpu,min_n,cluster_prefix,output_directory,cluster_directory,require_soi,mafft,phy_method)

def join_seqs(query,target,joint_fasta_fname,joint_ids_fname,blastp_outfile,verbose):
    cmd = f'cat {query} {target} > {joint_fasta_fname}_tmp; samtools faidx {joint_fasta_fname}_tmp'
    if verbose:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"cat {blastp_outfile} | awk '{{print $1\"\\n\"$2}}' | sort | uniq > {joint_ids_fname}"
    if verbose:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = f'xargs samtools faidx {joint_fasta_fname}_tmp < {joint_ids_fname} > {joint_fasta_fname};rm {joint_fasta_fname}_tmp'
    if verbose:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    logging.info('Done collecting BLASTP results')

def run_cluster(cl_id,cluster_directory,refnames_file,mafft_opt,phy_method,force,ncpu,verbose):
        input_file = os.path.join(cluster_directory,cl_id)
        cluster_fasta = os.path.join(cluster_directory,cl_id +  ".fasta")

        fname_aln = os.path.splitext(cluster_fasta)[0] + '.aln'
        tree_prefix = os.path.splitext(cluster_fasta)[0] + '.tree'
        fname_tree = tree_prefix + ".treefile"
        fname_possvm = fname_tree + ".ortholog_groups.csv"
        logfile = os.path.join(cluster_directory,cl_id + '.log')
        
        if os.path.isfile(fname_aln) and not force:
            print(f'Found alignment file: {fname_aln}! Skipping alignment')
        else:
            align_and_trim(input_file = cluster_fasta, output_file = fname_aln, ncpu = ncpu, mafft_opt = mafft_opt, logfile = logfile,verbose = verbose)
        if os.path.isfile(fname_tree) and not force:
            print(f'Found phylogeny file: {fname_tree}! Skipping alignment')
        else:
            #functions.phylogeny_fasttree(fasta_file = fname_aln, output_file = fname_tree)
            phylogeny(fasta_file = fname_aln, output_prefix = tree_prefix,ntmax = ncpu,method = phy_method)

        if os.path.isfile(fname_possvm) and not force:
            print(f'Found POSSVM file: {fname_tree}! Skipping')
        else:
            possvm(treefile = fname_tree,reference_names = refnames_file,ogprefix = "OG",min_support_transfer = 50)
            logging.info(f'Created {fname_possvm}')

def blastology_run(args,logging,verbose = False):
    query,target,temp_dir,prefix,soi,force,refnames_file,ncpu,min_n,cluster_prefix,output_directory,cluster_directory,require_soi,mafft,phy_method = parse_args(args)

    # Directories
    # check the temporary directory status:
    functions.check_dir(temp_dir,force = force)
    functions.check_dir(output_directory,force = force)
    functions.check_dir(cluster_directory,force = force)
    
    
    # Intermediate files 
    blastp_outfile = os.path.join(temp_dir,f'{prefix}.blastp.tsv')
    cluster_file = os.path.join(temp_dir,f'{prefix}_cluster.tsv')
    
    joint_fasta_fname = os.path.join(temp_dir, f"{prefix}.fasta")
    joint_ids_fname = os.path.join(temp_dir, f"{prefix}.hits.ids")
    
    blastp_log = os.path.join(temp_dir, 'blastp.log')
    cluster_log = os.path.join(temp_dir,'cluster.log')
    
    # BLASTP
    if not os.path.isfile(blastp_outfile) or force:
        if verbose:
            logging.info(f'BLASTP:\n Query: {args.query}\n Target: {args.target}\n Threads: {args.ncpu}\nOutput{blastp_outfile}')
        blastp(query = args.query, target = args.target,db = temp_dir + "/target", outfile = blastp_outfile,ncpu = args.ncpu,logfile = blastp_log)
    else:
        if verbose:
            logging.info(f'Found blastp output file {blastp_outfile}. Skipping')


    # Join queries and targets 
    join_seqs(query,target,joint_fasta_fname,joint_ids_fname,blastp_outfile,verbose)
    
    # Clustering 
    if os.path.isfile(cluster_file) and not force:
        logging.info(f'Found clustering file {cluster_file}. Skipping')
    else:
        clustering_method = 'diamond_mcl' 
        cluster(fasta_file = joint_fasta_fname,out_prefix = temp_dir + '/' + prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = ncpu,method = clustering_method,mcl_inflation = "1.1")
    
    # Cluster filtering 
    filter_clusters(query = query, temp_dir = temp_dir, cluster_file = cluster_file, soi = soi, require_soi = require_soi,min_n = min_n, refnames_file = refnames_file, cluster_prefix = cluster_prefix, cluster_directory = cluster_directory,output_directory = output_directory, prefix = prefix, joint_fasta_fname = joint_fasta_fname, verbose = True)
    cluster_fastas = [file for file in os.listdir(cluster_directory) if file.endswith('.fasta')]
    cluster_prefs = [x.replace('.fasta','') for x in cluster_fastas]
    
    # Now, for each cluster, run the easy-phylo
    for cl_id in cluster_prefs:
        run_cluster(cl_id = cl_id,cluster_directory=cluster_directory,refnames_file=refnames_file,mafft_opt=mafft,phy_method=phy_method,force=force,ncpu=ncpu,verbose=verbose)
 
    # Finally, retrieve the true orthologs!