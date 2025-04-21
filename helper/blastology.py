import os
import subprocess
import csv
import logging
from helper import functions 


from helper.functions import align_and_trim
from helper.functions import phylogeny
from helper.functions import possvm
from helper.functions import blastp
from helper.functions import cluster

# cluster function is imported from helper.functions 

def filter_clusters(query,temp_dir,cluster_file,soi,require_soi,min_n,refnames_file,cluster_prefix,output_directory,prefix,joint_fasta_fname,cluster_directory,verbose = False):
    logging.info('Cluster filtering')    
    query_ids_file = os.path.join(temp_dir,f'{prefix}_query.ids') 
    functions.get_fasta_names(fasta_file = query,out_file = query_ids_file,verbose = verbose)
    #functions.filter_clusters(cluster_file = cluster_file )
    with open(query_ids_file, "r") as file:
        query_ids = [line.strip() for line in file]

    clusters = {}
    with open(cluster_file, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for cluster_name, sequence_name in reader:
            clusters.setdefault(cluster_name, []).append(sequence_name)
    logging.info(f'{cluster_file}: {len(clusters)} clusters.')
    # report
    clusters_query = [k for k,v in clusters.items() if any(elem in query_ids for elem in v)]
    clusters_small = [k for k,v in clusters.items() if len(v) < min_n]
    clusters_soi = [k for k,v in clusters.items() if any(soi in elem for elem in v)]
    clusters_query_small = [x for x in clusters_query if len(clusters[x]) < min_n]
    clusters_soi_small = [x for x in clusters_soi if len(clusters[x]) < min_n]

    cluster_report = f'@N query sequences: {len(query_ids)}\n@Minimum cluster size allowed: {min_n}\n@Species of interest (SOI): {soi}\n\n@N clusters: {len(clusters)}\n@Clusters with query: {len(clusters_query)}\n@Clusters with SOI: {len(clusters_soi)}\n@Clusters small: {len(clusters_small)}\n@Clusters with query & small: {len(clusters_query_small)}\n@Clusters with SOI & small: {len(clusters_soi_small)}\n'
    cluster_report = cluster_report.replace('@',"\t\t\t     ")
    logging.info('Cluster filtering: \n' + cluster_report)


    # Filter by query sequences
    clusters_filt = [k for k,v in clusters.items() if any(elem in query_ids for elem in v) and len(v) >= min_n]
    cluters_filt = [x for x in clusters_filt if any(soi in elem for elem in clusters[x])]
    clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}

    logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters with queries passing size filtering.')
    clusters = clusters_filt_d
    
    if require_soi:
        if not soi == "":
            #logging.info(f'Filtering by species of interest {soi}')
            clusters_filt = [k for k,v in clusters.items() if any(soi in elem  for elem in v)]
            logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters with SOI ({soi}) sequences.')
            clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}
            clusters = clusters_filt_d
            if len(clusters_filt)==0:
                logging.error(f'No clusters containing the species of interest (--soi) {soi}. Aborting ...')          
                quit()

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

    cluster_fastas = []
    # For each cluster, create a separate file:
    for cl_id in clusters_renamed.keys():
        fasta_file = joint_fasta_fname
        ids_to_keep = clusters_renamed[cl_id]  # Replace with your list of IDs
        cluster_fasta = cluster_directory + "/" + cl_id +  ".fasta"
        cluster_fastas.append(cluster_fasta)
        functions.retrive_sequences(joint_fasta_fname, cluster_fasta, ids_to_keep)
    return(cluster_fastas)

def parse_args(args):
    
    print(args)
    
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
    mcl_inflation = args.mcl_inflation
    #globals().update(locals())
    #print(query,target,temp_dir,prefix,soi,force,refnames_file,ncpu,min_n,cluster_prefix,output_directory,cluster_directory,require_soi)
    return(query,target,temp_dir,prefix,soi,force,refnames_file,ncpu,min_n,cluster_prefix,output_directory,cluster_directory,require_soi,mafft,phy_method,mcl_inflation)

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

def get_results(cluster_directory,prefix,query_ids,soi = None,output_file = None,verbose = True):
    if not output_file:
        logging.error('Specify output file!!!')
        quit()
    cmd = f'cat  {cluster_directory}/{prefix}*groups.csv | grep -f <(cat {cluster_directory}/{prefix}*groups.csv | grep -f {query_ids} | cut -f 2 | sort | uniq)'
    if soi:
        cmd = cmd + f' | grep {soi} > {output_file}'
    else:
        cmd = cmd + f' > {output_file}'
    if verbose > 1:
        logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

def run_cluster(cl_id = None,cluster_directory = None,refnames_file = None,prefix = None,mafft_opt = None,phy_method = 'fasttree',force = True,ncpu = 1,verbose = False):
        input_file = os.path.join(cluster_directory,cl_id)
        cluster_fasta = os.path.join(cluster_directory,cl_id +  ".fasta")

        fname_aln = os.path.splitext(cluster_fasta)[0] + '.aln'
        #tree_prefix = os.path.splitext(cluster_fasta)[0] + '.tree'
        tree_prefix = os.path.splitext(cluster_fasta)[0]
        fname_tree = tree_prefix + ".tree"
        fname_possvm = fname_tree + ".ortholog_groups.csv"
        logfile = os.path.join(cluster_directory,cl_id + '.log')
        
        if os.path.isfile(fname_aln) and not force:
            logging.info(f'Found alignment file: {fname_aln}! Skipping alignment')
        else:
            align_and_trim(input_file = cluster_fasta, output_file = fname_aln, ncpu = ncpu, mafft_opt = mafft_opt, logfile = logfile,verbose = verbose)
        if os.path.isfile(fname_tree) and not force:
            logging.info(f'Found phylogeny file: {fname_tree}! Skipping alignment')
        else:
            #functions.phylogeny_fasttree(fasta_file = fname_aln, output_file = fname_tree)
            phylogeny(fasta_file = fname_aln, output_file = fname_tree,ntmax = ncpu,method = phy_method)

        if os.path.isfile(fname_possvm) and not force:
            logging.info(f'Found POSSVM file: {fname_tree}! Skipping')
        else:
            og_pref = f'{prefix}.OG' if prefix else 'OG'
            possvm(treefile = fname_tree,reference_names = refnames_file,ogprefix = og_pref,min_support_transfer = 50)
            logging.info(f'Created {fname_possvm}')

def blastology_run(args,logging,verbose = False):
    # this horrible argument parsing should be improved!!!
    query,target,temp_dir,prefix,soi,force,refnames_file,ncpu,min_n,cluster_prefix,output_directory,cluster_directory,require_soi,mafft,phy_method,mcl_inflation = parse_args(args)
    output_file = args.outputfile
    prefix = args.prefix
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
        #if verbose:
            #logging.info(f'BLASTP:\n Query: {args.query}\n Target: {args.target}\n Threads: {args.ncpu}\nOutput{blastp_outfile}')
        blastp(query = args.query, target = args.target,db = temp_dir + "/target", outfile = blastp_outfile,ncpu = args.ncpu,logfile = blastp_log, evalue=args.evalue,min_perc = args.min_perc,verbose = verbose)
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
        logging.info(f'Running MCL clustering of {joint_fasta_fname} with inflation {mcl_inflation} ...')
        cluster(fasta_file = joint_fasta_fname,out_prefix = temp_dir + '/' + prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = ncpu,method = clustering_method,mcl_inflation = mcl_inflation,verbose = verbose)
    # Cluster filtering 
    cluster_fastas = filter_clusters(query = query, temp_dir = temp_dir, cluster_file = cluster_file, soi = soi, require_soi = require_soi,min_n = min_n, refnames_file = refnames_file, cluster_prefix = cluster_prefix, cluster_directory = cluster_directory,output_directory = output_directory, prefix = prefix, joint_fasta_fname = joint_fasta_fname, verbose = True)
    cluster_fastas = [os.path.basename(x) for x in cluster_fastas]
    cluster_prefs = [x.replace('.fasta','') for x in cluster_fastas]
    
    # Now, for each cluster, run the easy-phylo
    logging.info(f'{len(cluster_prefs)} sequence clusters to process: {",".join(cluster_prefs)}')
    
    for cl_id in cluster_prefs:
        run_cluster(cl_id = cl_id,cluster_directory=cluster_directory,refnames_file=refnames_file,prefix = prefix,mafft_opt=mafft,phy_method=phy_method,force=force,ncpu=ncpu,verbose=verbose)

    # Finally, retrieve the true orthologs!
    #if not args.outputfile:
    #    logging.error('specify output file')
    #    if not args.soi:
    #        output_file = f'{prefix}_annotations.tsv'
    #    else:
    #        output_file = f'{prefix}_{soi}_annotations.tsv'

    get_results(cluster_directory=cluster_directory,prefix = prefix,query_ids =  os.path.join(temp_dir,f'{prefix}_query.ids') ,soi = args.soi,output_file = output_file,verbose = verbose)
    logging.info(f'Created {output_file}')
