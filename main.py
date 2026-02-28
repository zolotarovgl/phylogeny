import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
from Bio import SeqIO

from helper.hmmsearch import hmmsearch
#from helper.s02_cluster import cluster 
from helper.functions import align_and_trim
from helper.functions import phylogeny
from helper.functions import possvm
from helper.functions import blastp
from helper.functions import cluster

from helper import functions

# Ensure PyYAML is installed
try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is not installed. Install it using 'pip install pyyaml'.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
def load_config():
    config_file = "config.yaml"
    tool_directory = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(tool_directory, config_file)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file {config_path} not found.")
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

# Pipelines # 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Python wrapper around some useful commands
    """)
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # HMM Search 
    parser_search = subparsers.add_parser('hmmsearch', help='Search for a family using HMMER')
    parser_search.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_search.add_argument('-g', 
        '--gene_family_info', 
        required=True, 
        help="""
            Path to the gene family info file specifying HMMs and parameters.
            See example file data/genefam.tsv.
            """)
    parser_search.add_argument('gene_family_name', help='Name of the gene family to search')
    parser_search.add_argument('-o','--output_dir',required=True, help='Output directory')
    parser_search.add_argument('--hmm_dir',required=False, default = None, help='HMM directory')
    parser_search.add_argument('--pfam_db',required=False, default = None, help='Path to Pfam-A.hmm file. Default: None')
    parser_search.add_argument('--domain_expand',required=False, default = "50", help='Expand domain ranges to X aminoacids in both directions. Default: 50')
    parser_search.add_argument('-c', '--ncpu', required=False, default = int(1),  help='Number of CPU cores to use')
    parser_search.add_argument('--keep', required=False, default = False, action = 'store_true', help='Use this to keep temporary files')
    
    # Cluster
    parser_cluster = subparsers.add_parser('cluster', 
                                           help="Run clustering of the sequences",
                                           description = """Cluster sequences in provided fasta file.
                                           If the top cluster contains >= max_n sequences, recluster.""")
    parser_cluster.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_cluster.add_argument('--out_file', required=False, default = None, help='Output clustering file. CAVE: for legacy reasons, should be named PREFIX_cluster.tsv')
    parser_cluster.add_argument('--out_prefix', required=False,default = None, help='Output file prefix. Vis --out_file.')
    parser_cluster.add_argument('-c', '--ncpu', required=False, default = int(1), help='Number of CPU cores to use. Default: 1')
    parser_cluster.add_argument('-t', '--temp_dir', required=False, default = "tmp/", help='Temporary directory name. Default: tmp/')
    parser_cluster.add_argument('-i', '--inflation', default = float(1.1), help='Inflation parameter for MCL clustering')
    parser_cluster.add_argument('--inflation_step', default = float(0.1), help='Inflation parameter step for re-clustering')
    parser_cluster.add_argument('-m', '--maxn', default = int(1000), help='Maximum number of sequences in the cluster. Default: 1000')
    parser_cluster.add_argument('--method', default = "diamond_mcl", help='Clustering method. Default: diamond_mcl')
    parser_cluster.add_argument('--cluster_prefix', default = "HG", help='Cluster name prefix. Default: HG (i.e. HG1, HG2, ...)')
    

    # Alignment
    parser_align = subparsers.add_parser('align', help='Run alignment')
    parser_align.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_align.add_argument('-o', '--outfile', required=True, help='Output file')
    parser_align.add_argument('-c', '--ncpu', required=False, default = int(1), help='Number of CPU cores to use')
    parser_align.add_argument('-m', '--mafft', required=False, default ="", help='Mafft alignment options. Default  "", you can set it to --maxiterate 1000 --genafpair')
    parser_align.add_argument('--notrim', required=False, default = False, action = 'store_true', help='Use this to skip the alignment trimming step')

    # Trimming
    parser_trim = subparsers.add_parser('trim', help='Trim alignment')
    parser_trim.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_trim.add_argument('-o', '--outfile', required=True, help='Output file')
    parser_trim.add_argument('-l', '--logfile', required=False, default = None, help='Log file')
    parser_trim.add_argument('-m', '--method', required=False, default ="clipkit", help='Alignment trimming method. Default: clipkit')

    # Phylogeny
    parser_phylogeny = subparsers.add_parser('phylogeny', help='Run IQTREE2 for an alignment in --fasta')
    parser_phylogeny.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_phylogeny.add_argument('--outprefix', required=False, help='Output prefix for phylogeny files')
    parser_phylogeny.add_argument('--outfile', required=False, help='Output file name')
    parser_phylogeny.add_argument('-c', '--ncpu', required=True,  help='Number of CPU cores to use')
    parser_phylogeny.add_argument('--method', required=False, default = "fasttree",  help='Phylogeny method: [fasttree, iqtree2, iqtree3]. Default: fasttree')
    parser_phylogeny.add_argument('--iqtree2_model', required=False, default = "TEST",  help='IQTREE2 model. Default: TEST')

    # GeneRax
    parser_generax = subparsers.add_parser('generax', help='Run GeneRax')
    parser_generax.add_argument('--name', required=True, help='Family name')
    parser_generax.add_argument('--alignment', required=True, help='Path to alignment file')
    parser_generax.add_argument('--gene_tree', required=True, help='Path to gene tree file')
    parser_generax.add_argument('--species_tree', required=True, help='Path to species tree file')
    parser_generax.add_argument('--output_dir', required=True, help='Path to output directory')
    parser_generax.add_argument('--subs_model', help='Substitution model (e.g. LG+G)')
    parser_generax.add_argument('--iqtree_file', help='Optional IQ-TREE .iqtree file to extract model')
    parser_generax.add_argument('--per-family-rates', required=False, default = True, action = 'store_true', help='Whether to use per family rates')
    parser_generax.add_argument('--max_spr', required=False, default = int(5), help='Maximum SPR radius')
    parser_generax.add_argument('-c','--cpus', required=False, default = int(1), help='Number of CPU cores')
    parser_generax.add_argument('-o','--outfile', required=False, default = None, help='Name of the output tree file')
    parser_generax.add_argument('-l','--logfile', default = None, help='the log')


    # POSSVM
    parser_possvm = subparsers.add_parser('possvm', help='Run POSSVM')
    parser_possvm.add_argument('-t','--treefile', required = True, help='Input tree file')
    parser_possvm.add_argument('-r','--refnames', default = None, help='Reference gene names: gene \t name.')
    parser_possvm.add_argument('-o','--ogprefix', default = "OG", help='String. Prefix for ortholog clusters. Defaults to "OG".')
    parser_possvm.add_argument('-s','--refsps', default = None, help='POSSVM reference species')
    parser_possvm.add_argument('--sos', default = 0, help='POSSVM species overlap (--sos) param')
    parser_possvm.add_argument('--min_support_transfer', default = "50",dest = "possvm_minsupport", help='POSSVM Minimum support for label transfer')
    parser_possvm.add_argument('--itermidroot', default = "10", help='Number of rooting iterations')
    parser_possvm.add_argument('-l','--logfile', default = "/dev/null", help='the log')
    parser_possvm.add_argument('--outgroup', default = "", help='POSSVM: outgroup species file.')
   
    
    # EASY-PHYLO
    parser_easyphylo = subparsers.add_parser('easy-phylo',help = 'Build a phylogeny from a single fasta')
    parser_easyphylo.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_easyphylo.add_argument('-c','--ncpu', required=True, help='Number of CPU cores to use')
#    parser_easyphylo.add_argument('-m', '--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_easyphylo.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_easyphylo.add_argument('--outdir', default = None, help='Optional: output directory. By default, the output files are written to the directory of the input file')
    parser_easyphylo.add_argument('--ogprefix', default = "OG", help='POSSVM: String. Prefix for ortholog clusters. Defaults to "OG".')
    parser_easyphylo.add_argument('--sos', default = "0", help='POSSVM: --sos parameter. [0,1]. Default: 0')
    parser_easyphylo.add_argument('--force', required=False, default = False, action = 'store_true', help='Use this to rerun intermediate files (e.g. alignment)')
    parser_easyphylo.add_argument('--method', default = "iqtree3", help='Phylogeny method: fasttree, iqtree2, iqtree3. Default: iqtree3')
    parser_easyphylo.add_argument('--min_support_transfer', default = "50", dest = "easyphylo_minsupport", help='POSSVM Minimum support for label transfer')
    parser_easyphylo.add_argument('--mafft', required=False, default ="auto", help='Mafft alignment options. Default: auto - picks based on the number of sequences.\nAvailable options: auto, fast, linsi,einsi,ginsi')
    

    # BLASTOLOGY - phylosearch v2
    parser_blastology = subparsers.add_parser('blastology',
            help = 'Search query sequences in proteomes and annotate using phylogenies',
            description="""
    Annotate the proteins based on blastp search and focused phylogenies:

        1. Use BLASTP to search for the --query proteins in --target.
        2. Cluster the hits using diamond + MCL and filter the clusters 
            - should include any query sequence and, optionally a species of interest (--soi) sequence
        3. For each cluster run the phylogeny
        4. Concatenate the annotations
    """, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_blastology.add_argument('--query', required=True, help='Query fasta file containing sequences to search for.')
    parser_blastology.add_argument('--target', required=True, help='Target fasta file - usually a collection of proteomes')
    parser_blastology.add_argument('-c','--ncpu', required=True, help='Number of CPU cores to use.')

    parser_blastology.add_argument('-p','--prefix', default = "search", required = False, help='Prefix to use for output files, e.g. [PREFIX].HG0.treefile. Default: search')
    parser_blastology.add_argument('-o','--output_dir', default = "results", help='Output directory name.Default [results]')
    parser_blastology.add_argument('--outputfile', required=True, default = None, help='Name of output file with annotations')
    parser_blastology.add_argument('-s','--soi', default = "", required = False, help='Prefix of the species of interest - e.g. "Mlei"')
    parser_blastology.add_argument('--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_blastology.add_argument('--phymethod', required=False, default = "iqtree3",  help='Phylogeny method: fasttree, iqtree2, iqtree3. Default: iqtree3')
    parser_blastology.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_blastology.add_argument('--force', required=False, default = False, action = 'store_true', help='Use this to rerun intermediate files (e.g. alignment)')
    parser_blastology.add_argument('-t','--temp_dir', required=False, default = 'tmp/', help='Temporary directory name. Default: tmp/')
    parser_blastology.add_argument('--cluster_prefix', required=False, default = 'HG', help='Prefix to use with sequence clusters. Default: "HG"')
    parser_blastology.add_argument('--min_perc', required=False, default = 30, help='Minimum sequence percentage identity for BLASTP hit filtering. Default [30]')
    parser_blastology.add_argument('--evalue', required=False, default = "1e-5", help='BLAST E-value threshold. Default "1e-5"')
    parser_blastology.add_argument('--mcl_inflation', required=False, default = "1.1", help='MCL inflation parameter. Default "1.1"')

    args = parser.parse_args()

    verbose = True # for debugging 

    if args.command is None:
        parser.print_help()
        sys.exit(1)  # Exit with a non-zero code to indicate an error


    if args.command == 'hmmsearch':
        logging.info("Command: Search")
        if args.hmm_dir:
            logging.info(f'HMM directory specified: {args.hmm_dir}')
            hmm_dir = args.hmm_dir
        else:
            logging.error(f'No HMM directory specified. Setting to hmms/')
            hmm_dir = "hmms"
        
        domain_expand = int(args.domain_expand) 
        hmmsearch(fasta_file = args.fasta, gene_family_info = args.gene_family_info, gene_family_name=args.gene_family_name, output_dir=args.output_dir, pfam_db=args.pfam_db,hmm_dir=hmm_dir, ncpu = int(args.ncpu),
                   domain_expand = domain_expand, verbose = verbose,do_clean = args.keep)

    elif args.command == 'cluster':
        logging.info("Command: Cluster")
        #clustering_method = 'diamond_mcl'
        from helper import subcluster as subcl
        
        
        infasta = args.fasta
        temp_dir = 'tmp/'
        cluster_log = 'tmp/cluster.log'
        ncpu = args.ncpu 
        max_N = int(args.maxn) # maximum number of sequences in the biggest cluster
        clustering_method = args.method
        do_recluster = True # use for iterative reclusering  
        max_iterations = 100


        if not args.out_file and not args.out_prefix:
            logging.error("Provide either --out_file or --out_prefix for clustering command!")
            sys.exit(1) 
        elif not args.out_file and args.out_prefix:
            #print('out_prefix provided')
            out_prefix = args.out_prefix
        elif not args.out_prefix and args.out_file:
            # if out_file provided, try to guess the extension of the output file 
            if not args.out_file.endswith('_cluster.tsv'):
                logging.error('--out_file should end in "_cluster.tsv"!')
                sys.exit(1)
            out_prefix = args.out_file.replace('_cluster.tsv','')
        else:
            logging.error("Provide either --out_file or --out_prefix for clustering command!")
            sys.exit(1)
        cluster_file = out_prefix + '_cluster.tsv'
        
        # should create {out_prefix}_cluster.tsv 
        cluster(fasta_file = args.fasta,out_prefix = out_prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = args.ncpu,method = clustering_method, cluster_prefix = args.cluster_prefix, mcl_inflation = args.inflation, logging = logging)
        
        N_seqs = functions.count_seqs(args.fasta)
        N,max_N_obs = subcl.top_n(cluster_file)
        print(subcl.top_n(cluster_file))
        logging.info(f'{cluster_file}: {N_seqs} sequences => {N} clusters. Maximum cluster size: {max_N_obs}. Max allowed: {max_N}')
        

        # 2 re-clustering modes - global and local: 
        # global - triggers the re-clustering of the wholo group 
        # local - reclusters only the groups exceeding the threshold 
        if max_N_obs > max_N:
            recluster_mode = 'local'
            #recluster_mode = 'global'
            if do_recluster:
                logging.info(f'Reclustering mode: {recluster_mode}')
                if recluster_mode == 'global':
                    logging.error(f'N sequences in the biggest cluster is more ({max_N_obs}) than allowed ({max_N})!')
                    logging.info('Trying to recluster with higher inflation...')
                    N,max_N_obs = subcl.top_n(cluster_file)
                    inflation = float(args.inflation)
                    iteration = 0
                    

                    while max_N_obs > max_N and iteration < max_iterations:
                        inflation += 0.1
                        inflation = round(inflation,1)
                        infasta = args.fasta
                        cluster(fasta_file = infasta,out_prefix = out_prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = args.ncpu,method = clustering_method, cluster_prefix = "HG", mcl_inflation = inflation, verbose = False)
                        N,max_N_obs = subcl.top_n(cluster_file)
                        logging.info(f'Iteration: {iteration}\nFASTA: {infasta}\nInflation: {inflation}; N max: {max_N_obs}')
                        iteration += 1
                    if iteration >= max_iterations:
                        logging.error(f'Max iterations {max_iterations} reached and the max N seqs is still more ({max_N_obs}) than allowed ({max_N})!')
                        sys.exit(1)
                    else:
                        logging.info(f'Iterative clustering finished: Iteration: {iteration}; Inflation: {inflation}; N max: {max_N_obs}') 
                
                elif recluster_mode == 'local':
                    # Recluster the HGs exceeding the threshold while keeping the others 
                    clusters = subcl.cluster_dict(cluster_file)
                    counts = subcl.cluster_counts(cluster_file)                    
                    too_big = [k for k,v in counts.items() if v >= max_N]
                    logging.info(f'{cluster_file}: {len(too_big)} clusters exceeding max_N {max_N}: {",".join(too_big)}')

                    cluster_files = []
                    cluster_status = []
                    for hg_id in too_big:
                        ids = clusters[hg_id]
                        file,status = subcl.recluster_hg_local(args = args,hg_id = hg_id,sequence_ids = ids, fasta_file = args.fasta, 
                                                               temp_dir = temp_dir,ncpu = args.ncpu,max_N = args.maxn, inflation = args.inflation, inflation_step = float(args.inflation_step),
                                                               max_iterations = max_iterations,verbose = True)
                        cluster_files.append(file)
                        cluster_status.append(status)
                    
                    logging.info(f'{sum(cluster_status)} / {len(cluster_status)} clusters sucessfully reclustered.')
                    if sum(cluster_status) < len(cluster_status):
                        logging.error("ERROR: some clusters have failed to be sub-clustered!")
                        raise RuntimeError("Some clusters failed to be sub-clustered")
                        
                    
                    # Collect the subclusters  
                    clusters_subclust = {}

                    for sub_file in cluster_files:
                        d = subcl.cluster_dict(sub_file)
                        if not d:
                            logging.warning(f"{sub_file} produced no clusters.")
                            continue
                        clusters_subclust.update(d)

                    # ---- Load original clusters ----
                    clusters_init = subcl.cluster_dict(cluster_file)

                    # Keep clusters that were NOT reclustered
                    clusters_result = {
                        k: v for k, v in clusters_init.items()
                        if k not in too_big
                    }

                    # Add reclustered clusters
                    clusters_result.update(clusters_subclust)

                    logging.info(f"Total clusters after merging: {len(clusters_result)}")

                    sorted_clusters = sorted(
                        clusters_result.items(),
                        key=lambda x: len(x[1]),
                        reverse=True
                    )

                   
                    sorted_renamed = {}
                    for i, (old_id, genes) in enumerate(sorted_clusters, start=1):
                        sorted_renamed[f"HG{i}"] = genes

                    # ---- Backup original file ----
                    import shutil
                    backup_file = cluster_file + ".tmp"
                    shutil.copy(cluster_file, backup_file)
                    logging.info(f"Previous clustering saved as: {backup_file}")

                    # ---- Write final clustering ----
                    with open(cluster_file, "w") as f:
                        for hg, genes in sorted_renamed.items():
                            for gene in genes:
                                f.write(f"{hg}\t{gene}\n")

                    logging.info(f"Clusters after re-clustering written to: {cluster_file}")

                    # ---- Stats ----
                    N_seqs = functions.count_seqs(args.fasta)
                    N, max_N_obs = subcl.top_n(cluster_file)

                    logging.info(
                        f"{cluster_file}: {N_seqs} sequences => {N} clusters. "
                        f"Maximum cluster size: {max_N_obs}. Max allowed: {max_N}"
                    )
                    # CAVE: these results should be appended!
                else:
                    logging.error(f'Unknown clustering mode: {recluster_mode}')
        ################################################################################
        # Create homology group fastas: 
        from pathlib import Path
        from collections import defaultdict
        from Bio import SeqIO

        def split_clusters_to_fastas(cluster_tsv, input_fasta, out_prefix):
            """
            Create per-cluster FASTA files from cluster TSV output.

            Parameters
            ----------
            cluster_tsv : str or Path
                TSV file with: clusterID <tab> sequenceID
            domains_fasta : str or Path
                FASTA file containing all sequences
            out_prefix : str
                Prefix for output FASTA files
            """

            cluster_tsv = Path(cluster_tsv)
            domains_fasta = Path(input_fasta)

            # Read cluster assignments
            clusters = defaultdict(list)

            with cluster_tsv.open() as f:
                for line in f:
                    cluster_id, seq_id = line.strip().split("\t")
                    clusters[cluster_id].append(seq_id)

            if not clusters:
                print("No clusters found.")
                return

            # Load all sequences into memory (fast and simple)
            seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

            # Write one FASTA per cluster
            for cluster_id, seq_ids in clusters.items():
                output_file = Path(f"{out_prefix}.{cluster_id}.fasta")

                with output_file.open("w") as out_handle:
                    for sid in seq_ids:
                        if sid in seq_dict:
                            SeqIO.write(seq_dict[sid], out_handle, "fasta")

                print(f"Created {output_file}")
        split_clusters_to_fastas(cluster_tsv = cluster_file,input_fasta = args.fasta,out_prefix = out_prefix )




    elif args.command == 'align':
        logging.info("Command: Align")
        # --notrim - triggers skipping alignment trimming 
        functions.align_and_trim(input_file = args.fasta, output_file = args.outfile, ncpu = args.ncpu, mafft_opt = args.mafft, notrim = args.notrim)

    elif args.command == 'trim':
        logging.info("Command: trim")
        # trim(input_file,output_file,mode = "kpic-gappy", g = "0.7",logfile = '/dev/null 2>&1',verbose = True)
        if args.method == 'clipkit':
            logfile = args.logfile
            if not logfile:
                logfile = '/dev/null 2>&1'
            functions.clipkit_trim(input_file = args.fasta, output_file = args.outfile, mode = "kpic-gappy", g = "0.7",logfile = logfile, verbose = True)
        else:
            logging.error(f'trim: unknown alignment trimming method: {args.method}')
            sys.exit(1)
        logging.info(f'trim: Done. Output file: {args.outfile}. Log: {args.logfile}')

    elif args.command == 'phylogeny':
        logging.info("Command: Phylogeny")
        # Properly decide on the output prefix 
        method = args.method 
        outprefix = args.outprefix
        outfile = args.outfile 
        if not outfile and not outprefix:
            logging.error("Please, provide either --outprefix or --outfile option!")
            sys.exit(1)
        else:
            if outfile:
                outprefix = os.path.splitext(os.path.abspath(outfile))[0]
                logging.info(f'Phylogeny: Output file: {outfile} => Output prefix: {outprefix}')
            elif outprefix:
                outfile = outprefix + '.treefile'
                logging.info(f'Phylogeny: Output prefix: {outprefix} => Output prefix: {outfile}')

        # the phylogeny function should receive either the output file or the output prefix value
        phylogeny(fasta_file = args.fasta, output_file = outfile, output_prefix = outprefix,ntmax = args.ncpu, method = method,iqtree2_model = args.iqtree2_model)

    elif args.command == 'generax':
        logging.info("Command: GeneRax")
        print(args)

        from pathlib import Path
        script = Path(__file__).resolve().parent / "helper" / "generax.py"

        cmd = (
            f"python {script} "
            f"--name {args.name} "
            f"--alignment {args.alignment} "
            f"--gene_tree {args.gene_tree} "
            f"--species_tree {args.species_tree} "
            f"--output_dir {args.output_dir} "
            f"--subs_model {args.subs_model} "
            f"--max_spr {args.max_spr} "
            f"--cpus {args.cpus} "
            f"--logfile {args.logfile} "
            f"--outfile {args.outfile}"
        )

        logging.info(cmd)

        result = subprocess.run(cmd, shell=True)

        if args.logfile:
            logging.info(f"Log: {args.logfile}")

        # 🔥 Propagate exit code
        sys.exit(result.returncode)

    elif args.command == 'possvm':
        min_support_transfer = float(args.possvm_minsupport)
        #if not os.path.exists('submodules/possvm-orthology/possvm.py'):
        #    logging.error("Can't find submodules/possvm-orthology/possvm.py! Exiting ...")
        possvm(
            treefile  = args.treefile,   
            reference_names = args.refnames,
            ogprefix = args.ogprefix,
            refsps = args.refsps,
            min_support_transfer = min_support_transfer,
            logfile = args.logfile,
            itermidroot = int(args.itermidroot),
            sos = args.sos)

    elif args.command == 'easy-phylo':

        logging.info('Easy-phylo')
        inputfasta = args.fasta
        force = args.force
        method = args.method
        mafft = args.mafft
        
        mafft_opt = functions.pick_mafft(mafft,inputfasta, max_n = 500, maxiterate = 1000, logging = logging)
        print(mafft_opt)
        if args.outdir is None:
            fname_aln = os.path.splitext(inputfasta)[0] + '.aln'
            tree_prefix = os.path.splitext(inputfasta)[0] + '.tree'
        else:
            outdir = args.outdir
            if not os.path.exists(outdir):
                logging.info(f'Created output directory {outdir}/')
                os.makedirs(outdir)
            else:
                logging.info(f'Warning: Output directory {outdir}/ exists')
            basename = os.path.splitext(os.path.basename(inputfasta))[0]
            fname_aln = os.path.join(outdir, basename + '.aln')
            tree_prefix = os.path.join(outdir, basename + '.tree')

        print(tree_prefix)
        fname_tree = tree_prefix + ".treefile"
        
        log_aln = f'{os.path.dirname(fname_tree)}/alignment.log'
        log_phy = f'{os.path.dirname(fname_tree)}/phylogeny.log'
        log_possvm = f'{os.path.dirname(fname_tree)}/possvm.log'

        logging.info(f'Easy-phylo: {fname_aln} => {fname_tree}. Log: {log_phy}')
        
        if os.path.isfile(fname_aln) and not force:
            print(f'Found alignment file: {fname_aln}! Skipping alignment')
        else:
            align_and_trim(input_file = args.fasta, output_file = fname_aln, ncpu = args.ncpu, mafft_opt = mafft_opt,logfile = log_aln, notrim = False)
        if os.path.isfile(fname_tree) and not force:
            print(f'Found phylogeny file: {fname_tree}! Skipping alignment')
        else:
            phylogeny(fasta_file = fname_aln, output_file = fname_tree,output_prefix = tree_prefix,ntmax = args.ncpu, method = method, logfile = log_phy)
        min_support_transfer = float(args.easyphylo_minsupport)
        
        if args.refnames:
            possvm(
                    treefile = fname_tree,
                    reference_names = args.refnames,
                    ogprefix = args.ogprefix,
                    min_support_transfer = min_support_transfer, 
                    logfile = log_possvm,
                    sos = float(args.sos),
                    outgroup = args.outgroup)
            fname_out = f'{fname_tree}.ortholog_groups.csv'
        else: 
            fname_out = fname_tree
        # Check the output:
        if not os.path.exists(fname_out):
            raise FileNotFoundError(f'Expected output not found: {fname_out}')
        else:
            logging.info(f'Easy-phylo done: {fname_out}')

    elif args.command == 'phylo-search':
        raise(NotImplementedError("Use blastology command!"))
    
    elif args.command == 'blastology':
        logging.info('BLASTology')
        from helper import blastology 
        blastology.blastology_run(args,logging, verbose = True)
