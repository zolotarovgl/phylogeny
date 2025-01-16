import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
from Bio import SeqIO

from helper.hmmsearch import hmmsearch
from helper.s02_cluster import cluster
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
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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

def run_generax(hg_id):
    logging.info(f"Running GeneRax for homology group: {hg_id}")
    raise(NotImplementedError())
    pass

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
    
    # Cluster
    parser_cluster = subparsers.add_parser('cluster', help='Run clustering')
    parser_cluster.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_cluster.add_argument('--out_file', required=False, default = None, help='Output file. CAVE: should be named PREFIX_cluster.tsv')
    parser_cluster.add_argument('--out_prefix', required=False,default = None, help='Output file prefix. Vis --out_file.')
    parser_cluster.add_argument('-c', '--ncpu', required=False, default = int(1), help='Number of CPU cores to use. Default: 1')
    parser_cluster.add_argument('-t', '--temp_dir', required=False, default = "tmp/", help='Temporary directory name. Default: tmp/')
    parser_cluster.add_argument('-i', '--inflation', default = float(1.1), help='Inflation parameter for MCL clustering')
    parser_cluster.add_argument('-m', '--maxn', default = int(1000), help='Maximum number of sequences in the cluster. Default: 1000')

    # Alignment
    parser_align = subparsers.add_parser('align', help='Run alignment')
    parser_align.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_align.add_argument('-o', '--outfile', required=True, help='Output file')
    parser_align.add_argument('-c', '--ncpu', required=False, default = int(1), help='Number of CPU cores to use')
    parser_align.add_argument('-m', '--mafft', required=False, default ="", help='Mafft alignment options. Default  "", you can set it to --maxiterate 1000 --genafpair')
    
    # Phylogeny
    parser_phylogeny = subparsers.add_parser('phylogeny', help='Run IQTREE2 for an alignment in --fasta')
    parser_phylogeny.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_phylogeny.add_argument('-o', '--outprefix', required=True, help='Output prefix for IQTREE2 files')
    parser_phylogeny.add_argument('-c', '--ncpu', required=True,  help='Number of CPU cores to use')
    parser_phylogeny.add_argument('--method', required=False, default = "fasttree",  help='Phylogeny method: fasttree, iqtree2. Default: fasttree')

    # GeneRax
    parser_generax = subparsers.add_parser('generax', help='Run GeneRax [NOT IMPLEMENTED]')
    parser_generax.add_argument('-t','--tree', required = True, help='Species tree')
    parser_generax.add_argument('-f','--fasta', required = True, help='Alignment file')
    parser_generax.add_argument('-o','--outdir', required = True, help='Output folder')

    # POSSVM
    parser_possvm = subparsers.add_parser('possvm', help='Run POSSVM')
    parser_possvm.add_argument('-t','--treefile', required = True, help='Input tree file')
    parser_possvm.add_argument('-r','--refnames', default = None, help='Reference gene names: gene \t name.')
    parser_possvm.add_argument('-o','--ogprefix', default = "OG", help='String. Prefix for ortholog clusters. Defaults to "OG".')
    parser_possvm.add_argument('-s','--refsps', default = None, help='POSSVM reference species')
    parser_possvm.add_argument('--min_support_transfer', default = "50",dest = "possvm_minsupport", help='POSSVM Minimum support for label transfer')
    
    # EASY-PHYLO
    parser_easyphylo = subparsers.add_parser('easy-phylo',help = 'Build a phylogeny from a single fasta')
    parser_easyphylo.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_easyphylo.add_argument('-c','--ncpu', required=True, help='Number of CPU cores to use')
    parser_easyphylo.add_argument('-m', '--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_easyphylo.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_easyphylo.add_argument('-o','--ogprefix', default = "OG", help='POSSVM: String. Prefix for ortholog clusters. Defaults to "OG".')
    parser_easyphylo.add_argument('--force', required=False, help='Use this to rerun intermediate files (e.g. alignment)')
    parser_easyphylo.add_argument('--method', default = "fasttree", help='Phylogeny method: fasttree, iqtree2. Default: fasttree')
    parser_easyphylo.add_argument('--min_support_transfer', default = "50", dest = "easyphylo_minsupport", help='POSSVM Minimum support for label transfer')


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
    parser_blastology.add_argument('-o','--output_dir', default = "results", help='Output directory name.')
    parser_blastology.add_argument('-s','--soi', default = "", required = False, help='Prefix of the species of interest - e.g. "Mlei"')
    parser_blastology.add_argument('--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_blastology.add_argument('--phymethod', required=False, default = "iqtree2",  help='Phylogeny method: fasttree, iqtree2. Default: fasttree')
    parser_blastology.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_blastology.add_argument('--force', required=False, default = True, action = 'store_true', help='Use this to rerun intermediate files (e.g. alignment)')
    parser_blastology.add_argument('-t','--temp_dir', required=False, default = 'tmp/', help='Temporary directory name. Default: tmp/')
    parser_blastology.add_argument('--cluster_prefix', required=False, default = 'HG', help='Prefix to use with sequence clusters. Default: "HG"')
    parser_blastology.add_argument('--min_perc', required=False, default = 30, help='Minimum sequence percentage identity for BLASTP hit filtering. Default [30]')
    parser_blastology.add_argument('--evalue', required=False, default = "1e-5", help='BLAST E-value threshold. Default "1e-5"')

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
                   domain_expand = domain_expand, verbose = verbose)

    elif args.command == 'cluster':
        logging.info("Command: Cluster")
        clustering_method = 'diamond_mcl'
        
        functions.check_tool('diamond')
        functions.check_tool('mcl')
        
        infasta = args.fasta
        temp_dir = 'tmp/'
        cluster_log = 'tmp/cluster.log'
        ncpu = args.ncpu 
        max_N = int(args.maxn) # maximum number of sequences in the biggest cluster
 

        if not args.out_file and not args.out_prefix:
            print("Provide either --out_file or --out_prefix for clustering command!")
            sys.exit(1) 
        elif not args.out_file and args.out_prefix:
            print('out_prefix provided')
            out_prefix = args.out_prefix
        elif not args.out_prefix and args.out_file:
            print('out_file provided')
            out_prefix = args.out_file.replace('_cluster.tsv','')
        else:
            print("Provide either --out_file or --out_prefix for clustering command!")
        # should create {out_prefix}_cluster.tsv 
        cluster(fasta_file = args.fasta,out_prefix = out_prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = args.ncpu,method = clustering_method, cluster_prefix = "HG", mcl_inflation = args.inflation)
        cluster_file = out_prefix + '_cluster.tsv'

        def top_n(file_path):
            counts = {}
            with open(file_path, 'r') as file:
                for line in file:
                    first_col = line.split('\t')[0]  # Get the first column
                    counts[first_col] = counts.get(first_col, 0) + 1
                # Sort by counts (descending) and get the most common value
            return max(counts.values())

        max_N_obs = top_n(cluster_file)
        logging.info(f'{cluster_file}: max observed number of sequences: {max_N_obs}')
        
       
        if max_N_obs > max_N:
            logging.error(f'N sequences in the biggest cluster is more ({max_N_obs}) than allowed ({max_N})!')
            logging.info('Trying to recluster with higher inflation...')
            max_N_obs = top_n(cluster_file)
            inflation = args.inflation
            iteration = 0
            max_iterations = 100

            while max_N_obs > max_N and iteration < max_iterations:
                inflation += 0.1
                inflation = round(inflation,1)
                cluster(fasta_file = args.fasta,out_prefix = out_prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = args.ncpu,method = clustering_method, cluster_prefix = "HG", mcl_inflation = inflation, verbose = False)
                max_N_obs = top_n(cluster_file)
                logging.info(f'Iteration: {iteration}; Inflation: {inflation}; N max: {max_N_obs}')
                iteration += 1
            if iteration >= max_iterations:
                logging.error(f'Max iterations {max_iterations} reached and the max N seqs is still more ({max_N_obs}) than allowed ({max_N})!')
                sys.exit(1)
            else:
                logging.info(f'Iterative clustering finished: Iteration: {iteration}; Inflation: {inflation}; N max: {max_N_obs}') 


    elif args.command == 'align':
        logging.info("Command: Align")
        functions.check_tool('mafft')
        functions.check_tool('clipkit')
        align_and_trim(input_file = args.fasta, output_file = args.outfile, ncpu = args.ncpu, mafft_opt = args.mafft)

    elif args.command == 'phylogeny':
        logging.info("Command: Phylogeny")
        method = args.method 
        phylogeny(fasta_file = args.fasta, output_prefix = args.outprefix,ntmax = args.ncpu, method = method)

    elif args.command == 'generax':
        logging.info("Command: GeneRax")
        raise(NotImplementedError())
        #run_generax()

    elif args.command == 'possvm':
        min_support_transfer = float(args.possvm_minsupport)
        #if not os.path.exists('submodules/possvm-orthology/possvm.py'):
        #    logging.error("Can't find submodules/possvm-orthology/possvm.py! Exiting ...")
        possvm(treefile  = args.treefile,reference_names = args.refnames,ogprefix = args.ogprefix, refsps = args.refsps, min_support_transfer = min_support_transfer)

    elif args.command == 'easy-phylo':
        logging.info('Easy-phylo')
        fname_aln = os.path.splitext(args.fasta)[0] + '.aln'
        tree_prefix = os.path.splitext(args.fasta)[0] + '.tree'
        fname_tree = tree_prefix + ".treefile"
        force = args.force
        method = args.method 

        if os.path.isfile(fname_aln) and not force:
            print(f'Found alignment file: {fname_aln}! Skipping alignment')
        else:
            align_and_trim(input_file = args.fasta, output_file = fname_aln, ncpu = args.ncpu, mafft_opt = "")
        if os.path.isfile(fname_tree) and not force:
            print(f'Found phylogeny file: {fname_tree}! Skipping alignment')
        else:
            phylogeny(fasta_file = fname_aln, output_prefix = tree_prefix,ntmax = args.ncpu, method = method)
        min_support_transfer = float(args.easyphylo_minsupport)
        possvm(treefile = fname_tree,reference_names = args.refnames,ogprefix = args.ogprefix,min_support_transfer = min_support_transfer)
        print('Easy-phylo done!')

    elif args.command == 'phylo-search':
        raise(NotImplementedError("Use blastology command!"))
    
    elif args.command == 'blastology':
        logging.info('BLASTology')
        from helper import blastology 
        blastology.blastology_run(args,logging, verbose = True)
