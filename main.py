import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
from helper.s01_search import search
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


def check(tool_name):
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



#############################################################################################


def run_easyphylo(fasta_file,ncpu):
    logging.info(f"Easy-phylo: {fasta_file}")
    quit()









#############################################################################################

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

    # Search 
    parser_search = subparsers.add_parser('search', help='Search for a family using HMMER')
    parser_search.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_search.add_argument('-g', '--gene_family_info', required=True, help='Path to the gene family info file specifying HMMs and parameters')
    parser_search.add_argument('gene_family_name', help='Name of the gene family to search')

    # Cluster
    parser_cluster = subparsers.add_parser('cluster', help='Run clustering')
    parser_cluster.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_cluster.add_argument('-o', '--outfile', required=True, help='Output file')
    parser_cluster.add_argument('-c', '--ncpu', required=False, default = int(1), help='Number of CPU cores to use')
    parser_cluster.add_argument('-i', '--inflation', default = float(1.1), help='Inflation parameter for MCL clustering')

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

    # GeneRax
    parser_generax = subparsers.add_parser('generax', help='Run GeneRax')
    parser_generax.add_argument('-t','--tree', required = True, help='Species tree')
    parser_generax.add_argument('-f','--fasta', required = True, help='Alignment file')
    parser_generax.add_argument('-o','--outdir', required = True, help='Output folder')

    # POSSVM
    parser_possvm = subparsers.add_parser('possvm', help='Run POSSVM')
    parser_possvm.add_argument('-t','--treefile', required = True, help='Input tree file')
    parser_possvm.add_argument('-r','--refnames', default = None, help='Reference gene names: gene \t name.')
    parser_possvm.add_argument('-o','--ogprefix', default = "OG", help='String. Prefix for ortholog clusters. Defaults to "OG".')
    

    # EASY-PHYLO
    parser_easyphylo = subparsers.add_parser('easy-phylo',help = 'Build a phylogeny from a single fasta')
    parser_easyphylo.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_easyphylo.add_argument('-c','--ncpu', required=True, help='Number of CPU cores to use')
    parser_easyphylo.add_argument('-m', '--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_easyphylo.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_easyphylo.add_argument('-o','--ogprefix', default = "OG", help='POSSVM: String. Prefix for ortholog clusters. Defaults to "OG".')
    parser_easyphylo.add_argument('--force', required=False, help='Use this to rerun intermediate files (e.g. alignment)')

    # PHYLO-SEARCH
    parser_phylosearch = subparsers.add_parser('phylo-search',
            help = 'Annotate using phylogeny',
            description="""
    Annotate the proteins based on blastp search and focused phylogenies:

        1. Use BLASTP to search for the --query proteins in --target.
        2. Cluster the hits using MMSEQS2 and filter the clusters 
            - should include any query sequence and, optionally a species of interest (--soi) sequence
        3. For each cluster:
            - align using MAFFT 
            - trim the aligment using ClipKit 
            - run phylogeny using IQTREE2
            - parse the phylogeny using POSSVM 
        4. Concatenate the annotations
    """, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_phylosearch.add_argument('--query', required=True, help='Query fasta file containing sequences to search for.')
    parser_phylosearch.add_argument('--target', required=True, help='Target fasta file - usually a collection of proteomes')
    parser_phylosearch.add_argument('-c','--ncpu', required=True, help='Number of CPU cores to use.')
    parser_phylosearch.add_argument('-p','--prefix', default = None, required = True, help='Prefix to use for all files and orthogroups.')
    parser_phylosearch.add_argument('-s','--soi', default = "", required = False, help='Prefix of the species of interest - e.g. "Mlei"')
    parser_phylosearch.add_argument('-m', '--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_phylosearch.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_phylosearch.add_argument('--force', required=False, help='Use this to rerun intermediate files (e.g. alignment)')
    parser_phylosearch.add_argument('--temp_dir', required=False, default = 'tmp/', help='Temporary directory name. Default: tmp/')
    

    args = parser.parse_args()

    verbose = True # for debugging 

    if args.command is None:
        parser.print_help()
        sys.exit(1)  # Exit with a non-zero code to indicate an error

    config = load_config()

    if args.command == 'search':
        logging.info("Command: Search")
        search(args.fasta, args.gene_family_info, args.gene_family_name, config, verbose)

    elif args.command == 'cluster':
        logging.info("Command: Cluster")
        check('diamond')
        check('mcl')
        cluster(fasta_file = args.fasta, output_file = args.outfile ,inflation = args.inflation, ncpu = args.ncpu)

    elif args.command == 'align':
        logging.info("Command: Align")
        #check('mafft')
        check('clipkit')
        align_and_trim(input_file = args.fasta, output_file = args.outfile, ncpu = args.ncpu, mafft_opt = args.mafft)

    elif args.command == 'phylogeny':
        logging.info("Command: Phylogeny")
        check("iqtree2")
        phylogeny(fasta_file = args.fasta, output_prefix = args.outprefix,ntmax = args.ncpu)

    elif args.command == 'generax':
        logging.info("Command: GeneRax")
        raise(NotImplementedError())
        #run_generax()

    elif args.command == 'possvm':
        if not os.path.exists('submodules/possvm-orthology/possvm.py'):
            logging.error("Can't find submodules/possvm-orthology/possvm.py! Exiting ...")
        possvm(treefile  = args.treefile,reference_names = args.refnames,ogprefix = args.ogprefix)

    elif args.command == 'easy-phylo':
        logging.info('Easy-phylo')
        fname_aln = os.path.splitext(args.fasta)[0] + '.aln'
        tree_prefix = os.path.splitext(args.fasta)[0] + '.tree'
        fname_tree = tree_prefix + ".treefile"
        force = args.force 

        if os.path.isfile(fname_aln) and not force:
            print(f'Found alignment file: {fname_aln}! Skipping alignment')
        else:
            align_and_trim(input_file = args.fasta, output_file = fname_aln, ncpu = args.ncpu, mafft_opt = "")
        if os.path.isfile(fname_tree) and not force:
            print(f'Found phylogeny file: {fname_tree}! Skipping alignment')
        else:
            phylogeny(fasta_file = fname_aln, output_prefix = tree_prefix,ntmax = args.ncpu)
        possvm(treefile = fname_tree,reference_names = args.refnames,ogprefix = args.ogprefix)
        print('Easy-phylo done!')

    elif args.command == 'phylo-search':
        logging.info(f"Phylo-search\n Query: {args.query}\n Target: {args.target}\n Threads: {args.ncpu}\n Prefix: {args.prefix}\n Species of interest: {args.soi}\n Mafft: {args.mafft}\n Reference names: {args.refnames}")
        
        
        query = args.query
        target = args.target
        temp_dir = args.temp_dir
        prefix = args.prefix
        soi = args.soi

        # check the temporary directory status:
        force = False
        if force:
            functions.check_tempdir(args.temp_dir)
        
        blastp_outfile = temp_dir + "/blastp_result.tsv"
        cluster_file = temp_dir + '/' + prefix + '_cluster.tsv'
        
        if not os.path.isfile(blastp_outfile):
            logging.info(f'Phylo-sarch:BLASTP:\n Query: {args.query}\n Target: {args.target}\n Threads: {args.ncpu}\nOutput{blastp_outfile}')
            blastp(query = args.query, target = args.target,db = args.temp_dir + "/target", outfile = blastp_outfile,ncpu = args.ncpu)
        else:
            logging.info(f'Found blatp output file {blastp_outfile}. Skipping')

        # Gather the results and prepare for clustering 
        joint_fasta_fname = temp_dir + "/" + prefix + ".fasta"
        joint_ids_fname = temp_dir + "/" + prefix + ".hits.ids"
        cmd = f'cat {query} {target} > {joint_fasta_fname}_tmp; samtools faidx {joint_fasta_fname}_tmp'
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"cat {blastp_outfile} | awk '{{print $1\"\\n\"$2}}' | sort | uniq > {joint_ids_fname}"
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd = f'xargs samtools faidx {joint_fasta_fname}_tmp < {joint_ids_fname} > {joint_fasta_fname};rm {joint_fasta_fname}_tmp'
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)

        # Clustering and cluster filtering 
        if os.path.isfile(cluster_file):
            logging.info(f'Found clustering file {cluster_file}. Skipping')
        else:
            cluster(fasta_file = joint_fasta_fname,out_prefix = temp_dir + '/' + prefix,temp_dir = temp_dir)

        # Cluster filtering
        min_n = 0
        query_ids_file = temp_dir + "/query.ids" 
        functions.get_fasta_names(fasta_file = query,out_file = query_ids_file)
        #functions.filter_clusters(cluster_file = cluster_file )
        with open(query_ids_file, "r") as file:
            query_ids = [line.strip() for line in file]
        
        import csv
        clusters = {}
        with open(cluster_file, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            for cluster_name, sequence_name in reader:
                clusters.setdefault(cluster_name, []).append(sequence_name)
        print(f'# query sequences: {len(query_ids)}')
        print(f'# clusters: {len(clusters)}')

        # Filter by query sequences
        clusters_filt = [k for k,v in clusters.items() if any(elem in query_ids for elem in v) and len(v) >= min_n]
        logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters with query and >= {min_n} sequences.')
        clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}
        clusters = clusters_filt_d
        
        if not soi == "":
            logging.info('Filtering by species of interest {soi}')
            clusters_filt = [k for k,v in clusters.items() if any(soi in elem  for elem in v)]
            logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters with SOI ({soi}) sequences.')
            clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}
            clusters = clusters_filt_d

        # rename each cluster and store the sequences 
        clusters_renamed = {}
        for i,(k,v) in enumerate(clusters.items()):
            clusters_renamed.update({"c"+str(i):v})
        print(clusters_renamed)

        # this file handling takes sooooooo much time .... 
        # finally, get their sequences 
        # an rename? why am I wasting my time like this?


        from pyfaidx import Fasta

        def filter_fasta(input_fasta, output_fasta, ids_to_keep):
            fasta = Fasta(input_fasta)
            with open(output_fasta, "w") as outfile:
                for seq_id in ids_to_keep:
                    if seq_id in fasta:
                        outfile.write(f">{seq_id}\n{fasta[seq_id][:]}\n")

        # Usage example
        cl_id = 'c0'
        fasta_file = joint_fasta_fname
        ids_to_keep = clusters_renamed[cl_id]  # Replace with your list of IDs
        output_fasta = "suka.c0.fasta"
        filter_fasta(fasta_file, output_fasta, ids_to_keep)