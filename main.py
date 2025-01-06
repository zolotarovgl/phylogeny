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

    # PHYLO-SEARCH
    parser_phylosearch = subparsers.add_parser('phylo-search',
            help = 'Search query sequences in proteomes and annotate using phylogenies',
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
    parser_phylosearch.add_argument('-o','--output_dir', default = "results", help='Output directory name.')
    parser_phylosearch.add_argument('-s','--soi', default = "", required = False, help='Prefix of the species of interest - e.g. "Mlei"')
    parser_phylosearch.add_argument('-m', '--mafft', required=False, default ="--auto", help='MAFFT: Mafft alignment options. Default  --auto')
    parser_phylosearch.add_argument('-r','--refnames', default = None, help='POSSVM: Reference gene names: gene \t name')
    parser_phylosearch.add_argument('--force', required=False, default = False, action = 'store_true', help='Use this to rerun intermediate files (e.g. alignment)')
    parser_phylosearch.add_argument('-t','--temp_dir', required=False, default = 'tmp/', help='Temporary directory name. Default: tmp/')
    parser_phylosearch.add_argument('--cluster_prefix', required=False, default = 'HG', help='Prefix to use with sequence clusters. Default: "HG"')

    args = parser.parse_args()

    verbose = True # for debugging 

    if args.command is None:
        parser.print_help()
        sys.exit(1)  # Exit with a non-zero code to indicate an error

    config = load_config()

    if args.command == 'hmmsearch':
        logging.info("Command: Search")
        if args.hmm_dir:
            logging.info(f'HMM directory specified: {args.hmm_dir}')
            hmm_dir = args.hmm_dir
        else:
            logging.info(f'No HMM directory specified. Checking config.yaml ...')
            hmm_dir = config['hmm_dir']
            logging.info(f'HMM directory from config.yaml: {hmm_dir}')
        
        pfam_db = args.pfam_db
        domain_expand = int(args.domain_expand) 
        hmmsearch(args.fasta, args.gene_family_info, args.gene_family_name, args.output_dir, pfam_db,config = config, domain_expand = domain_expand, verbose = verbose)

    elif args.command == 'cluster':
        logging.info("Command: Cluster")
        clustering_method = 'diamond_mcl'
        
        #check('diamond')
        #check('mcl')
        
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
        #check('mafft')
        #check('clipkit')
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

    elif args.command == 'easy-phylo' or args.command == 'blastology':
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
        logging.info(f"Phylo-search\n Query: {args.query}\n Target: {args.target}\n Threads: {args.ncpu}\n Prefix: {args.prefix}\n Species of interest: {args.soi}\n Mafft: {args.mafft}\n Reference names: {args.refnames}")
        
        min_n = 3 # minimal number of sequences in the cluster
        cluster_prefix = args.prefix + '.' + args.cluster_prefix # a prefix to add to the cluster 
        output_directory = args.output_dir
        cluster_directory = os.path.join(output_directory,'clusters')
        
        query = args.query
        target = args.target
        temp_dir = args.temp_dir
        prefix = args.prefix
        soi = args.soi
        force = args.force
        refnames_file = args.refnames
        ncpu = args.ncpu

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
        # Intermediate files check
        if not os.path.isfile(blastp_outfile):
            logging.info(f'Phylo-search:BLASTP:\n Query: {args.query}\n Target: {args.target}\n Threads: {args.ncpu}\nOutput{blastp_outfile}')
            blastp(query = args.query, target = args.target,db = temp_dir + "/target", outfile = blastp_outfile,ncpu = args.ncpu,logfile = blastp_log)
        else:
            logging.info(f'Found blastp output file {blastp_outfile}. Skipping')
   
        # Gather the results and prepare for clustering 
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
            clustering_method = 'diamond_mcl' 
            cluster(fasta_file = joint_fasta_fname,out_prefix = temp_dir + '/' + prefix,temp_dir = temp_dir,logfile = cluster_log,ncpu = ncpu,method = clustering_method)

        # Cluster filtering
        query_ids_file = os.path.join(temp_dir,'query.ids') 
        functions.get_fasta_names(fasta_file = query,out_file = query_ids_file)
        #functions.filter_clusters(cluster_file = cluster_file )
        with open(query_ids_file, "r") as file:
            query_ids = [line.strip() for line in file]

#############################################################
# Cluster filtering 
#############################################################
        import csv
        clusters = {}
        with open(cluster_file, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            for cluster_name, sequence_name in reader:
                clusters.setdefault(cluster_name, []).append(sequence_name)

        logging.info(f'Cluster filtering: {cluster_file}')
        logging.info(f'N query sequences: {len(query_ids)}')
        logging.info(f'N clusters: {len(clusters)}')
        
        # report
        clusters_query = [k for k,v in clusters.items() if any(elem in query_ids for elem in v)]
        clusters_small = [k for k,v in clusters.items() if len(v) < min_n]
        clusters_soi = [k for k,v in clusters.items() if any(soi in elem for elem in v)]
        clusters_query_small = [x for x in clusters_query if len(clusters[x]) < min_n]
        clusters_soi_small = [x for x in clusters_soi if len(clusters[x]) < min_n]
        logging.info(f'Clusters with query: {len(clusters_query)}')
        logging.info(f'Clusters with soi: {len(clusters_soi)}')
        logging.info(f'Clusters small: {len(clusters_small)}')
        logging.info(f'Clusters with query & small: {len(clusters_query_small)}')
        logging.info(f'Clusters with SOI & small: {len(clusters_soi_small)}')
        
        # Filter by query sequences
        clusters_filt = [k for k,v in clusters.items() if any(elem in query_ids for elem in v) and len(v) >= min_n]
        cluters_filt = [x for x in clusters_filt if any(soi in elem for elem in clusters[x])]
        clusters_filt_d = {k:v for k,v in clusters.items() if k in clusters_filt}
#CAVE:   # any small clusters with SOI sequences?
        #clusters_small_soi = [k for k,v in clusters.items() if any(soi in elem for elem in v) and any(elem in query_ids for elem in v) and len(v) < min_n]
        #if len(clusters_small_soi) > 0:
        #    logging.info(f'Warning: found {len(clusters_small_soi)} clusters with {soi} and query sequences but small size (<= {min_n} sequences)!\nConsider chaning the clustering method') 

        logging.info(f'{len(clusters_filt)}/{len(clusters)} clusters passing filtering.')
        clusters = clusters_filt_d
        
        require_soi = True
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
# Add cluster information / metadata 

# CAVE: replace with easy-phylo call? 
        # Cluster alignment and phylogeny 
        for cl_id in clusters_renamed.keys():
            cluster_fasta = os.path.join(cluster_directory,cl_id +  ".fasta")

            fname_aln = os.path.splitext(cluster_fasta)[0] + '.aln'
            tree_prefix = os.path.splitext(cluster_fasta)[0] + '.tree'
            fname_tree = tree_prefix + ".treefile"
            logfile = os.path.join(cluster_directory,cl_id + '.log')

            if os.path.isfile(fname_aln):
                print(f'Found alignment file: {fname_aln}! Skipping alignment')
            else:
                align_and_trim(input_file = cluster_fasta, output_file = fname_aln, ncpu = ncpu, mafft_opt = "--maxiterate 1000 --localpair --quiet --reorder", logfile = logfile)
            if os.path.isfile(fname_tree):
                print(f'Found phylogeny file: {fname_tree}! Skipping alignment')
            else:
                #functions.phylogeny_fasttree(fasta_file = fname_aln, output_file = fname_tree)
                phylogeny(fasta_file = fname_aln, output_prefix = tree_prefix,ntmax = ncpu)


            ogprefix = cl_id + "."

            possvm(treefile = fname_tree,reference_names = refnames_file,ogprefix = ogprefix,logfile = logfile)
        # concatenate annotation files 
        os.makedirs('results',exist_ok = True)
        anno_files = [os.path.join(cluster_directory, file) for file in os.listdir(cluster_directory) if 'ortholog_groups.csv' in file]

        anno_outfile = '%s/%s.annotation.tsv' % (output_directory,prefix)
        cmd = 'cat %s | grep -v orthogroup > %s' % (" ".join(anno_files),anno_outfile)
        #logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f'Done! Annotations stored in {anno_outfile}')
         
