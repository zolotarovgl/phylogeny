import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
from helper.s01_search import search
from helper.s02_cluster import cluster

# Ensure PyYAML is installed
try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is not installed. Install it using 'pip install pyyaml'.")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config():
    config_file = "config.yaml"
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file {config_file} not found.")
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config


def check(tool_name):
    try:
        # Try to run the tool with a harmless argument like --help
        result = subprocess.run([tool_name, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # If return code is 0, it means the tool executed successfully
        if result.returncode == 0:
            print(f"{tool_name} is available and can be launched.")
        else:
            print(f"Error: {tool_name} is not functioning properly.")
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: {tool_name} is not available on your system.")
        sys.exit(1)

# Pipeline # 

def phylogeny(hg_id):
    logging.info(f"Running phylogeny for homology group: {hg_id}")
    pass

def run_generax(hg_id):
    logging.info(f"Running GeneRax for homology group: {hg_id}")
    pass

def run_possvm(hg_id):
    logging.info(f"Running POSSVM for homology group: {hg_id}")
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phylogeny tool")
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    parser_search = subparsers.add_parser('search', help='Search for a family using HMMER')
    parser_search.add_argument('-f', '--fasta', required=True, help='Path to the input fasta file')
    parser_search.add_argument('-g', '--gene_family_info', required=True, help='Path to the gene family info file specifying HMMs and parameters')
    parser_search.add_argument('gene_family_name', help='Name of the gene family to search')

    parser_cluster = subparsers.add_parser('cluster', help='Run clustering')
    parser_cluster.add_argument('-f', '--fasta', required=True, help='Path to the input fasta file')
    parser_cluster.add_argument('-o', '--outfile', required=True, help='Output file')
    parser_cluster.add_argument('-c', '--ncpu', required=False, default = int(1), help='Number of CPU cores to use')
    parser_cluster.add_argument('-i', '--inflation', default = float(1.1), help='Inflation parameter for MCL clustering')

    parser_phylogeny = subparsers.add_parser('phylogeny', help='Run phylogeny')
    parser_phylogeny.add_argument('hg_id', help='ID of the homology group')

    parser_generax = subparsers.add_parser('generax', help='Run GeneRax')
    parser_generax.add_argument('hg_id', help='ID of the homology group')

    parser_possvm = subparsers.add_parser('possvm', help='Run POSSVM')
    parser_possvm.add_argument('hg_id', help='ID of the homology group')

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

    elif args.command == 'phylogeny':
        logging.info("Command: Phylogeny")
        phylogeny(args.hg_id)

    elif args.command == 'generax':
        logging.info("Command: GeneRax")
        run_generax(args.hg_id)

    elif args.command == 'possvm':
        logging.info("Command: POSSVM")
        run_possvm(args.hg_id)

