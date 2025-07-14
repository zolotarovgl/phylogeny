# Function to run GeneRax for specified inputs 

import ete3   
from ete3 import Tree
import subprocess
import os
import logging
import sys
import argparse

def list_tips(treefile):
	t = Tree(treefile, format=1)
	return [leaf.name for leaf in t.iter_leaves()]

def prune_tree_to_tips(treefile, tips, outfile, strict = False):

	# if tips are stored in a file, load it 
	if isinstance(tips, str) and os.path.isfile(tips):

		tipsfile = tips
		with open(tipsfile) as f:
			tips = [line.strip() for line in f if line.strip()]

	t = Tree(treefile, format=1)
	missing = [tip for tip in tips if not t.search_nodes(name=tip)]
	if missing:
		if strict: 
			raise ValueError(f"Missing tips in tree: {', '.join(missing)}")
		else:
			print(f"Warning: missing tips in tree: {', '.join(missing)}")
			tips = [tip for tip in tips if t.search_nodes(name = tip)]
	t.prune(tips, preserve_branch_length=True)
	t.write(outfile=outfile)

# function to determine species prefixes from fasta 
def get_species_from_fasta(input_fasta, sep="_"):
	cmd = f"grep '>' {input_fasta} | cut -f1 -d'{sep}'"
	output = subprocess.check_output(cmd, shell=True, text=True)
	output = list(set(output.strip().splitlines()))
	output = [x.replace('>','') for x in output]
	return output

def extract_model(iqtree_file):
	with open(iqtree_file) as f:
		for line in f:
			if "Model of substitution:" in line:
				return line.strip().split()[-1]
	raise ValueError("Substitution model not found in IQ-TREE file.")

def create_generax_config(name,alignment_file,tree_file, output_file, subst_model =  None, iqtree_file = None):

	def extract_model(iqtree_file):
		with open(iqtree_file) as f:
			for line in f:
				if "Model of substitution:" in line:
					return line.strip().split()[-1]

	# if the substitution model has 
	if not subst_model and not iqtree_file:
		print(f'provide a substitution model or and iqtree_file')
		sys.exit(1)

	if not subst_model and iqtree_file:
		subst_model = extract_model(iqtree_file)

	with open(output_file, 'w') as f:
		f.write("[FAMILIES]\n")
		f.write(f"- {name}\n")
		f.write(f"alignment = {alignment_file}\n")
		f.write(f"starting_gene_tree = {tree_file}\n")
		f.write(f"subst_model = {subst_model}\n")

def check_binary(program,logging):
    import shutil
    if shutil.which(program) is None:
        logging.error(f'{program} not found in PATH. Please install or load the binary.')
        sys.exit(1)
    else:
        logging.info(f'Found {program}')


# Function to launch generax
def run_generax(config_file,species_tree,rec_model = 'UndatedDL',max_spr = 7,strategy = 'SPR', per_family_rates = True, ncpu = 1):

	check_binary('generax',logging = logging)
	check_binary('mpirun',logging = logging)

	if per_family_rates:
		per_family_rates = "--per-family-rates"
	else:
		per_family_rates = ""
	#generax -s species_tree.newick -f Tubulin.fam --per-family-rates -r UndatedDL -p results_generax --max-spr-radius 1 --strategy SPR
	cmd = f"mpirun --use-hwthread-cpus --oversubscribe -n {ncpu} generax -s {species_tree} -f {config_file} {per_family_rates} -r {rec_model} -p {output_dir} --max-spr-radius {max_spr} --strategy {strategy} > /dev/null 2>&1"
	logging.info(cmd)
	subprocess.run(cmd, shell=True, check=True)

####################

#python main.py align -f $INFASTA -o test/Tubulin.aln -c 10
#python main.py easy-phylo -f $INFASTA -c 12 --method iqtree2
#python main.py phylogeny -f test/Tubulin.aln -c 12 --outfile test/Tubulin.tree --method fasttree
#python helper/create_generax_config.py --name TEST --alignment test/Tubulin.aln --tree test/Tubulin.tree --output Tubulin.fam --model LG+G
#generax -s species_tree.newick -f Tubulin.fam --per-family-rates -r UndatedDL -p results_generax --max-spr-radius 1 --strategy SPR 


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Run GeneRax")
	parser.add_argument('--name', required=True, help='Family name')
	parser.add_argument('--alignment', required=True, help='Path to alignment file')
	parser.add_argument('--gene_tree', required=True, help='Path to gene tree file')
	parser.add_argument('--species_tree', required=True, help='Path to species tree file')
	parser.add_argument('--output_dir', required=True, help='Path to output directory')
	parser.add_argument('--subs_model', help='Substitution model (e.g. LG+G)')
	parser.add_argument('--iqtree_file', help='Optional IQ-TREE .iqtree file to extract model')
	parser.add_argument('--per-family-rates', required=False, default = True, action = 'store_true', help='Whether to use per family rates')
	parser.add_argument('--max-spr', required=False, default = int(5), help='Maximum SPR radius')
	parser.add_argument('-c','--cpus', required=False, default = int(1), help='Number of CPU cores')
	parser.add_argument('-o','--outfile', required=False, default = None, help='Name of the output tree file')

	args = parser.parse_args()

	logging.basicConfig(
		level=logging.INFO,
		format='%(asctime)s - %(levelname)s - %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S'
	)

	# Check files
	if not os.path.isfile(args.alignment):
		print(f"Error: alignment file doesn't exist! {args.alignment}")
		sys.exit(1)

	if not os.path.isfile(args.gene_tree):
		print(f"Error: gene tree file doesn't exist! {args.gene_tree}")
		sys.exit(1)

	if not os.path.isfile(args.species_tree):
		print(f"Error: gene tree file doesn't exist! {args.species_tree}")
		sys.exit(1)

	if not args.subs_model and not args.iqtree:
		print("Error: Provide either --subs_model or --iqtree_file")
		sys.exit(1)

	if not args.subs_model:
		if args.iqtree_file and not os.path.isfile(args.iqtree_file):
			print(f"Error: IQ-TREE file doesn't exist! {args.iqtree_file}")
			sys.exit(1)

	subs_model = args.subs_model or extract_model(args.iqtree_file)
	output_dir = args.output_dir
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	
	config_file = output_dir + '/' + 'generax.config'
	if not args.outfile:
		output_file = output_dir + '/' + 'gene_tree.newick'
	else:
		output_file = args.outfile

	create_generax_config(
		name=args.name,
		alignment_file=args.alignment,
		tree_file=args.gene_tree,
		output_file=config_file,
		subst_model=subs_model
	)

	logging.info(f'Created: {config_file}')

	run_generax(config_file = config_file, species_tree = args.species_tree, rec_model = 'UndatedDL', max_spr = args.max_spr, per_family_rates = args.per_family_rates, ncpu = args.cpus)
	# copy the results 
	result = f'{output_dir}/results/{args.name}/geneTree.newick'
	if os.path.isfile(result):
		cmd = f'cp {result} {output_file}'
		#logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)
		logging.info(f'GeneRax done: {output_file}')
	else:
		logging.error(f'ERROR: GeneRax output file is missing: {result}')
		sys.exit(1)
	
