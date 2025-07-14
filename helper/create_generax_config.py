import argparse
import os
import sys

def extract_model(iqtree_file):
	with open(iqtree_file) as f:
		for line in f:
			if "Model of substitution:" in line:
				return line.strip().split()[-1]
	raise ValueError("Substitution model not found in IQ-TREE file.")

def create_generax_config(name, alignment_file, tree_file, output_file, subst_model):
	with open(output_file, 'w') as f:
		f.write("[FAMILIES]\n")
		f.write(f"- {name}\n")
		f.write(f"alignment = {alignment_file}\n")
		f.write(f"starting_gene_tree = {tree_file}\n")
		f.write(f"subst_model = {subst_model}\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Generate GeneRax family config file.")
	parser.add_argument('--name', required=True, help='Family name')
	parser.add_argument('--alignment', required=True, help='Path to alignment file')
	parser.add_argument('--tree', required=True, help='Path to gene tree file')
	parser.add_argument('--output', required=True, help='Path to output config file')
	parser.add_argument('--model', help='Substitution model (e.g. LG+G)')
	parser.add_argument('--iqtree', help='Optional IQ-TREE .iqtree file to extract model')

	args = parser.parse_args()

	# Check files
	if not os.path.isfile(args.alignment):
		print(f"Error: alignment file doesn't exist! {args.alignment}")
		sys.exit(1)

	if not os.path.isfile(args.tree):
		print(f"Error: gene tree file doesn't exist! {args.tree}")
		sys.exit(1)

	if args.iqtree and not os.path.isfile(args.iqtree):
		print(f"Error: IQ-TREE file doesn't exist! {args.iqtree}")
		sys.exit(1)

	if not args.model and not args.iqtree:
		print("Error: Provide either --model or --iqtree")
		sys.exit(1)

	model = args.model or extract_model(args.iqtree)

	create_generax_config(
		name=args.name,
		alignment_file=args.alignment,
		tree_file=args.tree,
		output_file=args.output,
		subst_model=model
	)

	print(f'Created: {args.output}')
