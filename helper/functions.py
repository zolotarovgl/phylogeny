import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml
import ete3


#### Tool functions #####
def align(fasta_file, output_file, ncpu, mafft_opt, verbose = True):
	if verbose:
		n  = count_seqs(fasta_file,verbose)
		logging.info(f"Aligning {n} sequences with mafft")

	cmd = (f"mafft --reorder --quiet --thread {ncpu} {mafft_opt} {fasta_file} > {output_file}")
	if verbose:
		logging.info(cmd)
	subprocess.run(cmd, shell=True, check=True)
	if verbose:
		logging.info(f"Alignment done: {output_file}")

def clipkit_trim(input_file,output_file,mode = "kpic-gappy", g = "0.7",logfile = '/dev/null 2>&1',verbose = True):
	cmd = f"clipkit {input_file} -m {mode} -o {output_file} -g {g} > {logfile}"
	if verbose:
		logging.info(cmd)
	subprocess.run(cmd, shell=True, check=True)


def align_and_trim(input_file,output_file,ncpu = 1,mafft_opt = "",clipkit_mode = "kpic-gappy",clipkit_g = 0.7, clean = True, notrim = False, logfile = '/dev/null',verbose = True):
	if not os.path.exists(input_file):
		logging.error(f"{input_file} doesn't exist")
		sys.exit(1)
	tmpfile = input_file + '.tmp'
	if mafft_opt == 'fast':
		mafft_opt = ""
	do_trim = not notrim
	if do_trim:
		align(input_file,tmpfile,ncpu = ncpu,mafft_opt = mafft_opt,verbose = verbose)
		clipkit_trim(tmpfile,output_file,mode = clipkit_mode,g = clipkit_g,logfile = logfile,verbose = verbose)
		drop_gaponly(output_file, logging = logging)
		if clean:
			cmd = f"rm {tmpfile}"
			#logging.info(cmd)
			subprocess.run(cmd, shell=True, check=True)
	else:
		align(input_file,output_file,ncpu = ncpu,mafft_opt = mafft_opt,verbose = verbose)
	drop_gaponly(output_file, logging = logging)


def drop_gaponly(aln_file, logging = None):
	"""Remove sequences that are ALL gaps after alignment/trimming.

	IQ-TREE aborts outright on them -- "Sequence X contains only gaps or missing data /
	ERROR: Some sequences (see above) are problematic" -- which kills the whole family.
	Measured 2026-08-19 on Calponin: one Capitella fragment
	(Ctel_gnl_WGS_AMQN_CAPTEDRAFT_mRNA92005) survived MAFFT, was trimmed to nothing by
	clipkit -m kpic-gappy, and took the run down with it. Short fragments reach the
	alignment more often now that recruitment is no longer truncated by -max_target_seqs.

	phylohpc already does exactly this between its align and phylogeny rules
	(workflow/remove_gaponly.py, called from step2.smk); blastology had no equivalent.
	"""
	if not os.path.exists(aln_file):
		return
	kept, dropped = [], []
	name, seq = None, []
	def _flush():
		if name is None:
			return
		s = "".join(seq)
		(kept if s.replace("-", "").replace(".", "") != "" else dropped).append((name, s))
	with open(aln_file) as fh:
		for line in fh:
			line = line.rstrip("\n")
			if line.startswith(">"):
				_flush(); name, seq = line, []
			else:
				seq.append(line)
	_flush()
	if not dropped:
		return
	with open(aln_file, "w") as fh:
		for n, s in kept:
			fh.write(f"{n}\n{s}\n")
	msg = (f'ALIGNMENT: dropped {len(dropped)} gap-only sequence(s) after trimming '
		   f'({len(kept)} kept): {",".join(n[1:] for n, _ in dropped[:5])}'
		   + (" ..." if len(dropped) > 5 else ""))
	if logging is not None:
		logging.warning(msg)
	else:
		print(msg)


# Phylogeny wrappers
def get_node_support_range(treefile):
	# Load the tree
	#tree = Tree(treefile, format=1)  # format=1 assumes Newick with bootstrap values
	tree = ete3.PhyloTree("%s" % (treefile))
	# Collect support values from internal nodes

	support_values = [node.support for node in tree.traverse() if not node.is_leaf() and node.support is not None]
	min_support = min(support_values)
	max_support = max(support_values)
	return(min_support,max_support)


def phylogeny_get_prefix(output_file = None,output_prefix = None,verbose = False):
	# given on the output file name, get the file prefix - needed for the iqtree!
	if not output_file and not output_prefix or output_file and output_prefix:
		logging.error('Please, provide either --output_file or --output_prefix!')
		sys.exit(1)
	else:
		if output_file:
			logging.info(f'Output file provided. Guessing file prefix from the name:')
		if output_file.endswith('.treefile'):
			output_prefix = output_file.replace('.treefile','')
			if verbose:
				logging.info(f"IQTREE3: outputfile ends with .treefile => {output_prefix}")
		elif output_file.endswith('.tree'):
			output_prefix = output_file.replace('.tree','')
			if verbose:
				logging.info(f"IQTREE3: outputfile ends with .treefile => {output_prefix}")
			else:
				logging.info(f'Output prefix provided: {output_prefix}')

	return(output_prefix)

def check_binary(program,logging):
	import shutil
	if shutil.which(program) is None:
		logging.error(f'{program} not found in PATH. Please install or load the binary.')
		sys.exit(1)
	else:
		logging.info(f'Found {program}')



def phylogeny(fasta_file, output_file, output_prefix=None, ntmax=1, method='iqtree2',iqtree2_model = 'TEST', logfile='/dev/null'):

	logging.info(f'Phylogeny: {method}')
	check_binary(method, logging)

	if not output_prefix:
		logging.error('ERROR: phylogeny: specify the output prefix!')
		sys.exit(1)

	if method == 'iqtree2':
		phylogeny_iqtree2(
				fasta_file=fasta_file, output_file=output_file,
				output_prefix=output_prefix,
				model = iqtree2_model,
				cptime = 1000, nstop = 200,nm = 1000, bb = 1000, quiet = '',
				iqtree2 = 'iqtree2',
				ntmax=ntmax, logfile=logfile)

	elif method == 'iqtree3':
		phylogeny_iqtree3(fasta_file=fasta_file, output_file=output_file,
						  output_prefix=output_prefix, ntmax=ntmax, logfile=logfile)

	elif method == 'fasttree':
		phylogeny_fasttree(fasta_file, output_file, logfile=logfile)

	if not os.path.isfile(output_file):
		logging.error(f'Phylogeny has failed. Cannot find {output_file}! Aborting ...')
		sys.exit(1)

def phylogeny_iqtree2(fasta_file, output_file=None, output_prefix=None,
					  model='TEST', cptime=1000, nstop=200, nm=1000,
					  ntmax=15, bb=1000, quiet="",
					  iqtree2="iqtree2", logfile='/dev/null', verbose=True):

	logging.info(f"Phylogeny iqtree2: {fasta_file} {output_file}")

	if not output_prefix:
		logging.error('ERROR: specify output prefix!')
		sys.exit(1)

	cmd = f"{iqtree2} -pre {output_prefix} -s {fasta_file} -m {model} -mset LG,WAG,JTT -nt AUTO -ntmax {ntmax} -bb {bb} -nm {nm} -nstop {nstop} -cptime {cptime} {quiet} --redo > {logfile} 2>&1"
	logging.info(cmd)

	ret = subprocess.run(cmd, shell=True).returncode
	if ret != 0:
		sys.exit(ret)

	if output_file != f"{output_prefix}.treefile":
		mv_cmd = f"mv {output_prefix}.treefile {output_file}"
		ret = subprocess.run(mv_cmd, shell=True).returncode
		if ret != 0:
			sys.exit(ret)


def phylogeny_iqtree3(fasta_file, output_file = None, output_prefix = None, model = "MFP", cptime = 1000, nstop = 200, nm = 1000, ntmax = 15, bb = 1000, quiet = "",iqtree3 = "iqtree3",logfile = '/dev/null',verbose = True):
	# iqtree creates the files given a prefix {PREFIX}.treefile 
	# If the output file name provided and ends in .tree - use as a prefix 
	logging.info(f"Phylogeny: {fasta_file} {output_file}")

	cmd = f"{iqtree3} -s {fasta_file} -m {model} -T AUTO -ntmax {ntmax} -bb {bb} -nm {nm} -nstop {nstop} -cptime {cptime} {quiet} --prefix {output_prefix} --redo > {logfile} 2>&1"
	logging.info(cmd)
	subprocess.run(cmd, shell=True, check=True)
	#logging.info(f'IQTREE2: Created {output_prefix}.treefile')
	if output_file != f"{output_prefix}.treefile":
		cmd = f"mv {output_prefix}.treefile {output_file}"
		subprocess.run(cmd, shell=True, check=True)

def phylogeny_fasttree(fasta_file, output_file, logfile=''):

	logging.info(f"Phylogeny fasttree: {fasta_file} {output_file}")

	if logfile == '/dev/null':
		logfile = ''
	if logfile:
		logfile = f'-log {logfile}'

	cmd = f"fasttree {logfile} -quiet -gtr {fasta_file} | grep -v Ignoring > {output_file} 2>&1"
	logging.info(cmd)

	ret = subprocess.run(cmd, shell=True).returncode
	if ret != 0:
		sys.exit(ret)

	logging.info(f'Created {output_file}')


def possvm(treefile,
			output_prefix = None,
			reference_names = None, 
			ogprefix = "OG", 
			possvm = 'submodules/possvm-orthology/possvm.py',
			logfile = False,
			refsps = None,
			min_support_transfer = 50,
			itermidroot = 10,
			sos = 0,
			outgroup = "",
			phy = ""):

	if logfile:
		logging.info(f"Possvm: {treefile}\nLog: {logfile}")
	else:
		logfile = '/dev/null'
	# Adjust min_support according to the value range in provided tree 
	if not os.path.isfile(treefile):
		logging.error(f'ERROR: {treefile} does not exist!')
		sys.exit(1)
	nsr = get_node_support_range(treefile)
	if min_support_transfer:
		if min_support_transfer > nsr[1]:
			logging.info(f"Minimum node support ({min_support_transfer}) is bigger than the maximum observed support value ({nsr[1]}); Adjusting the threshold to {round(min_support_transfer/100,2)}") 
			min_support_transfer = min_support_transfer / 100
	# get the location of the possvm submodule 
	scriptdir = os.path.dirname(os.path.abspath(__file__))
	possvm = scriptdir + '/../' + possvm
	
	if reference_names:
		reference_names = f"-r {reference_names}"
	else:
		reference_names = ""
	
	if refsps:
		reference_species = f"-refsps {refsps}" 
	else:
		reference_species = ""
	
	if outgroup != '':
		outgroup  = f'--outgroup {outgroup}'
	
	if phy != '':
		phy =  f'--phy {phy}'
	# NOTE: -skipprint is intentionally NOT passed here -- POSSVM prints the annotated
	# phylogeny as a PDF by default, and the pipeline should always produce that visualization.
	cmd = f"python {possvm} --sos {sos} -ogprefix {ogprefix} -method lpa -itermidroot {itermidroot} -min_support_transfer {min_support_transfer}  -i {treefile} {reference_names} {reference_species} {outgroup}  {phy} >> {logfile} 2>&1"
	#print(cmd)
	logging.info(cmd)
	os.system(f'echo "{cmd}" > {logfile}')    
	subprocess.run(cmd, shell=True, check=True)

# Phylo-search functions 

def blastp(query,target,db,outfile,ncpu=1,evalue = "1e-5",min_perc = None,outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",logfile = '/dev/null',verbose = False,max_target_seqs = 5000):
	if not os.path.isfile(f'{db}.pdb'):
		cmd = f"makeblastdb -in {target} -dbtype prot -out {db} > {logfile}"
		if verbose:
			logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)
	else:
		logging.info(f'Found db files {db}. Skipping db building')
	# -max_target_seqs: BLAST's default is 500 subjects PER QUERY. Against a large target
	# (e.g. ~1e6 proteins) that cap binds by similarity rank, so the slots are consumed by
	# whichever lineage has the most near-identical paralogs and divergent homologs are
	# silently dropped. Measured 2026-08-19 on a 1,032,337-protein target, 9 myosin queries:
	# every query saturated at exactly 500, and only 3 of the 26 Aurelia myosins were
	# recruited; at 5000 the same search recovers 22 of 26. See
	# 2022_Mlei/muscle/docs/blastology_recruitment_bias.md and Issues.md I36.
	# NB the parameter is not a simple top-N filter -- it interacts with BLAST's heuristics.
	cmd = f'blastp -query {query} -out {outfile} -db {db} -evalue {evalue} -num_threads {ncpu} -max_target_seqs {max_target_seqs} -outfmt "{outfmt}" >> {logfile} 2>&1'
	logging.info(cmd)
	subprocess.run(cmd, shell=True, check=True)
	# Saturation warning: if any query returns exactly max_target_seqs distinct subjects the
	# cap is binding and the recruitment is truncated -- do NOT read absences from such a run.
	try:
		import collections
		per_q = collections.defaultdict(set)
		with open(outfile) as fh:
			for line in fh:
				f = line.split('\t')
				if len(f) > 1:
					per_q[f[0]].add(f[1])
		sat = [q for q, subs in per_q.items() if len(subs) >= int(max_target_seqs)]
		if sat:
			logging.warning(
				f'BLASTP: {len(sat)}/{len(per_q)} queries SATURATED -max_target_seqs '
				f'({max_target_seqs}); recruitment is truncated and absences are unreliable. '
				f'Raise it. Queries: {",".join(sorted(sat)[:5])}'
				+ (' ...' if len(sat) > 5 else ''))
	except Exception as e:
		logging.warning(f'BLASTP: could not check max_target_seqs saturation: {e}')
	if min_perc:
		logging.info(f'Minimum BLASTP hit percetage is set to {min_perc}. Filtering')
		cmd = f"cat {outfile} | awk '$3>={min_perc}' > {outfile}.filtered; mv {outfile}.filtered {outfile}"
		if verbose:
			logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)


def cluster(fasta_file,out_prefix,temp_dir,logfile = '/dev/null',method = 'mmseqs2',ncpu = 1,mcl_inflation = "1.1",cluster_prefix = "HG",verbose = True, logging = None,per_species_n = 6,graph_max_target_seqs = None,mmseqs_cov = 0.3):
	if method == 'mmseqs2':
		# ------------------------------------------------------------------------------
		# Connected-component clustering. --cluster-mode 1 (connected component) is the ONLY
		# mode that recovers divergent lineages, because it propagates TRANSITIVE homology:
		# A-B and B-C edges pull A and C together even when A-C similarity is below any
		# threshold. Modes 0 (set cover) and 2 (greedy incremental) ask instead "is this
		# similar enough to a representative?", which is a fundamentally weaker criterion.
		#
		# Measured 2026-08-19 on Tropomyosin (204 seqs, 13 queries), "does the Mnemiopsis
		# tropomyosin share a cluster with a query?":
		#     --cluster-mode 1  ->  13/13   <- this
		#     --cluster-mode 0  ->   4/13
		#     --cluster-mode 2  ->   0/13   <- what this branch used to run
		# Mode 2 fails at EVERY coverage setting tested (-c 0.8 / 0.5 / 0.3 / 0.2 / 0.0).
		#
		# Advantage over the diamond_mcl route: needs neither a neighbour cap nor an MCL
		# inflation value, so there is no threshold to justify per family.
		# ------------------------------------------------------------------------------
		cmd = f"mmseqs easy-cluster -s 7.5 -c {mmseqs_cov} --cov-mode 0 --cluster-mode 1 --min-seq-id 0.0 {fasta_file} {out_prefix} {temp_dir} --cluster-reassign >> {logfile} 2>&1"
		if verbose:
			 logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)
		# ADAPTER: mmseqs writes "<representative_seqid>\t<member_seqid>", but the rest of the
		# pipeline (blastology.filter_clusters) expects "<cluster_prefix><n>\t<seqid>" exactly
		# as the diamond_mcl branch emits. Without this the run still completes, but every
		# cluster is named after a representative SEQUENCE -- which then propagates into the
		# per-cluster fastas, tree prefixes and orthogroup names. Renumber by size, descending,
		# to match MCL's convention.
		import collections as _c
		_cl = _c.OrderedDict()
		with open(f'{out_prefix}_cluster.tsv') as _fh:
			for _line in _fh:
				_f = _line.rstrip('\n').split('\t')
				if len(_f) < 2:
					continue
				_cl.setdefault(_f[0], []).append(_f[1])
		_order = sorted(_cl.items(), key=lambda kv: -len(kv[1]))
		with open(f'{out_prefix}_cluster.tsv', 'w') as _fh:
			for _i, (_rep, _members) in enumerate(_order, start=1):
				for _m in _members:
					_fh.write(f'{cluster_prefix}{_i}\t{_m}\n')
		logging.info(f'CLUSTER (mmseqs cluster-mode 1): {len(_order)} clusters, '
					 f'largest {len(_order[0][1]) if _order else 0}; '
					 f'ids renumbered to {cluster_prefix}<n>')
	elif method == 'diamond_mcl':
		cmd = f"diamond makedb --in {fasta_file} -d {fasta_file} --quiet"
		if verbose: 
		   logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)

		# ------------------------------------------------------------------------------
		# All-vs-all graph for MCL.
		#
		# The cap here used to be a hard-coded GLOBAL 30 neighbours per sequence, which
		# fragments families along lineage lines: a sequence whose 30 best hits are all its
		# own lineage's paralogs gets NO edge to the rest of the family, MCL separates it,
		# and (in blastology) the query-free component is discarded before alignment.
		# Measured 2026-08-19 on Tropomyosin: the 30 neighbours of Mlei_v05_G010695 were all
		# ctenophore, and Mlei+Nvec+Aaur+Spis (148 seqs) were dropped.
		#
		# Fix, after Broccoli (Derelle et al. 2020, MBE 37:3389, "the N best hits per species
		# are reported, N=6 by default"): keep the N best hits PER SPECIES instead. One
		# lineage's paralogs then cannot monopolise the budget. Measured on the same family,
		# per-species N=3 succeeds with 4,817 edges where global top-30 FAILS with 4,170 --
		# i.e. it is not about graph size but about which edges are kept.
		#
		# The global diamond cap is derived, not guessed: to guarantee N per species you need
		# at least N * n_species slots.
		# Set per_species_n = None (or 0) to restore the old global-cap behaviour.
		# ------------------------------------------------------------------------------
		species = set()
		with open(fasta_file) as _fh:
			for _line in _fh:
				if _line.startswith('>'):
					species.add(_line[1:].split()[0].split('_')[0])
		n_species = len(species)
		if per_species_n and n_species > 1:
			diamond_max_target_seqs = int(graph_max_target_seqs or max(int(per_species_n) * n_species, 30))
		else:
			diamond_max_target_seqs = int(graph_max_target_seqs or 30)
			if per_species_n and n_species <= 1:
				logging.warning(
					f'CLUSTER: only {n_species} species prefix parsed from {fasta_file}; '
					f'per-species capping DISABLED (ids must look like <Species>_<rest>).')
				per_species_n = None
		logging.info(f'CLUSTER: {n_species} species, per_species_n={per_species_n}, '
					 f'diamond --max-target-seqs {diamond_max_target_seqs}')
		cmd = f"diamond blastp --more-sensitive --max-target-seqs {diamond_max_target_seqs} -d {fasta_file} -q  {fasta_file} -o {out_prefix}_diamond.csv --quiet --threads {ncpu}"
		if verbose:
			logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)
		cmd = f"awk '{{ print $1,$2,$12 }}' {out_prefix}_diamond.csv > {out_prefix}_diamond.abc"
		if verbose:
			 logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True) 
		# Broccoli-style per-species filter: keep the N best hits per (query, subject species).
		if per_species_n:
			import collections as _c
			_best = _c.defaultdict(list)
			_n_in = 0
			with open(f'{out_prefix}_diamond.abc') as _fh:
				for _line in _fh:
					_f = _line.split()
					if len(_f) < 3:
						continue
					_n_in += 1
					_best[(_f[0], _f[1].split('_')[0])].append((float(_f[2]), _line))
			_n_out = 0
			with open(f'{out_prefix}_diamond.abc', 'w') as _fh:
				for _k, _v in _best.items():
					_v.sort(key=lambda x: -x[0])
					for _sc, _line in _v[:int(per_species_n)]:
						_fh.write(_line); _n_out += 1
			logging.info(f'CLUSTER: per-species filter (N={per_species_n}): '
						 f'{_n_in} -> {_n_out} edges')
			if _n_out == 0:
				logging.error('CLUSTER: per-species filter produced an EMPTY graph -- aborting')
				sys.exit(1)
		cmd = f"mcl {out_prefix}_diamond.abc --abc -I {mcl_inflation} -o {out_prefix}_mcl.tsv 2> /dev/null"
		
		if verbose:
			logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)  
		cmd = f"""
		cat {out_prefix}_mcl.tsv | awk '{{ for (i = 1; i <= NF; i++) print "{cluster_prefix}"NR"\\t"$i }}' > {out_prefix}_cluster.tsv
		"""
		if verbose:
			 logging.info(cmd)
		subprocess.run(cmd, shell=True, check=True)  

	else:
		logging.info(f'Unknown clustering method {method}!')
		sys.exit(1)

def count_seqs(fasta_file, verbose=False):
	cmd = f"grep -c '>' {fasta_file}"
	if verbose > 1:
		logging.info(cmd)
	result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
	return int(result.stdout.strip())


def pick_mafft(mafft,fasta,max_n = 500, maxiterate = 1000, logging = None):
	# pick mafft options  based on the number of sequences
	logging.info(f'MAFFT: {mafft}')
	if mafft == 'auto':
			logging.info(f'MAFFT: picking MAFFT option based on the number of sequences.')
			n_seq = count_seqs(fasta, verbose=False)
			if n_seq <= max_n:
				mafft_mode = 'linsi'
			else:
				mafft_mode = 'fast'
			logging.info(f'MAFFT: # {n_seq} input sequences (max. {max_n}) => {mafft_mode}')
	else:
		mafft_mode = mafft

	if mafft_mode == 'fast':
		mafft_opt = ''
	elif mafft_mode == 'linsi':
		mafft_opt = f'--maxiterate {maxiterate} --localpair'
	else:
		logging.error(f"ERROR: unknown mafft mode: {mafft_mode}")
		sys.exit(1)
	return(mafft_opt)

def get_fasta_names(fasta_file,out_file,verbose = False):
	cmd = f"grep '>' {fasta_file} | sed 's/>//g' | sort | uniq > {out_file}"
	if verbose > 1:
		logging.info(cmd)
	subprocess.run(cmd, shell=True, check=True)

def check_tool(tool_name):
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

def check_dir(dirpath,force = False):
	if not os.path.exists(dirpath):
		os.makedirs(dirpath)
		logging.info(f'Directory created: {dirpath}')
	else:
		if not force:
			print(f'Directory {dirpath} exists, but --force is not set! Continuing ...')
			


from Bio import SeqIO

def retrieve_sequences(input_fasta, output_fasta, ids_to_keep):
	with open(output_fasta, "w") as outfile:
		for record in SeqIO.parse(input_fasta, "fasta"):
			if record.id in ids_to_keep:
				SeqIO.write(record, outfile, "fasta")

	logging.info(f"Created {output_fasta}")


def filter_clusters(cluster_file,fasta_file,query_ids_file,soi = None):
	# fasta_file should contain all the sequences 

	# Filtering :
	# The clusters should contain: the SOI and the query sequences
	logging.info('Filtering clusters:')
	print(cluster_file)
	print(fasta_file)
	print(query_ids_file)
	print(soi)

	print('filtering done')
