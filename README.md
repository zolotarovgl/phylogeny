# Installation:
```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylogeny.git
```

To use GeneRax, install it from source (`https://github.com/BenoitMorel/GeneRax`) or from conda (`https://anaconda.org/bioconda/generax`).


# Commands   


```bash
usage: main.py [-h] {hmmsearch,pfamscan,cluster,align,phylogeny,generax,possvm,easy-phylo,blastology} ...

Python wrapper around some useful commands

positional arguments:
  {hmmsearch,pfamscan,cluster,align,phylogeny,generax,possvm,easy-phylo,blastology}
                        Sub-command help
    hmmsearch           Search for a family using HMMER
    pfamscan            Scan sequences against a pressed HMM database using hmmscan
    cluster             Run clustering
    align               Run alignment
    phylogeny           Run phylogeny 
    generax             Run GeneRax 
    possvm              Run POSSVM
    easy-phylo          Build a phylogeny from a single fasta
    blastology          Search query sequences in proteomes and annotate using phylogenies
```




# Command Examples  

### HMMSEARCH
Use the gene family info in `data/genefam.tsv` to search for proteins and to group them into defined gene families.  

```bash
PFAM_DB=~/ant/xgraubove/data/pfam/Pfam-A.hmm #location of PFAM database for .hmm fetching 
python main.py hmmsearch -f data/sample.fasta -g data/genefam.tsv Insulin -o results --pfam_db $PFAM_DB --domain_expand 50 
```

- `--domain_expand` option controls the number of aminoacids added from left and right of the extracted domain range   


Main outputs: 
- `bet.Insulin.domains.fasta` - domain ranges sequences  
- `bet.Insulin.domains.csv` - domain ranges .bed file   
- `bet.Insulin.seqs.fasta` - full protein sequences  

### PFAMSCAN
Use `hmmscan` against a pressed HMM database such as `Pfam-A.hmm` without going through the family-specific `hmmsearch` workflow.

First, press the database once:

```bash
hmmpress ~/Documents/db/Pfam/Pfam-A.hmm
```

Then run the scan:

```bash
python main.py pfamscan \
  -f tmp/query.fasta \
  --pfam_db ~/Documents/db/Pfam/Pfam-A.hmm \
  -o results/pfamscan \
  -c 8
```

Use `--outprefix` if you want to override the default output filename prefix derived from the input FASTA basename.
Use `--arch_sep` if you want a different separator in `*.pfamscan_archs.csv`; the default is a comma.

Main outputs:
- `query.pfamscan.domtblout` - raw `hmmscan --domtblout` result
- `query.pfamscan.out` - full `hmmscan` stdout/stderr capture
- `query.pfamscan.tsv` - parsed tab-separated hit table
- `query.pfamscan.domains.csv` - bed-like per-domain hit table: `sequence_id  start  end  pfam_name`
- `query.pfamscan_archs.csv` - one merged record per protein with comma-separated domains ordered left-to-right

For non-Pfam HMM databases that do not define gathering thresholds, disable `--cut_ga` and set a domain E-value explicitly:

```bash
python main.py pfamscan \
  -f tmp/query.fasta \
  --pfam_db tmp/toy.hmm \
  -o results/pfamscan \
  --no-cut_ga \
  --domE 1e-3
```

### Clustering   


An example of an iterative re-clustering for a case where the biggest cluster exceeds the maximum number of sequences (`-m`) staring from a 1.1 inflation & with an inflation step of 1 (`--inflation_step`):    


```bash
python main.py cluster -f data/cluster_test.fasta --out_prefix boo -c 4 -i 1.1 -m 50 --inflation_step 1
```

__TODOs:__ add inflation step, add a possibility of keeping a cluster if fails, add letter suffixes instead of re-numbering  



### BLASTOLOGY 
The tool will use BLASTP to search, then it will cluster the sequences with MCL&diamond and run the phylogeny for each cluster. 
Finally, it will parse the resulting trees and will output the best orthologs for the queries.  

Use the reference sequences for search. For instance, given the reference sequences stored in `data/BCL2.fasta`, search `data/sample.fasta` for homologs in species `Owefus` using 10 cores:  
```bash
python main.py blastology --query data/BCL2.fasta --refnames data/BCL2.names --target data/sample.fasta -c 5 --force --soi Owefus --outputfile Owenia_bcl2.tab --phymethod fasttree
```
Not specifying `--soi` (species of interest) will return the annotation file with all sequences in the orthogroups with queries.  



#### QUERY FORMATTING  
It's important that the query names are formatted as following:  
```
>Hsap_QUERY_P10415
>Mmus_QUERY_P10417
>Dmel_QUERY_Q7KM33
```
First prefix __MUST__ correspond to prefixes in Xavi's database to not to mess up the species reconciliation. It's also good to include the word `QUERY` in the name to be able to search for queries later on.  
The `.names` file should be map the names in the fasta to the protein names (these will be used for name transfer):   
```
Hsap_QUERY_P10415       hsapBCL2
Mmus_QUERY_P10417       mmusBCL2
Dmel_QUERY_Q7KM33       dmelBCL2
```
The names (2nd column) can be any, but try to keep them distinguishable!  

**TODO**: add species prefix check - strict and not. Request that all the reference species are present in the target file.


The results should look like the following:  
```bash
cat Owenia_bcl2.tab 
Owefus_OFUSG13935.2     search.OG0:dmelBCL2     0.997   dmelBCL2        0.84
Owefus_OFUSG16636.1     search.OG3:hsapBCL2/mmusBCL2    0.956   hsapBCL2/mmusBCL2       0.956/0.956
Owefus_OFUSG14207.1     search.OG3:hsapBCL2/mmusBCL2    0.956   hsapBCL2/mmusBCL2       0.956/0.956
```

### Example
Below is an example of how to fetch the sequences from the uniprot, rename them and run the blastology to get Owefus orthologs: 
```bash
python helper/fetch_uniprot.py data/PAX9.uniprot.txt  > test.fasta
# to the names file (from the fasta headers)
cat test.fasta | grep '>' | sed 's/>//g' | tr '.' '\t' | awk '{print $1"."$2"\t"$2}' > test.names
python main.py blastology --query test.fasta --refnames test.names --target data/sample.fasta -c 5 --force --soi Owefus --outputfile test.tab --phymethod fasttree --mafft ""
```




# Tests   


```bash 
# alignment and trimming 
INFASTA=data/arp23.fasta
python main.py align -f $INFASTA -o test.aln --notrim
python main.py trim -f test.aln -o test.alt --logfile test.log

# keep the pre-ClipKIT alignment as test.alt.untrimmed
python main.py align -f $INFASTA -o test.alt --keep

# easy-phylo
python main.py easy-phylo --method fasttree -f data/arp23.fasta -c 16 --outdir results
# blastology
python main.py blastology --query data/BCL2.fasta --refnames data/BCL2.names --target data/sample.fasta -c 5 --force --soi Owefus --outputfile Owenia_bcl2.tab --phymethod fasttree

# blastology smoke test
bash tests/smoke_blastology.sh

# pfamscan smoke test
bash tests/smoke_pfamscan.sh
```



# Generax   

GeneRax Example:  

```bash
module load OpenMPI
python main.py generax \
  --name Tubulin \
  --alignment generax_test/Tubulin.aln \
  --gene_tree generax_test/Tubulin.tree \
  --species_tree generax_test/species_tree.newick \
  --output_dir results_generax \
  --subs_model LG \
  --max_spr 1 \
  --cpus 1 \
  --logfile generax.log \
  --outfile updated.tree
```
