# Installation:
```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylogeny.git
```

To use GeneRax, install it from source (`https://github.com/BenoitMorel/GeneRax`) or from conda (`https://anaconda.org/bioconda/generax`).


# Commands   


```bash
usage: main.py [-h] {hmmsearch,cluster,align,phylogeny,generax,possvm,easy-phylo,blastology} ...

Python wrapper around some useful commands

positional arguments:
  {hmmsearch,cluster,align,phylogeny,generax,possvm,easy-phylo,blastology}
                        Sub-command help
    hmmsearch           Search for a family using HMMER
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

### Clustering   


An example of an iterative re-clustering for a case where the biggest cluster exceeds the maximum number of sequences (`-m`):  


```bash
python main.py cluster -f data/cluster_test.fasta --out_prefix boo -c 4 -i 2 -m 70
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

# Tests   


```bash 
# alignment and trimming 
INFASTA=data/arp23.fasta
python main.py align -f $INFASTA -o test.aln --notrim
python main.py trim -f test.aln -o test.alt --logfile test.log

# easy-phylo
python main.py easy-phylo --method fasttree -f data/arp23.fasta -c 16 --outdir results
# blastology
python main.py blastology --query data/BCL2.fasta --refnames data/BCL2.names --target data/sample.fasta -c 5 --force --soi Owefus --outputfile Owenia_bcl2.tab --phymethod fasttree
```



# Generax   

GeneRax Example:  

```bash
python main.py generax \
  --name Tubulin \
  --alignment generax_test/Tubulin.aln \
  --gene_tree generax_test/Tubulin.tree \
  --species_tree generax_test/species_tree.newick \
  --output_dir results_generax \
  --subs_model LG \
  --max-spr 3 \
  --cpus 1 \
  --outfile updated.tree
```