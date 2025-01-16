# Installation:
```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylogeny.git
```

# Examples  

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


### BLASTOLOGY 
The tool will use BLASTP to search, then it will cluster the sequences with MCL&diamond and run the phylogeny for each cluster. 
Finally, it will parse the resulting trees and will output the best orthologs for the queries.  

Use the reference sequences for search. For instance, given the reference sequences stored in `data/BCL2.fasta`, search `data/sample.fasta` for homologs in species `Owefus` using 10 cores:  
```bash
python main.py blastology --query data/BCL2.fasta -r data/BCL2.names  --target data/sample.fasta -c 10 --force --soi Owefus --outputfile Owenia_bcl2.tab
```
Not specifying `--soi` will return the annotation file with all sequences in the orthogroups with queries.  

__QUERY FORMATTING:__  
It's important that the query names are formatted as following:  
```
>Hsap_QUERY_P10415
>Mmus_QUERY_P10417
>Dmel_QUERY_Q7KM33
```
First prefix MUST correspond to prefixes in Xavi's databas to not to mess up the species reconciliation. It's also good to include the word `QUERY` in the name. 
The `.names` file should be map the names in the fasta to the protein names (these will be used for name transfer):   
```
Hsap_QUERY_P10415       hsapBCL2
Mmus_QUERY_P10417       mmusBCL2
Dmel_QUERY_Q7KM33       dmelBCL2
```
The names (2nd column) can be any, but try to keep them distinguishable!  

**TODO**: add species prefix check - strict and not. Request that all the reference species are present in the target file.
