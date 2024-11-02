- add more possvm options   
- make phylogeny run for less time 
- cluster stats output 
- separate log files 
- move cluster filtering from the main 
- make nicer temporary and output directories  
- manage --force and re-running  
- Add fasttree 
- Add HMM search 
- Clustering guided by representative sequence alignment and phylogeny  
- Sequence reclustering after initial clusters have been built - to include distant homologs in the species of interest  
- MMSeqs clustering should output a canonical clustering file - akin to broccoli and mcl (each cluster per line) => update parsing functions!  
- Think carefully how would you like to handle the clustering  
- Make it species-focused / Ctenphora focused. Such that you can flag the lineage-specific duplicates  

Installation:
```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylogeny.git
```

Examples: 
```bash
python main.py search -f data/sample.fasta -g data/genefam.tsv Myosin # creates results_annotation/searches/myo.Myosin.domains.fasta

mkdir -p results_annotation/alignments
python main.py align -f results_annotation/searches/myo.Myosin.domains.fasta -o results_annotation/alignments/test.aln -c 10
mkdir -p results_annotation/gene_trees
python main.py phylogeny -f results_annotation/alignments/test.aln -o results_annotation/gene_trees/test -c 15


# search the proteomes and annotate the hits
python main.py phylo-search --query example/rho_cdc42.fasta --target example/Mlei.fasta -c 3 -o test -p x


python main.py phylo-search --query Hsap_Myosin.fasta --target example/proteomes.fasta -c 15 -o results -s Mlei -r Hsap_gene_names.csv -p boo
```
