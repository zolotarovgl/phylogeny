- make the search outputs obvious from the command line name  
- mmseqs2 clustering  
- easy-phylo: implement intermediate file checks   
- add more possvm options   
- make phylogeny run for less time 

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
```
