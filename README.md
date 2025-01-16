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

Use the reference sequences for search. For instance, given the reference sequences stored in `data/arp23.fasta`, search `data/sample.fasta`

```bash
python main.py blastology --query data/arp23.fasta --target data/sample.fasta -c 10 -p boo --force -r data/arp23.names
```

*TODO*: add reference sequence renaming!
