import sys, requests, re

if len(sys.argv) < 2:
	print("Usage: python fetch_uniprot.py <ids_file> [output_fasta]")
	sys.exit(1)


def fetch_uniprot_desc(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    r = requests.get(url)
    if r.status_code != 200:
        return None, None

    data = r.json()

    return(data["proteinDescription"]["recommendedName"])




def fetch_uniprot_names(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    r = requests.get(url)
    if r.status_code != 200:
        return None, None

    data = r.json()

    long_name = None
    short_name = None
    
    desc = fetch_uniprot_desc(acc)
    

    try:
        long_name = desc["fullName"]["value"]
    except KeyError:
        pass

    try:
        short_name = desc["shortNames"][0]["value"]
    except KeyError:
        pass
    return long_name, short_name


ids_file = sys.argv[1]
out_file = sys.argv[2] if len(sys.argv) > 2 else None

with open(ids_file) as f:
	ids = [line.strip() for line in f if line.strip()]

d = {}
for id in ids:
    d[id] = fetch_uniprot_names(id)


for k,v in d.items():
    long = v[0]
    short = v[1]
    if not short:
        short = ''
    if not long:
        long = ''
    name = short
    if not name:
        if not long: 
            name = ""
        else:
            name = long

    print(f'{k}\t' + long)
