# givn a file with Uniprot ids, fetch the sequences 
import sys, requests, re

if len(sys.argv) < 2:
	print("Usage: python fetch_uniprot.py <ids_file> [output_fasta]")
	sys.exit(1)
	

do_rename = True

ids_file = sys.argv[1]
out_file = sys.argv[2] if len(sys.argv) > 2 else None

with open(ids_file) as f:
	ids = [line.strip() for line in f if line.strip()]

query = " OR ".join(f"accession:{i}" for i in ids)
url = "https://rest.uniprot.org/uniprotkb/stream"
params = {"format": "fasta", "query": query}

species_dict = {'HUMAN' : 'Hsap', 'MOUSE': 'Mmus', 'DROME' : 'Dmel', 'RAT' : 'Ratnor', 'CELE' : 'Cele'}

r = requests.get(url, params=params)
r.raise_for_status()

renamed = []
for line in r.text.splitlines():
	if line.startswith(">"):
		m = re.match(r">(\S+)", line)
		if m:
			h = m.group(1)
			parts = h.split("_")
			db_acc = parts[0].split("|")
			species = parts[-1]
			
			if species in species_dict and do_rename:
				species = species_dict[species]
			if len(db_acc) == 3:
				acc = db_acc[1]
				protein_name = db_acc[2]
				if protein_name == acc:
					print(f'ERROR: could not find a proper protein name for {acc}! Remove it from the list')
					quit()
				new_id = f"{species}_QUERY_{acc}.{protein_name}"
			else:
				print("ERROR: Don't know how to process the header:")
				print(db_acc)
				quit()
			line = f">{new_id}"
	renamed.append(line)

if out_file:
	with open(out_file, "w") as out:
		out.write("\n".join(renamed) + "\n")
else:
	print("\n".join(renamed))


