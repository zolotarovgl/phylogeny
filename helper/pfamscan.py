import csv
import logging
import os
import subprocess
import sys
from collections import defaultdict


def _check_pressed_hmm_db(hmm_db):
    required_suffixes = [".h3f", ".h3i", ".h3m", ".h3p"]
    missing = [hmm_db + suffix for suffix in required_suffixes if not os.path.exists(hmm_db + suffix)]
    if missing:
        logging.error(
            "Pressed HMM database files are missing for %s. Run `hmmpress %s` first.",
            hmm_db,
            hmm_db,
        )
        logging.error("Missing files: %s", ", ".join(missing))
        sys.exit(1)


def _parse_domtblout(domtblout_file, output_tsv):
    fieldnames = [
        "sequence_id",
        "pfam_name",
        "pfam_accession",
        "domain_index",
        "domain_count",
        "sequence_evalue",
        "sequence_score",
        "domain_cevalue",
        "domain_ievalue",
        "domain_score",
        "hmm_from",
        "hmm_to",
        "ali_from",
        "ali_to",
        "env_from",
        "env_to",
        "domain_acc",
        "description",
    ]

    hit_count = 0
    sequence_ids = set()
    rows = []
    with open(domtblout_file, "r") as infile, open(output_tsv, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames)
        for line in infile:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split(maxsplit=22)
            if len(columns) < 22:
                logging.warning("Skipping malformed hmmscan line: %s", line.rstrip("\n"))
                continue

            sequence_id = columns[3]
            description = columns[22] if len(columns) > 22 else ""
            row = {
                "sequence_id": sequence_id,
                "pfam_name": columns[0],
                "pfam_accession": columns[1],
                "domain_index": columns[9],
                "domain_count": columns[10],
                "sequence_evalue": columns[6],
                "sequence_score": columns[7],
                "domain_cevalue": columns[11],
                "domain_ievalue": columns[12],
                "domain_score": columns[13],
                "hmm_from": columns[15],
                "hmm_to": columns[16],
                "ali_from": columns[17],
                "ali_to": columns[18],
                "env_from": columns[19],
                "env_to": columns[20],
                "domain_acc": columns[21],
                "description": description,
            }
            writer.writerow([row[fieldname] for fieldname in fieldnames])
            rows.append(row)
            hit_count += 1
            sequence_ids.add(sequence_id)

    return rows, hit_count, len(sequence_ids)


def _write_domains_bed(rows, output_file):
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in sorted(rows, key=lambda x: (x["sequence_id"], int(x["ali_from"]), int(x["ali_to"]), x["pfam_name"])):
            writer.writerow(
                [
                    row["sequence_id"],
                    row["ali_from"],
                    row["ali_to"],
                    row["pfam_name"],
                ]
            )

def _write_architecture_file(rows, output_file, arch_sep=" "):
    grouped_rows = defaultdict(list)
    for row in rows:
        grouped_rows[row["sequence_id"]].append(row)

    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for sequence_id in sorted(grouped_rows):
            # Preserve domain order along the protein from left to right.
            protein_rows = sorted(
                grouped_rows[sequence_id],
                key=lambda x: (int(x["ali_from"]), int(x["ali_to"]), x["pfam_name"]),
            )
            start = min(int(row["ali_from"]) for row in protein_rows)
            end = max(int(row["ali_to"]) for row in protein_rows)
            architecture = arch_sep.join(row["pfam_name"] for row in protein_rows)
            writer.writerow([sequence_id, start, end, architecture])


def pfamscan(
    fasta_file,
    pfam_db,
    output_dir,
    prefix=None,
    ncpu=1,
    cut_ga=True,
    dom_evalue=None,
    arch_sep=",",
    force=False,
    verbose=True,
):
    if not os.path.exists(fasta_file):
        logging.error("Input FASTA file does not exist: %s", fasta_file)
        sys.exit(1)
    if not os.path.exists(pfam_db):
        logging.error("HMM database file does not exist: %s", pfam_db)
        sys.exit(1)
    if cut_ga and dom_evalue is not None:
        logging.error("Use either --cut_ga or --domE, not both.")
        sys.exit(1)

    _check_pressed_hmm_db(pfam_db)

    os.makedirs(output_dir, exist_ok=True)
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(fasta_file))[0]

    domtblout_file = os.path.join(output_dir, f"{prefix}.pfamscan.domtblout")
    stdout_file = os.path.join(output_dir, f"{prefix}.pfamscan.out")
    hits_file = os.path.join(output_dir, f"{prefix}.pfamscan.tsv")
    domains_file = os.path.join(output_dir, f"{prefix}.pfamscan.domains.csv")
    archs_file = os.path.join(output_dir, f"{prefix}.pfamscan_archs.csv")

    if os.path.exists(domtblout_file) and not force:
        logging.info("Found existing hmmscan result %s. Skipping hmmscan.", domtblout_file)
    else:
        cmd = ["hmmscan", "--cpu", str(ncpu), "--domtblout", domtblout_file]
        if cut_ga:
            cmd.append("--cut_ga")
        elif dom_evalue is not None:
            cmd.extend(["--domE", str(dom_evalue)])
        cmd.extend([pfam_db, fasta_file])

        if verbose:
            logging.info("Running hmmscan: %s", " ".join(cmd))
        with open(stdout_file, "w") as stdout_handle:
            subprocess.run(cmd, stdout=stdout_handle, stderr=subprocess.STDOUT, check=True)

    rows, hit_count, sequence_count = _parse_domtblout(domtblout_file, hits_file)
    _write_domains_bed(rows, domains_file)
    _write_architecture_file(rows, archs_file, arch_sep=arch_sep)
    logging.info(
        "Pfam scan done: %s hits across %s sequences. Parsed results: %s",
        hit_count,
        sequence_count,
        hits_file,
    )

    return {
        "domtblout": domtblout_file,
        "stdout": stdout_file,
        "tsv": hits_file,
        "domains": domains_file,
        "archs": archs_file,
        "hit_count": hit_count,
        "sequence_count": sequence_count,
    }
