#!/usr/bin/env python3
"""
PARP1 Automodification Domain (AD) Comparison Script
=====================================================
Fetches ~100 representative PARP1 sequences (>900 aa) from selected phyla
(Chordata, Echinodermata, Cnidaria, Porifera and related — excluding
Arthropoda, Nematoda, Annelida, Mollusca), aligns them with MAFFT, extracts
the region homologous to human PARP1 residues 481-526, and generates a
sequence logo.

Requirements:
    pip install biopython logomaker matplotlib numpy pandas
    apt install mafft   (or conda install -c bioconda mafft)

Usage:
    python parp1_ad_comparison.py [output_dir]
"""

import os, sys, time, re, textwrap, subprocess, tempfile
from io import StringIO
from collections import Counter

from Bio import Entrez, SeqIO, SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import MafftCommandline

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker

# ── Configuration ──────────────────────────────────────────────────────────
Entrez.email = "user@example.com"          # NCBI requires an email
Entrez.api_key = None                      # Optional: set your NCBI API key for faster queries

HUMAN_PARP1_ACCESSION = "NP_001609.2"      # RefSeq accession for human PARP1
AD_START = 481                             # 1-based start of AD region in human PARP1
AD_END   = 526                             # 1-based end   of AD region in human PARP1
TARGET_COUNT = 100                         # Approximate number of sequences desired
MIN_LENGTH = 900                           # Minimum protein length (amino acids)

# Taxonomy IDs
INCLUDE_TAXA = [
    7711,   # Chordata
    7586,   # Echinodermata
    6073,   # Cnidaria
    6040,   # Porifera
    10197,  # Ctenophora        (related basal metazoan)
    7568,   # Hemichordata      (related deuterostome)
    7735,   # Cephalochordata   (related chordate)
    7712,   # Urochordata       (related chordate)
]

EXCLUDE_TAXA = [
    6656,   # Arthropoda
    6231,   # Nematoda
    6340,   # Annelida
    6447,   # Mollusca
]

# ── Helper functions ───────────────────────────────────────────────────────

def fetch_human_parp1():
    """Download the human PARP1 protein sequence from NCBI."""
    print("Fetching human PARP1 sequence...")
    handle = Entrez.efetch(db="protein", id=HUMAN_PARP1_ACCESSION,
                           rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print(f"  {record.id}: {len(record.seq)} aa")
    return record


def build_entrez_query():
    """Build an NCBI Entrez query for PARP1 orthologs in target taxa."""
    gene_terms = (
        '("poly ADP-ribose polymerase 1"[Protein Name] OR '
        'PARP1[Gene Name] OR PARP-1[Gene Name] OR '
        '"poly(ADP-ribose) polymerase 1"[Protein Name] OR '
        'ADPRT[Gene Name] OR ADPRT1[Gene Name])'
    )
    include = " OR ".join(f"txid{t}[Organism:exp]" for t in INCLUDE_TAXA)
    exclude = " AND ".join(f"NOT txid{t}[Organism:exp]" for t in EXCLUDE_TAXA)
    length_filter = f'{MIN_LENGTH}:100000[Sequence Length]'
    query = f'{gene_terms} AND ({include}) AND {exclude} AND {length_filter}'
    return query


def search_sequences():
    """Search NCBI Protein for PARP1 orthologs matching our criteria."""
    query = build_entrez_query()
    print(f"Searching NCBI Protein database...")
    print(f"  Query: {query[:120]}...")

    handle = Entrez.esearch(db="protein", term=query, retmax=500, usehistory="y")
    results = Entrez.read(handle)
    handle.close()

    count = int(results["Count"])
    print(f"  Found {count} sequences")
    return results


def download_sequences(search_results):
    """Download sequences from NCBI in batches."""
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    count = int(search_results["Count"])

    records = []
    batch_size = 100
    for start in range(0, count, batch_size):
        print(f"  Downloading batch {start+1}-{min(start+batch_size, count)} of {count}...")
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text",
                                       retstart=start, retmax=batch_size,
                                       webenv=webenv, query_key=query_key)
                batch = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                records.extend(batch)
                break
            except Exception as e:
                print(f"    Retry {attempt+1}/3: {e}")
                time.sleep(2 ** (attempt + 1))
        time.sleep(0.4)  # Rate limiting
    return records


def blast_for_more(human_seq, existing_ids):
    """Use BLASTP to find additional PARP1 orthologs if Entrez search was insufficient."""
    print("Running BLASTP to find additional orthologs (this may take several minutes)...")

    include_str = " OR ".join(f"txid{t}[ORGN]" for t in INCLUDE_TAXA)
    exclude_str = " AND ".join(f"NOT txid{t}[ORGN]" for t in EXCLUDE_TAXA)
    entrez_query = f"({include_str}) AND {exclude_str}"

    result_handle = NCBIWWW.qblast(
        "blastp", "nr", human_seq.seq,
        entrez_query=entrez_query,
        hitlist_size=300,
        expect=1e-10,
    )
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    result_handle.close()

    new_ids = []
    for alignment in blast_record.alignments:
        acc = alignment.accession
        if acc not in existing_ids:
            # Check length from hit info
            for hsp in alignment.hsps:
                if hsp.sbjct_end - hsp.sbjct_start + 1 > MIN_LENGTH * 0.5:
                    new_ids.append(acc)
                    break

    print(f"  BLAST found {len(new_ids)} additional candidate accessions")
    return new_ids[:200]


def fetch_by_ids(id_list):
    """Fetch protein sequences by accession IDs."""
    records = []
    batch_size = 50
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i+batch_size]
        ids = ",".join(batch)
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="protein", id=ids,
                                       rettype="fasta", retmode="text")
                recs = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                records.extend(recs)
                break
            except Exception as e:
                print(f"    Retry {attempt+1}/3: {e}")
                time.sleep(2 ** (attempt + 1))
        time.sleep(0.4)
    return records


def filter_and_deduplicate(records, human_record):
    """Filter by length, remove exact duplicates, and pick one per species."""
    # Filter by minimum length
    long_enough = [r for r in records if len(r.seq) >= MIN_LENGTH]
    print(f"  After length filter (>={MIN_LENGTH} aa): {len(long_enough)} sequences")

    # Remove sequences with too many X residues
    clean = [r for r in long_enough
             if str(r.seq).count("X") / len(r.seq) < 0.05]
    print(f"  After removing X-heavy sequences: {len(clean)} sequences")

    # Extract organism from description and keep one per species
    species_map = {}
    for rec in clean:
        # Parse species from FASTA description, e.g. "[Homo sapiens]"
        m = re.search(r'\[([^\]]+)\]', rec.description)
        species = m.group(1) if m else rec.id
        if species not in species_map:
            species_map[species] = rec
        else:
            # Keep the longer sequence for each species
            if len(rec.seq) > len(species_map[species].seq):
                species_map[species] = rec

    unique = list(species_map.values())
    print(f"  After one-per-species dedup: {len(unique)} sequences")

    # Make sure human is included
    human_species = "Homo sapiens"
    if human_species not in species_map:
        unique.insert(0, human_record)

    return unique


def select_representatives(records, target_n):
    """Select approximately target_n representative sequences, balancing phyla."""
    if len(records) <= target_n:
        print(f"  Using all {len(records)} sequences (below target of {target_n})")
        return records

    # Just take a diverse subset — simple approach: sort by description
    # and take evenly-spaced entries
    records_sorted = sorted(records, key=lambda r: r.description)
    step = len(records_sorted) / target_n
    selected = [records_sorted[int(i * step)] for i in range(target_n)]

    # Ensure human is included
    human_ids = [r for r in selected
                 if "Homo sapiens" in r.description]
    if not human_ids:
        for r in records:
            if "Homo sapiens" in r.description:
                selected[0] = r
                break

    print(f"  Selected {len(selected)} representative sequences")
    return selected


def run_mafft(records, outdir):
    """Run MAFFT multiple sequence alignment."""
    input_fasta = os.path.join(outdir, "parp1_sequences.fasta")
    output_aln  = os.path.join(outdir, "parp1_alignment.fasta")

    # Make sequence IDs short and unique for MAFFT
    clean_records = []
    seen_ids = set()
    for i, rec in enumerate(records):
        m = re.search(r'\[([^\]]+)\]', rec.description)
        species = m.group(1) if m else f"seq_{i}"
        # Create a short, unique ID
        parts = species.split()
        if len(parts) >= 2:
            short_id = f"{parts[0][:3]}_{parts[1][:6]}".replace(".", "")
        else:
            short_id = species[:10].replace(" ", "_")
        # Ensure uniqueness
        base_id = short_id
        counter = 1
        while short_id in seen_ids:
            short_id = f"{base_id}_{counter}"
            counter += 1
        seen_ids.add(short_id)
        clean_records.append(SeqRecord(rec.seq, id=short_id, description=""))

    SeqIO.write(clean_records, input_fasta, "fasta")
    print(f"Wrote {len(clean_records)} sequences to {input_fasta}")

    print("Running MAFFT alignment (this may take a few minutes)...")
    cmd = ["mafft", "--auto", "--thread", "-1", input_fasta]
    with open(output_aln, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"  MAFFT stderr: {result.stderr[:500]}")
        raise RuntimeError("MAFFT alignment failed")
    print(f"  Alignment saved to {output_aln}")
    return output_aln


def find_human_in_alignment(aln_records):
    """Find the human PARP1 record in the alignment."""
    for rec in aln_records:
        rid = rec.id.lower()
        if "hom_sapien" in rid or "homo" in rid or "hom_sap" in rid:
            return rec
    # Fallback: first record (should be human if we inserted it first)
    return aln_records[0]


def extract_ad_region(aln_file, ad_start, ad_end):
    """
    Extract the alignment columns corresponding to human PARP1 positions
    ad_start to ad_end (1-based, ungapped numbering).
    """
    aln_records = list(SeqIO.parse(aln_file, "fasta"))
    human_rec = find_human_in_alignment(aln_records)
    human_aln_seq = str(human_rec.seq)

    # Map ungapped positions to alignment columns
    ungapped_pos = 0
    col_start = None
    col_end = None
    for col_idx, aa in enumerate(human_aln_seq):
        if aa != "-":
            ungapped_pos += 1
            if ungapped_pos == ad_start:
                col_start = col_idx
            if ungapped_pos == ad_end:
                col_end = col_idx
                break

    if col_start is None or col_end is None:
        raise ValueError(f"Could not map human positions {ad_start}-{ad_end} "
                         f"to alignment columns (reached pos {ungapped_pos})")

    print(f"  Human PARP1 {ad_start}-{ad_end} maps to alignment columns "
          f"{col_start+1}-{col_end+1}")

    # Extract the sub-alignment
    ad_records = []
    for rec in aln_records:
        sub_seq = str(rec.seq)[col_start:col_end+1]
        # Skip sequences that are entirely gaps in this region
        if sub_seq.replace("-", ""):
            ad_records.append(SeqRecord(Seq(sub_seq), id=rec.id, description=""))

    print(f"  Extracted AD region from {len(ad_records)} sequences "
          f"(length {col_end - col_start + 1} columns)")
    return ad_records, col_end - col_start + 1


def make_sequence_logo(ad_records, region_length, outdir):
    """Generate a sequence logo from the AD region alignment."""
    output_png = os.path.join(outdir, "parp1_ad_logo.png")
    output_pdf = os.path.join(outdir, "parp1_ad_logo.pdf")
    output_csv = os.path.join(outdir, "parp1_ad_matrix.csv")

    # Build a counts matrix
    standard_aa = list("ACDEFGHIKLMNPQRSTVWY")
    n_cols = region_length
    counts = {aa: [0] * n_cols for aa in standard_aa}
    counts["-"] = [0] * n_cols  # track gaps

    for rec in ad_records:
        seq = str(rec.seq)
        for pos, aa in enumerate(seq):
            aa_upper = aa.upper()
            if aa_upper in counts:
                counts[aa_upper][pos] += 1

    # Convert to DataFrame (exclude gaps for the logo)
    df = pd.DataFrame({aa: counts[aa] for aa in standard_aa})
    df.index = list(range(AD_START, AD_START + n_cols))

    # Normalize to frequency
    row_sums = df.sum(axis=1)
    df_freq = df.div(row_sums, axis=0).fillna(0)

    # Compute information content
    df_info = logomaker.transform_matrix(df_freq, from_type="probability",
                                         to_type="information")

    # Save matrix
    df_info.to_csv(output_csv)

    # Create logo
    fig, ax = plt.subplots(1, 1, figsize=(max(16, n_cols * 0.4), 3.5))

    logo = logomaker.Logo(df_info, ax=ax, shade_below=0.5,
                          fade_below=0.5, color_scheme="chemistry")
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)

    ax.set_ylabel("Information (bits)", fontsize=12)
    ax.set_xlabel(f"Human PARP1 position (residues {AD_START}-{AD_END})", fontsize=12)
    ax.set_title("PARP1 Automodification Domain — Sequence Logo", fontsize=14)

    plt.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  Sequence logo saved to {output_png} and {output_pdf}")
    print(f"  Information matrix saved to {output_csv}")


def save_ad_alignment(ad_records, outdir):
    """Save the extracted AD region alignment to a FASTA file."""
    outpath = os.path.join(outdir, "parp1_ad_region.fasta")
    SeqIO.write(ad_records, outpath, "fasta")
    print(f"  AD region alignment saved to {outpath}")


# ── Main pipeline ──────────────────────────────────────────────────────────

def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "parp1_ad_results"
    os.makedirs(outdir, exist_ok=True)
    print(f"Output directory: {outdir}\n")

    # Step 1: Fetch human PARP1
    human_rec = fetch_human_parp1()

    # Step 2: Search for orthologs via Entrez
    search_results = search_sequences()
    total_found = int(search_results["Count"])

    # Step 3: Download sequences
    all_records = download_sequences(search_results)
    print(f"  Downloaded {len(all_records)} sequences total")

    # Step 4: If we didn't get enough, try BLAST
    if len(all_records) < TARGET_COUNT:
        existing_ids = {r.id for r in all_records}
        blast_ids = blast_for_more(human_rec, existing_ids)
        if blast_ids:
            extra = fetch_by_ids(blast_ids)
            all_records.extend(extra)
            print(f"  Total after BLAST supplement: {len(all_records)} sequences")

    # Step 5: Filter and deduplicate
    print("\nFiltering sequences...")
    filtered = filter_and_deduplicate(all_records, human_rec)

    # Step 6: Select representatives
    print("\nSelecting representatives...")
    representatives = select_representatives(filtered, TARGET_COUNT)

    # Step 7: Multiple sequence alignment
    print(f"\nAligning {len(representatives)} sequences...")
    aln_file = run_mafft(representatives, outdir)

    # Step 8: Extract AD region (human PARP1 481-526)
    print(f"\nExtracting AD region (human PARP1 positions {AD_START}-{AD_END})...")
    ad_records, region_length = extract_ad_region(aln_file, AD_START, AD_END)
    save_ad_alignment(ad_records, outdir)

    # Step 9: Generate sequence logo
    print("\nGenerating sequence logo...")
    make_sequence_logo(ad_records, region_length, outdir)

    print(f"\n{'='*60}")
    print("Pipeline complete! Output files:")
    for f in sorted(os.listdir(outdir)):
        fpath = os.path.join(outdir, f)
        size = os.path.getsize(fpath)
        print(f"  {f} ({size:,} bytes)")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
