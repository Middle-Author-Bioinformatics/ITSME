#!/usr/bin/env python3
"""Find and classify eukaryotic rRNA/ITS loci in an assembly.

Starting from an assembly, the script screens contigs against the targeted NCBI
SSU, LSU, and ITS databases, extracts matching contigs, runs ITSx, and applies
ITSME's marker-specific BLAST and near-top-hit consensus. SSU and LSU determine
the primary locus taxonomy; ITS is reported as independent supporting evidence.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


RANKS = ("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
RAW_FIELDS = (
    "query", "query_length", "accession", "percent_identity", "aligned_bp",
    "query_coverage_hsp", "evalue", "bitscore", "taxids", "title",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Screen an assembly for eukaryotic rRNA/ITS contigs, run ITSx, and "
            "generate an ITSME-style taxonomy summary."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-i", "--input",
        help="Original SPAdes contigs FASTA; performs screening and ITSx automatically",
    )
    source.add_argument(
        "--itsx-prefix",
        help="Reuse existing ITSx output instead of starting from an assembly",
    )
    parser.add_argument(
        "--db-dir", required=True,
        help="ITSME database directory containing NCBI_rRNA_BLAST/ and taxdump/",
    )
    parser.add_argument(
        "--taxdump-dir",
        help="Override the taxdump directory inferred from --db-dir",
    )
    parser.add_argument("-o", "--output-dir", default="itsx_taxonomy")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("--max-hits", type=int, default=50)
    parser.add_argument("--report-hits", type=int, default=10)
    parser.add_argument(
        "--screen-task", choices=("megablast", "dc-megablast", "blastn"),
        default="megablast",
        help="Initial assembly-screening sensitivity (default: megablast)",
    )
    parser.add_argument(
        "--screen-evalue", type=float, default=1e-10,
        help="E-value threshold for the initial contig screen (default: 1e-10)",
    )
    parser.add_argument(
        "--screen-min-aligned", type=int, default=80,
        help="Minimum aligned bases for retaining a candidate contig (default: 80)",
    )
    parser.add_argument(
        "--itsx-taxa", default="M",
        help="ITSx taxon code(s); M is Metazoa (default: M)",
    )
    parser.add_argument(
        "--partial-min-length", type=int, default=100,
        help="Minimum length for ITSx partial-region reporting (default: 100)",
    )
    parser.add_argument(
        "--near-top-fraction", type=float, default=0.95,
        help="Include hits with bitscore >= this fraction of the top hit (default: 0.95)",
    )
    parser.add_argument(
        "--skip-blast", action="store_true",
        help="Reuse raw BLAST tables already present in the output directory",
    )
    args = parser.parse_args()
    if not 0 < args.near_top_fraction <= 1:
        parser.error("--near-top-fraction must be >0 and <=1")
    if (args.threads < 1 or args.max_hits < 1 or args.report_hits < 1 or
            args.screen_min_aligned < 1 or args.partial_min_length < 1):
        parser.error("thread and hit counts must be positive integers")
    return args


def require_program(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"ERROR: required program not found in PATH: {name}")


def prefix_file(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}")


def choose_region_file(prefix: Path, region: str) -> Path | None:
    for suffix in (f"{region}.full_and_partial.fasta", f"{region}.fasta"):
        path = prefix_file(prefix, suffix)
        if path.exists() and path.stat().st_size:
            return path
    return None


def read_fasta(path: Path | None) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    if path is None:
        return {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records.setdefault(current, [])
            elif current is not None:
                records[current].append(line)
    return {name: "".join(sequence) for name, sequence in records.items()}


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for name in sorted(records):
            sequence = records[name]
            if not sequence:
                continue
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def extract_fasta_records(source: Path, wanted: set[str], output: Path) -> int:
    """Stream selected records without loading a potentially huge assembly."""
    written = 0
    keep = False
    with source.open() as incoming, output.open("w") as outgoing:
        for line in incoming:
            if line.startswith(">"):
                sequence_id = line[1:].split()[0]
                keep = sequence_id in wanted
                if keep:
                    written += 1
            if keep:
                outgoing.write(line)
    return written


def parse_positions(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"ERROR: ITSx positions file not found: {path}")
    catalog: dict[str, dict[str, str]] = {}
    labels = {"SSU": "SSU", "ITS1": "ITS1", "5.8S": "5.8S", "ITS2": "ITS2", "LSU": "LSU"}
    with path.open() as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            contig = fields[0].split()[0]
            length_match = re.search(r"(\d+)\s*bp", fields[1])
            row = {"contig": contig, "length_bp": length_match.group(1) if length_match else "NA"}
            for label in labels:
                row[label] = "NA"
            notes = []
            for field in fields[2:]:
                match = re.match(r"^(SSU|ITS1|5\.8S|ITS2|LSU):\s*(.*)$", field)
                if match:
                    value = match.group(2).strip()
                    row[match.group(1)] = value if value else "NA"
                elif field.strip():
                    notes.append(field.strip())
            row["itsx_note"] = " ".join(notes) if notes else "NA"
            catalog[contig] = row
    return catalog


def fasta_ids(path: Path) -> set[str]:
    return set(read_fasta(path)) if path.exists() and path.stat().st_size else set()


def find_databases(root: Path) -> tuple[Path, Path]:
    blast_dir = root / "NCBI_rRNA_BLAST" if (root / "NCBI_rRNA_BLAST").is_dir() else root
    taxdump = root / "taxdump"
    return blast_dir, taxdump


def build_blast_alias(blast_dir: Path, alias: Path) -> None:
    databases = [
        blast_dir / "SSU_eukaryote_rRNA",
        blast_dir / "LSU_eukaryote_rRNA",
        blast_dir / "ITS_eukaryote_sequences",
    ]
    subprocess.run(
        [
            "blastdb_aliastool", "-dblist", " ".join(map(str, databases)),
            "-dbtype", "nucl", "-out", str(alias),
            "-title", "Eukaryotic SSU LSU ITS",
        ],
        check=True,
    )


def screen_assembly(assembly: Path, blast_dir: Path, outdir: Path,
                    task: str, evalue: float, minimum_aligned: int,
                    threads: int) -> Path:
    alias = outdir / "euk_rrna"
    if not prefix_file(alias, "nal").exists():
        print("Building combined SSU/LSU/ITS BLAST alias...", file=sys.stderr)
        build_blast_alias(blast_dir, alias)

    hits = outdir / "screen_hits.tsv"
    print(f"Screening assembly with {task}: {assembly}", file=sys.stderr)
    subprocess.run(
        [
            "blastn", "-query", str(assembly), "-db", str(alias),
            "-task", task, "-evalue", str(evalue),
            "-max_target_seqs", "1", "-max_hsps", "1",
            "-num_threads", str(threads),
            "-outfmt", "6 qseqid saccver pident length qcovhsp evalue bitscore",
            "-out", str(hits),
        ],
        check=True,
    )

    identifiers: set[str] = set()
    with hits.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                try:
                    if int(fields[3]) >= minimum_aligned:
                        identifiers.add(fields[0])
                except ValueError:
                    continue
    ids_path = outdir / "candidate_contig_ids.txt"
    ids_path.write_text("".join(f"{name}\n" for name in sorted(identifiers)))
    if not identifiers:
        raise SystemExit(
            "ERROR: no candidate contigs passed the screen. Try "
            "--screen-task dc-megablast or reduce --screen-min-aligned."
        )

    candidates = outdir / "rrna_candidate_contigs.fasta"
    extracted = extract_fasta_records(assembly, identifiers, candidates)
    if extracted != len(identifiers):
        print(
            f"WARNING: requested {len(identifiers)} contigs but extracted {extracted}; "
            "check for truncated or duplicate FASTA identifiers.",
            file=sys.stderr,
        )
    print(f"Extracted {extracted} candidate contigs: {candidates}", file=sys.stderr)
    return candidates


def run_itsx(candidates: Path, prefix: Path, taxa: str,
             partial_minimum: int, threads: int) -> None:
    print(f"Running ITSx on candidate contigs: {candidates}", file=sys.stderr)
    subprocess.run(
        [
            "ITSx", "-i", str(candidates), "-o", str(prefix),
            "-t", taxa, "--save_regions", "all",
            "--partial", str(partial_minimum), "--positions", "T",
            "--detailed_results", "T", "--cpu", str(threads),
        ],
        check=True,
    )


def run_blast(query: Path, database: Path, output: Path, threads: int,
              max_hits: int, marker: str) -> None:
    output.write_text("")
    if not query.exists() or not query.stat().st_size:
        return
    command = [
        "blastn", "-query", str(query), "-db", str(database),
        "-task", "blastn", "-max_target_seqs", str(max_hits),
        "-max_hsps", "1", "-num_threads", str(threads),
        "-outfmt", "6 qseqid qlen saccver pident length qcovhsp evalue bitscore staxids stitle",
        "-out", str(output),
    ]
    if marker == "ITS":
        command += ["-word_size", "7", "-dust", "no", "-evalue", "1e-5"]
    else:
        command += ["-evalue", "1e-20"]
    subprocess.run(command, check=True)


def load_taxdump(path: Path):
    nodes = path / "nodes.dmp"
    names_file = path / "names.dmp"
    if not nodes.exists() or not names_file.exists():
        raise SystemExit(f"ERROR: taxdump requires nodes.dmp and names.dmp: {path}")
    parents: dict[str, str] = {}
    tax_ranks: dict[str, str] = {}
    names: dict[str, str] = {}
    merged: dict[str, str] = {}
    with nodes.open() as handle:
        for line in handle:
            fields = [value.strip() for value in line.split("|")]
            parents[fields[0]], tax_ranks[fields[0]] = fields[1], fields[2]
    with names_file.open() as handle:
        for line in handle:
            fields = [value.strip() for value in line.split("|")]
            if len(fields) > 3 and fields[3] == "scientific name":
                names[fields[0]] = fields[1]
    merged_file = path / "merged.dmp"
    if merged_file.exists():
        with merged_file.open() as handle:
            for line in handle:
                fields = [value.strip() for value in line.split("|")]
                merged[fields[0]] = fields[1]
    return parents, tax_ranks, names, merged


def lineage_for(taxids: str, taxonomy) -> tuple[dict[str, str], str]:
    parents, tax_ranks, names, merged = taxonomy
    taxid = next((value for value in re.split(r"[;,]", taxids) if value.isdigit()), "")
    while taxid in merged:
        taxid = merged[taxid]
    aliases = {"superkingdom": "domain", **{rank: rank for rank in RANKS}}
    result: dict[str, str] = {rank: "NA" for rank in RANKS}
    current, seen = taxid, set()
    while current and current not in seen and current in parents:
        seen.add(current)
        rank = aliases.get(tax_ranks.get(current, ""))
        if rank and result[rank] == "NA":
            result[rank] = names.get(current, "NA")
        parent = parents[current]
        if parent == current:
            break
        current = parent
    return result, names.get(taxid, "NA")


def parse_raw_hits(path: Path, marker: str, taxonomy) -> list[dict[str, str]]:
    rows = []
    if not path.exists() or not path.stat().st_size:
        return rows
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 9)
            if len(fields) != 10:
                continue
            row = dict(zip(RAW_FIELDS, fields))
            ranks, scientific_name = lineage_for(row["taxids"], taxonomy)
            row.update(ranks)
            row["scientific_names"] = scientific_name
            row["marker"] = marker
            rows.append(row)
    return rows


def resolved(value: str | None) -> bool:
    return bool(value and value not in {"NA", "N/A", "Unclassified", "unclassified"})


def marker_consensus(rows: list[dict[str, str]], fraction: float) -> dict[str, str]:
    consensus = {rank: "NA" for rank in RANKS}
    if not rows:
        return consensus
    top = max(float(row["bitscore"]) for row in rows)
    selected = [row for row in rows if float(row["bitscore"]) >= top * fraction]
    for rank in RANKS:
        values = {row[rank] for row in selected if resolved(row.get(rank))}
        if not values:
            continue
        if len(values) != 1:
            break
        consensus[rank] = values.pop()
    return consensus


def lineage_text(lineage: dict[str, str]) -> str:
    values = [lineage[rank] for rank in RANKS if resolved(lineage.get(rank))]
    return "; ".join(values) if values else "NA"


def intersect_lineages(ssu: dict[str, str], lsu: dict[str, str]) -> dict[str, str]:
    result = {rank: "NA" for rank in RANKS}
    for rank in RANKS:
        left, right = ssu[rank], lsu[rank]
        if not resolved(left) or not resolved(right):
            continue
        if left.casefold() != right.casefold():
            break
        result[rank] = left
    return result


def top_hit(rows: list[dict[str, str]]) -> dict[str, str]:
    return max(rows, key=lambda row: float(row["bitscore"])) if rows else {}


def spades_depth(contig: str) -> str:
    match = re.search(r"_cov_([0-9.]+)(?:\s|$)", contig)
    return match.group(1) if match else "NA"


def main() -> None:
    args = parse_args()
    require_program("blastn")
    root = Path(args.db_dir).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    blast_dir, inferred_taxdump = find_databases(root)
    taxdump = Path(args.taxdump_dir).expanduser().resolve() if args.taxdump_dir else inferred_taxdump

    if args.input:
        require_program("blastdb_aliastool")
        require_program("ITSx")
        assembly = Path(args.input).expanduser().resolve()
        if not assembly.exists() or not assembly.stat().st_size:
            raise SystemExit(f"ERROR: input assembly is missing or empty: {assembly}")
        candidates = screen_assembly(
            assembly, blast_dir, outdir, args.screen_task,
            args.screen_evalue, args.screen_min_aligned, args.threads,
        )
        prefix = outdir / "itsx"
        run_itsx(candidates, prefix, args.itsx_taxa,
                 args.partial_min_length, args.threads)
    else:
        prefix = Path(args.itsx_prefix).expanduser().resolve()

    positions = prefix_file(prefix, "positions.txt")
    catalog = parse_positions(positions)
    complete_ids = fasta_ids(prefix_file(prefix, "full.fasta"))

    queries = {
        "SSU": choose_region_file(prefix, "SSU"),
        "LSU": choose_region_file(prefix, "LSU"),
    }
    its_parts = {
        region: read_fasta(choose_region_file(prefix, region))
        for region in ("ITS1", "5_8S", "ITS2")
    }
    its_records = {
        contig: "".join(its_parts[region].get(contig, "") for region in ("ITS1", "5_8S", "ITS2"))
        for contig in catalog
    }
    its_query = outdir / "ITS_queries.fasta"
    write_fasta({name: seq for name, seq in its_records.items() if seq}, its_query)
    queries["ITS"] = its_query

    raw_paths = {marker: outdir / f"ncbi_{marker}_hits.tsv" for marker in ("SSU", "ITS", "LSU")}
    databases = {
        "SSU": blast_dir / "SSU_eukaryote_rRNA",
        "ITS": blast_dir / "ITS_eukaryote_sequences",
        "LSU": blast_dir / "LSU_eukaryote_rRNA",
    }
    if not args.skip_blast:
        for marker in ("SSU", "ITS", "LSU"):
            query = queries[marker]
            if query is None:
                raw_paths[marker].write_text("")
                continue
            print(f"BLASTing {marker}: {query}", file=sys.stderr)
            run_blast(query, databases[marker], raw_paths[marker], args.threads, args.max_hits, marker)

    print("Loading NCBI taxonomy...", file=sys.stderr)
    taxonomy = load_taxdump(taxdump)
    all_rows = []
    for marker in ("SSU", "ITS", "LSU"):
        marker_rows = parse_raw_hits(raw_paths[marker], marker, taxonomy)
        grouped_marker = defaultdict(list)
        for row in marker_rows:
            grouped_marker[row["query"]].append(row)
        for query in grouped_marker:
            grouped_marker[query].sort(key=lambda row: float(row["bitscore"]), reverse=True)
            all_rows.extend(grouped_marker[query][:args.report_hits])

    top_fields = [
        "marker", "query", "rank", "query_length", "accession", "percent_identity",
        "aligned_bp", "query_coverage_hsp", "evalue", "bitscore", "taxids",
        "scientific_names", *RANKS, "title",
    ]
    top_path = outdir / "ncbi_blast_top_hits.tsv"
    grouped = defaultdict(list)
    with top_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=top_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for marker in ("SSU", "ITS", "LSU"):
            marker_rows = [row for row in all_rows if row["marker"] == marker]
            by_query = defaultdict(list)
            for row in marker_rows:
                by_query[row["query"]].append(row)
            for query in sorted(by_query):
                for hit_rank, row in enumerate(by_query[query], 1):
                    output = {field: row.get(field, "NA") for field in top_fields}
                    output["rank"] = hit_rank
                    writer.writerow(output)
                    grouped[(query, marker)].append(row)

    summary_fields = [
        "contig", "locus_type", "length_bp", "SPAdes_mean_depth",
        "SSU_coordinates", "ITS1_coordinates", "5.8S_coordinates",
        "ITS2_coordinates", "LSU_coordinates", "taxonomy_status",
        "consensus_taxonomy", "SSU_consensus_taxonomy", "ITS_consensus_taxonomy",
        "LSU_consensus_taxonomy", "SSU_top_hit", "SSU_identity", "SSU_aligned_bp",
        "ITS_top_hit", "ITS_identity", "ITS_aligned_bp", "LSU_top_hit",
        "LSU_identity", "LSU_aligned_bp", "itsx_note",
    ]
    summary_path = outdir / "master_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        for contig in catalog:
            marker_tax = {
                marker: marker_consensus(grouped.get((contig, marker), []), args.near_top_fraction)
                for marker in ("SSU", "ITS", "LSU")
            }
            ssu, its, lsu = marker_tax["SSU"], marker_tax["ITS"], marker_tax["LSU"]
            ssu_phylum, lsu_phylum = ssu["phylum"], lsu["phylum"]
            if resolved(ssu_phylum) and resolved(lsu_phylum):
                if ssu_phylum.casefold() != lsu_phylum.casefold():
                    status = "CHIMERA"
                    combined = intersect_lineages(ssu, lsu)
                else:
                    status = "PASS_taxonomically_coherent"
                    combined = intersect_lineages(ssu, lsu)
            elif resolved(ssu_phylum):
                status, combined = "PARTIAL_taxonomically_assigned", ssu
            elif resolved(lsu_phylum):
                status, combined = "PARTIAL_taxonomically_assigned", lsu
            else:
                status, combined = "UNRESOLVED", {rank: "NA" for rank in RANKS}

            tops = {marker: top_hit(grouped.get((contig, marker), [])) for marker in ("SSU", "ITS", "LSU")}
            record = catalog[contig]
            output = {
                "contig": contig,
                "locus_type": "COMPLETE_ITS" if contig in complete_ids else "PARTIAL_ITS",
                "length_bp": record["length_bp"],
                "SPAdes_mean_depth": spades_depth(contig),
                "SSU_coordinates": record["SSU"],
                "ITS1_coordinates": record["ITS1"],
                "5.8S_coordinates": record["5.8S"],
                "ITS2_coordinates": record["ITS2"],
                "LSU_coordinates": record["LSU"],
                "taxonomy_status": status,
                "consensus_taxonomy": lineage_text(combined),
                "SSU_consensus_taxonomy": lineage_text(ssu),
                "ITS_consensus_taxonomy": lineage_text(its),
                "LSU_consensus_taxonomy": lineage_text(lsu),
                "itsx_note": record["itsx_note"],
            }
            for marker in ("SSU", "ITS", "LSU"):
                hit = tops[marker]
                output[f"{marker}_top_hit"] = hit.get("scientific_names", "NA")
                output[f"{marker}_identity"] = hit.get("percent_identity", "NA")
                output[f"{marker}_aligned_bp"] = hit.get("aligned_bp", "NA")
            writer.writerow(output)

    print(f"Wrote {len(catalog)} loci to {summary_path}", file=sys.stderr)
    print(f"Detailed BLAST hits: {top_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
