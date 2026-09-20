#!/usr/bin/env python3
"""Compare per-sample ITSME and assembly-wide ITSx taxonomy summaries.

Only samples represented in the ITSME full-summary directory are included.
ITSx-only samples are deliberately ignored. Multiple loci and taxa are retained
as separate rows.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


RANKS = ("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
COORDINATE_FIELDS = (
    "SSU_coordinates", "ITS1_coordinates", "5.8S_coordinates",
    "ITS2_coordinates", "LSU_coordinates",
)
MISSING = {"", "na", "n/a", "none", "not found"}
PARTIAL_COORDINATES = {"no start", "no end"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare ITSME full/partial results with assembly-wide ITSx summaries. "
            "ITSx-only samples are ignored."
        )
    )
    parser.add_argument(
        "--itsme-summaries", required=True,
        help="Directory containing ITSME full/reconstructed summary CSV files",
    )
    parser.add_argument(
        "--itsme-partials", required=True,
        help="Directory containing ITSME partial-locus CSV files",
    )
    parser.add_argument(
        "--itsx-summaries", required=True,
        help="Directory containing assembly-wide ITSx master-summary CSV files",
    )
    parser.add_argument(
        "-o", "--output-dir", default="itsme_itsx_comparison",
        help="Output directory (default: itsme_itsx_comparison)",
    )
    return parser.parse_args()


def canonical_sample(path: Path, source: str) -> str:
    """Convert method-specific result filenames to a shared sample identifier."""
    name = path.stem
    if source == "ITSME":
        name = re.sub(r"^itsme_", "", name, flags=re.IGNORECASE)
        name = re.sub(r"-2$", "", name)
    else:
        name = re.sub(r"(?:_spades|-spades)$", "", name, flags=re.IGNORECASE)
    return name


def require_directory(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise SystemExit(f"ERROR: {label} directory not found: {path}")
    return path


def csv_files_by_sample(directory: Path, source: str) -> dict[str, list[Path]]:
    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in sorted(directory.glob("*.csv")):
        grouped[canonical_sample(path, source)].append(path)
    return dict(grouped)


def read_csv_rows(paths: list[Path]) -> list[tuple[Path, dict[str, str]]]:
    records: list[tuple[Path, dict[str, str]]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                print(f"WARNING: empty CSV skipped: {path}", file=sys.stderr)
                continue
            for row in reader:
                if not any((value or "").strip() for value in row.values()):
                    continue
                records.append((path, {key: (value or "").strip() for key, value in row.items()}))
    return records


def coordinate_state(value: str) -> str:
    normalized = value.strip().casefold()
    if normalized in MISSING:
        return "ABSENT"
    if normalized in PARTIAL_COORDINATES:
        return "PARTIAL"
    if re.fullmatch(r"\d+\s*-\s*\d+", value.strip()):
        return "PRESENT"
    return "PARTIAL"


def regions_present(row: dict[str, str]) -> str:
    labels = {
        "SSU_coordinates": "SSU",
        "ITS1_coordinates": "ITS1",
        "5.8S_coordinates": "5.8S",
        "ITS2_coordinates": "ITS2",
        "LSU_coordinates": "LSU",
    }
    regions = []
    for field in COORDINATE_FIELDS:
        state = coordinate_state(row.get(field, "NA"))
        if state == "PRESENT":
            regions.append(labels[field])
        elif state == "PARTIAL":
            regions.append(f"{labels[field]}(partial)")
    if regions:
        return "; ".join(regions)
    locus_type = row.get("locus_type", "").upper()
    if "18S" in locus_type:
        return "SSU"
    if "28S" in locus_type:
        return "LSU"
    return "NA"


def sequence_class(row: dict[str, str], table_kind: str) -> str:
    """Use numeric ITS1+5.8S+ITS2 coordinates as the strict full-ITS criterion."""
    internal = (
        coordinate_state(row.get("ITS1_coordinates", "NA")),
        coordinate_state(row.get("5.8S_coordinates", "NA")),
        coordinate_state(row.get("ITS2_coordinates", "NA")),
    )
    if all(state == "PRESENT" for state in internal):
        return "FULL"
    if table_kind == "ITSME_PARTIAL":
        return "PARTIAL"
    return "PARTIAL"


def split_lineage(value: str) -> dict[str, str]:
    ranks = {rank: "NA" for rank in RANKS}
    if value.strip().casefold() in MISSING:
        return ranks
    parts = [part.strip() for part in value.split(";") if part.strip()]
    for rank, taxon in zip(RANKS, parts):
        ranks[rank] = taxon
    return ranks


def normalize_record(sample: str, method: str, table_kind: str,
                     source_file: Path, row: dict[str, str]) -> dict[str, str]:
    lineage = row.get("consensus_taxonomy", "NA") or "NA"
    ranks = split_lineage(lineage)
    normalized = {
        "sample": sample,
        "method": method,
        "sequence_class": sequence_class(row, table_kind),
        "contig": row.get("contig", "NA") or "NA",
        "locus_type": row.get("locus_type", "NA") or "NA",
        "length_bp": row.get("length_bp", "NA") or "NA",
        "mean_depth": (
            row.get("mean_depth") or row.get("SPAdes_mean_depth") or "NA"
        ),
        "regions_present": regions_present(row),
        "taxonomy_status": row.get("taxonomy_status", "NA") or "NA",
        "consensus_taxonomy": lineage,
        "SSU_consensus_taxonomy": row.get("SSU_consensus_taxonomy", "NA") or "NA",
        "ITS_consensus_taxonomy": row.get("ITS_consensus_taxonomy", "NA") or "NA",
        "LSU_consensus_taxonomy": row.get("LSU_consensus_taxonomy", "NA") or "NA",
        "itsx_note": row.get("itsx_note", "NA") or "NA",
        "source_file": source_file.name,
    }
    for field in COORDINATE_FIELDS:
        normalized[field] = row.get(field, "NA") or "NA"
    normalized.update(ranks)
    return normalized


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})


def taxa_at_rank(records: list[dict[str, str]], method: str, rank: str,
                 sequence_class_filter: str = "FULL") -> set[str]:
    return {
        row[rank] for row in records
        if row["method"] == method
        and row["sequence_class"] == sequence_class_filter
        and row.get(rank, "NA").casefold() not in MISSING
    }


def join_taxa(values: set[str]) -> str:
    return " | ".join(sorted(values, key=str.casefold)) if values else "NA"


def comparison_status(itsme: set[str], itsx: set[str]) -> str:
    if not itsme and not itsx:
        return "NO_ASSIGNED_FULL_LOCI"
    if not itsme:
        return "ITSX_ONLY"
    if not itsx:
        return "ITSME_ONLY"
    overlap = itsme & itsx
    if itsme == itsx:
        return "SAME_PHYLA"
    if overlap:
        return "PARTIAL_PHYLA_OVERLAP"
    return "DISCORDANT_PHYLA"


def main() -> None:
    args = parse_args()
    itsme_summary_dir = require_directory(args.itsme_summaries, "ITSME summary")
    itsme_partial_dir = require_directory(args.itsme_partials, "ITSME partial")
    itsx_summary_dir = require_directory(args.itsx_summaries, "ITSx summary")
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    itsme_summary_files = csv_files_by_sample(itsme_summary_dir, "ITSME")
    itsme_partial_files = csv_files_by_sample(itsme_partial_dir, "ITSME")
    itsx_files = csv_files_by_sample(itsx_summary_dir, "ITSX")

    # An ITSME summary file, even if header-only, is the completion marker.
    included_samples = set(itsme_summary_files)
    if not included_samples:
        raise SystemExit(f"ERROR: no ITSME summary CSV files found in {itsme_summary_dir}")

    ignored_itsx = sorted(set(itsx_files) - included_samples)
    if ignored_itsx:
        print(
            f"Ignoring {len(ignored_itsx)} ITSx-only sample(s) not completed by ITSME: "
            + ", ".join(ignored_itsx),
            file=sys.stderr,
        )

    records: list[dict[str, str]] = []
    for sample in sorted(included_samples):
        for path, row in read_csv_rows(itsme_summary_files.get(sample, [])):
            records.append(normalize_record(sample, "ITSME", "ITSME_SUMMARY", path, row))
        for path, row in read_csv_rows(itsme_partial_files.get(sample, [])):
            records.append(normalize_record(sample, "ITSME", "ITSME_PARTIAL", path, row))
        for path, row in read_csv_rows(itsx_files.get(sample, [])):
            records.append(normalize_record(sample, "ITSX", "ITSX_SUMMARY", path, row))

    records.sort(
        key=lambda row: (
            row["sample"], row["method"],
            0 if row["sequence_class"] == "FULL" else 1,
            -(int(row["length_bp"]) if str(row["length_bp"]).isdigit() else 0),
            row["contig"],
        )
    )

    locus_fields = [
        "sample", "method", "sequence_class", "contig", "locus_type",
        "length_bp", "mean_depth", "regions_present", *COORDINATE_FIELDS,
        "taxonomy_status", "consensus_taxonomy", *RANKS,
        "SSU_consensus_taxonomy", "ITS_consensus_taxonomy",
        "LSU_consensus_taxonomy", "itsx_note", "source_file",
    ]
    write_csv(output_dir / "locus_comparison_long.csv", records, locus_fields)

    taxon_counts: dict[tuple[str, str, str], dict[str, int]] = defaultdict(
        lambda: {
            "itsme_full_loci": 0, "itsme_partial_loci": 0,
            "itsx_full_loci": 0, "itsx_partial_loci": 0,
        }
    )
    for row in records:
        for rank in RANKS:
            taxon = row.get(rank, "NA")
            if taxon.casefold() in MISSING:
                continue
            count_field = f"{row['method'].lower()}_{row['sequence_class'].lower()}_loci"
            taxon_counts[(row["sample"], rank, taxon)][count_field] += 1

    rank_order = {rank: index for index, rank in enumerate(RANKS)}
    taxon_rows = []
    for (sample, rank, taxon), counts in sorted(
        taxon_counts.items(),
        key=lambda item: (item[0][0], rank_order[item[0][1]], item[0][2].casefold()),
    ):
        taxon_rows.append({"sample": sample, "rank": rank, "taxon": taxon, **counts})
    taxon_fields = [
        "sample", "rank", "taxon", "itsme_full_loci", "itsme_partial_loci",
        "itsx_full_loci", "itsx_partial_loci",
    ]
    write_csv(output_dir / "taxon_comparison_long.csv", taxon_rows, taxon_fields)

    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in records:
        by_sample[row["sample"]].append(row)
    sample_rows = []
    for sample in sorted(included_samples):
        sample_records = by_sample.get(sample, [])
        itsme_phyla = taxa_at_rank(sample_records, "ITSME", "phylum")
        itsx_phyla = taxa_at_rank(sample_records, "ITSX", "phylum")
        row: dict[str, object] = {
            "sample": sample,
            "itsme_full_loci": sum(
                item["method"] == "ITSME" and item["sequence_class"] == "FULL"
                for item in sample_records
            ),
            "itsme_partial_loci": sum(
                item["method"] == "ITSME" and item["sequence_class"] == "PARTIAL"
                for item in sample_records
            ),
            "itsx_full_loci": sum(
                item["method"] == "ITSX" and item["sequence_class"] == "FULL"
                for item in sample_records
            ),
            "itsx_partial_loci": sum(
                item["method"] == "ITSX" and item["sequence_class"] == "PARTIAL"
                for item in sample_records
            ),
            "itsme_full_phyla": join_taxa(itsme_phyla),
            "itsx_full_phyla": join_taxa(itsx_phyla),
            "shared_full_phyla": join_taxa(itsme_phyla & itsx_phyla),
            "itsme_only_full_phyla": join_taxa(itsme_phyla - itsx_phyla),
            "itsx_only_full_phyla": join_taxa(itsx_phyla - itsme_phyla),
            "comparison_status": comparison_status(itsme_phyla, itsx_phyla),
            "itsx_summary_available": 1 if sample in itsx_files else 0,
        }
        sample_rows.append(row)

    sample_fields = [
        "sample", "itsme_full_loci", "itsme_partial_loci", "itsx_full_loci",
        "itsx_partial_loci", "itsme_full_phyla", "itsx_full_phyla",
        "shared_full_phyla", "itsme_only_full_phyla", "itsx_only_full_phyla",
        "comparison_status", "itsx_summary_available",
    ]
    write_csv(output_dir / "sample_comparison.csv", sample_rows, sample_fields)

    print(f"Included ITSME-completed samples: {len(included_samples)}")
    print(f"Combined locus rows: {len(records)}")
    print(f"Wrote: {output_dir / 'sample_comparison.csv'}")
    print(f"Wrote: {output_dir / 'locus_comparison_long.csv'}")
    print(f"Wrote: {output_dir / 'taxon_comparison_long.csv'}")


if __name__ == "__main__":
    main()
