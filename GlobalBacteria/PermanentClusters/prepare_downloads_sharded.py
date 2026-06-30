#!/usr/bin/env python3
"""
Streaming preparation of sharded FASTA.gz download files.

Designed for very large inputs:
  variants:  hash, samples, abundances, taxonomy_id, sequence
  taxonomy:  id, cluster, species, genus, ...

The program never loads all taxonomy IDs or variants into RAM.

Pipeline:
  1. Stream variants once:
       - compute per-sample totals
       - partition annotated rows by taxonomy-ID range
  2. Stream taxonomy in ascending ID order, one range at a time:
       - externally sort only the matching variant partition by taxonomy ID
       - write cluster files directly
       - spool species/genus records into 256 hash buckets
  3. Externally sort each species/genus spool by target filename and
     write every final gzip file exactly once.

Requires GNU sort.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import BinaryIO, Iterator


BUFFER_SIZE = 4 * 1024 * 1024
FASTA_SUFFIX = b".fasta.gz"
_CREATED_OUTPUT_DIRS: set[Path] = set()


def open_binary_input(path: Path) -> BinaryIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rb")
    return path.open("rb", buffering=BUFFER_SIZE)


def ensure_empty_directory(path: Path, label: str) -> None:
    if path.exists():
        if not path.is_dir():
            raise RuntimeError(f"{label} is not a directory: {path}")
        if any(path.iterdir()):
            raise RuntimeError(f"{label} must be empty: {path}")
    else:
        path.mkdir(parents=True)


def validate_filename(filename: bytes) -> None:
    for invalid in (b"/", b"\\", b"\x00", b"\t", b"\r", b"\n"):
        if invalid in filename:
            raise ValueError(
                f"Unsafe character in generated filename: {filename!r}"
            )


def output_kind(filename: bytes) -> str:
    if filename.startswith(b"GB_"):
        return "GB"
    if filename.startswith(b"species_"):
        return "species"
    if filename.startswith(b"genus_"):
        return "genus"
    return "other"


def sharded_output_path(output_root: Path, filename: bytes) -> Path:
    """
    Deterministic path derived only from the exact final filename.

    Example:
        species_Foo_bar.fasta.gz
        -> output/species/ab/cd/species_Foo_bar.fasta.gz
    """
    validate_filename(filename)
    digest = hashlib.sha256(filename).hexdigest()
    decoded = filename.decode("utf-8")
    return output_root / output_kind(filename) / digest[:2] / digest[2:4] / decoded


def group_bucket(filename: bytes) -> int:
    # First SHA-256 byte: 0..255.
    return hashlib.sha256(filename).digest()[0]


class BucketWriters:
    def __init__(
        self,
        directory: Path,
        prefix: str,
        max_open: int,
    ) -> None:
        self.directory = directory
        self.prefix = prefix
        self.max_open = max_open
        self.directory.mkdir(parents=True, exist_ok=True)
        self.handles: OrderedDict[int, BinaryIO] = OrderedDict()

    def path_for(self, bucket: int) -> Path:
        return self.directory / f"{self.prefix}_{bucket:06d}.tsv"

    def write(self, bucket: int, data: bytes) -> None:
        fp = self.handles.pop(bucket, None)

        if fp is None:
            if len(self.handles) >= self.max_open:
                _, old_fp = self.handles.popitem(last=False)
                old_fp.close()

            fp = self.path_for(bucket).open("ab", buffering=BUFFER_SIZE)

        self.handles[bucket] = fp
        fp.write(data)

    def close(self) -> None:
        while self.handles:
            _, fp = self.handles.popitem(last=False)
            fp.close()


def parse_variant_line(
    line: bytes,
    line_number: int,
) -> tuple[bytes, bytes, bytes, bytes, bytes]:
    fields = line.rstrip(b"\r\n").split(b"\t", 4)
    if len(fields) != 5:
        raise ValueError(
            f"Variants line {line_number}: expected 5 tab-separated fields, "
            f"found {len(fields)}"
        )
    return fields[0], fields[1], fields[2], fields[3], fields[4]


def sample_pairs(
    samples_field: bytes,
    abundances_field: bytes,
    context: str,
) -> Iterator[tuple[bytes, bytes]]:
    samples = samples_field.split(b";")
    abundances = abundances_field.split(b";")

    if len(samples) != len(abundances):
        raise ValueError(
            f"{context}: {len(samples)} sample IDs but "
            f"{len(abundances)} abundances"
        )

    return zip(samples, abundances)


def partition_variants(
    variants_path: Path,
    variants_dir: Path,
    range_size: int,
    max_open_partitions: int,
    progress_every: int,
) -> tuple[dict[bytes, int], list[int], int, int]:
    totals: dict[bytes, int] = defaultdict(int)
    writers = BucketWriters(
        variants_dir,
        prefix="variants",
        max_open=max_open_partitions,
    )
    used_buckets: set[int] = set()
    all_rows = 0
    annotated_rows = 0

    try:
        with open_binary_input(variants_path) as src:
            for line_number, line in enumerate(src, 1):
                if not line.strip():
                    continue

                _, samples, abundances, tax_id_raw, _ = parse_variant_line(
                    line, line_number
                )

                for sample_id, abundance in sample_pairs(
                    samples, abundances, f"variants line {line_number}"
                ):
                    totals[sample_id] += int(abundance)

                all_rows += 1

                if tax_id_raw != b"-":
                    tax_id = int(tax_id_raw)
                    if tax_id < 1:
                        raise ValueError(
                            f"Variants line {line_number}: invalid taxonomy ID "
                            f"{tax_id}"
                        )

                    bucket = (tax_id - 1) // range_size
                    writers.write(bucket, line)
                    used_buckets.add(bucket)
                    annotated_rows += 1

                if progress_every and all_rows % progress_every == 0:
                    print(
                        f"Partitioned {all_rows:,} variant rows; "
                        f"{annotated_rows:,} annotated",
                        flush=True,
                    )
    finally:
        writers.close()

    return totals, sorted(used_buckets), all_rows, annotated_rows


def save_sample_totals(path: Path, totals: dict[bytes, int]) -> None:
    with path.open("wb", buffering=BUFFER_SIZE) as fp:
        for sample_id in sorted(totals, key=lambda value: int(value)):
            fp.write(sample_id + b"\t" + str(totals[sample_id]).encode() + b"\n")


def taxonomy_rows(
    taxonomy_path: Path,
) -> Iterator[tuple[int, bytes, bytes, bytes]]:
    previous_id = 0

    with open_binary_input(taxonomy_path) as src:
        for line_number, line in enumerate(src, 1):
            if not line.strip():
                continue

            fields = line.rstrip(b"\r\n").split(b"\t", 4)
            if len(fields) < 4:
                raise ValueError(
                    f"Taxonomy line {line_number}: expected at least "
                    f"4 tab-separated fields"
                )

            tax_id = int(fields[0])
            if tax_id <= previous_id:
                raise ValueError(
                    "Taxonomy table must be strictly sorted by numeric ID; "
                    f"line {line_number} has {tax_id} after {previous_id}"
                )
            previous_id = tax_id

            cluster = fields[1]
            species = fields[2].replace(b" ", b"_")
            genus = fields[3].replace(b" ", b"_")

            yield tax_id, cluster, species, genus


def load_taxonomy_range(
    tax_iter: Iterator[tuple[int, bytes, bytes, bytes]],
    pending: tuple[int, bytes, bytes, bytes] | None,
    start_id: int,
    end_id: int,
    range_size: int,
) -> tuple[
    list[tuple[bytes, bytes, bytes] | None],
    tuple[int, bytes, bytes, bytes] | None,
]:
    while pending is not None and pending[0] < start_id:
        pending = next(tax_iter, None)

    mapping: list[tuple[bytes, bytes, bytes] | None] = [None] * range_size

    while pending is not None and pending[0] <= end_id:
        tax_id, cluster, species, genus = pending
        mapping[tax_id - start_id] = (cluster, species, genus)
        pending = next(tax_iter, None)

    return mapping, pending


def sort_process(
    input_path: Path,
    key: str,
    sort_memory: str,
    sort_parallel: int,
    sort_tmp: Path,
) -> subprocess.Popen[bytes]:
    env = os.environ.copy()
    env["LC_ALL"] = "C"

    command = [
        "sort",
        "-t",
        "\t",
        key,
        "-S",
        sort_memory,
        "--parallel",
        str(sort_parallel),
        "-T",
        str(sort_tmp),
        str(input_path),
    ]

    return subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        env=env,
        bufsize=BUFFER_SIZE,
    )


def open_gzip_output(
    output_root: Path,
    filename: bytes,
    gzip_level: int,
) -> BinaryIO:
    path = sharded_output_path(output_root, filename)
    parent = path.parent
    if parent not in _CREATED_OUTPUT_DIRS:
        parent.mkdir(parents=True, exist_ok=True)
        _CREATED_OUTPUT_DIRS.add(parent)
    return gzip.open(path, "wb", compresslevel=gzip_level)


def fasta_record(
    variant_hash: bytes,
    sample_id: bytes,
    taxonomy_label: bytes,
    abundance: bytes,
    total: int,
    sequence: bytes,
) -> bytes:
    return (
        b">"
        + variant_hash
        + b"|SampleID_"
        + sample_id
        + b"|"
        + taxonomy_label
        + b"|marker_16S|abund_"
        + abundance
        + b"_total_"
        + str(total).encode()
        + b"\n"
        + sequence
        + b"\n"
    )


def group_spool_record(
    filename: bytes,
    variant_hash: bytes,
    samples: bytes,
    abundances: bytes,
    sequence: bytes,
) -> bytes:
    validate_filename(filename)
    return (
        filename
        + b"\t"
        + variant_hash
        + b"\t"
        + samples
        + b"\t"
        + abundances
        + b"\t"
        + sequence
        + b"\n"
    )


def process_variant_partition(
    partition_path: Path,
    bucket: int,
    range_size: int,
    taxonomy: list[tuple[bytes, bytes, bytes] | None],
    sample_totals: dict[bytes, int],
    output_root: Path,
    group_writers: BucketWriters,
    sort_memory: str,
    sort_parallel: int,
    sort_tmp: Path,
    gzip_level: int,
) -> tuple[int, int]:
    start_id = bucket * range_size + 1
    end_id = start_id + range_size - 1

    process = sort_process(
        partition_path,
        key="-k4,4n",
        sort_memory=sort_memory,
        sort_parallel=sort_parallel,
        sort_tmp=sort_tmp,
    )
    assert process.stdout is not None

    current_tax_id: int | None = None
    cluster_fp: BinaryIO | None = None
    cluster_name: bytes | None = None
    rows = 0
    fasta_records = 0

    try:
        for line_number, line in enumerate(process.stdout, 1):
            (
                variant_hash,
                samples,
                abundances,
                tax_id_raw,
                sequence,
            ) = parse_variant_line(line, line_number)

            tax_id = int(tax_id_raw)
            if not start_id <= tax_id <= end_id:
                raise RuntimeError(
                    f"Taxonomy ID {tax_id} is outside partition "
                    f"{start_id}..{end_id}"
                )

            tax = taxonomy[tax_id - start_id]
            if tax is None:
                raise KeyError(
                    f"Taxonomy ID {tax_id} from variants was not found "
                    f"in the taxonomy table"
                )

            row_cluster, species, genus = tax

            if tax_id != current_tax_id:
                if cluster_fp is not None:
                    cluster_fp.close()

                cluster_name = row_cluster
                cluster_filename = b"GB_" + cluster_name + FASTA_SUFFIX
                cluster_fp = open_gzip_output(
                    output_root, cluster_filename, gzip_level
                )
                current_tax_id = tax_id

            assert cluster_fp is not None
            assert cluster_name is not None

            pairs = list(
                sample_pairs(
                    samples,
                    abundances,
                    f"partition {bucket}, sorted line {line_number}",
                )
            )

            for sample_id, abundance in pairs:
                try:
                    total = sample_totals[sample_id]
                except KeyError as exc:
                    raise KeyError(
                        f"Missing total for sample {sample_id!r}"
                    ) from exc

                cluster_fp.write(
                    fasta_record(
                        variant_hash,
                        sample_id,
                        cluster_name,
                        abundance,
                        total,
                        sequence,
                    )
                )
                fasta_records += 1

            if species != b"-" and b"_sp." not in species:
                filename = b"species_" + species + FASTA_SUFFIX
                group_writers.write(
                    group_bucket(filename),
                    group_spool_record(
                        filename,
                        variant_hash,
                        samples,
                        abundances,
                        sequence,
                    ),
                )

            if genus != b"-":
                filename = b"genus_" + genus + FASTA_SUFFIX
                group_writers.write(
                    group_bucket(filename),
                    group_spool_record(
                        filename,
                        variant_hash,
                        samples,
                        abundances,
                        sequence,
                    ),
                )

            rows += 1
    finally:
        if cluster_fp is not None:
            cluster_fp.close()
        process.stdout.close()

    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(
            f"GNU sort failed for {partition_path} with status {return_code}"
        )

    return rows, fasta_records


def finalize_group_bucket(
    spool_path: Path,
    sample_totals: dict[bytes, int],
    output_root: Path,
    sort_memory: str,
    sort_parallel: int,
    sort_tmp: Path,
    gzip_level: int,
) -> tuple[int, int]:
    process = sort_process(
        spool_path,
        key="-k1,1",
        sort_memory=sort_memory,
        sort_parallel=sort_parallel,
        sort_tmp=sort_tmp,
    )
    assert process.stdout is not None

    current_filename: bytes | None = None
    current_label: bytes | None = None
    fp: BinaryIO | None = None
    rows = 0
    records = 0

    try:
        for line_number, line in enumerate(process.stdout, 1):
            fields = line.rstrip(b"\r\n").split(b"\t", 4)
            if len(fields) != 5:
                raise ValueError(
                    f"Group spool {spool_path}, sorted line {line_number}: "
                    f"expected 5 fields"
                )

            filename, variant_hash, samples, abundances, sequence = fields

            if filename != current_filename:
                if fp is not None:
                    fp.close()

                if not filename.endswith(FASTA_SUFFIX):
                    raise ValueError(f"Unexpected filename: {filename!r}")

                current_filename = filename
                current_label = filename[: -len(FASTA_SUFFIX)]
                fp = open_gzip_output(output_root, filename, gzip_level)

            assert fp is not None
            assert current_label is not None

            for sample_id, abundance in sample_pairs(
                samples,
                abundances,
                f"group spool {spool_path}, sorted line {line_number}",
            ):
                fp.write(
                    fasta_record(
                        variant_hash,
                        sample_id,
                        current_label,
                        abundance,
                        sample_totals[sample_id],
                        sequence,
                    )
                )
                records += 1

            rows += 1
    finally:
        if fp is not None:
            fp.close()
        process.stdout.close()

    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(
            f"GNU sort failed for {spool_path} with status {return_code}"
        )

    return rows, records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare sharded cluster/species/genus FASTA.gz files"
    )
    parser.add_argument("variants", type=Path)
    parser.add_argument("taxonomy", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--range-size", type=int, default=500_000)
    parser.add_argument("--group-buckets", type=int, default=256)
    parser.add_argument("--max-open-partitions", type=int, default=400)
    parser.add_argument("--sort-memory", default="8G")
    parser.add_argument("--sort-parallel", type=int, default=4)
    parser.add_argument("--gzip-level", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=5_000_000)
    parser.add_argument("--keep-work", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.range_size < 1:
        raise ValueError("--range-size must be positive")
    if args.group_buckets != 256:
        raise ValueError(
            "This implementation uses exactly 256 group buckets "
            "(one per first SHA-256 byte)"
        )
    if not 0 <= args.gzip_level <= 9:
        raise ValueError("--gzip-level must be between 0 and 9")

    if shutil.which("sort") is None:
        raise RuntimeError("GNU sort was not found in PATH")

    ensure_empty_directory(args.output_dir, "Output directory")
    ensure_empty_directory(args.work_dir, "Work directory")

    variants_dir = args.work_dir / "variant_partitions"
    groups_dir = args.work_dir / "group_spools"
    sort_tmp = args.work_dir / "sort_tmp"
    sort_tmp.mkdir(parents=True)

    print("Stage 1/3: partition variants and compute sample totals", flush=True)

    (
        sample_totals,
        used_buckets,
        all_rows,
        annotated_rows,
    ) = partition_variants(
        args.variants,
        variants_dir,
        args.range_size,
        args.max_open_partitions,
        args.progress_every,
    )

    save_sample_totals(args.work_dir / "sample_totals.tsv", sample_totals)

    print(
        f"Variant rows: {all_rows:,}; annotated: {annotated_rows:,}; "
        f"samples: {len(sample_totals):,}; range partitions: "
        f"{len(used_buckets):,}",
        flush=True,
    )

    print("Stage 2/3: range join and cluster output", flush=True)

    group_writers = BucketWriters(
        groups_dir,
        prefix="groups",
        max_open=args.group_buckets,
    )

    tax_iter = taxonomy_rows(args.taxonomy)
    pending = next(tax_iter, None)
    processed_rows = 0

    try:
        for position, bucket in enumerate(used_buckets, 1):
            start_id = bucket * args.range_size + 1
            end_id = start_id + args.range_size - 1

            taxonomy, pending = load_taxonomy_range(
                tax_iter,
                pending,
                start_id,
                end_id,
                args.range_size,
            )

            partition_path = variants_dir / f"variants_{bucket:06d}.tsv"

            rows, records = process_variant_partition(
                partition_path,
                bucket,
                args.range_size,
                taxonomy,
                sample_totals,
                args.output_dir,
                group_writers,
                args.sort_memory,
                args.sort_parallel,
                sort_tmp,
                args.gzip_level,
            )
            processed_rows += rows

            if not args.keep_work:
                partition_path.unlink()

            print(
                f"Range {position:,}/{len(used_buckets):,} "
                f"({start_id:,}..{end_id:,}): "
                f"{rows:,} variants, {records:,} cluster FASTA records",
                flush=True,
            )
    finally:
        group_writers.close()

    if processed_rows != annotated_rows:
        raise RuntimeError(
            f"Processed {processed_rows:,} annotated rows, "
            f"expected {annotated_rows:,}"
        )

    print("Stage 3/3: species/genus grouping and output", flush=True)

    total_group_rows = 0
    total_group_records = 0

    for bucket in range(args.group_buckets):
        spool_path = groups_dir / f"groups_{bucket:06d}.tsv"
        if not spool_path.exists() or spool_path.stat().st_size == 0:
            continue

        rows, records = finalize_group_bucket(
            spool_path,
            sample_totals,
            args.output_dir,
            args.sort_memory,
            args.sort_parallel,
            sort_tmp,
            args.gzip_level,
        )
        total_group_rows += rows
        total_group_records += records

        if not args.keep_work:
            spool_path.unlink()

        print(
            f"Group bucket {bucket:03d}: "
            f"{rows:,} variant-target rows, {records:,} FASTA records",
            flush=True,
        )

    print(
        f"Done. Group rows: {total_group_rows:,}; "
        f"group FASTA records: {total_group_records:,}",
        flush=True,
    )

    if not args.keep_work:
        shutil.rmtree(args.work_dir)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(1)
