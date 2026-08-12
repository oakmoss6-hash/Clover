"""Input normalization for Clover CLI."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator


def iter_clover_records(path: str | Path) -> Iterator[str]:
    """
    Yield Clover-native ``read_id sequence`` records.

    Supported inputs:
        FASTQ / FQ
        FASTA / FA
        Clover text: read_id sequence
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".fastq", ".fq"}:
        with path.open(
            "r",
            encoding="utf-8",
            errors="replace",
        ) as source:
            while True:
                header = source.readline()

                if not header:
                    return

                sequence = source.readline()
                plus = source.readline()
                quality = source.readline()

                if (
                    not sequence
                    or not plus
                    or not quality
                ):
                    raise ValueError(
                        f"truncated FASTQ record in {path}"
                    )

                if not header.startswith("@"):
                    raise ValueError(
                        f"invalid FASTQ header: "
                        f"{header[:80]!r}"
                    )

                read_id = (
                    header[1:]
                    .strip()
                    .split()[0]
                )

                sequence = (
                    sequence.strip().upper()
                )

                yield f"{read_id} {sequence}"

        return

    if suffix in {".fasta", ".fa", ".fna"}:
        with path.open(
            "r",
            encoding="utf-8",
            errors="replace",
        ) as source:
            read_id = None
            parts = []

            for raw in source:
                line = raw.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if read_id is not None:
                        yield (
                            f"{read_id} "
                            f"{''.join(parts).upper()}"
                        )

                    read_id = (
                        line[1:].split()[0]
                    )
                    parts = []
                else:
                    if read_id is None:
                        raise ValueError(
                            "FASTA sequence before header"
                        )

                    parts.append(line)

            if read_id is not None:
                yield (
                    f"{read_id} "
                    f"{''.join(parts).upper()}"
                )

        return

    # Historical Clover text format.
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as source:
        for lineno, raw in enumerate(source, 1):
            line = raw.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected "
                    "'read_id sequence'"
                )

            yield (
                f"{fields[0]} "
                f"{fields[1].upper()}"
            )
