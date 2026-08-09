"""Adapter from Clover membership and passive core export to ClusterRecord."""
from __future__ import annotations
import ast
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
from .clover_io import read_clover_input
from .models import ClusterRecord,CoreRecord,ReadRecord
class MalformedMembershipError(ValueError): pass
class DuplicateMembershipError(ValueError): pass
class ConflictingMembershipError(ValueError): pass
class CoreSequenceUnavailable(ValueError): pass
class MalformedCoreExportError(ValueError): pass
class DuplicateCoreRecordError(ValueError): pass
class CoreSequenceMismatchError(ValueError): pass
@dataclass(frozen=True)
class AdapterStats:
    total_input_reads:int; assigned_reads:int; unassigned_reads:int; cluster_count:int
@dataclass(frozen=True)
class ClusterBuildResult:
    clusters:list[ClusterRecord]; stats:AdapterStats
def read_membership(path:str|Path)->dict[str,str]:
    path=Path(path); text=path.read_text(encoding="utf-8",errors="replace")
    if not text.strip(): return {}
    try: parsed=ast.literal_eval(text.strip())
    except (SyntaxError,ValueError) as exc: raise MalformedMembershipError(f"{path}: malformed membership literal") from exc
    if not isinstance(parsed,list): raise MalformedMembershipError(f"{path}: membership output must be a list")
    out={}
    for i,item in enumerate(parsed):
        if not isinstance(item,tuple) or len(item)!=2: raise MalformedMembershipError(f"{path}: membership item {i} must be a 2-tuple")
        rid,cid=str(item[0]),str(item[1])
        if not rid: raise MalformedMembershipError(f"{path}: membership item {i} has empty read_id")
        if not cid: raise MalformedMembershipError(f"{path}: membership item {i} has empty cluster_id")
        if rid in out:
            if out[rid]!=cid: raise ConflictingMembershipError(f"read_id {rid} assigned to both {out[rid]} and {cid}")
            raise DuplicateMembershipError(f"duplicate membership for read_id {rid}")
        out[rid]=cid
    return out
def read_core_sequences(path: str | Path) -> dict[str, CoreRecord]:
    path = Path(path)
    cores: dict[str, CoreRecord] = {}

    with path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        expected = ["cluster_id", "core_read_id", "core_sequence"]

        if header != expected:
            raise MalformedCoreExportError(
                f"{path}: expected TSV header {expected}, got {header}"
            )

        for line_number, raw in enumerate(f, start=2):
            line = raw.rstrip("\r\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != 3:
                raise MalformedCoreExportError(
                    f"{path}:{line_number}: expected 3 tab-delimited fields"
                )

            try:
                record = CoreRecord(
                    cluster_id=fields[0],
                    core_read_id=fields[1],
                    sequence=fields[2],
                )
            except ValueError as exc:
                raise MalformedCoreExportError(
                    f"{path}:{line_number}: {exc}"
                ) from exc

            if record.cluster_id in cores:
                raise DuplicateCoreRecordError(
                    f"{path}:{line_number}: duplicate core record "
                    f"for cluster {record.cluster_id}"
                )

            cores[record.cluster_id] = record

    return cores


def _normalize(
    core_sequences: Mapping[str, CoreRecord | str] | None,
) -> dict[str, CoreRecord]:
    if core_sequences is None:
        return {}

    normalized: dict[str, CoreRecord] = {}

    for cluster_id, value in core_sequences.items():
        if isinstance(value, CoreRecord):
            normalized[cluster_id] = value
        else:
            normalized[cluster_id] = CoreRecord(
                cluster_id=cluster_id,
                core_read_id=cluster_id,
                sequence=value,
            )

    return normalized


def build_cluster_records(input_path:str|Path,membership_path:str|Path,core_sequences:Mapping[str,CoreRecord|str]|None=None)->ClusterBuildResult:
    reads=read_clover_input(input_path); membership=read_membership(membership_path); cores=_normalize(core_sequences)
    grouped:dict[str,list[ReadRecord]]=defaultdict(list)
    for rid,cid in membership.items():
        if rid not in reads: raise ValueError(f"membership references unknown read_id: {rid}")
        grouped[cid].append(reads[rid])
    missing=sorted((set(grouped)|set(cores))-set(grouped).intersection(cores))
    if missing: raise CoreSequenceUnavailable("CORE_SEQUENCE_UNAVAILABLE: missing reliable core_sequence or membership for cluster(s): "+", ".join(missing))
    clusters=[]; assigned_ids=set(membership)
    for cid,members in sorted(grouped.items()):
        core=cores[cid]
        if core.core_read_id not in reads: raise ValueError(f"core_read_id not found in original input: {core.core_read_id}")
        core_read=reads[core.core_read_id]
        if core_read.sequence!=core.sequence: raise CoreSequenceMismatchError(f"core sequence mismatch for {cid}: export does not match original input")
        if core.core_read_id not in assigned_ids: members=[core_read,*members]
        clusters.append(ClusterRecord(cid,core.sequence,members,core.core_read_id))
    stats=AdapterStats(len(reads),len(membership),len(reads)-len(membership),len(clusters))
    return ClusterBuildResult(clusters,stats)
