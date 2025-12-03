#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mutation_caller.py

Per-UMI mutation calling against a reference with robust handling of 'N' reference bases.
- Skips sites where the reference base is 'N' (does not call SNPs there; does not count them in the denominator)
- Skips deletions whose deleted span contains any 'N'
- (Conservative) Skips insertions anchored at reference-N sites
- Computes per-UMI callable length and mutation rate over callable sites only

Inputs:
  - Reference FASTA (indexed or not; we read it fully)
  - Consensus FASTA (multi-FASTA of per-UMI consensus sequences)
Outputs:
  - Summary CSV with per-UMI metrics
  - Optional per-mutation TSV

Dependencies:
  - Standard library only (alignment fallback included)
  - Optional speedups:
      * parasail (fast alignment)    pip install parasail
      * biopython (pairwise2 align)  pip install biopython
"""

from __future__ import annotations
import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# ---------- Optional accelerators ----------
_HAS_PARASAIL = False
_HAS_BIO = False
try:
    import parasail  # type: ignore
    _HAS_PARASAIL = True
except Exception:
    pass

try:
    from Bio import SeqIO  # type: ignore
    from Bio import pairwise2  # type: ignore
    _HAS_BIO = True
except Exception:
    pass

# ---------- I/O ----------

def read_fasta(path: str | Path) -> Dict[str, str]:
    """
    Read a FASTA (single or multi) -> dict: {header_without_spaces: UPPERCASE_SEQ}
    Falls back to a tiny parser if Biopython isn't installed.
    """
    p = Path(path)
    if _HAS_BIO:
        d: Dict[str, str] = {}
        for rec in SeqIO.parse(str(p), "fasta"):
            name = str(rec.id).split()[0]
            d[name] = str(rec.seq).upper()
        if not d:
            raise ValueError(f"No sequences found in {p}")
        return d

    # Minimal FASTA parser
    d: Dict[str, str] = {}
    name: Optional[str] = None
    buf: List[str] = []
    with open(p, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf).upper()
                    buf = []
                name = line[1:].split()[0]
            else:
                buf.append(line)
        if name is not None:
            d[name] = "".join(buf).upper()
    if not d:
        raise ValueError(f"No sequences found in {p}")
    return d

# ---------- Alignment ----------

@dataclass
class Alignment:
    ref_aln: str
    qry_aln: str
    score: float

def _align_parasail(ref: str, qry: str) -> Alignment:
    # Global alignment with simple scoring
    # match +2, mismatch -3, gap open -5, gap extend -1 (tunable)
    match = 2
    mismatch = -3
    open_pen = 5
    extend_pen = 1
    matrix = parasail.matrix_create("ACGTN", match, mismatch)
    result = parasail.nw_trace_scan_16(qry, ref, open_pen, extend_pen, matrix)
    t = parasail.traceback(result, qry, ref, matrix)
    # parasail returns aligned query and reference
    return Alignment(ref_aln=t.ref, qry_aln=t.query, score=result.score)

def _align_biopython(ref: str, qry: str) -> Alignment:
    # Needleman-Wunsch global alignment (pairwise2)
    # Scores tuned for long ONT consensus vs ref; adjust if needed
    # match +2, mismatch -3, gap open -5, gap extend -1
    aln = pairwise2.align.globalms(ref, qry, 2, -3, -5, -1, one_alignment_only=True)
    if not aln:
        return Alignment(ref_aln="", qry_aln="", score=float("-inf"))
    a = aln[0]
    return Alignment(ref_aln=a.seqA, qry_aln=a.seqB, score=a.score)

def _align_builtin(ref: str, qry: str) -> Alignment:
    # Simple Needleman-Wunsch fallback; O(n*m) but fine for amplicons
    match = 2
    mismatch = -3
    go = -5
    ge = -1

    n = len(ref)
    m = len(qry)
    # DP matrices
    # M = match/mismatch; Ix = gap in ref (insert in qry); Iy = gap in qry (del in qry)
    INF = -10**9
    M = [[INF]*(m+1) for _ in range(n+1)]
    Ix = [[INF]*(m+1) for _ in range(n+1)]
    Iy = [[INF]*(m+1) for _ in range(n+1)]
    TB = [[0]*(m+1) for _ in range(n+1)]  # 0:M,1:Ix,2:Iy

    M[0][0] = 0
    for i in range(1, n+1):
        Iy[i][0] = go + (i-1)*ge
        M[i][0] = Iy[i][0]
        TB[i][0] = 2
    for j in range(1, m+1):
        Ix[0][j] = go + (j-1)*ge
        M[0][j] = Ix[0][j]
        TB[0][j] = 1

    def s(a: str, b: str) -> int:
        if a == b:
            return match
        return mismatch

    for i in range(1, n+1):
        ri = ref[i-1]
        for j in range(1, m+1):
            qj = qry[j-1]
            # Extend/open gaps
            Ix[i][j] = max(M[i][j-1] + go, Ix[i][j-1] + ge)
            Iy[i][j] = max(M[i-1][j] + go, Iy[i-1][j] + ge)
            # Match/Mismatch
            M[i][j] = max(M[i-1][j-1] + s(ri, qj), Ix[i-1][j-1] + s(ri, qj), Iy[i-1][j-1] + s(ri, qj))
            # track simple TB on best of M/Ix/Iy
            if M[i][j] >= Ix[i][j] and M[i][j] >= Iy[i][j]:
                TB[i][j] = 0
            elif Ix[i][j] >= Iy[i][j]:
                TB[i][j] = 1
            else:
                TB[i][j] = 2

    # Traceback (coarse; follow largest among M/Ix/Iy)
    i, j = n, m
    ref_aln = []
    qry_aln = []
    cur = 0  # assume M at end
    while i > 0 or j > 0:
        # choose the best layer at (i,j)
        best = M[i][j]; cur = 0
        if Ix[i][j] > best:
            best = Ix[i][j]; cur = 1
        if Iy[i][j] > best:
            best = Iy[i][j]; cur = 2

        if cur == 0:
            # match/mismatch
            ref_aln.append(ref[i-1] if i > 0 else "-")
            qry_aln.append(qry[j-1] if j > 0 else "-")
            i -= 1; j -= 1
        elif cur == 1:
            # gap in ref: insert in query
            ref_aln.append("-")
            qry_aln.append(qry[j-1] if j > 0 else "-")
            j -= 1
        else:
            # gap in query: deletion in query (relative to ref)
            ref_aln.append(ref[i-1] if i > 0 else "-")
            qry_aln.append("-")
            i -= 1

    ref_aln.reverse()
    qry_aln.reverse()
    return Alignment(ref_aln="".join(ref_aln), qry_aln="".join(qry_aln), score=best)

def align_sequences(ref: str, qry: str) -> Alignment:
    """
    Try parasail -> Biopython -> builtin in that order.
    """
    if _HAS_PARASAIL:
        try:
            return _align_parasail(ref, qry)
        except Exception:
            pass
    if _HAS_BIO:
        try:
            return _align_biopython(ref, qry)
        except Exception:
            pass
    return _align_builtin(ref, qry)

# ---------- Mutation calling with ref-N skipping ----------

@dataclass
class Mutation:
    type: str            # "substitution" | "insertion" | "deletion"
    ref_position: int    # 0-based coordinate in the unaligned reference
    ref_base: str
    cons_base: str
    size: int = 1
    extra: Optional[str] = None

def _find_mutations_skip_ref_N(
    ref_aln: str,
    cons_aln: str,
    ref_raw: str,
    skip_insertions_at_refN: bool = True
) -> List[Mutation]:
    """
    Identify differences between aligned ref/cons sequences while:
      - skipping any position where the original reference base is 'N'/'n'
      - skipping deletions if any deleted reference base is 'N'
      - (optionally) skipping insertions anchored at 'N' (conservative)

    Coordinates reported as 0-based positions on the *raw* reference.
    """
    muts: List[Mutation] = []
    i = 0                # index in aligned strings
    ref_pos = 0          # 0-based position in raw reference

    def ref_is_N(pos: int) -> bool:
        return 0 <= pos < len(ref_raw) and (ref_raw[pos] == "N" or ref_raw[pos] == "n")
    
    # Pre-compute all continuous N regions in the reference
    # This helps us efficiently check if insertions are within UMI regions
    n_regions: List[Tuple[int, int]] = []  # List of (start, end) positions for each N region
    in_n_region = False
    n_start = 0
    for idx, base in enumerate(ref_raw):
        if ref_is_N(idx):
            if not in_n_region:
                in_n_region = True
                n_start = idx
        else:
            if in_n_region:
                in_n_region = False
                n_regions.append((n_start, idx - 1))
    if in_n_region:
        n_regions.append((n_start, len(ref_raw) - 1))
    
    def is_in_N_region(pos: int) -> bool:
        """Check if position is within or near any N region (UMI regions)"""
        if pos < 0 or pos >= len(ref_raw):
            return False
        # Check current position
        if ref_is_N(pos):
            return True
        # Check if position is within any N region (with adaptive tolerance for boundaries)
        # Insertions can be anchored at positions near but not exactly at N regions
        # due to alignment artifacts. For each N region, we use a window that's
        # the size of the region plus 5 bases on each side
        for n_start, n_end in n_regions:
            n_length = n_end - n_start + 1
            # Use a window that's the N region size plus 5 bases on each side
            # This ensures we catch insertions near UMI regions regardless of their size
            window_size = max(15, n_length + 10)  # At least 15 bases, or N region size + 10
            if n_start - window_size <= pos <= n_end + window_size:
                return True
        return False

    L = len(ref_aln)
    while i < L:
        r = ref_aln[i]
        q = cons_aln[i]

        # INSERTION relative to reference (gap in reference)
        if r == "-" and q != "-":
            # insertion is anchored *between* ref_pos-1 and ref_pos
            # Check both the anchor position and the previous position (where we just were)
            anchor = ref_pos  # left-anchored: consider anchor at current ref_pos
            prev_anchor = ref_pos - 1 if ref_pos > 0 else -1
            
            # Skip insertion if it's anchored at or near an N position
            # Check both anchor positions because insertions can be anchored between positions
            if skip_insertions_at_refN:
                skip_insertion = False
                # Check current anchor position
                if is_in_N_region(anchor):
                    skip_insertion = True
                # Check previous anchor position (insertion might be anchored between positions)
                elif prev_anchor >= 0 and is_in_N_region(prev_anchor):
                    skip_insertion = True
                
                if skip_insertion:
                    # consume insertion but do not record
                    while i < L and ref_aln[i] == "-" and cons_aln[i] != "-":
                        i += 1
                    continue

            ins_bases: List[str] = []
            while i < L and ref_aln[i] == "-" and cons_aln[i] != "-":
                ins_bases.append(cons_aln[i])
                i += 1
            if ins_bases:
                seq = "".join(ins_bases)
                muts.append(Mutation(
                    type="insertion",
                    ref_position=anchor,
                    ref_base="-",
                    cons_base=seq,
                    size=len(seq),
                    extra=seq
                ))
            continue

        # DELETION relative to reference (gap in consensus)
        if r != "-" and q == "-":
            del_start = ref_pos
            del_bases: List[str] = []
            contains_N = False
            while i < L and ref_aln[i] != "-" and cons_aln[i] == "-":
                del_bases.append(ref_aln[i])
                if ref_is_N(ref_pos):
                    contains_N = True
                ref_pos += 1
                i += 1
            if contains_N:
                # skip deletions spanning unknown reference
                continue
            seq = "".join(del_bases)
            muts.append(Mutation(
                type="deletion",
                ref_position=del_start,
                ref_base=seq,
                cons_base="-",
                size=len(seq),
                extra=seq
            ))
            continue

        # SUBSTITUTION
        if r != "-" and q != "-" and r != q:
            if not ref_is_N(ref_pos):
                muts.append(Mutation(
                    type="substitution",
                    ref_position=ref_pos,
                    ref_base=r,
                    cons_base=q,
                    size=1
                ))
            ref_pos += 1
            i += 1
            continue

        # MATCH (counts toward callable if ref is not N and qry not gap)
        if r != "-" and q != "-":
            ref_pos += 1

        i += 1

    return muts

def _callable_sites(ref_aln: str, cons_aln: str, ref_raw: str) -> int:
    """
    Count positions that are valid for calling:
      - reference base in A/C/G/T (not gap, not N) in the raw reference
      - consensus not a gap
    """
    callable_ = 0
    ref_pos = 0  # 0-based position in raw reference
    for r, q in zip(ref_aln, cons_aln):
        # Only count if: aligned ref is not a gap, raw ref at this position is ACGT (not N), and consensus is not a gap
        if r != "-" and q != "-":
            if 0 <= ref_pos < len(ref_raw) and ref_raw[ref_pos] in "ACGTacgt":
                callable_ += 1
            ref_pos += 1
        elif r != "-":
            # Reference base (could be N) but consensus is gap - don't count
            ref_pos += 1
    return callable_

# ---------- Main runner ----------

def run(
    reference_fa: str | Path,
    consensus_fa: str | Path,
    out_csv: str | Path,
    per_mutation_tsv: Optional[str | Path] = None,
    ref_contig: Optional[str] = None,
    skip_insertions_at_refN: bool = True
) -> None:
    # Load reference(s)
    ref_dict = read_fasta(reference_fa)
    if not ref_contig:
        # default to the first contig
        ref_name = next(iter(ref_dict.keys()))
    else:
        if ref_contig not in ref_dict:
            raise ValueError(f"Requested contig {ref_contig!r} not in reference FASTA.")
        ref_name = ref_contig
    ref_seq = ref_dict[ref_name]

    # Load consensus sequences (multi-FASTA)
    cons_dict = read_fasta(consensus_fa)

    # Prepare outputs
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    per_mutation_path = Path(per_mutation_tsv) if per_mutation_tsv else None
    if per_mutation_path:
        per_mutation_path.parent.mkdir(parents=True, exist_ok=True)

    # Write headers
    with open(out_csv, "w", newline="") as csvfh:
        writer = csv.writer(csvfh, delimiter=",")
        writer.writerow([
            "consensus_id",
            "ref_contig",
            "consensus_length",
            "callable_bases",
            "mutation_count",
            "mutation_rate",
            "snp_count",
            "ins_count",
            "del_count",
            "alignment_score"
        ])

    if per_mutation_path:
        with open(per_mutation_path, "w", newline="") as tsvfh:
            w = csv.writer(tsvfh, delimiter="\t")
            w.writerow(["consensus_id","type","ref_position","ref_base","cons_base","size","extra"])

    # Process each consensus
    for cons_id, cons_seq in cons_dict.items():
        # Align
        aln = align_sequences(ref_seq, cons_seq)
        if not aln.ref_aln or not aln.qry_aln:
            # write a row with zeros
            with open(out_csv, "a", newline="") as csvfh:
                writer = csv.writer(csvfh, delimiter=",")
                writer.writerow([cons_id, ref_name, len(cons_seq), 0, 0, 0.0, 0, 0, 0, float("-inf")])
            continue

        # Call mutations with ref-N skipping
        muts = _find_mutations_skip_ref_N(
            aln.ref_aln, aln.qry_aln, ref_seq, skip_insertions_at_refN=skip_insertions_at_refN
        )

        # Callable denominator (excludes N positions in reference)
        callable_len = _callable_sites(aln.ref_aln, aln.qry_aln, ref_seq)
        denom = callable_len if callable_len > 0 else len(cons_seq)
        mut_rate = (len(muts) / denom) if denom > 0 else 0.0

        # Type counts
        snp = sum(1 for m in muts if m.type == "substitution")
        ins = sum(1 for m in muts if m.type == "insertion")
        dele = sum(1 for m in muts if m.type == "deletion")

        # Write summary row
        with open(out_csv, "a", newline="") as csvfh:
            writer = csv.writer(csvfh, delimiter=",")
            writer.writerow([
                cons_id, ref_name, len(cons_seq),
                callable_len, len(muts), f"{mut_rate:.6f}",
                snp, ins, dele, f"{aln.score:.2f}"
            ])

        # Optional: per-mutation TSV
        if per_mutation_path:
            with open(per_mutation_path, "a", newline="") as tsvfh:
                w = csv.writer(tsvfh, delimiter="\t")
                for m in muts:
                    w.writerow([
                        cons_id, m.type, m.ref_position, m.ref_base, m.cons_base, m.size, (m.extra or "")
                    ])

# ---------- Wrapper for belt CLI ----------

def call_mutations_for_sample(
    consensus_file: str | Path,
    reference_file: str | Path,
    sample_name: str,
    output_dir: str | Path
) -> None:
    """
    Wrapper function for belt mutation-call CLI command.
    
    This function generates:
    - {sample_name}_mutations.json: Summary JSON with mutation statistics
    - samples/{sample_name}/{sample_name}_per_read.csv: Per-UMI mutation summary CSV
    
    Args:
        consensus_file: Path to compiled consensus FASTA file
        reference_file: Path to reference FASTA file
        sample_name: Name of the sample
        output_dir: Output directory for results
    """
    import json
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output paths
    per_read_csv = output_dir / "samples" / sample_name / f"{sample_name}_per_read.csv"
    json_output = output_dir / f"{sample_name}_mutations.json"
    
    # Run mutation calling
    run(
        reference_fa=reference_file,
        consensus_fa=consensus_file,
        out_csv=per_read_csv,
        per_mutation_tsv=None,  # Don't generate per-mutation TSV by default
        ref_contig=None,
        skip_insertions_at_refN=True
    )
    
    # Read the CSV and create a JSON summary
    mutations_data = []
    total_mutations = 0
    total_callable = 0
    total_snps = 0
    total_ins = 0
    total_dels = 0
    
    if per_read_csv.exists():
        with open(per_read_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mutations_data.append({
                    "consensus_id": row["consensus_id"],
                    "consensus_length": int(row["consensus_length"]),
                    "callable_bases": int(row["callable_bases"]),
                    "mutation_count": int(row["mutation_count"]),
                    "mutation_rate": float(row["mutation_rate"]),
                    "snp_count": int(row["snp_count"]),
                    "ins_count": int(row["ins_count"]),
                    "del_count": int(row["del_count"]),
                    "alignment_score": float(row["alignment_score"])
                })
                total_mutations += int(row["mutation_count"])
                total_callable += int(row["callable_bases"])
                total_snps += int(row["snp_count"])
                total_ins += int(row["ins_count"])
                total_dels += int(row["del_count"])
    
    # Create JSON summary
    summary = {
        "sample_name": sample_name,
        "reference_file": str(reference_file),
        "consensus_file": str(consensus_file),
        "total_umis": len(mutations_data),
        "total_mutations": total_mutations,
        "total_callable_bases": total_callable,
        "overall_mutation_rate": (total_mutations / total_callable) if total_callable > 0 else 0.0,
        "total_snps": total_snps,
        "total_insertions": total_ins,
        "total_deletions": total_dels,
        "per_umi_mutations": mutations_data
    }
    
    # Write JSON output
    with open(json_output, 'w') as f:
        json.dump(summary, f, indent=2)

# ---------- CLI ----------

def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Per-UMI mutation caller with reference-N masking (no changes needed to consensus step)."
    )
    ap.add_argument("--reference", "-r", required=True, help="Reference FASTA")
    ap.add_argument("--consensus", "-c", required=True, help="Consensus multi-FASTA (per-UMI sequences)")
    ap.add_argument("--out-csv", "-o", required=True, help="Summary CSV path")
    ap.add_argument("--per-mutation-tsv", "-m", default=None, help="(Optional) path for per-mutation TSV")
    ap.add_argument("--contig", default=None, help="Reference contig name to use (default: first in FASTA)")
    ap.add_argument("--allow-insertions-at-refN", action="store_true",
                    help="By default we skip insertions anchored at reference-N sites. Set this to ALLOW them.")
    args = ap.parse_args(argv)

    try:
        run(
            reference_fa=args.reference,
            consensus_fa=args.consensus,
            out_csv=args.out_csv,
            per_mutation_tsv=args.per_mutation_tsv,
            ref_contig=args.contig,
            skip_insertions_at_refN=(not args.allow_insertions_at_refN),
        )
        return 0
    except Exception as e:
        sys.stderr.write(f"[mutation_caller] ERROR: {e}\n")
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
