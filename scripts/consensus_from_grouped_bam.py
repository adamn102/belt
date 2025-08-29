#!/usr/bin/env python3
import argparse, os, re, sys, itertools, subprocess
from pathlib import Path
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
import functools
import tempfile


SAFE = re.compile(r"[^A-Za-z0-9._-]")

def sanitize(x: str) -> str:
    return SAFE.sub("_", x or "NA")

def write_group_fastq(umi, group_iter, fq_path):
    with open(fq_path, "w", buffering=1024*1024) as fq:
        for r in group_iter:
            rn = r.query_name
            seq = r.query_sequence or ""
            qual = r.qual if r.qual is not None else ("I" * len(seq))
            fq.write(f"@{rn}\n{seq}\n+\n{qual}\n")


def run_one(reads_fq, ref_fa, meta_dir, con_dir, bam_dir, umi_threads, model, log_file, umi_safe):
    """Run consensus.sh for one UMI and return the consensus FASTA path."""
    logfh = open(log_file, "a") if log_file else subprocess.DEVNULL
    try:
        subprocess.run(
            [
                "./scripts/consensus.sh",
                str(reads_fq),
                str(ref_fa),
                str(meta_dir),
                str(con_dir),
                str(umi_threads),
                str(model),
                str(bam_dir),  # optional 7th arg handled by consensus.sh
            ],
            stdout=logfh,
            stderr=logfh,
            check=True,
        )
    finally:
        if logfh is not subprocess.DEVNULL:
            logfh.close()
    # remove temp FASTQ
    try:
        os.remove(reads_fq)
    except OSError:
        pass
    # consensus.sh writes <UMI>_consensus.fasta
    return str(Path(con_dir) / f"{umi_safe}_consensus.fasta")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grouped-bam", required=True)
    ap.add_argument("--tag", default="ZU")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--umi-threads", type=int, default=1)
    ap.add_argument("--model", default="r1041_e82_400bps_sup_variant_v5.0.0")
    ap.add_argument("--tmp-root", default=None)
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--log-file", default=None)

    # ✅ parse BEFORE using args
    args = ap.parse_args()

    sample = args.sample
    out_root = Path(args.outdir) / sample
    con_dir = out_root / "consensus"
    bam_dir = out_root / "bam"
    tmp_root = Path(args.tmp_root) if args.tmp_root else (out_root / "tmp")
    con_dir.mkdir(parents=True, exist_ok=True)
    bam_dir.mkdir(parents=True, exist_ok=True)
    tmp_root.mkdir(parents=True, exist_ok=True)
    log_file = args.log_file or str(con_dir / "runner.log")

    produced = []
    futures = []

    def umi_paths(umi_safe: str):
        umi_root = tmp_root / umi_safe
        reads_fq = umi_root / f"{umi_safe}_reads.fastq"  # matches consensus.sh UMI parsing
        meta_dir = umi_root / "meta"                     # safe for cleanup inside consensus.sh
        meta_dir.mkdir(parents=True, exist_ok=True)
        return reads_fq, meta_dir

    with pysam.AlignmentFile(args.grouped_bam, "rb") as bam, \
         ProcessPoolExecutor(max_workers=args.procs) as pool:

        key = lambda r: r.get_tag(args.tag)
        for umi, group_iter in itertools.groupby(bam.fetch(until_eof=True), key=key):
            umi_safe = sanitize(umi)
            reads_fq, meta_dir = umi_paths(umi_safe)
            reads_fq.parent.mkdir(parents=True, exist_ok=True)
            write_group_fastq(umi, group_iter, reads_fq)

            fut = pool.submit(
                run_one,
                str(reads_fq), str(args.reference), str(meta_dir), str(con_dir),
                str(bam_dir), int(args.umi_threads), str(args.model), str(log_file),
                umi_safe
            )
            futures.append(fut)

            # throttle queued tasks
            if len(futures) >= args.procs * 4:
                for done in as_completed(futures[:args.procs]):
                    produced.append(done.result())
                futures = futures[args.procs:]

        for done in as_completed(futures):
            produced.append(done.result())

    if args.manifest:
        with open(args.manifest, "w") as fh:
            for p in produced:
                fh.write(p + "\n")
    for p in produced:
        print(p)



#def run_one(reads_fq, ref_fa, meta_dir, con_dir, umi_threads, model, log_file):
#    # call consensus.sh and delete the reads file afterwards
#    logfh = open(log_file, "a") if log_file else subprocess.DEVNULL
#    try:
#        subprocess.run(
#            [
#                "./scripts/consensus.sh",
#                str(reads_fq),
#                str(ref_fa),
#                str(meta_dir),
#                str(con_dir),
#                str(umi_threads),
#                str(model),
#            ],
#            stdout=logfh, stderr=logfh, check=True,
#        )
#    finally:
#        if logfh is not subprocess.DEVNULL:
#            logfh.close()
#    # cleanup temp fastq
#    try:
#        os.remove(reads_fq)
#    except OSError:
#        pass
#    # return the expected fasta path (consensus.sh names it <UMI>_consensus.fasta)
#    return str(Path(con_dir) / (Path(reads_fq).stem + "_consensus.fasta"))
#
#def main():
#    ap = argparse.ArgumentParser()
#    ap.add_argument("--grouped-bam", required=True)
#    ap.add_argument("--tag", default="ZU")
#    ap.add_argument("--reference", required=True)
#    ap.add_argument("--outdir", required=True)
#    ap.add_argument("--sample", required=True)
#    ap.add_argument("--threads", type=int, default=2)      # total threads (not used directly)
#    ap.add_argument("--procs", type=int, default=1)        # parallel UMIs
#    ap.add_argument("--umi-threads", type=int, default=1)  # threads per UMI for consensus.sh
#    ap.add_argument("--model", default="r1041_e82_400bps_sup_variant_v5.0.0")
#    ap.add_argument("--tmp-root", default=None)
#    ap.add_argument("--manifest", default=None)
#    ap.add_argument("--log-file", default=None)
#    args = ap.parse_args()
#
#    sample = args.sample
#    out_root = Path(args.outdir) / sample
#    con_dir = out_root / "consensus"
#    tmp_root = Path(args.tmp_root) if args.tmp_root else (out_root / "tmp")
#    con_dir.mkdir(parents=True, exist_ok=True)
#    tmp_root.mkdir(parents=True, exist_ok=True)
#    log_file = args.log_file or str(con_dir / "runner.log")
#
#    produced = []
#    futures = []
#
#    def umi_paths(umi_safe: str):
#        umi_root = tmp_root / umi_safe
#        reads_fq = umi_root / f"{umi_safe}.fastq"
#        meta_dir = umi_root / "meta"  # safe: consensus.sh cleans this dir
#        meta_dir.mkdir(parents=True, exist_ok=True)
#        return reads_fq, meta_dir
#
#    with pysam.AlignmentFile(args.grouped_bam, "rb") as bam, \
#         ProcessPoolExecutor(max_workers=args.procs) as pool:
#
#        key = lambda r: r.get_tag(args.tag)
#        for umi, group_iter in itertools.groupby(bam.fetch(until_eof=True), key=key):
#            umi_safe = re.sub(r"[^A-Za-z0-9._-]", "_", umi or "NA")
#            reads_fq, meta_dir = umi_paths(umi_safe)
#
#            # dump this UMI's reads (FASTQ lives OUTSIDE meta_dir)
#            reads_fq.parent.mkdir(parents=True, exist_ok=True)
#            write_group_fastq(umi, group_iter, reads_fq)
#
#            # submit job
#            fut = pool.submit(
#                run_one,
#                str(reads_fq), str(args.reference), str(meta_dir), str(con_dir),
#                int(args.umi_threads), str(args.model), str(log_file)
#            )
#            futures.append(fut)
#
#            # throttle queued tasks so memory doesn't balloon
#            if len(futures) >= args.procs * 4:
#                for done in as_completed(futures[:args.procs]):
#                    produced.append(done.result())
#                futures = futures[args.procs:]
#
#        # drain the rest
#        for done in as_completed(futures):
#            produced.append(done.result())
#
#    if args.manifest:
#        with open(args.manifest, "w") as fh:
#            for p in produced:
#                fh.write(p + "\n")
#    for p in produced:
#        print(p)

if __name__ == "__main__":
    sys.exit(main())
