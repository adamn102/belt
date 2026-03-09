# Feasibility: Writing UMI Data Directly into the BAM File

## Summary

**Yes, it is feasible.** The pipeline can be extended so that corrected UMIs are written into the BAM as auxiliary tags. `umi_group` can then read UMIs from the BAM instead of (or in addition to) the CSV, removing the CSV as a hard dependency and improving performance for large files.

---

## Current Pipeline Flow

```
┌─────────────────┐     ┌──────────────────────────────────────────────────┐     ┌─────────────────┐
│  Aligned BAM    │────▶│  umi_map (umi_extraction_map)                     │────▶│  CSV            │
│  (minimap2 etc) │     │  • findUMI(BAM) → raw UMIs                        │     │  Paired_UMIs_   │
└─────────────────┘     │  • whitelist/cluster correction                  │     │  read_names.csv │
                         │  • Writes CSV only                               │     └────────┬────────┘
                         └──────────────────────────────────────────────────┘              │
                                                                                            │
┌─────────────────┐     ┌──────────────────────────────────────────────────┐               │
│  Same BAM       │────▶│  umi_group (long_read_umi_grouping_fast)          │◀──────────────┘
└─────────────────┘     │  • Reads CSV (Read Name, UMI1, UMI2, …)          │
                         │  • Filters by threshold / UMI2 majority         │
                         │  • Tags kept reads with ZU, sorts → grouped.bam  │
                         └──────────────────────────────────────────────────┘
                                                        │
                                                        ▼
                         ┌──────────────────────────────────────────────────┐
                         │  consensus (process_grouped_bam)                  │
                         │  • Reads grouped BAM, uses ZU tag only           │
                         └──────────────────────────────────────────────────┘
```

- **umi_map**: Reads BAM, computes corrected UMIs, writes **only** CSV (`*_Paired_UMIs_read_names.csv`).
- **umi_group**: Reads the **same** BAM and the CSV; uses CSV for all UMI values and filtering; writes a new BAM with **ZU** tag and sorted by ZU.
- **consensus**: Uses only the **grouped BAM** and the **ZU** tag; no CSV.

So the only consumer of the UMI CSV is `umi_group`. The rest of the pipeline (consensus) only needs the BAM with ZU.

---

## Proposed Flow with BAM Tags

```
┌─────────────────┐     ┌──────────────────────────────────────────────────┐     ┌─────────────────────┐
│  Aligned BAM    │────▶│  umi_map                                           │────▶│  Tagged BAM          │
└─────────────────┘     │  • findUMI(BAM) → raw UMIs, correct                 │     │  (U1, U2, … tags)   │
                         │  • Write BAM with UMI tags (U1, U2, …)             │     │  + optional CSV      │
                         └──────────────────────────────────────────────────┘     └──────────┬──────────┘
                                                                                              │
                         ┌──────────────────────────────────────────────────┐               │
                         │  umi_group                                       │◀──────────────┘
                         │  • If BAM has UMI tags: build table from tags    │
                         │  • Else: load CSV (backward compatible)          │
                         │  • Filter, tag ZU, sort → grouped.bam            │
                         └──────────────────────────────────────────────────┘
```

- **umi_map**: Still does extraction + correction. **Additionally** (or optionally) writes a **tagged BAM**: same alignments plus one tag per UMI (e.g. U1=UMI1, U2=UMI2, …).
- **umi_group**: Prefer **BAM as source of UMI data** when tags are present; otherwise keep using CSV. No change to output (grouped BAM with ZU) or to consensus.

---

## BAM Tag Conventions (SAM spec)

- Auxiliary tags: `TAG:TYPE:VALUE`, with TAG = two characters matching `[A-Za-z][A-Za-z0-9]`.
- Suggested:
  - **U1, U2, U3, …** (letter `U` + digit): corrected UMI1, UMI2, UMI3, …  
    - One tag per UMI column; missing UMI = omit tag or use empty string (to be defined).
  - **ZU**: used by `umi_group` for the **group** UMI in the **output** BAM (already the case).  
    - Input BAM from umi_map would **not** set ZU; umi_group sets ZU after filtering.

So:
- **umi_map output BAM**: tags **U1, U2, …** (corrected UMIs only).
- **umi_group output BAM**: **ZU** (group UMI) as today; consensus continues to use ZU only.

---

## What Needs to Change

### 1. `umi_map.py` (umi_extraction_map)

- After building `output_df` (Read Name, UMI1, UMI2, …):
  - Open the **alignment BAM** and a **new BAM** (e.g. `{sampleName}_umi_tagged.bam`) with same header.
  - Build a **single** lookup: `read_name → (UMI1, UMI2, …)` from `output_df` (e.g. dict or pandas index).
  - For each read in the alignment BAM:
    - Look up read name in the dict.
    - Set tags **U1**, **U2**, … from the corrected values (handle missing: omit tag or set to `""`).
    - Write the read to the new BAM.
  - Close BAMs; optionally keep writing the CSV for debugging / compatibility.
- Add a **parameter** (e.g. `write_tagged_bam=True`) and **output path** (e.g. `outputDir/{sampleName}_umi_tagged.bam`) so the pipeline can request the tagged BAM.

### 2. `umi_group.py` (long_read_umi_grouping_fast)

- **`correctedUMIs`** becomes **optional** when reading from BAM is enabled.
- **Logic**:
  - If **BAM has UMI tags** (e.g. at least one of U1, U2 present) **and** we are allowed to use them (e.g. `corrected_umis` is None or a flag `use_bam_tags=True`):
    - **Single pass** over the BAM: for each read, read tags U1, U2, … and build an in-memory structure equivalent to the current `df`: `read_name → {UMI1, UMI2, …}`.
    - Run **existing** filtering (threshold, UMI2 majority, 1000-read cap) on this structure.
    - **Second pass**: same BAM (or the same tagged BAM), write only kept reads with **ZU** set, then sort by ZU.  
  - Else:
    - Keep current behavior: load CSV, build `df`, then same filtering and BAM write as today.
- **Grouping column**: When reading from BAM, “UMI1” maps to tag **U1**, “UMI2” to **U2**, etc., so grouping and UMI2-majority logic stay the same conceptually.
- **FASTQ / “all UMIs”**: Today built from `df` columns; when using BAM, build from the same in-memory map (UMI1=U1, UMI2=U2, …).

### 3. `belt.py` (CLI and run_pipeline)

- **umi-map**:
  - Add option, e.g. `--write-tagged-bam` and/or `--tagged-bam PATH`, and pass through to `umi_extraction_map`.
- **umi-group**:
  - Make `--corrected_umis` optional when `--alignment` is the tagged BAM and BAM has UMI tags.
  - Add something like `--use-bam-tags` to force “read UMIs from BAM only”.
- **run_pipeline (config)**:
  - Add config keys, e.g. `write_tagged_bam: true`, `tagged_bam_path: ...`, `use_bam_tags: true`.
  - When `use_bam_tags` is true: umi_map writes tagged BAM; umi_group is called with `alignment=<tagged_bam>`, `corrected_umis=None` (or omitted), and uses BAM tags.

---

## Backward Compatibility

- **Keep CSV path**: If `corrected_umis` is provided and `use_bam_tags` is false (or not set), `umi_group` behaves exactly as today (CSV only).
- **Optional tagged BAM**: If `write_tagged_bam` is false or not set, umi_map only writes CSV; pipeline stays as today.
- **Consensus**: No change; it only sees the grouped BAM and ZU.

So “write UMI data into the BAM” is an **add-on**: same pipeline as today, with an extra option to carry UMIs in the BAM and have umi_group use them instead of the CSV.

---

## Benefits

1. **Single file**: UMI data lives in the BAM; no separate CSV for the main path.
2. **Performance**: umi_group can build the read→UMI map in one streaming pass over the BAM (no CSV load, no repeated DataFrame lookups by read name).
3. **Scalability**: Avoids loading millions of rows into pandas; memory can be “one read at a time” plus the final keep-set and count structures.
4. **Consistency**: Read order and identity come from the BAM; no risk of BAM/CSV ordering or name mismatches.

---

## Edge Cases / Design Choices

| Topic | Recommendation |
|-------|----------------|
| **Missing UMI** | Omit the tag for that UMI (e.g. no U2) or set to `""`; umi_group treats missing tag as NA, same as empty CSV cell. |
| **Many motifs (e.g. >9)** | Use U1–U9, then two-letter tags (e.g. Ua, Ub, …) or a single delimited tag (e.g. `UX:Z:val1\|val2\|...`) and parse in umi_group. |
| **Raw vs corrected** | Only corrected UMIs need to be in the BAM for grouping. Raw can stay in CSV if needed for analysis. |
| **Config `corrected_umis`** | When using BAM tags, config can omit `corrected_umis` or set it to the tagged BAM path for logging; umi_group detects tags and does not load CSV. |

---

## Conclusion

Writing UMI data directly into the BAM is **feasible** and fits the current pipeline:

- **umi_map** already has the corrected UMI table and the alignment BAM; adding a second output that is the same BAM with U1, U2, … tags is straightforward.
- **umi_group** can be extended to build the equivalent of the current `df` from BAM tags in one pass, then run the same filtering and ZU-tagging logic; CSV remains supported for backward compatibility.
- **consensus** already uses only the grouped BAM and ZU; no change.

Recommended implementation order: (1) add tagged BAM output to umi_map, (2) add “read from BAM tags” path in umi_group with fallback to CSV, (3) wire CLI and config so the pipeline can run in “BAM-only” mode.
