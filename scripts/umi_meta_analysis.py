#!/usr/bin/env python3
"""
BELT UMI Meta-Analysis Module
Analyzes UMI patterns and statistics from UMI mapping data
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional
import logging
import pysam
from collections import Counter
import re
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def analyze_umi_patterns(umi_file: Path) -> Dict:
    """
    Analyze UMI patterns for a single sample
    
    Args:
        umi_file: Path to UMI mapping CSV file
    
    Returns:
        Dictionary containing UMI analysis results
    """
    logger.info(f"Analyzing UMI patterns for: {umi_file}")
    
    try:
        # Read the UMI mapping file
        df = pd.read_csv(umi_file)
        
        # Extract UMI sequences
        umi_sequences = []
        if 'UMI1' in df.columns:
            umi_sequences.extend(df['UMI1'].dropna().tolist())
        if 'UMI2' in df.columns:
            umi_sequences.extend(df['UMI2'].dropna().tolist())
        
        if not umi_sequences:
            logger.warning(f"No UMI sequences found in {umi_file}")
            return {
                'sample_name': umi_file.stem.replace('_Paired_UMIs_read_names', ''),
                'total_umis': 0,
                'unique_umis': 0,
                'umi_length_stats': {},
                'umi_composition': {},
                'umi_diversity': 0,
                'read_counts': {}
            }
        
        # Basic statistics
        total_umis = len(umi_sequences)
        unique_umis = len(set(umi_sequences))
        
        # UMI length statistics
        umi_lengths = [len(umi) for umi in umi_sequences]
        length_stats = {
            'mean': np.mean(umi_lengths),
            'std': np.std(umi_lengths),
            'min': min(umi_lengths),
            'max': max(umi_lengths),
            'median': np.median(umi_lengths)
        }
        
        # UMI composition analysis
        all_bases = ''.join(umi_sequences)
        base_counts = Counter(all_bases)
        total_bases = len(all_bases)
        composition = {base: count/total_bases for base, count in base_counts.items()}
        
        # UMI diversity (Shannon entropy)
        umi_counts = Counter(umi_sequences)
        total_unique = len(umi_counts)
        if total_unique > 1:
            proportions = [count/total_umis for count in umi_counts.values()]
            diversity = -sum(p * np.log2(p) for p in proportions if p > 0)
        else:
            diversity = 0
        
        # Read count analysis
        read_counts = {}
        if 'read_count' in df.columns:
            read_counts = {
                'mean': df['read_count'].mean(),
                'std': df['read_count'].std(),
                'min': df['read_count'].min(),
                'max': df['read_count'].max(),
                'median': df['read_count'].median()
            }
        
        return {
            'sample_name': umi_file.stem.replace('_Paired_UMIs_read_names', ''),
            'total_umis': total_umis,
            'unique_umis': unique_umis,
            'umi_length_stats': length_stats,
            'umi_composition': composition,
            'umi_diversity': diversity,
            'read_counts': read_counts
        }
        
    except Exception as e:
        logger.error(f"Error analyzing {umi_file}: {e}")
        return {
            'sample_name': umi_file.stem.replace('_Paired_UMIs_read_names', ''),
            'total_umis': 0,
            'unique_umis': 0,
            'umi_length_stats': {},
            'umi_composition': {},
            'umi_diversity': 0,
            'read_counts': {},
            'error': str(e)
        }


def umi_meta_analysis(umi_dir: str, output_dir: str, samples: Optional[List[str]] = None):
    """
    Perform meta-analysis on UMI mapping data
    
    Args:
        umi_dir: Directory containing UMI mapping CSV files
        output_dir: Output directory for results
        samples: Optional list of sample names
    """
    logger.info("Starting UMI meta-analysis")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find UMI files
    umi_path = Path(umi_dir)
    umi_files = list(umi_path.glob("*_Paired_UMIs_read_names.csv"))
    
    if not umi_files:
        raise ValueError(f"No UMI mapping files found in {umi_dir}")
    
    # Filter by samples if provided
    if samples:
        umi_files = [f for f in umi_files if any(sample in f.name for sample in samples)]
    
    logger.info(f"Found {len(umi_files)} UMI files to analyze")
    
    # Analyze each sample
    sample_results = []
    umi1_umi2_summaries = []
    for umi_file in umi_files:
        result = analyze_umi_patterns(umi_file)
        sample_results.append(result)
        # New analysis: UMI1→UMI2 structure per sample
        try:
            per_group_df, per_sample_summary = analyze_umi1_to_umi2_structure(umi_file)
            # Write per-sample detailed table
            per_sample_dir = output_path / "umi1_umi2"
            per_sample_dir.mkdir(parents=True, exist_ok=True)
            out_csv = per_sample_dir / f"{per_sample_summary['sample_name']}_UMI1_to_UMI2.csv"
            per_group_df.to_csv(out_csv, index=False)
            umi1_umi2_summaries.append(per_sample_summary)
        except Exception as e:
            logger.warning(f"UMI1→UMI2 analysis failed for {umi_file}: {e}")

    # Generate summary
    generate_umi_summary(sample_results, output_path)

    # Generate pipeline-level metrics by parsing logs
    try:
        metrics_df = collect_pipeline_metrics(umi_dir=umi_dir, sample_results=sample_results)
        # Merge compact UMI1→UMI2 summaries if present
        if umi1_umi2_summaries:
            extra_df = pd.DataFrame(umi1_umi2_summaries)
            metrics_df = metrics_df.merge(
                extra_df,
                left_on='Sample',
                right_on='sample_name',
                how='left'
            ).drop(columns=['sample_name'])
        metrics_csv = output_path / "umi_pipeline_metrics.csv"
        metrics_df.to_csv(metrics_csv, index=False)
        logger.info(f"UMI pipeline metrics written to {metrics_csv}")
    except Exception as e:
        logger.warning(f"Failed to collect pipeline metrics: {e}")
    
    logger.info(f"UMI meta-analysis complete. Results saved to {output_dir}")


def generate_umi_summary(sample_results: List[Dict], output_path: Path):
    """
    Generate summary files from UMI analysis results
    
    Args:
        sample_results: List of UMI analysis results
        output_path: Path to output directory
    """
    logger.info("Generating UMI summary")
    
    # Create summary DataFrame
    summary_data = []
    for result in sample_results:
        summary_data.append({
            'Sample': result['sample_name'],
            'Total_UMIs': result['total_umis'],
            'Unique_UMIs': result['unique_umis'],
            'UMI_Diversity': result['umi_diversity'],
            'Mean_Length': result['umi_length_stats'].get('mean', 0),
            'Length_Std': result['umi_length_stats'].get('std', 0),
            'Min_Length': result['umi_length_stats'].get('min', 0),
            'Max_Length': result['umi_length_stats'].get('max', 0),
            'A_Composition': result['umi_composition'].get('A', 0),
            'T_Composition': result['umi_composition'].get('T', 0),
            'G_Composition': result['umi_composition'].get('G', 0),
            'C_Composition': result['umi_composition'].get('C', 0),
            'N_Composition': result['umi_composition'].get('N', 0)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_path / "umi_summary.csv", index=False)
    
    # Create detailed report
    report = {
        "analysis_metadata": {
            "total_samples": len(sample_results),
            "total_umis": sum(r['total_umis'] for r in sample_results),
            "total_unique_umis": sum(r['unique_umis'] for r in sample_results),
            "overall_diversity": np.mean([r['umi_diversity'] for r in sample_results if r['umi_diversity'] > 0])
        },
        "results": sample_results
    }
    
    with open(output_path / "umi_analysis_report.json", 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info("UMI summary files generated successfully")


def collect_pipeline_metrics(umi_dir: str, sample_results: List[Dict]) -> pd.DataFrame:
    """
    Aggregate per-sample metrics from pipeline logs into a single DataFrame.

    Expected locations (relative to pipeline root):
      - logs/1-fastq_filter/json/{sample}.json
      - logs/2-mapping_logs/{sample}.log
      - logs/5-umi_map/json/{sample}_umi_extraction_map.log.json
      - logs/6-umi_group/{sample}_umi_group.log.json

    The pipeline root is inferred from umi_dir (usually data/5-umi_maps): root = Path(umi_dir).parents[1]
    """
    umi_path = Path(umi_dir).resolve()
    # infer pipeline root (…/pipelines/ERRSeq2)
    try:
        pipeline_root = umi_path.parents[1]
    except Exception:
        pipeline_root = umi_path

    logs_root = pipeline_root / "logs"
    fastq_json_dir = logs_root / "1-fastq_filter" / "json"
    mapping_logs_dir = logs_root / "2-mapping_logs"
    umi_map_json_dir = logs_root / "5-umi_map" / "json"
    umi_group_json_dir = logs_root / "6-umi_group"

    # index results by sample name for convenience
    samples = [r['sample_name'] for r in sample_results]
    rows = []

    for sample in samples:
        row = {
            'Sample': sample,
            # Step 1: fastq filter
            'S1_PreFilter_Reads': None,
            'S1_PostFilter_Reads': None,
            'S1_Filter_Pass_Pct': None,
            # Step 2: mapping
            'S2_Mapped_Reads': None,
            'S2_Map_Pct_of_PostFilter': None,
            # Step 5: umi-map
            'S5_Total_Reads': None,
            'S5_AllUMIs_Reads': None,
            # Step 6: umi-group
            'S6_Total_Groups': None,
            'S6_Groups_Passing': None,
            'S6_Fraction_Passing': None,
        }

        # 1) fastq filter JSON
        fq_json = fastq_json_dir / f"{sample}.json"
        if fq_json.exists():
            try:
                with open(fq_json) as f:
                    fj = json.load(f)
                # fastp-like schema
                before = fj.get('summary', {}).get('before_filtering', {}).get('total_reads')
                after = fj.get('summary', {}).get('after_filtering', {}).get('total_reads')
                # fallback keys
                if before is None and 'before_filtering' in fj:
                    before = fj['before_filtering'].get('total_reads')
                if after is None and 'after_filtering' in fj:
                    after = fj['after_filtering'].get('total_reads')
                row['S1_PreFilter_Reads'] = int(before) if before is not None else None
                row['S1_PostFilter_Reads'] = int(after) if after is not None else None
                if row['S1_PreFilter_Reads'] is not None and row['S1_PostFilter_Reads'] is not None:
                    row['S1_Filter_Pass_Pct'] = round((row['S1_PostFilter_Reads'] / row['S1_PreFilter_Reads']), 6) if row['S1_PreFilter_Reads'] > 0 else 0.0
            except Exception:
                pass

        # 2) mapping logs
        map_log = mapping_logs_dir / f"{sample}.log"
        if map_log.exists():
            try:
                txt = map_log.read_text()
                # try several regex patterns
                patterns = [
                    r"mapped\s*[:=]\s*(\d+)",
                    r"mapped reads\s*[:=]\s*(\d+)",
                    r"primary mapped\s*[:=]\s*(\d+)",
                    r"reads mapped\s*[:=]\s*(\d+)"
                ]
                mapped = None
                for pat in patterns:
                    m = re.search(pat, txt, flags=re.IGNORECASE)
                    if m:
                        mapped = int(m.group(1))
                        break
                if mapped is not None:
                    row['S2_Mapped_Reads'] = mapped
                    if row['S1_PostFilter_Reads'] is not None:
                        row['S2_Map_Pct_of_PostFilter'] = round((mapped / row['S1_PostFilter_Reads']), 6) if row['S1_PostFilter_Reads'] > 0 else 0.0
            except Exception:
                pass

        # Fallback: derive mapped reads directly from sorted BAM when logs are ambiguous
        if row['S2_Mapped_Reads'] is None:
            bam_path = pipeline_root / "data" / "4-sorted_reads" / f"{sample}.bam"
            if bam_path.exists():
                try:
                    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                        mapped = 0
                        for r in bam.fetch(until_eof=True):
                            if r.is_unmapped:
                                continue
                            if r.is_secondary or r.is_supplementary:
                                continue
                            mapped += 1
                    row['S2_Mapped_Reads'] = mapped
                    if row['S1_PostFilter_Reads'] is not None:
                        row['S2_Map_Pct_of_PostFilter'] = round((mapped / row['S1_PostFilter_Reads']), 6) if row['S1_PostFilter_Reads'] > 0 else 0.0
                except Exception:
                    pass

        # 3) umi-map JSON
        umi_map_json = umi_map_json_dir / f"{sample}_umi_extraction_map.log.json"
        if umi_map_json.exists():
            try:
                with open(umi_map_json) as f:
                    umj = json.load(f)
                total_reads = umj.get('total_reads')
                all_umis_reads = umj.get('reads_with_all_umis')
                reads_with_only = umj.get('reads_with_only', {}) or {}
                row['S5_Total_Reads'] = int(total_reads) if total_reads is not None else None
                row['S5_AllUMIs_Reads'] = int(all_umis_reads) if all_umis_reads is not None else None
                # add per-UMI only fractions
                denom = row['S5_Total_Reads'] if row['S5_Total_Reads'] else 0
                for umi_key, count in reads_with_only.items():
                    try:
                        cnt = int(count)
                    except Exception:
                        continue
                    frac = round((cnt / denom), 6) if denom > 0 else 0.0
                    row[f"S5_Only_{umi_key}_Reads"] = cnt
                    row[f"S5_Only_{umi_key}_Frac"] = frac
                # add whitelist mapping fractions (overall and per-column) if available
                wms = umj.get('whitelist_mapping_stats') or {}
                overall = wms.get('overall') or {}
                overall_frac = overall.get('fraction_mapped')
                if isinstance(overall_frac, (int, float)):
                    row['S5_UmiWhitelist_FracMapped_Overall'] = float(round(overall_frac, 6))
                per_col = wms.get('per_column') or {}
                for col_name, stats in per_col.items():
                    frac_col = stats.get('fraction_mapped')
                    if isinstance(frac_col, (int, float)):
                        safe_col = str(col_name).replace(' ', '_')
                        row[f'S5_UmiWhitelist_FracMapped_{safe_col}'] = float(round(frac_col, 6))
            except Exception:
                pass

        # 4) umi-group JSON
        umi_group_json = umi_group_json_dir / f"{sample}_umi_group.log.json"
        if umi_group_json.exists():
            try:
                with open(umi_group_json) as f:
                    ug = json.load(f)
                total_groups = ug.get('total_unique_umis')
                passing_groups = ug.get('groups_passing_threshold')
                frac = ug.get('fraction_passing_threshold')
                row['S6_Total_Groups'] = int(total_groups) if total_groups is not None else None
                row['S6_Groups_Passing'] = int(passing_groups) if passing_groups is not None else None
                row['S6_Fraction_Passing'] = float(round(frac, 6)) if isinstance(frac, (int, float)) else (
                    round((row['S6_Groups_Passing'] / row['S6_Total_Groups']), 6) if (row.get('S6_Total_Groups') not in (None, 0) and row.get('S6_Groups_Passing') is not None) else None
                )

                # New: UMI2-majority filter metrics (if applied)
                if ug.get('umi2_majority_filter_applied'):
                    gp2 = ug.get('groups_passing_umi2_majority')
                    fp2 = ug.get('fraction_passing_umi2_majority')
                    row['S6_Groups_Passing_UMI2Maj'] = int(gp2) if gp2 is not None else None
                    # Prefer logged fraction; otherwise compute relative to threshold-passing groups
                    if isinstance(fp2, (int, float)):
                        row['S6_Fraction_Passing_UMI2Maj'] = float(round(fp2, 6))
                    elif row.get('S6_Groups_Passing') not in (None, 0) and row.get('S6_Groups_Passing_UMI2Maj') is not None:
                        row['S6_Fraction_Passing_UMI2Maj'] = round(
                            row['S6_Groups_Passing_UMI2Maj'] / row['S6_Groups_Passing'], 6
                        )
            except Exception:
                pass

        rows.append(row)

    # Normalize into DataFrame
    df = pd.DataFrame(rows)
    # Order columns for readability
    ordered_cols = [
        'Sample',
        'S1_PreFilter_Reads','S1_PostFilter_Reads','S1_Filter_Pass_Pct',
        'S2_Mapped_Reads','S2_Map_Pct_of_PostFilter',
        'S5_Total_Reads','S5_AllUMIs_Reads'
    ]
    # Add dynamic S5 Only_* columns in sorted order
    s5_only_cols = sorted([c for c in df.columns if c.startswith('S5_Only_')])
    ordered_cols.extend(s5_only_cols)
    # Add whitelist mapping fraction columns
    whitelist_cols = [c for c in df.columns if c.startswith('S5_UmiWhitelist_FracMapped_')]
    # Ensure Overall first if present
    overall_first = [c for c in whitelist_cols if c.endswith('_Overall')]
    percol_rest = sorted([c for c in whitelist_cols if not c.endswith('_Overall')])
    ordered_cols.extend(overall_first + percol_rest)
    ordered_cols.extend(['S6_Total_Groups','S6_Groups_Passing','S6_Fraction_Passing','S6_Groups_Passing_UMI2Maj','S6_Fraction_Passing_UMI2Maj'])
    # Keep any extras at the end
    ordered_cols += [c for c in df.columns if c not in ordered_cols]
    df = df.reindex(columns=ordered_cols)
    return df


def analyze_umi1_to_umi2_structure(umi_file: Path) -> (pd.DataFrame, Dict):
    """
    For each UMI1 group, count distinct UMI2 species, compute entropy of UMI2 counts,
    and identify the majority UMI2. Returns a per-group DataFrame and a compact per-sample
    summary with aggregate stats.
    """
    df = pd.read_csv(umi_file)
    sample_name = umi_file.stem.replace('_Paired_UMIs_read_names', '')

    # Ensure expected columns
    if 'UMI1' not in df.columns:
        raise ValueError('UMI1 column not found')
    if 'UMI2' not in df.columns:
        # If only a single UMI present, synthesize UMI2 as NA to keep schema consistent
        df['UMI2'] = np.nan

    # Work only on rows with UMI1 present
    df_valid = df[df['UMI1'].notna()].copy()

    # Group by UMI1 and analyze UMI2 distributions
    records = []
    distinct_counts = []
    entropies = []
    for umi1, grp in df_valid.groupby('UMI1', dropna=False):
        umi2_vals = grp['UMI2'].dropna()
        umi1_count = int(len(grp))
        if umi2_vals.empty:
            num_distinct = 0
            entropy = 0.0
            majority_umi2 = None
            majority_count = 0
            majority_fraction = 0.0
        else:
            counts = Counter(umi2_vals)
            num_distinct = len(counts)
            total = sum(counts.values())
            probs = [c / total for c in counts.values()]
            # Shannon entropy (base e)
            entropy = float(-sum(p * np.log(p) for p in probs if p > 0))
            majority_umi2, majority_count = max(counts.items(), key=lambda kv: (kv[1], kv[0]))
            majority_fraction = float(majority_count / total) if total > 0 else 0.0

        records.append({
            'Sample': sample_name,
            'UMI1': umi1,
            'UMI1_Count': umi1_count,
            'UMI2_Distinct_Count': int(num_distinct),
            'UMI2_Entropy': float(entropy),
            'UMI2_Majority': majority_umi2,
            'UMI2_Majority_Count': int(majority_count),
            'UMI2_Majority_Fraction': float(majority_fraction)
        })
        distinct_counts.append(num_distinct)
        entropies.append(entropy)

    per_group_df = pd.DataFrame.from_records(records)

    # Compact per-sample summary (aggregates across UMI1 groups)
    if records:
        summary = {
            'sample_name': sample_name,
            'U12_Mean_UMI2_Distinct_Per_UMI1': float(np.mean(distinct_counts)) if distinct_counts else 0.0,
            'U12_Median_UMI2_Distinct_Per_UMI1': float(np.median(distinct_counts)) if distinct_counts else 0.0,
            'U12_Mean_Entropy_Per_UMI1': float(np.mean(entropies)) if entropies else 0.0,
            'U12_Median_Entropy_Per_UMI1': float(np.median(entropies)) if entropies else 0.0,
            'U12_Total_UMI1_Groups': int(len(per_group_df)) if not per_group_df.empty else 0,
            'U12_Mean_UMI1_Count': float(per_group_df['UMI1_Count'].mean()) if not per_group_df.empty else 0.0
        }
    else:
        summary = {
            'sample_name': sample_name,
            'U12_Mean_UMI2_Distinct_Per_UMI1': 0.0,
            'U12_Median_UMI2_Distinct_Per_UMI1': 0.0,
            'U12_Mean_Entropy_Per_UMI1': 0.0,
            'U12_Median_Entropy_Per_UMI1': 0.0,
            'U12_Total_UMI1_Groups': 0,
            'U12_Mean_UMI1_Count': 0.0
        }

    return per_group_df, summary


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze UMI patterns and statistics")
    parser.add_argument("--umi-dir", required=True, help="Directory containing UMI mapping CSV files")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--samples", help="Comma-separated list of sample names")
    
    args = parser.parse_args()
    
    samples = None
    if args.samples:
        samples = [s.strip() for s in args.samples.split(',')]
    
    umi_meta_analysis(
        umi_dir=args.umi_dir,
        output_dir=args.output_dir,
        samples=samples
    )
