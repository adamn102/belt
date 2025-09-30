#!/usr/bin/env python3
"""
BELT Barcode Statistics Analyzer

This script analyzes barcode detection statistics from the demultiplexing CSV output.
It provides detailed statistics on:
- Reads with only 5' barcode found
- Reads with only 3' barcode found  
- Reads with no barcodes found
- Reads with both barcodes found
"""

import pandas as pd
import argparse
import os
import sys
from pathlib import Path


def analyze_barcode_stats(csv_file, output_file=None, verbose=True):
    """
    Analyze barcode detection statistics from demultiplexing CSV output.
    
    Args:
        csv_file (str): Path to the motif_extraction_results.csv file
        output_file (str): Optional path to save detailed statistics
        verbose (bool): Print detailed statistics to console
        
    Returns:
        dict: Statistics summary
    """
    
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    
    # Load the CSV data
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        raise ValueError(f"Error reading CSV file: {e}")
    
    # Check required columns
    required_cols = ['start_barcode_name', 'end_barcode_name', 'assignment_reason']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    total_reads = len(df)
    
    # Analyze barcode detection patterns
    # Create boolean masks for different barcode detection scenarios
    has_5p = (df['start_barcode_name'].notna()) & (df['start_barcode_name'] != '') & (df['start_barcode_name'] != 'unknown')
    has_3p = (df['end_barcode_name'].notna()) & (df['end_barcode_name'] != '') & (df['end_barcode_name'] != 'unknown')
    
    # Categorize reads based on barcode detection
    only_5p = has_5p & ~has_3p
    only_3p = ~has_5p & has_3p
    both_barcodes = has_5p & has_3p
    no_barcodes = ~has_5p & ~has_3p
    
    # Calculate statistics
    stats = {
        'total_reads': total_reads,
        'only_5p_barcode': only_5p.sum(),
        'only_3p_barcode': only_3p.sum(),
        'both_barcodes': both_barcodes.sum(),
        'no_barcodes': no_barcodes.sum(),
        'successful_assignments': (df['assignment_reason'] == 'success').sum(),
        'motif_successful_assignments': df['assignment_reason'].str.startswith('motif_success').sum(),
        'sample_not_found': (df['assignment_reason'] == 'sample_not_found').sum(),
        'barcode_not_found': (df['assignment_reason'] == 'barcode_not_found').sum()
    }
    
    # Calculate percentages
    for key in ['only_5p_barcode', 'only_3p_barcode', 'both_barcodes', 'no_barcodes']:
        stats[f'{key}_percent'] = (stats[key] / total_reads * 100) if total_reads > 0 else 0
    
    # Print statistics
    if verbose:
        print(f"\n" + "="*80)
        print(f"📊 BARCODE DETECTION STATISTICS")
        print(f"="*80)
        print(f"   Total reads analyzed: {total_reads:,}")
        print()
        
        print(f"🔍 BARCODE DETECTION BREAKDOWN:")
        print(f"-" * 80)
        print(f"   Both 5' and 3' barcodes found:     {stats['both_barcodes']:6,} reads ({stats['both_barcodes_percent']:5.1f}%)")
        print(f"   Only 5' barcode found:            {stats['only_5p_barcode']:6,} reads ({stats['only_5p_barcode_percent']:5.1f}%)")
        print(f"   Only 3' barcode found:            {stats['only_3p_barcode']:6,} reads ({stats['only_3p_barcode_percent']:5.1f}%)")
        print(f"   No barcodes found:                {stats['no_barcodes']:6,} reads ({stats['no_barcodes_percent']:5.1f}%)")
        print()
        
        print(f"📋 ASSIGNMENT SUCCESS BREAKDOWN:")
        print(f"-" * 80)
        print(f"   Successful assignments (soft-clip): {stats['successful_assignments']:6,} reads ({stats['successful_assignments']/total_reads*100:5.1f}%)")
        print(f"   Successful assignments (motif):    {stats['motif_successful_assignments']:6,} reads ({stats['motif_successful_assignments']/total_reads*100:5.1f}%)")
        print(f"   Sample not found:                   {stats['sample_not_found']:6,} reads ({stats['sample_not_found']/total_reads*100:5.1f}%)")
        print(f"   Barcode not found:                  {stats['barcode_not_found']:6,} reads ({stats['barcode_not_found']/total_reads*100:5.1f}%)")
        print()
        
        # Show detailed breakdown for reads with only one barcode
        if stats['only_5p_barcode'] > 0:
            print(f"🔬 READS WITH ONLY 5' BARCODE:")
            print(f"-" * 80)
            only_5p_df = df[only_5p]
            assignment_reasons = only_5p_df['assignment_reason'].value_counts()
            for reason, count in assignment_reasons.items():
                print(f"   {reason:30} {count:6,} reads ({count/stats['only_5p_barcode']*100:5.1f}%)")
            
            # Show top 5' barcode names found
            if 'start_barcode_name' in only_5p_df.columns:
                top_5p_barcodes = only_5p_df['start_barcode_name'].value_counts().head(5)
                if len(top_5p_barcodes) > 0:
                    print(f"   Top 5' barcodes found:")
                    for barcode, count in top_5p_barcodes.items():
                        if pd.notna(barcode) and barcode != '':
                            print(f"     {barcode:25} {count:6,} reads")
            print()
        
        if stats['only_3p_barcode'] > 0:
            print(f"🔬 READS WITH ONLY 3' BARCODE:")
            print(f"-" * 80)
            only_3p_df = df[only_3p]
            assignment_reasons = only_3p_df['assignment_reason'].value_counts()
            for reason, count in assignment_reasons.items():
                print(f"   {reason:30} {count:6,} reads ({count/stats['only_3p_barcode']*100:5.1f}%)")
            
            # Show top 3' barcode names found
            if 'end_barcode_name' in only_3p_df.columns:
                top_3p_barcodes = only_3p_df['end_barcode_name'].value_counts().head(5)
                if len(top_3p_barcodes) > 0:
                    print(f"   Top 3' barcodes found:")
                    for barcode, count in top_3p_barcodes.items():
                        if pd.notna(barcode) and barcode != '':
                            print(f"     {barcode:25} {count:6,} reads")
            print()
        
        if stats['no_barcodes'] > 0:
            print(f"🔬 READS WITH NO BARCODES:")
            print(f"-" * 80)
            no_barcodes_df = df[no_barcodes]
            assignment_reasons = no_barcodes_df['assignment_reason'].value_counts()
            for reason, count in assignment_reasons.items():
                print(f"   {reason:30} {count:6,} reads ({count/stats['no_barcodes']*100:5.1f}%)")
            print()
        
        # Show summary of barcode detection quality
        print(f"📈 BARCODE DETECTION QUALITY SUMMARY:")
        print(f"-" * 80)
        detection_rate = (stats['both_barcodes'] + stats['only_5p_barcode'] + stats['only_3p_barcode']) / total_reads * 100
        both_barcode_rate = stats['both_barcodes'] / total_reads * 100
        print(f"   Overall barcode detection rate: {detection_rate:.1f}%")
        print(f"   Both barcodes detection rate:  {both_barcode_rate:.1f}%")
        print(f"   Single barcode detection rate: {(stats['only_5p_barcode'] + stats['only_3p_barcode']) / total_reads * 100:.1f}%")
        print(f"   No barcode detection rate:     {stats['no_barcodes_percent']:.1f}%")
        
        print(f"="*80)
    
    # Save detailed statistics to file if requested
    if output_file:
        detailed_stats = {
            'metric': [
                'Total reads',
                'Both 5\' and 3\' barcodes found',
                'Only 5\' barcode found', 
                'Only 3\' barcode found',
                'No barcodes found',
                'Successful assignments (soft-clip)',
                'Successful assignments (motif)',
                'Sample not found',
                'Barcode not found'
            ],
            'count': [
                stats['total_reads'],
                stats['both_barcodes'],
                stats['only_5p_barcode'],
                stats['only_3p_barcode'],
                stats['no_barcodes'],
                stats['successful_assignments'],
                stats['motif_successful_assignments'],
                stats['sample_not_found'],
                stats['barcode_not_found']
            ],
            'percentage': [
                100.0,
                stats['both_barcodes_percent'],
                stats['only_5p_barcode_percent'],
                stats['only_3p_barcode_percent'],
                stats['no_barcodes_percent'],
                stats['successful_assignments']/total_reads*100,
                stats['motif_successful_assignments']/total_reads*100,
                stats['sample_not_found']/total_reads*100,
                stats['barcode_not_found']/total_reads*100
            ]
        }
        
        stats_df = pd.DataFrame(detailed_stats)
        stats_df.to_csv(output_file, index=False)
        
        if verbose:
            print(f"📄 Detailed statistics saved to: {output_file}")
    
    return stats


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Analyze barcode detection statistics from BELT demultiplexing output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze statistics from CSV file
  python barcode_stats.py motif_extraction_results.csv
  
  # Save detailed statistics to file
  python barcode_stats.py motif_extraction_results.csv --output barcode_stats.csv
  
  # Analyze from demultiplexed directory
  python barcode_stats.py demultiplexed_fastq/motif_extraction_results.csv
        """
    )
    
    parser.add_argument(
        "csv_file", 
        help="Path to motif_extraction_results.csv file"
    )
    parser.add_argument(
        "--output", "-o",
        help="Path to save detailed statistics CSV file"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress console output"
    )
    
    args = parser.parse_args()
    
    try:
        stats = analyze_barcode_stats(
            csv_file=args.csv_file,
            output_file=args.output,
            verbose=not args.quiet
        )
        
        if not args.quiet:
            print("✅ Barcode statistics analysis completed successfully!")
        
        # Exit with success code
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    exit(main())
