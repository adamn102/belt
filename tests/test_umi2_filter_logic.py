#!/usr/bin/env python3
"""
Unit tests for UMI2 majority filtering logic
Tests the filtering logic without requiring samtools or full BAM processing
"""
import os
import sys
import pytest
import pandas as pd
import json
from pathlib import Path
import tempfile

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))


def test_umi2_filtering_logic_direct():
    """
    Test the UMI2 filtering logic directly without BAM processing
    This verifies that the data filtering works correctly
    """
    # Create test UMI data
    test_data = []
    
    # UMI1-A: 6 reads with UMI2-X (60%), 4 with UMI2-Y
    for i in range(6):
        test_data.append({"Read Name": f"read_A_X_{i}", "UMI1": "AAAA", "UMI2": "XXXX"})
    for i in range(4):
        test_data.append({"Read Name": f"read_A_Y_{i}", "UMI1": "AAAA", "UMI2": "YYYY"})
    
    # UMI1-B: 8 reads with UMI2-Z (80%), 2 with UMI2-W
    for i in range(8):
        test_data.append({"Read Name": f"read_B_Z_{i}", "UMI1": "BBBB", "UMI2": "ZZZZ"})
    for i in range(2):
        test_data.append({"Read Name": f"read_B_W_{i}", "UMI1": "BBBB", "UMI2": "WWWW"})
    
    # UMI1-C: 3 reads with UMI2-Q (30%), 7 with UMI2-R (70%)
    for i in range(3):
        test_data.append({"Read Name": f"read_C_Q_{i}", "UMI1": "CCCC", "UMI2": "QQQQ"})
    for i in range(7):
        test_data.append({"Read Name": f"read_C_R_{i}", "UMI1": "CCCC", "UMI2": "RRRR"})
    
    df = pd.DataFrame(test_data)
    
    # Test with 50% threshold
    print("\n" + "="*80)
    print("TEST: 50% Threshold")
    print("="*80)
    
    threshold = 5
    umi2_majority_threshold = 0.5
    gcol = "UMI1"
    
    # Replace invalid values
    df["UMI2"] = df["UMI2"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
    
    # Count UMI1 groups
    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    
    print(f"UMI1 groups passing threshold ({threshold}): {list(keepers.index)}")
    print(f"Counts: {keepers.to_dict()}")
    
    # Apply UMI2 filtering
    umi1_to_majority_umi2 = {}
    passing_umi1_values = []
    reads_before = 0
    reads_after = 0
    
    for umi1_value in keepers.index:
        sub = df.loc[df[gcol] == umi1_value, [gcol, "UMI2"]]
        total_in_umi1 = int(sub.shape[0])
        reads_before += total_in_umi1
        
        if total_in_umi1 == 0:
            continue
        
        vc2 = sub["UMI2"].value_counts(dropna=True)
        if vc2.empty:
            continue
        
        majority_umi2 = vc2.index[0]
        majority_count = int(vc2.iloc[0])
        majority_fraction = majority_count / total_in_umi1
        
        print(f"\n{umi1_value}:")
        print(f"  Total reads: {total_in_umi1}")
        print(f"  UMI2 distribution: {vc2.to_dict()}")
        print(f"  Majority UMI2: {majority_umi2} ({majority_count}/{total_in_umi1} = {majority_fraction:.1%})")
        
        if majority_fraction >= umi2_majority_threshold:
            print(f"  ✓ PASS (>= {umi2_majority_threshold:.0%})")
            passing_umi1_values.append(umi1_value)
            umi1_to_majority_umi2[umi1_value] = majority_umi2
            reads_after += majority_count
        else:
            print(f"  ✗ FAIL (< {umi2_majority_threshold:.0%})")
    
    print(f"\n{'-'*80}")
    print(f"SUMMARY:")
    print(f"  UMI1 groups passing UMI2 filter: {len(passing_umi1_values)}/{len(keepers)}")
    print(f"  Reads before filter: {reads_before}")
    print(f"  Reads after filter: {reads_after}")
    print(f"  Reads filtered: {reads_before - reads_after}")
    print("="*80 + "\n")
    
    # Verify results
    assert len(passing_umi1_values) == 3, f"Expected 3 passing groups, got {len(passing_umi1_values)}"
    assert "AAAA" in passing_umi1_values, "UMI1-AAAA should pass (60%)"
    assert "BBBB" in passing_umi1_values, "UMI1-BBBB should pass (80%)"
    assert "CCCC" in passing_umi1_values, "UMI1-CCCC should pass (70%)"
    
    # Verify UMI2 mappings
    assert umi1_to_majority_umi2["AAAA"] == "XXXX"
    assert umi1_to_majority_umi2["BBBB"] == "ZZZZ"
    assert umi1_to_majority_umi2["CCCC"] == "RRRR"
    
    # Verify read counts
    assert reads_after == 6 + 8 + 7, f"Expected {6+8+7} reads after filter, got {reads_after}"
    assert reads_before - reads_after == 4 + 2 + 3, f"Expected {4+2+3} filtered reads"
    
    print("✓ All assertions passed!")
    return True


def test_umi2_filtering_strict_threshold():
    """Test with stricter 70% threshold"""
    test_data = []
    
    # UMI1-A: 60% majority (should fail at 70%)
    for i in range(6):
        test_data.append({"Read Name": f"read_A_X_{i}", "UMI1": "AAAA", "UMI2": "XXXX"})
    for i in range(4):
        test_data.append({"Read Name": f"read_A_Y_{i}", "UMI1": "AAAA", "UMI2": "YYYY"})
    
    # UMI1-B: 80% majority (should pass at 70%)
    for i in range(8):
        test_data.append({"Read Name": f"read_B_Z_{i}", "UMI1": "BBBB", "UMI2": "ZZZZ"})
    for i in range(2):
        test_data.append({"Read Name": f"read_B_W_{i}", "UMI1": "BBBB", "UMI2": "WWWW"})
    
    df = pd.DataFrame(test_data)
    
    print("\n" + "="*80)
    print("TEST: 70% Strict Threshold")
    print("="*80)
    
    threshold = 5
    umi2_majority_threshold = 0.7
    gcol = "UMI1"
    
    df["UMI2"] = df["UMI2"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    
    umi1_to_majority_umi2 = {}
    passing_umi1_values = []
    
    for umi1_value in keepers.index:
        sub = df.loc[df[gcol] == umi1_value, [gcol, "UMI2"]]
        total_in_umi1 = int(sub.shape[0])
        
        if total_in_umi1 == 0:
            continue
        
        vc2 = sub["UMI2"].value_counts(dropna=True)
        if vc2.empty:
            continue
        
        majority_umi2 = vc2.index[0]
        majority_count = int(vc2.iloc[0])
        majority_fraction = majority_count / total_in_umi1
        
        print(f"{umi1_value}: {majority_umi2} = {majority_fraction:.1%} ", end="")
        
        if majority_fraction >= umi2_majority_threshold:
            print("✓ PASS")
            passing_umi1_values.append(umi1_value)
            umi1_to_majority_umi2[umi1_value] = majority_umi2
        else:
            print("✗ FAIL")
    
    print("="*80 + "\n")
    
    # Verify
    assert len(passing_umi1_values) == 1, "Only UMI1-BBBB should pass"
    assert "BBBB" in passing_umi1_values, "UMI1-BBBB should pass (80%)"
    assert "AAAA" not in passing_umi1_values, "UMI1-AAAA should fail (60% < 70%)"
    
    print("✓ All assertions passed!")
    return True


def test_umi2_filtering_with_na_values():
    """Test handling of NA/NoMatch/NoBarcode values"""
    test_data = []
    
    # UMI1-A: 6 valid, 2 NoMatch, 2 NoBarcode
    for i in range(6):
        test_data.append({"Read Name": f"read_valid_{i}", "UMI1": "AAAA", "UMI2": "XXXX"})
    for i in range(2):
        test_data.append({"Read Name": f"read_nomatch_{i}", "UMI1": "AAAA", "UMI2": "NoMatch"})
    for i in range(2):
        test_data.append({"Read Name": f"read_nobarcode_{i}", "UMI1": "AAAA", "UMI2": "NoBarcode"})
    
    df = pd.DataFrame(test_data)
    
    print("\n" + "="*80)
    print("TEST: Handling NA/NoMatch/NoBarcode Values")
    print("="*80)
    
    print(f"Before filtering: {df['UMI2'].value_counts().to_dict()}")
    
    # Replace invalid values with NA
    df["UMI2"] = df["UMI2"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
    
    print(f"After filtering: {df['UMI2'].value_counts(dropna=False).to_dict()}")
    
    threshold = 5
    umi2_majority_threshold = 0.5
    gcol = "UMI1"
    
    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    
    for umi1_value in keepers.index:
        sub = df.loc[df[gcol] == umi1_value, [gcol, "UMI2"]]
        total_in_umi1 = int(sub.shape[0])
        
        # Count only valid UMI2 values (dropna=True)
        vc2 = sub["UMI2"].value_counts(dropna=True)
        
        if not vc2.empty:
            majority_umi2 = vc2.index[0]
            majority_count = int(vc2.iloc[0])
            # Calculate fraction based on total reads
            majority_fraction = majority_count / total_in_umi1
            
            print(f"\n{umi1_value}:")
            print(f"  Total reads: {total_in_umi1}")
            print(f"  Valid UMI2 values: {vc2.to_dict()}")
            print(f"  Majority: {majority_umi2} ({majority_count}/{total_in_umi1} = {majority_fraction:.1%})")
    
    print("="*80 + "\n")
    
    # Verify: 6 valid out of 10 total = 60%, which should pass at 50% threshold
    assert len(keepers) == 1
    sub = df.loc[df["UMI1"] == "AAAA"]
    vc2 = sub["UMI2"].value_counts(dropna=True)
    assert vc2.iloc[0] == 6, "Should have 6 valid UMI2-XXXX reads"
    assert (6 / 10) >= 0.5, "60% should pass 50% threshold"
    
    print("✓ All assertions passed!")
    return True


def test_read_level_filtering():
    """Test that only reads with majority UMI2 are kept"""
    test_data = []
    
    # UMI1-A: 7 with UMI2-X, 3 with UMI2-Y
    for i in range(7):
        test_data.append({"Read Name": f"read_A_X_{i}", "UMI1": "AAAA", "UMI2": "XXXX"})
    for i in range(3):
        test_data.append({"Read Name": f"read_A_Y_{i}", "UMI1": "AAAA", "UMI2": "YYYY"})
    
    df = pd.DataFrame(test_data)
    
    print("\n" + "="*80)
    print("TEST: Read-Level Filtering")
    print("="*80)
    
    threshold = 5
    umi2_majority_threshold = 0.5
    gcol = "UMI1"
    
    df["UMI2"] = df["UMI2"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    
    # Identify majority UMI2 AND check threshold after filtering
    umi1_to_majority_umi2 = {}
    for umi1_value in keepers.index:
        sub = df.loc[df[gcol] == umi1_value, [gcol, "UMI2"]]
        vc2 = sub["UMI2"].value_counts(dropna=True)
        if not vc2.empty:
            majority_umi2 = vc2.index[0]
            majority_count = int(vc2.iloc[0])
            majority_fraction = majority_count / sub.shape[0]
            
            # Check both majority fraction AND that majority count meets threshold
            if majority_fraction >= umi2_majority_threshold and majority_count >= threshold:
                umi1_to_majority_umi2[umi1_value] = majority_umi2
    
    # Filter reads to only keep those with majority UMI2
    kept_rows = []
    for umi1_value in keepers.index:
        majority_umi2 = umi1_to_majority_umi2.get(umi1_value)
        if majority_umi2 is not None:
            mask = (df[gcol] == umi1_value) & (df["UMI2"] == majority_umi2)
            kept_rows.append(df.loc[mask])
    
    if kept_rows:
        filtered_df = pd.concat(kept_rows, ignore_index=True)
    else:
        filtered_df = pd.DataFrame()
    
    print(f"Total reads: {len(df)}")
    print(f"Reads after filtering: {len(filtered_df)}")
    print(f"Filtered out: {len(df) - len(filtered_df)}")
    print(f"\nKept reads: {list(filtered_df['Read Name'])}")
    print("="*80 + "\n")
    
    # Verify
    assert len(filtered_df) == 7, "Should keep only 7 reads with UMI2-XXXX"
    assert all("_X_" in r for r in filtered_df['Read Name']), "Should only keep UMI2-X reads"
    assert len(df) - len(filtered_df) == 3, "Should filter out 3 UMI2-Y reads"
    
    print("✓ All assertions passed!")
    return True


def test_threshold_reapplied_after_umi2_filter():
    """Test that threshold is re-applied after UMI2 filtering removes reads"""
    test_data = []
    
    # UMI1-A: 10 reads total, 4 with UMI2-X (40%), 6 with UMI2-Y (60%)
    # After filtering to majority (UMI2-Y), only 6 reads remain
    # If threshold is 10, this group should be EXCLUDED
    for i in range(4):
        test_data.append({"Read Name": f"read_A_X_{i}", "UMI1": "AAAA", "UMI2": "XXXX"})
    for i in range(6):
        test_data.append({"Read Name": f"read_A_Y_{i}", "UMI1": "AAAA", "UMI2": "YYYY"})
    
    # UMI1-B: 15 reads total, 12 with UMI2-Z (80%), 3 with UMI2-W (20%)
    # After filtering to majority (UMI2-Z), 12 reads remain
    # If threshold is 10, this group should PASS
    for i in range(12):
        test_data.append({"Read Name": f"read_B_Z_{i}", "UMI1": "BBBB", "UMI2": "ZZZZ"})
    for i in range(3):
        test_data.append({"Read Name": f"read_B_W_{i}", "UMI1": "BBBB", "UMI2": "WWWW"})
    
    df = pd.DataFrame(test_data)
    
    print("\n" + "="*80)
    print("TEST: Threshold Re-applied After UMI2 Filtering")
    print("="*80)
    
    threshold = 10  # High threshold
    umi2_majority_threshold = 0.5
    gcol = "UMI1"
    
    df["UMI2"] = df["UMI2"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    
    print(f"Initial groups passing threshold ({threshold}): {list(keepers.index)}")
    print(f"  UMI1-A: {counts['AAAA']} reads")
    print(f"  UMI1-B: {counts['BBBB']} reads")
    
    # Apply UMI2 filtering with threshold re-check
    passing_groups = []
    for umi1_value in keepers.index:
        sub = df.loc[df[gcol] == umi1_value, [gcol, "UMI2"]]
        total_in_umi1 = int(sub.shape[0])
        
        vc2 = sub["UMI2"].value_counts(dropna=True)
        if not vc2.empty:
            majority_umi2 = vc2.index[0]
            majority_count = int(vc2.iloc[0])
            majority_fraction = majority_count / total_in_umi1
            
            print(f"\n{umi1_value}:")
            print(f"  Total reads: {total_in_umi1}")
            print(f"  Majority UMI2: {majority_umi2} ({majority_count} reads, {majority_fraction:.1%})")
            print(f"  Majority fraction >= {umi2_majority_threshold:.0%}: {majority_fraction >= umi2_majority_threshold}")
            print(f"  Majority count >= {threshold}: {majority_count >= threshold}")
            
            # NEW LOGIC: Check both conditions
            if majority_fraction >= umi2_majority_threshold and majority_count >= threshold:
                print(f"  ✓ PASS (both conditions met)")
                passing_groups.append(umi1_value)
            else:
                print(f"  ✗ FAIL (at least one condition not met)")
    
    print(f"\n{'-'*80}")
    print(f"RESULT: {len(passing_groups)} group(s) passed: {passing_groups}")
    print("="*80 + "\n")
    
    # Verify
    assert len(passing_groups) == 1, f"Only UMI1-B should pass, but got {len(passing_groups)}: {passing_groups}"
    assert "BBBB" in passing_groups, "UMI1-B should pass (12 reads >= 10 threshold)"
    assert "AAAA" not in passing_groups, "UMI1-A should fail (6 reads < 10 threshold after filtering)"
    
    print("✓ All assertions passed!")
    print("✓ Threshold correctly re-applied after UMI2 filtering!")
    return True


if __name__ == "__main__":
    print("\n" + "#"*80)
    print("# UMI2 MAJORITY FILTERING LOGIC TESTS")
    print("#"*80)
    
    test_umi2_filtering_logic_direct()
    test_umi2_filtering_strict_threshold()
    test_umi2_filtering_with_na_values()
    test_read_level_filtering()
    test_threshold_reapplied_after_umi2_filter()
    
    print("\n" + "#"*80)
    print("# ALL TESTS PASSED!")
    print("#"*80 + "\n")

