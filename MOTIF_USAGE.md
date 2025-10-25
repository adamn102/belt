# BELT Motif Parameter Usage

All BELT commands now handle motifs consistently, supporting both single motifs and multiple motifs.

## Updated Commands

The following commands have been updated to handle motifs consistently:
1. `extract` - Extract UMIs from BAM files
2. `umi-map` - Map UMIs to whitelist
3. `whitelist` - Generate UMI whitelist

## How Motif Parameters Work

### Single Motif (String)
```bash
# Just pass a single motif string
belt extract --motif "NNNNNNNNNNNN" ...
belt umi-map --motif "NNNNNNNNNNNN" ...
belt whitelist --motif "NNNNNNNNNNNN" ...
```

### Multiple Motifs (List)
```bash
# Pass as a Python list string
belt extract --motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"]' ...
belt umi-map --motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"]' ...
belt whitelist --motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"]' ...
```

### Default Behavior
If no `--motif` parameter is provided, all commands default to a single 12bp N-motif:
```python
default = ["NNNNNNNNNNNN"]
```

## Examples

### Example 1: Single UMI Extraction
```bash
belt extract \
    --reference reference.fa \
    --alignment reads.bam \
    --output output/ \
    --motif "NNNNNNNNNNNN"
```

### Example 2: Five UMI Extraction
```bash
belt extract \
    --reference reference.fa \
    --alignment reads.bam \
    --output output/ \
    --motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"]'
```

### Example 3: UMI Mapping with Multiple Motifs
```bash
belt umi-map \
    --reference reference.fa \
    --alignment reads.bam \
    --whitelist whitelist.csv \
    --output_dir output/ \
    --motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"]'
```

### Example 4: Whitelist Generation with Multiple Motifs
```bash
belt whitelist \
    --reference reference.fa \
    --alignment reads.bam \
    --output_dir output/ \
    --motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN"]'
```

## Technical Details

### Parsing Logic
All commands use the same parsing logic:

```python
# Parse motif if it's a string representation of a list
if isinstance(motif, str):
    import ast
    try:
        motif = ast.literal_eval(motif)  # Parse list string
    except (ValueError, SyntaxError):
        motif = [motif]  # Treat as single motif
elif not isinstance(motif, list):
    motif = [motif]  # Convert to list
```

### Output Format
When multiple motifs are used, the output CSV will have columns:
- `Read Name`
- `UMI1`, `UMI2`, `UMI3`, ... (one column per motif)

FASTQ headers will include all UMIs:
```
@read_name|UMI1:GTAGCCTTAACG|UMI2:TCCGATGATCGA|UMI3:GGCATTAAGCTT|...
```

## Migration from Old Format

### Legacy 5pUMI/3pUMI Format
Old commands that used:
```bash
--motif "NHHHHHHHHHHHHHM"  # Single motif for 5' and 3' ends
```

### New Multi-UMI Format
Now use:
```bash
--motif '["NNNNNNNNNNNN", "NNNNNNNNNNNN"]'  # Separate motifs for each position
```

## Benefits of Consistent Handling

✅ **Flexible**: Supports 1 to N motifs per read
✅ **Consistent**: All commands work the same way
✅ **Backwards Compatible**: Single motif strings still work
✅ **Type Safe**: Automatic conversion to list format
✅ **User Friendly**: Clear error messages if parsing fails

## Notes

- Motif lists must use Python list syntax with square brackets and quotes
- Use single quotes around the list to avoid shell interpretation issues
- The number of UMIs extracted will match the number of motifs provided
- All motifs are searched sequentially in the reference sequence

