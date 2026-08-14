#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict

def distribute_counts(start_pos, length, window_size):
    """
    Yields (window_index, overlap_length) for a CIGAR operation 
    that might span across one or more window boundaries.
    """
    end_pos = start_pos + length
    current_pos = start_pos
    
    while current_pos < end_pos:
        win_idx = current_pos // window_size
        win_end = (win_idx + 1) * window_size
        overlap = min(end_pos, win_end) - current_pos
        yield win_idx, overlap
        current_pos += overlap

def calculate_windowed_metrics(bam_path, output_path, window_size):
    """
    Parses an --eqx BAM file and calculates windowed variant metrics.
    """
    # Dictionary structure: windows[ref_name][window_index] = counts_dict
    windows = defaultdict(lambda: defaultdict(lambda: {'=': 0, 'X': 0, 'I': 0, 'D': 0, 'M': 0}))
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # Get reference lengths to write empty windows at the end if needed
        ref_lengths = {bam.get_reference_name(i): bam.get_reference_length(i) for i in range(bam.nreferences)}
        
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                continue
                
            ref_name = read.reference_name
            ref_pos = read.reference_start
            
            for op, length in read.cigartuples:
                if op == 7:  # '=' Sequence Match
                    for w_idx, count in distribute_counts(ref_pos, length, window_size):
                        windows[ref_name][w_idx]['='] += count
                    ref_pos += length
                    
                elif op == 8:  # 'X' Sequence Mismatch
                    for w_idx, count in distribute_counts(ref_pos, length, window_size):
                        windows[ref_name][w_idx]['X'] += count
                    ref_pos += length
                    
                elif op == 2:  # 'D' Deletion (Query is missing bases the reference has)
                    for w_idx, count in distribute_counts(ref_pos, length, window_size):
                        windows[ref_name][w_idx]['D'] += count
                    ref_pos += length
                    
                elif op == 1:  # 'I' Insertion (Query has bases the reference lacks)
                    # Insertions don't consume reference space, so they fall into the current window
                    w_idx = ref_pos // window_size
                    windows[ref_name][w_idx]['I'] += length
                    
                elif op == 0:  # 'M' Standard Match/Mismatch (Fallback if --eqx wasn't used)
                    for w_idx, count in distribute_counts(ref_pos, length, window_size):
                        windows[ref_name][w_idx]['M'] += count
                    ref_pos += length
                    
                elif op == 3:  # 'N' Skipped region from the reference
                    ref_pos += length
                    
                # Ops 4 (S), 5 (H), 6 (P) do not consume reference bases

    # Write output TSV
    with open(output_path, 'w') as out:
        header = [
            "chrom", "start", "end", "matches", "mismatches", 
            "insertions", "deletions", "unclassified_M", "base_identity"
        ]
        out.write("\t".join(header) + "\n")
        
        for ref_name, ref_len in ref_lengths.items():
            if ref_name not in windows:
                continue
                
            max_win_idx = (ref_len // window_size) + 1
            
            for w_idx in range(max_win_idx):
                start = w_idx * window_size
                end = min(start + window_size, ref_len)
                
                # If window is beyond the reference length, stop
                if start >= ref_len:
                    break
                    
                counts = windows[ref_name].get(w_idx, {'=': 0, 'X': 0, 'I': 0, 'D': 0, 'M': 0})
                
                exact_matches = counts['=']
                mismatches = counts['X']
                insertions = counts['I']
                deletions = counts['D']
                unclass_m = counts['M']
                
                total_errors = mismatches + insertions + deletions
                aligned_bases = exact_matches + unclass_m + total_errors
                
                if aligned_bases > 0:
                    # Identity = Matches / Total Aligned sequence (including indels)
                    identity = (exact_matches + unclass_m) / aligned_bases
                else:
                    identity = 0.0
                    
                row = [
                    ref_name,
                    str(start),
                    str(end),
                    str(exact_matches),
                    str(mismatches),
                    str(insertions),
                    str(deletions),
                    str(unclass_m),
                    f"{identity:.6f}"
                ]
                out.write("\t".join(row) + "\n")

    print(f"Windowed metrics written to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate windowed metrics from an extended CIGAR BAM.")
    parser.add_argument("--bam", required=True, help="Path to the filtered BAM file.")
    parser.add_argument("--out", required=True, help="Path to the output TSV file.")
    parser.add_argument("--window", type=int, default=100000, help="Window size in base pairs (default: 100,000).")
    args = parser.parse_args()
    
    calculate_windowed_metrics(args.bam, args.out, args.window)