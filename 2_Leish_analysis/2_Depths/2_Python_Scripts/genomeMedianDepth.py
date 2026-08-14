#!/usr/bin/env python

import os
import sys
import numpy as np

def calc_median_depth(input_file):
    sample_id = os.path.basename(input_file).split(".")[0]
    with open(input_file, "r") as f:
        depths = []
        for line in f:
            cols = line.strip().split()
            if len(cols) != 2:
                raise ValueError("Input file is not in the expected format!")

            depth = float(cols[1])
            depths.append(depth)
        
        median = np.median(depths)
        return sample_id, median

def main(input_dir, output_file):
    files = [f for f in os.listdir(input_dir) if f.endswith(".cds.median.depth")]
    sample_medians = {}
    for f in files:
        file_path = os.path.join(input_dir, f)
        sample_id, median_depth = calc_median_depth(file_path)
        sample_medians[sample_id] = median_depth

    sorted_dict = dict(sorted(sample_medians.items()))

    with open(output_file, 'w') as out:
        for sample_id, median_depth in sorted_dict.items():
            if median_depth < 10:
                out.write(sample_id + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python genomeMedianDepth.py <input_directory> <output_file_name>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.isdir(input_dir):
        print("Input directory does not exist.")
        sys.exit(1)

    main(input_dir, output_file)
