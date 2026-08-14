#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns

CDS_depths = {}

def extract_CDS_depths(input_file):
    global CDS_depths
    
    with open(input_file, "r") as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) != 2:
                raise ValueError("Input file is not in the expected format!")
            
            info_field, depth = cols
            CDS = info_field.split(";")[-1]
            depth = float(depth)
            
            if CDS not in CDS_depths:
                CDS_depths[CDS] = [depth]
            else:
                CDS_depths[CDS].append(depth)
                

def  main(input_dir, output_dir_name):
    files = [f for f in os.listdir(input_dir) if f.endswith(".cds.median.depth")]
    for file in files:
        file_path = os.path.join(input_dir, file)
        extract_CDS_depths(file_path)
    
    output_dir = os.path.join(os.path.dirname(os.path.dirname(input_dir)), output_dir_name)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    low_depth_cds = []  # CDS keys with at least one depth less than 5
    low_median_cds = []  # CDS keys with median depth less than 5
    low_mean_cds = []  # CDS keys with mean depth less than 5
    
      
    for CDS, depths in CDS_depths.items():
        #plt.figure(figsize=(8, 6))
        #sns.violinplot(depths)
        #plt.title(f'Depth distribution for CDS: {CDS}')
        #plt.xlabel('CDS')
        #plt.ylabel('Depth')
        #plt.savefig(os.path.join(output_dir, f'{CDS}_depth_violinplot.pdf'))
        #plt.close()
    
        if any(depth < 10 for depth in depths):
            low_depth_cds.append(CDS)
            
        median_depth = np.median(depths)
        if median_depth < 10:
            low_median_cds.append(CDS)
            
        mean_depth = np.mean(depths)
        if mean_depth < 10:
            low_mean_cds.append(CDS)

    sys.stderr.write('\n==> Wrote violinplots for CDS depths to respective files.')

    with open(os.path.join(output_dir, 'low_depth_cds.txt'), 'w') as file:
        for cds in sorted(low_depth_cds):
            file.write(cds + '\n')
    sys.stderr.write('\n==> Wrote CDS with at least 1 depth < 10 to: %s' % (os.path.abspath(os.path.join(output_dir, 'low_depth_cds.txt'))))

    with open(os.path.join(output_dir, 'low_median_cds.txt'), 'w') as file:
        for cds in sorted(low_median_cds):
            file.write(cds + '\n')
    sys.stderr.write('\n==> Wrote CDS with median depth < 10 to: %s' % (os.path.abspath(os.path.join(output_dir, 'low_median_cds.txt'))))
    
    with open(os.path.join(output_dir, 'low_mean_cds.txt'), 'w') as file:
        for cds in sorted(low_mean_cds):
            file.write(cds + '\n')
    sys.stderr.write('\n==> Wrote CDS with mean depth < 10 to: %s' % (os.path.abspath(os.path.join(output_dir, 'low_mean_cds.txt'))))

    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python CDS_violinplots.py <input_directory> <output_directory>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir_name = sys.argv[2]
    
    if not os.path.isdir(input_dir):
        print("Input directory does not exist.")
        sys.exit(1)

    main(input_dir, output_dir_name)
