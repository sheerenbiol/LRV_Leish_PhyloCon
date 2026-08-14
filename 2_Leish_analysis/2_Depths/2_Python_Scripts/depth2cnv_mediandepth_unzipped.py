#!/usr/bin/env python

'''
INPUT = samtools depth file, bed file
OUTPUTs:
  * haploid cds depth (out.cds.depth)
  * cds length (out.cds.length)
  * somy (out.somy)
  * histogram of depths (out.hist.pdf)

Also, output something you can see in igv: YES. If cds depth = 1, then gray; else red = << 1 and blue >> 1
'''

from __future__ import division
from collections import defaultdict
from collections import Counter

import os
import sys
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt


def my_mode(sample):
    c = Counter(sample)
    return [k for k, v in c.items() if v == c.most_common(1)[0][1]]


def _make_bed_allpos(bedfile):
  '''
  For internal use only: used to quickly parse depth files, avoiding for loops
  This script may need to be adapted depending on the information in the bedfile!
  '''
  newbed=defaultdict(list)
  ngene=0
  with open(bedfile, "r") as f:
    for line in f:
      ngene = ngene+1
      gene = line.split()
      gene_chr = str(gene[0])
      gene_start = int(gene[1])
      gene_end = int(gene[2])
      gene_strand = gene[3]
      gene_id = str(gene[5])
      new_gene_id="chr=" + str(gene_chr) + ";start=" + str(gene_start) + ";end=" + str(gene_end) + ";strand=" + str(gene_strand) + ";" + str(gene_id)
      
      for n in range(gene_start,(gene_end+1)):
        newloc=gene_chr+":"+str(n)
        newbed[newloc]=new_gene_id
  
  sys.stderr.write('\n==> Found and processed %i genes.\n' % (ngene))
  return newbed


def main(depthfile, bedfile, output):
  sys.stderr.write('\n==> Creating a bed dictionary.')
  bed = _make_bed_allpos(bedfile)

  sys.stderr.write('==> Reading and storing depth file. This might take a while.')
  chrdepths = defaultdict(list)
  genedepths = defaultdict(list)
  with open(depthfile, 'rb') as f:
    for line in f:
      #if line.startswith(('LbrM_0','LbrM_1','LbrM_2','LbrM_3')):
        line=bytes.decode(line)
        loc=line.split('\t')
        loc_chr=str(loc[0])
        loc_pos=int(loc[1])
        loc_depth=float(loc[2])
        loc_chr_pos=loc_chr + ":" + str(loc_pos)
        
        if loc_chr_pos == '[]':
          print(loc)
          
        ## Storing chromosomal depths
        chrdepths[loc_chr].append(loc_depth)
        
        ## Storing CDS depths
        if loc_chr_pos in bed:
          genedepths[str(bed[loc_chr_pos])].append(loc_depth)
      
  sys.stderr.write('\n\n==> Estimating median and haploid chromosomal depths.')
  
  ## concatenating genome-wide list of integer depth values
  genome_depths = []
  chromosomes = ('LbrM.01_v4_pilon','LbrM.02_v4_pilon','LbrM.03_v4_pilon','LbrM.04_v4_pilon','LbrM.05_v4_pilon','LbrM.06_v4_pilon','LbrM.07_v4_pilon','LbrM.08_v4_pilon','LbrM.09_v4_pilon','LbrM.10_v4_pilon','LbrM.11_v4_pilon','LbrM.12_v4_pilon','LbrM.13_v4_pilon','LbrM.14_v4_pilon','LbrM.15_v4_pilon','LbrM.16_v4_pilon','LbrM.17_v4_pilon','LbrM.18_v4_pilon','LbrM.19_v4_pilon','LbrM.20_v4_pilon','LbrM.21_v4_pilon','LbrM.22_v4_pilon','LbrM.23_v4_pilon','LbrM.24_v4_pilon','LbrM.25_v4_pilon','LbrM.26_v4_pilon','LbrM.27_v4_pilon','LbrM.28_v4_pilon','LbrM.29_v4_pilon','LbrM.30_v4_pilon','LbrM.31_v4_pilon','LbrM.32_v4_pilon','LbrM.33_v4_pilon','LbrM.34_v4_pilon','LbrM.35_v4_pilon',)
  chrdepths_chromosomes={ key:value for key,value in chrdepths.items() if key in chromosomes}
  for items in chrdepths_chromosomes.values():
    genome_depths.extend(items)
  

  #########################################################
  ## calculating genome-wide median, mean and mode depth ##
  #########################################################
  genome_depths=np.array(genome_depths)
  genomemedian=np.median(genome_depths)
  genomemean=round(np.mean(genome_depths),1)
  genomemode=my_mode(genome_depths[genome_depths>10])
  genomemode=float("".join([str(integer) for integer in genomemode]))
  
  ## plotting histogram
  xcoords = [genomemode, genomemean, genomemedian]
  colors = ['r','k','b']
  leg = ['mode', 'mean', 'median']
  
  y=range(0,300,5)
  plt.ioff()
  fig = plt.figure()
  plt.xlim(0,300)
  plt.hist(genome_depths, bins = y)
  for xc,c,t in zip(xcoords,colors,leg):
    plt.axvline(x=xc, label='{}'.format(t), c=c)
  plt.title('Read depth distribution ' + str(output))
  plt.xlabel('Read Depth')
  plt.ylabel('Frequency')
  plt.legend()
  plt.savefig(str(output) + '.depth.hist.pdf')
  plt.close(fig)

  sys.stderr.write('\n==> Wrote depth histogram to: %s' % (os.path.abspath(output+'.depth.hist.pdf')))

  
  #########################################################
  ## calculating chromosomal median, mean and mode depth ##
  #########################################################
  medianchrdepths=dict()
  meanchrdepths=dict()
  modechrdepths=dict()
  somy=dict()

  for c in chrdepths_chromosomes:
    chrvalues=np.array(chrdepths_chromosomes[c])
    medianchrdepths[c]=np.median(chrvalues)
    meanchrdepths[c]=round(np.mean(chrvalues),1)
    modechrdepths[c]=my_mode(chrvalues[chrvalues>10])[0]
    #modechrdepths[c]=float("".join([str(integer) for integer in mode]))
    somy[c]=round(2*modechrdepths[c]/genomemode, 3)
    
    xcoords = [modechrdepths[c], meanchrdepths[c], medianchrdepths[c]]
    colors = ['r','k','b']
    leg = ['mode', 'mean', 'median']
  
    y=range(0,300,5)
    fig = plt.figure()
    plt.xlim(0,300)
    plt.hist(chrvalues, bins = y)
    for xc,xcol,t in zip(xcoords,colors,leg):
        plt.axvline(x=xc, label='{}'.format(t), c=xcol)
    plt.title('Read depth distribution ' + str(output) + c)
    plt.xlabel('Read Depth')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(str(output) + '_' + c + '_' + '.depth.hist.pdf')
    plt.close(fig)
    

  with open(output + '.somy', "w") as file:
    file.write("Genome median depth: %s\n" % (genomemedian))
    file.write("Genome mean depth: %s\n" % (genomemean))
    file.write("Genome mode depth: %s\n\n" % (genomemode))
    file.write("chromosome mean median mode somy\n")
    for item in somy:
      file.write("%s %s %s %s %s\n" % (item, meanchrdepths[item], medianchrdepths[item], modechrdepths[item], round(somy[item],1)))
  sys.stderr.write('\n==> Wrote somy estimates to: %s\n' % (os.path.abspath(output+'.somy')))


  sys.stderr.write('\n==> Estimating Median CDS depth, haploid CDS depth and CDS lengths. \n')
  cdsmediandepth=dict()
  for d in genedepths:
    try:
      dchr=str(d).split(';')[0].split('=')[1]
      ddepths=np.array(genedepths[d])
      try:
        cdsmediandepth[d]=round(np.median(ddepths), 3)
      except KeyError:
        cdsmediandepth[d]=-1
    except IndexError: 
      print(d)

  with open(output + '.cds.median.depth', "w") as file:
    for item in cdsmediandepth:
      file.write("%s %s\n" % (item, cdsmediandepth[item]))
  sys.stderr.write('==> Wrote median CDS depth to: %s\n' % (os.path.abspath(output+'.cds.median.depth')))
  
  cdsdepth=dict()
  for d in genedepths:
    try:
      dchr=str(d).split(';')[0].split('=')[1]
      ddepths=np.array(genedepths[d])
      try:
        cdsdepth[d]=round(np.median(ddepths)/modechrdepths[dchr], 3)
      except KeyError:
        cdsdepth[d]=-1
    except IndexError: 
      print(d)

  with open(output + '.cds.depth', "w") as file:
    for item in cdsdepth:
      file.write("%s %s\n" % (item, cdsdepth[item]))
  sys.stderr.write('==> Wrote haploid CDS depth to: %s\n' % (os.path.abspath(output+'.cds.depth')))


  cdslength=dict()
  for d in genedepths:
    cdslength[d]=np.sum(np.array(genedepths[d])>=5)/len(genedepths[d]) #Calculate the percentage covered by a minimum of 5 reads for each gene.
 
  with open(output + '.cds.length', "w") as file:
    for item in cdslength:
      file.write("%s %s\n" % (item, cdslength[item]))
  sys.stderr.write('==> Wrote CDS length to: %s\n' % (os.path.abspath(output+'.cds.length')))

    

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='Calculates haploid chromosomal and cds read depths.', \
	  usage = 'depth2cnv [options] <depthfile> <bedfile> <out>')
  parser.add_argument('depthfile', help='Read depths as outputted from samtools', metavar='depth')
  parser.add_argument('bedfile', help='BED file', metavar='bed')
  parser.add_argument('out', help='Prefix used for labeling output files', metavar='out')
  options = parser.parse_args()
  
  main(options.depthfile, options.bedfile, options.out)
