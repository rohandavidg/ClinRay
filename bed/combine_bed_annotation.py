#!/usr/bin/env python

"""
This script combines
the annotation from the stratification
files into a single
dataframe
"""

import pysam
import argparse
import collections
import pandas as pd
import glob
import numpy as np
from functools import reduce

def main():
    args = parse_args()
    run(args.input_dir,
        args.querry_bed,
        args.outname)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest='input_dir',
                       help='dir to chunks',
                       required=True)
    parser.add_argument('-q', dest='querry_bed', 
                        help="querry bed file used for analysis", 
                        required=True)
    parser.add_argument('-o', dest='outname', 
                        help="name of the output file name", 
                        required=True)    
    args = parser.parse_args()
    return args


def run(input_dir,querry_bed, outname):
    merge_files = merge_chunk_tsv(input_dir,
                                  querry_bed,
                                  outname)

    
def merge_chunk_tsv(input_dir, querry_bed, outname):
    all_files = glob.glob(input_dir + "/*.tsv")
    df = pd.read_csv(querry_bed, sep='\t', names = ['chrom', 'start', 'end',
                                                    'gene', 'number', 'strand'])
    df['chrom'] = df['chrom'].astype(str)
    df['start'] =  df['start'].astype(str)
    df['end'] = df['end'].astype(str)
    list_df = [df]
    for filename in all_files:
        new_df = pd.read_csv(filename, index_col=None, sep='\t')
        new_df['chrom'] = new_df['chrom'].astype(str)
        new_df['start'] = new_df['start'].astype(str)
        new_df['end'] = new_df['end'].astype(str)
        list_df.append(new_df)
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['chrom', 'start', 'end'],
                                                    how='outer'), list_df)
    df_merged.fillna(0, inplace=True)
    pd.DataFrame.to_csv(df_merged, outname, sep='\t',index=False)


if __name__ == "__main__":
    main()
