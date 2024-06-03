#!/usr/bin/env python


import pysam
import argparse
import collections
import pandas as pd
import glob
import numpy as np

def main():
    args = parse_args()
    run(args.input_dir, args.out_dir)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest='input_dir',
                       help='dir to chunks',
                       required=True)
    parser.add_argument('-o', dest='out_dir', 
                        help="out dir", 
                        required=True)
    args = parser.parse_args()
    return args


def run(input_dir, out_dir):
    outname='outfile.txt'
    merge_files = merge_chunk_csv(input_dir, out_dir, outname)



def merge_chunk_csv(input_dir, out_dir, outname):
    all_files = glob.glob(input_dir + "/*.csv")
    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)
    frame = pd.concat(li, axis=0, ignore_index=True)
    out_name = out_dir + "/" + outname
    pd.DataFrame.to_csv(frame, out_name, index=False)


if __name__ == "__main__":
    main()
