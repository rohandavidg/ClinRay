#!/usr/bin/env python


"""
This script annotates a given bed file with
overlap and reciprocal overlap

__author__ = "Rohan Gnanaolivu"
__github__ = "ADO repo"
"""

import argparse
import pybedtools
import csv
import collections
import pandas as pd
import os
import gzip
import tempfile
import zlib
import random
import logging
import time
import datetime

def main():
    args = parse_args()
    run(args.target_bed,
        args.anno_bed,
        args.anno_name,
        args.outdir)

    
def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i',dest='target_bed',
                        help="bed file to annotate",
                        required=True)
    parser.add_argument('-a',dest='anno_bed',
                        help="annotation bed file track",
                        required=True)
    parser.add_argument('-n',dest='anno_name',
                        help="name of the annnotation to use",
                        required=True)    
    parser.add_argument('-o',dest='outdir',
                        help="name of the annnotation to use",
                        required=True)
    args = parser.parse_args()
    return args


def run(target_bed,
        anno_bed,
        anno_name,
        outdir):
    assert os.path.isfile(target_bed)
    assert os.path.isfile(anno_bed)
    random_int = random.randint(0, 99999)
    logger = configure_logger(random_int)
    querry_filename = outdir + "/" + "tmp" + str(random_int) + ".querry.bed"
    anno_filename = outdir + "/" + "tmp." + anno_name + ".bed"
    tmp_querry_bed = create_tmp_bed(target_bed,
                                    "querry",
                                    querry_filename,
                                    logger)
    tmp_anno_bed = create_tmp_bed(anno_bed,
                                  anno_name,
                                  anno_filename,
                                  logger)
    overlap_count_bed = bedtools_overlap(querry_filename,
                                         anno_filename,
                                         anno_name,
                                         outdir,
                                         logger)
    os.remove(querry_filename)
    os.remove(anno_filename)

def configure_logger(value):
    """
    setting up logging
    """
    logger = logging.getLogger('Stratification_Add')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime("Stratification_Add" + str(value) + "-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

    
class BedFile(object):

    def __init__(self, chrom, start, stop, anno):
        self.chrom = chrom
        self.start =  start
        self.stop  = stop
        self.anno = anno
        
    def __iter__(self):
        return self
    
    def _chrom(self):
        if self.chrom.startswith("chr"):
            return self.chrom.strip()
        else:
            new_chrom = "chr" + self.chrom.strip()
            return new_chrom

    def _start(self):
        if self.start.isnumeric():
            return(self.start.strip())
        else:
            if isinstance(int(self.start), int):
                return(self.stop.strip())
            else:
                print("start position is not an integer {0}".format(self.start))
                exit(1)

    def _stop(self):
        if self.stop.isnumeric():
            return self.stop
        else:
            if isinstance(int(self.stop), int):
                return self.stop
            else:
                print("stop position is not an integer {0}".format(self.stop))
                exit(1)        
    
    def _anno(self):
        if self.anno:
            return self.anno
    
    def _compute_length(self):
        length = int(self.stop) - int(self.start) 
        return length

    
def create_tmp_bed(bed_file,
                   anno,
                   tmp_filename,
                   logger):
    if bed_file.endswith("gz"):
        with gzip.open(bed_file, 'rt') as fh:
            logger.info("bed file in gzip format {0}".format(bed_file))
            BedLine = fh.readlines()
            parse_bed(BedLine,
                      anno,
                      tmp_filename,
                      logger)
    else:
        with open(bed_file, 'rb') as fh:
            lines = [x.decode('utf8').strip() for x in fh.readlines()]
            parse_bed(lines,
                      anno,
                      tmp_filename,
                      logger)

def parse_bed(bed,
              anno,
              tmp_filename,
              logger):
    with open(tmp_filename, 'w') as fout:
        for line in bed:
            try:
                chrom = line.split('\t')[0].strip()
                start = line.split('\t')[1].strip()
                stop = line.split('\t')[2].strip()
                get_bed = BedFile(chrom,
                                  start,
                                  stop,
                                  anno)
                tmp_chrom = get_bed._chrom()
                tmp_start = get_bed._start()
                tmp_stop = get_bed._stop()
                tmp_anno = get_bed._anno()
                reg_length = get_bed._compute_length()
                out = [tmp_chrom,
                       tmp_start,
                       tmp_stop,
                       tmp_anno,
                       reg_length]
                fout.write('\t'.join(str(i) for i in out) + '\n')
            except IndexError:
                logger.debug("bed file is not in valid bed format {0}".format(bed))
                exit(1)
            
            
def bedtools_overlap(querry_tmp,
                     anno_tmp,
                     anno_name,
                     outdir,
                     logger):
    assert os.path.isfile(querry_tmp)
    assert os.path.isfile(anno_tmp)
    querry_tmp_bed = pybedtools.BedTool(querry_tmp)
    a  = querry_tmp_bed.intersect(anno_tmp, wo=True)
    if a:
        b =  a.to_dataframe(names=['chrom',
                                   'start',
                                   'end',
                                   'querry',
                                   'length',
                                   'anno_chrom',
                                   'anno_start',
                                   'anno_end',
                                   'anno',
                                   'anno_length',
                                   'overlap'])
        b[anno_name +'.PO'] = b['overlap']/b['length']
        b[anno_name + '.RPO'] = b['overlap']/b['anno_length']
        req_cols = ['chrom',
                    'start',
                    'end',
                    anno_name +'.PO',
                    anno_name + '.RPO']
        out_df = b[req_cols]
        out_df  = out_df.groupby(['chrom', 'start', 'end'], as_index=False).mean()
        out_df = out_df.drop_duplicates()
        out_df[anno_name +'.PO'] = round(out_df[anno_name +'.PO'],3)
        out_df[anno_name +'.RPO'] = round(out_df[anno_name +'.RPO'],3)    
        outfile = outdir + "/" + anno_name+".PO.tsv"
        out_df.to_csv(outfile, sep='\t', index=False)
    else:
        logger.debug("bedtools intersect failed {0}".format(anno_tmp))

        
if __name__ == "__main__":
    main()
