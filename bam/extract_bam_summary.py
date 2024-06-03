#!/usr/bin/env python

import pysam
import argparse
import collections
import pandas as pd
from collections import defaultdict
import numpy as np

def main():
    args = parse_args()
    run(args.bam_file, args.bed_file, args.sample_name, args.out_dir)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest='bam_file',
                       help='bam file',
                       required=True)
    parser.add_argument('-s', dest='sample_name',
                       help='sample name',
                       required=True)
    parser.add_argument('-b', dest='bed_file', type=argparse.FileType('r'),
                        help="bed file", 
                        required=True)
    parser.add_argument('-o', dest='out_dir', 
                        help="out dir", 
                        required=True)
    args = parser.parse_args()
    return args


def run(bam_file, bed_file, sample_name, out_dir):
    bed_list = read_bed(bed_file)
    bam_bed_dict, raw_count_dict = index_out_bam(bam_file, bed_list)
    create_file = create_dataframe(bam_bed_dict, 
                                   sample_name, 
                                   out_dir, raw_count_dict)


def read_bed(bed_file):
    bed_file_list = []
    f = bed_file.readlines()
    for raw_line in f:
        value = raw_line.strip().split('\t')
        chrom = value[0].strip()
        start =  value[1].strip()
        stop = value[2].strip()
        string = [chrom, start, stop]
        bed_file_list.append(string)
    return bed_file_list

def generate_list(bed_list):
    for i in bed_list:
        yield i


def index_out_bam(bam_file, bed_list):
    bed_metrics_dict = defaultdict(list)
    read_depth_dict = {}
    for i in generate_list(bed_list):
        string = i[0] + "_" + i[1] + "_"+  i[2]
        samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
        count=0
        chrom = str(i[0]).strip()
        start = int(i[1])
        stop = int(i[2])
        for read in samfile.fetch(str(chrom), start, stop):
            count+=1
            mapq = read.mapping_quality
            if read.is_proper_pair:
                isize = read.template_length
                isize_dict = {'isize' : isize, 'proper_pair': 1, 'MAPQ' : mapq}
                tag_dict = dict(read.get_tags())
                try:
                    new_tag = {k: tag_dict[k] for k in ('AS', 'XS', 'MQ')}
                    new_tag_update = dict(list(isize_dict.items()) + list(new_tag.items()))
                    bed_metrics_dict[string].append(new_tag_update)
                except KeyError:
                    try:
                        new_tag = {k: tag_dict[k] for k in ('AS', 'XS')}
                        new_tag_update = dict(list(isize_dict.items()) + list(new_tag.items()))
                        bed_metrics_dict[string].append(new_tag_update)
                    except KeyError:
                        new_tag = {k: tag_dict[k] for k in ('AS')}
                        new_tag_update = dict(list(isize_dict.items()) + list(new_tag.items()))
                        bed_metrics_dict[string].append(new_tag_update)
            else:
                isize = read.template_length
                isize_dict = {'isize' : isize, 'proper_pair': 0, 'MAPQ' : mapq}
                tag_dict = dict(read.get_tags())
                try:
                    new_tag = {k: tag_dict[k] for k in ('AS', 'XS', 'MQ')}
                    new_tag_update = dict(list(isize_dict.items()) + list(new_tag.items()))
                    bed_metrics_dict[string].append(new_tag_update)
                except KeyError:
                    try:
                        new_tag = {k: tag_dict[k] for k in ('AS', 'XS')}
                        new_tag_update = dict(list(isize_dict.items()) + list(new_tag.items()))
                        bed_metrics_dict[string].append(new_tag_update)
                    except KeyError:
                        new_tag = {k: tag_dict[k] for k in ('AS')}
                        new_tag_update = dict(list(isize_dict.items()) + list(new_tag.items()))
                        bed_metrics_dict[string].appned(new_tag_update)
        read_depth_dict[string] = count
    return bed_metrics_dict, read_depth_dict


def modify_dict(bed_metrics_dict, read_depth_dict):
    new_merged_dict =  {}
    for k, v in bed_metrics_dict.items():
        new_v = defaultdict(list)
        for x in v:
            for p, q in x.items():
                new_v[p].append(q)
        less_than10 = len([x for x in new_v['MAPQ'] if x < 10])
        total_items = len(new_v['MAPQ'])
        percentage_less_than_10 = (less_than10 / total_items) * 100 if total_items != 0 else 0
        new_pct_lt10 = {'pct_count_mapq_Lt10':percentage_less_than_10}
        new_cnt_lt10 = {'count_mapq_Lt10':less_than10}
        new_median_v = {r +"_median": np.median(s)  for r, s in new_v.items()}
        new_mean_v = {r +"_mean": np.mean(s)  for r, s in new_v.items()}
        new_std_v = {r + "_std": np.std(s)  for r, s in new_v.items()}
        new_min_v = {r +"_min": np.min(s)  for r, s in new_v.items()}
        dp_dict = {'raw_dp': read_depth_dict[k]  if read_depth_dict[k] else True}
        result_dict = {**new_median_v, **new_mean_v, **new_std_v, **new_min_v, **new_pct_lt10, **new_cnt_lt10,**dp_dict}
        new_merged_dict[k] = result_dict
    print(new_merged_dict)
    return new_merged_dict


def normalize_list(input_list):
    min_val = min(input_list)
    max_val = max(input_list)
    if min_val == max_val:
        return [0.5] * len(input_list)  # Set all values to 0.5 (midpoint)

    normalized_list = [(x - min_val) / (max_val - min_val) for x in input_list]
    return normalized_list


#def calculate_mq_count(


def create_dataframe(bed_metrics_dict, 
                     sample_name, 
                     out_dir, 
                     read_dict):
    df = pd.DataFrame.from_dict((modify_dict(bed_metrics_dict, read_dict)), orient='index')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'probe'})
    sample =  sample_name.split('.')[0]
    df['SAMPLE_NAME'] = sample
    filename = out_dir + "/" + sample_name + '.metrics.csv'
    pd.DataFrame.to_csv(df, filename, index=False)



if __name__ == "__main__":
    main()
