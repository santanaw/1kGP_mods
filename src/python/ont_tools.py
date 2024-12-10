#!/usr/bin/env python

"""
@author: Walter Santana-Garcia
Tools to assist in the processing of files related to ONT datasets.
"""

import re
import os
import sys
import copy
import pysam
import random
import textwrap
import argparse
import subprocess
import numpy as np
import pandas as pd

VERSION = "0.1.2"
PROGRAM = "ont_tools"
AUTHOR  = "Walter Santana-Garcia"
CONTACT = "wsantana@ebi.ac.uk"



def filter_bam(input_f, output_f, min_len = None, min_qual = None, rm_dup = True):
    '''
    Filter reads in bam file
    
    Parameter
    ---------
    
    input_f : str 
        an input bam file
    output_f : str
        an output file name for filtered bam
        (*.pass.bam and *.fail.bam will be created)
    min_len : int
        minimum read legth to keep
    min_qual : float
        minimum mean read quality to keep
    rm_dup : bool
        remove reads with duplicated name

    Returns
    ---------
    : list
        the file name of the pass and fail bam
    '''


    prefix = os.path.splitext(output_f)[0]
    output_f = prefix + ".pass.bam"
    failed_f = prefix + ".fail.bam"
    
    read_ids = set()
    read_stats = {"total":0, "pass": 0, "fail_short": 0, "fail_lowqual": 0, "fail_duplicate": 0}
    
    sys.stdout.write("; Filtering bam file {0}\n".format(input_f))
    sys.stdout.write("; Using filters:\n")
    sys.stdout.write("; Minimal read length:{0}\n".format(min_len))
    sys.stdout.write("; Minimal read quality:{0}\n".format(min_qual))    
    sys.stdout.write("; Remove reads with same ID:{0}\n\n\n".format(rm_dup))

    # Fix for index error
    save = pysam.set_verbosity(0)
    
    with pysam.AlignmentFile(input_f, check_sq=False) as bam_fh:
        with pysam.AlignmentFile(output_f, 'wb', template=bam_fh, add_sam_header=True) as output_fh:
            with pysam.AlignmentFile(failed_f, 'wb', template=bam_fh, add_sam_header=True) as failed_fh:

                for read in bam_fh.fetch(until_eof=True):

                    if min_len and len(read.query_sequence) < min_len:
                        read_stats["fail_short"] += 1
                        failed_fh.write(read)

                    elif min_qual and np.mean(read.query_qualities) < min_qual:
                        read_stats["fail_lowqual"] += 1
                        failed_fh.write(read)
                    
                    elif rm_dup and read.query_name in read_ids:
                        read_stats["fail_duplicate"] += 1
                        failed_fh.write(read)
                        
                    else:
                        read_stats["pass"] += 1
                        output_fh.write(read)
                    
                    read_ids.add(read.query_name)
                    read_stats["total"]+=1

    # Set verbosity back
    pysam.set_verbosity(save)

    sys.stdout.write("; Filtering completed\n")
    sys.stdout.write("; Read statistics:\n")
    sys.stdout.write("; {0} total reads processed\n".format(read_stats["total"]))
    sys.stdout.write("; {0} reads retained\n".format(read_stats["pass"]))
    sys.stdout.write("; {0} short reads discarded\n".format(read_stats["fail_short"]))
    sys.stdout.write("; {0} low quality reads discarded\n".format(read_stats["fail_lowqual"]))    
    sys.stdout.write("; {0} duplicated reads discarded\n".format(read_stats["fail_duplicate"]))

    return [output_f,failed_f]

def main():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''
        Utils for processing ONT reads  
        -----------------------------------------------------------        
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='ont_tools ' + VERSION)
    subparsers = parser.add_subparsers(dest='util_tool', description=textwrap.dedent('''
        ont_tools.py is a toolkit for processing nanopore related files.
        
        ont_tools.py { filter_bam } -h
        -------------------------------------------------------
        '''))
    
    parser_filterbam = subparsers.add_parser('filter_bam', help="Filter reads in bam file")
    parser_filterbam.add_argument('-i', '--input_file', help='Input aligned/unaligned bam with read qualities.', required=True)
    parser_filterbam.add_argument('-o', '--output_file', help='Output bam file to save filtered reads.', required=True)
    parser_filterbam.add_argument('-ml', '--min_length', type=int, help='Minimum read legth to retain.', default=None)
    parser_filterbam.add_argument('-mq', '--min_quality', type=float,help='Minimum mean read quality to retain.', default=None)
    parser_filterbam.add_argument('-rd', '--rm_duplicates', action= "store_true", help='Remove duplicated reads in bam file.', default=False)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.util_tool == "filter_bam":
        input_f = args.input_file
        output_f = args.output_file
        min_len = args.min_length
        min_qual = args.min_quality
        rm_dup = args.rm_duplicates

        filter_bam(input_f, output_f, min_len, min_qual, rm_dup)

    else:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
if __name__ == "__main__":
    main()
    
