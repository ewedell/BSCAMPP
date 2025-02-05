import sys, os, utils, shutil
import json
import time
import argparse
import treeswift
import copy
import multiprocessing as mp
import threading
import itertools
from collections import Counter

def parseArgs():
    # example usage
    example_usage = '''Example usages:
> default
    %(prog)s -i raxml.info 
'''

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--info", type=str,
                        help="Path to model parameters",
                        required=True, default=None)
    parser.add_argument("-t", "--tree", type=str,
                        help="Path to reference tree with estimated branch lengths",
                        required=True, default=None)
    parser.add_argument("-d", "--outdir", type=str,
                        help="Directory path for output",
                        required=True, default=None)
    parser.add_argument("-a", "--alignment", type=str,
                        help="Path for query and reference sequence alignment "
                        "in fasta format", required=True, default=None)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name",
                        required=False, default="EPA-ng-BSCAMPP")
    parser.add_argument("-m", "--model", type=str,
                        help="Model used for edge distances",
                        required=False, default="GTR")
    parser.add_argument("-b", "--subtreesize", type=int,
                        help="Integer size of the subtree",
                        required=False, default=2000)
    parser.add_argument("-V", "--votes", type=int,
                        help="Integer number of votes per query sequence",
                        required=False, default=5)
    parser.add_argument("-s", "--similarityflag", type=str2bool,
                        help="boolean, True if maximizing sequence similarity "
                        "instead of simple Hamming distance (ignoring gap "
                        "sites in the query)",
                        required=False, default=True)
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
    parser.add_argument("-q", "--qalignment", type=str,
                        help="Path to query sequence alignment in fasta format (ref alignment separate)",
                        required=False, default="")
    parser.add_argument("-f", "--fragmentflag", type=str2bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=True)
    parser.add_argument('--threads', type=int,
            help='number of threads for EPA-ng, default: all',
            required=False, default=-1)
    parser.add_argument('--placement-method', type=str,
            help='The base placement method to run, default: epa-ng',
            choices=['epa-ng', 'pplacer'], default='epa-ng',
            required=False)

def str2bool(b):
    if isinstance(b, bool):
       return b
    if b.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif b.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')        
