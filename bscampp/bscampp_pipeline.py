import json, time, sys, os, shutil
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter

from bscampp import get_logger, __version__
from bscampp.configs import *
import bscampp.utils as utils

from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor

_LOG = get_logger(__name__)

# process pool initializer
def initial_pool(parser, cmdline_args):
    # avoid redundant logging for child process
    buildConfigs(parser, cmdline_args, child_process=True)

# main pipeline for BSCAMPP
def bscampp_pipeline(*args, **kwargs):
    s1 = time.time()
    m = Manager(); lock = m.Lock()

    # parse command line arguments and build configurations
    parser, cmdline_args = parseArguments()

    # initialize multiprocessing (if needed)
    _LOG.warning('Initializing ProcessPoolExecutor...')
    pool = ProcessPoolExecutor(Configs.num_cpus, initializer=initial_pool,
            initargs=(parser, cmdline_args,))

    # (1) TODO

    # shutdown pool
    _LOG.warning('Shutting down ProcessPoolExecutor...')
    pool.shutdown()
    _LOG.warning('ProcessPoolExecutor shut down.')
    
    # clean up temp files if not keeping
    if not Configs.keeptemp:
        _LOG.info('Removing temporary files...')
        clean_temp_files()

    # stop BSCAMPP 
    send = time.time()
    _LOG.info('BSCAMPP completed in {} seconds...'.format(send - s1))

'''
Parse arguments from commandline and config file
'''
def parseArguments():
    global _root_dir, main_config_path
    parser = _init_parser()
    cmdline_args = sys.argv[1:]
    
    # build config
    buildConfigs(parser, cmdline_args)
    _LOG.info('BSCAMPP is running with: {}'.format(
        ' '.join(cmdline_args)))

    return parser, cmdline_args

def _init_parser():
    # example usage
    example_usages = '''Example usages:
> default
    %(prog)s -i raxml.info 
'''

    parser = ArgumentParser(
            description=(
                "This program runs BSCAMPP, a scalable phylogenetic "
                "placement framework that scales EPA-ng/pplacer "
                "to very large tree placement."
                ),
            conflict_handler='resolve',
            epilog=example_usages,
            formatter_class=utils.SmartHelpFormatter,
            )
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
    return parser

def str2bool(b):
    if isinstance(b, bool):
       return b
    if b.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif b.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')        
