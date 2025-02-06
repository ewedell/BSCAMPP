import json, time, sys, os, shutil
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter

from bscampp import get_logger, log_exception, __version__
from bscampp.configs import *
from bscampp.functions import *
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
    t0 = time.perf_counter()
    m = Manager(); lock = m.Lock()

    # parse command line arguments and build configurations
    parser, cmdline_args = parseArguments()

    # initialize multiprocessing (if needed)
    _LOG.warning('Initializing ProcessPoolExecutor...')
    pool = ProcessPoolExecutor(Configs.num_cpus, initializer=initial_pool,
            initargs=(parser, cmdline_args,))

    # (0) temporary files wrote to here
    workdir = os.path.join(Configs.outdir, f'TEMP{Configs.tmpfilenbr}')
    try:
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
    except OSError:
        log_exception(_LOG)

    # (1) read in tree, alignment, and separate reference sequences from
    # query sequences
    tree, aln_path, aln, qaln_path, qaln = readData(workdir)

    # (2) compute closest leaves for all query sequences
    query_votes_dict, query_top_vote_dict = getClosestLeaves(
            aln_path, qaln_path, aln, qaln, workdir)

    # shutdown pool
    _LOG.warning('Shutting down ProcessPoolExecutor...')
    pool.shutdown()
    _LOG.warning('ProcessPoolExecutor shut down.')
    
    # clean up temp files if not keeping
    #if not Configs.keeptemp:
    #    _LOG.info('Removing temporary files...')
    #    clean_temp_files()

    # stop BSCAMPP 
    send = time.perf_counter()
    _LOG.info('BSCAMPP completed in {} seconds...'.format(send - t0))

def clean_temp_files():
    # all temporary files/directories to remove
    temp_items = [f'TEMP{Configs.tmpfilenbr}']
    for temp in temp_items:
        temp_path = os.path.join(Configs.outdir, temp)
        if os.path.isfile(temp_path):
            os.remove(temp_path)
        elif os.path.isdir(temp_path):
            shutil.rmtree(temp_path)
        else:
            continue
        _LOG.info(f'\tRemoved {temp}')

def parseArguments():
    global _root_dir, main_config_path
    parser = _init_parser()
    cmdline_args = sys.argv[1:]
    
    # build config
    buildConfigs(parser, cmdline_args)
    _LOG.info('BSCAMPP is running with: {}'.format(
        ' '.join(cmdline_args)))
    getConfigs()

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
    parser.add_argument('-v', '--version', action='version',
            version="%(prog)s " + __version__)
    parser.groups = dict()

    # basic group
    basic_group = parser.add_argument_group(
            "Basic parameters".upper(),
            "These are the basic parameters for BSCAMPP.")
    parser.groups['basic_group'] = basic_group

    basic_group.add_argument('--placement-method', type=str,
                  help='The base placement method to use. Default: epa-ng',
                  choices=['epa-ng', 'pplacer'], default='epa-ng',
                  required=False)
    basic_group.add_argument("-i", "--info", "--info-path", type=str,
                  dest="info_path",
                  help=("Path to model parameters. E.g., .bestModel "
                  "from RAxML/RAxML-ng"),
                  required=True, default=None)
    basic_group.add_argument("-t", "--tree", "--tree-path", type=str,
                  dest="tree_path",
                  help="Path to reference tree with estimated branch lengths",
                  required=True, default=None)
    basic_group.add_argument("-a", "--alignment", "--aln-path", type=str,
                  dest="aln_path",
                  help=("Path for reference sequence alignment in "
                  "FASTA format. Optionally with query sequences. "
                  "Query alignment can be specified with --qaln-path"), 
                  required=True, default=None)
    basic_group.add_argument("-q", "--qalignment", "--qaln-path", type=str,
                  dest="qaln_path",
                  help=("Optionally provide path to query sequence alignment "
                  "in FASTA format. Default: None"),
                  required=False, default=None)
    basic_group.add_argument("-d", "--outdir", type=str,
                  help="Directory path for output. Default: bscampp_output/",
                  required=False, default="bscampp_output")
    basic_group.add_argument("-o", "--output", type=str, dest="outname",
                  help="Output file name. Default: bscampp_result.jplace",
                  required=False, default="bscampp_result.jplace")
    basic_group.add_argument("--threads", "--num-cpus", type=int,
                  dest="num_cpus",
                  help="Number of cores for parallelization, default: -1 (all)",
                  required=False, default=-1)

    # advanced parameter settings
    advance_group = parser.add_argument_group(
            "Advance parameters".upper(),
            ("These parameters control how BSCAMPP is run. "
             "The default values are set based on experiments."
             ))
    parser.groups['advance_group'] = advance_group

    advance_group.add_argument("-m", "--model", type=str,
                  help="Model used for edge distances. Default: GTR",
                  required=False, default="GTR")
    advance_group.add_argument("-b", "--subtreesize", type=int,
                  help="Integer size of the subtree. Default: 2000",
                  required=False, default=2000)
    advance_group.add_argument("-V", "--votes", type=int,
                  help="Number of votes per query sequence. Default: 5",
                  required=False, default=5)
    advance_group.add_argument("-s", "--similarityflag", type=str2bool,
                  help="boolean, True if maximizing sequence similarity "
                  "instead of simple Hamming distance (ignoring gap "
                  "sites in the query). Default: True",
                  required=False, default=True)
    
    # miscellaneous group
    misc_group = parser.add_argument_group(
            "Miscellaneous parameters".upper(),)
    parser.groups['misc_group'] = misc_group

    misc_group.add_argument("-n","--tmpfilenbr", type=int,
                  help="temporary file indexing. Default: 0",
                  required=False, default=0)
    misc_group.add_argument("-f", "--fragmentflag", type=str2bool,
                  help="If queries contains fragments. Default: True",
                  required=False, default=True)
    #misc_group.add_argument('--filter-min', type=int,
    #        help=('(EPA-ng setting) Minimum number of placements to report. '
    #              'Default: 1'),
    #        required=False, default=1)
    #misc_group.add_argument('--filter-max', type=int,
    #        help=('(EPA-ng setting) Maximum number of placements to report. '
    #              'Default: 7'),
    #        required=False, default=7)
    #misc_group.add_argument('--filter-acc-lwr', type=float,
    #        help=('(EPA-ng setting) Accumulated threshold to stop reporting '
    #        'placements. Default: None'),
    #        required=False, default=None)
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
