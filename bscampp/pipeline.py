import json, time, sys, os, shutil
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter
import argparse

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

    # set up a dry run if specified
    dry_run = False
    if 'dry_run' in kwargs and isinstance(kwargs['dry_run'], bool):
        dry_run = kwargs['dry_run']

    # parse command line arguments and build configurations
    parser, cmdline_args = parseArguments(dry_run=dry_run)

    # initialize multiprocessing (if needed)
    _LOG.warning('Initializing ProcessPoolExecutor...')
    # maximally concurrently run Configs.num_cpus // 2 jobs, each job
    # can use 2 threads
    pool = ProcessPoolExecutor(Configs.max_workers, initializer=initial_pool,
            initargs=(parser, cmdline_args,))

    # (0) temporary files wrote to here
    if not dry_run:
        workdir = os.path.join(Configs.outdir, f'tmp{Configs.tmpfilenbr}')
        try:
            if not os.path.isdir(workdir):
                os.makedirs(workdir)
        except OSError:
            log_exception(_LOG)
    else:
        workdir = os.getcwd()

    # (1) read in tree, alignment, and separate reference sequences from
    # query sequences
    tree, leaf_dict, aln_path, aln, qaln_path, qaln, qname_map, qname_map_rev \
            = readData(workdir, dry_run=dry_run)

    # (2) compute closest leaves for all query sequences
    query_votes_dict, query_top_vote_dict = getClosestLeaves(
            aln_path, qaln_path, aln, qaln, workdir, dry_run=dry_run)

    # (3) first assign all queries to their closest-leaf subtrees,
    # then do reassignment to minimize distance between each's top vote
    # and the subtree's seed leaf
    new_subtree_dict, placed_query_list = assignQueriesToSubtrees(
            query_votes_dict, query_top_vote_dict, tree, leaf_dict,
            dry_run=dry_run)

    # (4) perform placement for each subtree
    output_jplace = placeQueriesToSubtrees(tree, leaf_dict, new_subtree_dict,
            placed_query_list, aln, qaln, cmdline_args, workdir,
            qname_map, qname_map_rev,
            pool, lock, dry_run=dry_run)

    # (5) write the output jplace to local
    writeOutputJplace(output_jplace, dry_run=dry_run)

    # shutdown pool
    _LOG.warning('Shutting down ProcessPoolExecutor...')
    pool.shutdown()
    _LOG.warning('ProcessPoolExecutor shut down.')
    
    # clean up temp files if not keeping
    if not Configs.keeptemp:
        _LOG.info('Removing temporary files...')
        clean_temp_files()

    # stop BSCAMPP 
    send = time.perf_counter()
    _LOG.info('BSCAMPP completed in {} seconds...'.format(send - t0))

    if dry_run:
        return True
    else:
        return False


# main pipeline for SCAMPP
def scampp_pipeline(*args, **kwargs):
    t0 = time.perf_counter()
    m = Manager(); lock = m.Lock()

    # set up a dry run if specified
    dry_run = False
    if 'dry_run' in kwargs and isinstance(kwargs['dry_run'], bool):
        dry_run = kwargs['dry_run']

    # parse command line arguments and build configurations
    parser, cmdline_args = parseArguments(dry_run=dry_run, method="SCAMPP")

    # initialize multiprocessing (if needed)
    _LOG.warning('Initializing ProcessPoolExecutor...')
    pool = ProcessPoolExecutor(Configs.num_cpus, initializer=initial_pool,
            initargs=(parser, cmdline_args,))

    # (0) temporary files wrote to here
    if not dry_run:
        workdir = os.path.join(Configs.outdir, f'tmp{Configs.tmpfilenbr}')
        try:
            if not os.path.isdir(workdir):
                os.makedirs(workdir)
        except OSError:
            log_exception(_LOG)
    else:
        workdir = os.getcwd()

    # (1) read in tree, alignment, and separate reference sequences from
    # query sequences
    tree, leaf_dict, aln_path, aln, qaln_path, qaln, qname_map, qname_map_rev \
            = readData(workdir, dry_run=dry_run)

    # (2) compute closest leaves for all query sequences
    query_votes_dict, query_top_vote_dict = getClosestLeaves(
            aln_path, qaln_path, aln, qaln, workdir, dry_run=dry_run)

    # (3) first assign each query to the subtree built using the closest 
    # leaf as the seed sequence 
    new_subtree_dict, placed_query_list = buildQuerySubtrees(
            query_votes_dict, query_top_vote_dict, tree, leaf_dict,
            dry_run=dry_run)

    # (4) perform placement for each subtree
    output_jplace = placeQueriesToSubtrees(tree, leaf_dict, new_subtree_dict,
            placed_query_list, aln, qaln, cmdline_args, workdir,
            qname_map, qname_map_rev,
            pool, lock, dry_run=dry_run)

    # (5) write the output jplace to local
    writeOutputJplace(output_jplace, dry_run=dry_run)

    # shutdown pool
    _LOG.warning('Shutting down ProcessPoolExecutor...')
    pool.shutdown()
    _LOG.warning('ProcessPoolExecutor shut down.')
    
    # clean up temp files if not keeping
    if not Configs.keeptemp:
        _LOG.info('Removing temporary files...')
        clean_temp_files()

    # stop SCAMPP 
    send = time.perf_counter()
    _LOG.info('SCAMPP completed in {} seconds...'.format(send - t0))

    if dry_run:
        return True
    else:
        return False


def clean_temp_files():
    # all temporary files/directories to remove
    temp_items = [f'tmp{Configs.tmpfilenbr}']
    for temp in temp_items:
        temp_path = os.path.join(Configs.outdir, temp)
        if os.path.isfile(temp_path):
            os.remove(temp_path)
        elif os.path.isdir(temp_path):
            shutil.rmtree(temp_path)
        else:
            continue
        _LOG.info(f'- Removed {temp}')

def parseArguments(dry_run=False, method="BSCAMPP"):
    default_outdir = f"{method.lower()}_output"
    default_outname = f"{method.lower()}_result"

    parser = _init_parser(default_outdir=default_outdir,
            default_outname=default_outname)
    cmdline_args = sys.argv[1:]

    if dry_run:
        cmdline_args = ['-i', 'dummy.info', '-t', 'dummy.tre',
                '-a', 'dummy.fa'] 
    
    # build config
    buildConfigs(parser, cmdline_args)
    _LOG.info('{} is running with: {}'.format(method,
        ' '.join(cmdline_args)))
    getConfigs()

    return parser, cmdline_args

def _init_parser(default_outdir="bscampp_output",
        default_outname="bscampp_result"):
    # example usage
    example_usages = '''Example usages:
> (1) Default
    %(prog)s -i raxml.bestModel -t reference.tre -a alignment.fa  
> (2) Separate alignment file for query sequences
    %(prog)s -i raxml.bestModel -t reference.tre -a reference.fa -q query.fa
> (3) Use pplacer instead of EPA-ng as base method (need RAxML-ng info or FastTree log file)
    %(prog)s -i fasttree.log -t reference.tre -a alignment.fa --placement-method pplacer
'''

    parser = ArgumentParser(
            description=(
                "This program runs BSCAMPP/SCAMPP, a scalable phylogenetic "
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
    required = True

    ## add a subcommand for updating configuration file without running
    ## the BSCAMPP pipeline
    #subparsers = parser.add_subparsers(dest='command',
    #        help='Subcommands for BSCAMPP')
    #update_parser = subparsers.add_parser('update-configs',
    #        help='Update the configuration file without running BSCAMPP.')

    ## try update args requirement if subcommand(s) are used
    #if 'update-configs' in sys.argv:
    #    required = False

    # basic group
    basic_group = parser.add_argument_group(
            "Basic parameters".upper(),
            "These are the basic parameters for BSCAMPP/SCAMPP.")
    parser.groups['basic_group'] = basic_group

    basic_group.add_argument('--placement-method', type=str,
                  help='The base placement method to use. Default: epa-ng',
                  choices=['epa-ng', 'pplacer'], default='epa-ng',
                  required=False)
    basic_group.add_argument("-i", "--info", "--info-path", type=str,
                  dest="info_path",
                  help=("Path to model parameters. E.g., .bestModel "
                  "from RAxML/RAxML-ng"),
                  required=required, default=None)
    basic_group.add_argument("-t", "--tree", "--tree-path", type=str,
                  dest="tree_path",
                  help="Path to reference tree with estimated branch lengths",
                  required=required, default=None)
    basic_group.add_argument("-a", "--alignment", "--aln-path", type=str,
                  dest="aln_path",
                  help=("Path for reference sequence alignment in "
                  "FASTA format (can be a .gz file). "
                  "Optionally with query sequences. "
                  "Query alignment can be specified with --qaln-path"), 
                  required=required, default=None)
    basic_group.add_argument("-q", "--qalignment", "--qaln-path", type=str,
                  dest="qaln_path",
                  help=("Optionally provide path to query sequence alignment "
                  "in FASTA format (can be a .gz file). Default: None"),
                  required=False, default=None)
    #basic_group.add_argument("--molecule", type=str,
    #              choices=['nucl', 'nucleotide', 'prot', 'protein'],
    #              help=("Specify nucleotide or protein sequences. "
    #              "Default: infer datatype"),
    #              required=False, default=None)
    basic_group.add_argument("-d", "--outdir", type=str,
                  help="Directory path for output. Default: bscampp_output/",
                  required=False, default=default_outdir)
    basic_group.add_argument("-o", "--output", type=str, dest="outname",
                  help="Output file name. Default: bscampp_result.jplace",
                  required=False, default=f"{default_outname}.jplace")
    basic_group.add_argument("--threads", "--num-cpus", type=int,
                  dest="num_cpus",
                  help="Number of cores for parallelization, default: -1 (all)",
                  required=False, default=-1)
    basic_group.add_argument("--cpus-per-job", type=int,
                  help="Number of cores to use for each job, default: 2",
                  required=False, default=2)

    # advanced parameter settings
    advance_group = parser.add_argument_group(
            "Advance parameters".upper(),
            ("These parameters control how BSCAMPP is run. "
             "The default values are set based on experiments."
             ))
    parser.groups['advance_group'] = advance_group

    #advance_group.add_argument("-m", "--model", type=str,
    #              help=("Model used for edge distances. EPA-ng will use the "
    #              "provided info_path (*.bestModel) for model. "
    #              "Default: GTR for nucleotide, LG for protein"),
    #              required=False, default=None)
    advance_group.add_argument("-b", "--subtreesize", type=int,
                  help="Integer size of the subtree. Default: 2000",
                  required=False, default=2000)
    advance_group.add_argument("-V", "--votes", type=int,
                  help="(BSCAMPP only) Number of votes per query sequence. "
                  "Default: 5",
                  required=False, default=5)
    advance_group.add_argument("--subtreetype", type=str,
                  help="(SCAMPP only) Options for collecting "
                  "nodes for the subtree - d for edge weighted "
                  "distances, n for node distances, h for Hamming "
                  "distances. Default: d",
                  required=False, default='d')
    advance_group.add_argument("--similarityflag", type=str2bool,
                  help="Boolean, True if maximizing sequence similarity "
                  "instead of simple Hamming distance (ignoring gap "
                  "sites in the query). Default: True",
                  required=False, default=True)
    
    # miscellaneous group
    misc_group = parser.add_argument_group(
            "Miscellaneous parameters".upper(),)
    parser.groups['misc_group'] = misc_group

    misc_group.add_argument("-n","--tmpfilenbr", type=int,
                  help="Temporary file indexing (e.g., tmp0/). Default: 0",
                  required=False, default=0)
    misc_group.add_argument("--fragmentflag", type=str2bool,
                  help="If queries contains fragments. Does not do anything "
                  "if similarity flag is set to True. Default: True",
                  required=False, default=True)
    misc_group.add_argument("--keeptemp", type=str2bool,
                  help="Boolean, True to keep all temporary files. "
                  "Default: False",
                  required=False, default=False)
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
