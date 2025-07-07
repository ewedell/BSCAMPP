import os, time
try:
    import configparser
except ImportError:
    from ConfigParser import configparser
from argparse import ArgumentParser, Namespace
import subprocess

from bscampp.init_configs import init_config_file
from bscampp import get_logger, log_exception
#from bscampp.utils import inferDataType

# detect home.path or create if missing
homepath = os.path.dirname(__file__) + '/home.path'
_root_dir, main_config_path = init_config_file(homepath)

# set valid configparse section names
valid_config_sections = []
logging_levels = set(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])

_LOG = get_logger(__name__)

'''
Configuration defined by users and by default values
'''
class Configs:

    # basic input paths
    info_path = None        # info file for pplacer or EPA-ng
    tree_path = None        # placement tree path
    aln_path = None         # alignment for backbone. Optinally with queries
    qaln_path = None        # (optional) alignment for query.
    outdir = None           # output directory
    outname = None          # output name for the final jplace file
    keeptemp = False        # whether to keep all temporary files
    verbose = 'INFO'        # default verbose level to print
    num_cpus = 1            # number of cores to use for parallelization
    cpus_per_job = 2        # number of cores to use per job
    max_workers = 1         # max_workers for ProcessPoolExecutor
                            # ... = max(1, num_cpus // cpus_per_job)

    # binaries
    pplacer_path = None
    epang_path = None
    taxit_path = None
    hamming_distance_dir = None

    # placement settings
    placement_method = 'epa-ng'
    subtreesize = 2000
    votes = 5
    similarityflag = True

    # miscellaneous
    tmpfilenbr = 0
    fragmentflag = True
    subtreetype = 'd'

# check if the given configuration is valid to add
def set_valid_configuration(name, conf):
    if not isinstance(conf, Namespace):
        _LOG.warning(
            "Looking for Namespace object from \'{}\' but find {}".format(
                name, type(conf)))
        return

    # basic section defined in main.config
    if name == 'basic':
        for k in conf.__dict__.keys():
            k_attr = getattr(conf, k)
            if not k_attr:
                continue
            if k in Configs.__dict__:
                setattr(Configs, k, k_attr)
    else:
        pass

# valid attribute check for print out
def valid_attribute(k, v):
    if not isinstance(k, str):
        return False
    if k.startswith('_'):
        return False
    return True

# print out current configuration
def getConfigs():
    msg = '\n************ Configurations ************\n' + \
            f'\thome.path: {homepath}\n' + \
            f'\tmain.config: {main_config_path}\n\n'
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            msg += f'\tConfigs.{k}: {v}\n'
    print(msg, flush=True)

# read in config file if it exists
def _read_config_file(filename, cparser, opts,
        child_process=False, expand=None):
    config_defaults = []
    with open(filename, 'r') as f:
        cparser.read_file(f)
        if cparser.has_section('commandline'):
            for k, v in cparser.items('commandline'):
                config_defaults.append(f'--{k}')
                config_defaults.append(v)

        for section in cparser.sections():
            if section == 'commandline':
                continue
            if getattr(opts, section, None):
                section_name_space = getattr(opts, section)
            else:
                section_name_space = Namespace()
            for k, v in cparser.items(section):
                if expand and k == 'path':
                    v = os.path.join(expand, v)
                setattr(section_name_space, k, v)
            setattr(opts, section, section_name_space)
    return config_defaults

# validate a given binary file is executable
def validate_binary_executable(name, path):
    # 1. make sure path exists
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} is not found!")

    # 2. make sure path is executable with no IO error
    try:
        p = subprocess.Popen([path], text=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        _, __ = p.communicate()
    except OSError:
        raise AssertionError(' '.join([
            f"{path} is not executable! Please do the followings:",
            f"\n\n\t1. obtain a working executable for the software: {name}",
            f"\n\t2. modify {main_config_path} to replace the old executable path\n",
            ]))

'''
Build Config class
'''
def buildConfigs(parser, cmdline_args, child_process=False, rerun=False):
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    args = parser.parse_args(cmdline_args)

    # Check if only updating config files, if so, re-initialize the 
    # configuration file at ~/.bscampp/main.config and exit
    #if args.command == 'update-configs':
    #    _ = init_config_file(homepath, rerun=True)
    #    _LOG.warning('Finished re-initializing the configuration file '
    #            f'at {main_config_path}, exiting...')
    #    exit(0)

    # first load arguments from main.configs
    main_args = Namespace()
    cmdline_main = _read_config_file(main_config_path,
            cparser, main_args, child_process=child_process)

    # merge arguments, in the correct order so things are overridden correctly
    args = parser.parse_args(cmdline_main + cmdline_args,
            namespace=main_args)

    # directly add all arguments that's defined in the Configs class
    for k in args.__dict__.keys():
        k_attr = getattr(args, k)
        if k in Configs.__dict__:
            # valid argument that's defined in the Configs class
            setattr(Configs, k, k_attr)
        else:
            # check if the argument is valid
            set_valid_configuration(k, k_attr)

    # create outdir
    if not os.path.isdir(Configs.outdir):
        os.makedirs(Configs.outdir)

    # modify outname if it does not have a .jplace suffix
    if Configs.outname.split('.')[-1].lower() != 'jplace':
        Configs.outname += '.jplace'

    # modify num_cpus if it is the default value
    if Configs.num_cpus > 0:
        Configs.num_cpus = min(os.cpu_count(), Configs.num_cpus)
    else:
        Configs.num_cpus = os.cpu_count()
    # compute max_workers based on num_cpus and cpus_per_job
    Configs.max_workers = max(1, Configs.num_cpus // Configs.cpus_per_job)

    # sanity check for existence of base placement binary path
    if Configs.placement_method == 'epa-ng':
        validate_binary_executable(Configs.placement_method, Configs.epang_path)
        #assert os.path.exists(Configs.epang_path), 'epa-ng not detected!'
    elif Configs.placement_method == 'pplacer':
        validate_binary_executable(Configs.placement_method, Configs.pplacer_path)
        #assert os.path.exists(Configs.pplacer_path), 'pplacer not detected!'
