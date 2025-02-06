import os, time
try:
    import configparser
except ImportError:
    from ConfigParser import configparser
from argparse import ArgumentParser, Namespace
from bscampp.init_configs import init_config_file
from bscampp import get_logger, log_exception

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
    global _root_dir

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

    # binaries
    pplacer_path = None
    epang_path = None

    # placement settings
    placement_method = 'epa-ng'
    model = 'GTR'
    subtreesize = 2000
    votes = 5
    similarityflag = True

    # miscellaneous
    tmpfilenbr = 0
    fragmentflag = True

# check if the given configuration is valid to add
def set_valid_configuration(name, conf):
    if not isinstance(conf, Namespace):
        _LOG.warning(
            "Looking for Namespace object from \'{}\' but find {}".format(
                name, type(conf)))
        log_exception(_LOG)
        exit()

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

'''
Build Config class
'''
def buildConfigs(parser, cmdline_args, child_process=False, rerun=False):
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    args = parser.parse_args(cmdline_args)

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

    # modify num_cpus if it is the default value
    if Configs.num_cpus > 0:
        Configs.num_cpus = min(os.cpu_count(), Configs.num_cpus)
    else:
        Configs.num_cpus = os.cpu_count()
