import os, time
try:
    import configparser
except ImportError:
    from ConfigParser import configparser
from argparse import ArgumentParser, Namespace
from bscampp.init_configs import init_config_file
from bscampp import get_logger

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
    info_path = None        # info file for pplacer
    tree_path = None
    aln_path = None
    qaln_path = None
    outdir = None
    outname = None
    keeptemp = False
    verbose = 'INFO'
