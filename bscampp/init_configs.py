import os, sys, shutil
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from argparse import ArgumentParser, Namespace
from platform import platform

def find_main_config(homepath):
    with open(homepath, 'r') as f:
        _root_dir = f.read().strip()
        main_config_path = os.path.join(_root_dir, 'main.config')
        if os.path.exists(main_config_path):
            return _root_dir, main_config_path
        else:
            return None, None

'''
Initialize the config file at ~/.bscampp/main.config
By default will prioritize existing software from the user environment
'''
def init_config_file(homepath, rerun=False, prioritize_user_software=True):
    if not rerun:
        # make sure home.path exists
        if not os.path.exists(homepath):
            print(f'Cannot find home.path: {homepath}, regenerating...')
        else:
            _root_dir, main_config_path = find_main_config(homepath)
            if _root_dir is not None:
                return _root_dir, main_config_path
            else:
                print(f'Cannot find main.config, regenerating...')

    _root_dir = os.path.expanduser('~/.bscampp')
    main_config_path = os.path.join(_root_dir, 'main.config')
    print(f'Initializing the config file at: {main_config_path}')

    # write to local home.path and _root_dir
    if not os.path.isdir(_root_dir):
        os.mkdir(_root_dir)
    with open(homepath, 'w') as f:
        f.write(_root_dir)

    # create main.config based on the default.config at this file's location
    _config_path = os.path.join(os.path.dirname(__file__), 'default.config')
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    assert os.path.exists(_config_path), \
            'default config file missing! Please redownload from GitHub\n'

    if os.path.exists(main_config_path):
        print(f'Main configuration file {main_config_path} exists...')
        print('Overwriting the existing config file...')
    print('\n')

    with open(_config_path, 'r') as f:
        cparser.read_file(f)

    # check platform, e.g., macOS or linux, etc.
    platform_name = platform()
    tools_dir = os.path.join(os.path.dirname(__file__), 'tools')

    # TODO
    return _root_dir, main_config_path
