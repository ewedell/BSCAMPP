import os, shutil, subprocess, stat, re, traceback, shlex
from subprocess import Popen
from abc import abstractmethod

from bscampp import get_logger, log_exception
#from bscampp.configs import Configs

_LOG = get_logger(__name__)

'''
Template class Job for running external software/jobs
'''
class Job(object):
    def __init__(self):
        self.job_type = ""
        self.errors = []
        self.b_ignore_error = False
        self.pid = -1
        self.returncode = 0

    def __call__(self):
        return self.run()

    def get_pid(self):
        return self.pid

    # run the job with given invocation and raise errors when encountered
    def run(self, stdin="", lock=None, logging=None, shell=False):
        try:
            cmd, outpath = self.get_invocation()
            _LOG.debug(f'Running job_type: {self.job_type}, output: {outpath}')

            # failsafe for NotImplemented jobs
            if len(cmd) == 0:
                raise ValueError(
                        f'{self.job_type} does not have a valid run command. '
                        'It might be due to (invalid input type, etc.).')

            # identify binaries as the first field
            binpath = cmd[0]
            # deal with special cases, e.g., python, java
            if binpath == 'java':
                binpath = cmd[2]
            elif binpath == 'python' or binpath == 'python3':
                binpath = cmd[1]
            assert os.path.exists(binpath) or binpath == 'gzip', \
                    ('executable for %s does not exist: %s' %
                        (self.job_type, binpath))
            assert \
                (binpath.count('/')== 0 or os.path.exists(binpath)), \
                ('path for %s does not exist (%s)' %
                    (self.job_type, binpath))

            _LOG.debug('Arguments: %s', ' '.join(
                (str(x) if x is not None else '?NoneType?' for x in cmd)))

            # logging to local or to PIPE
            stderr, stdout = '', ''
            scmd = ' '.join(cmd)
            if logging != None:
                logpath = os.path.join(
                        os.path.dirname(outpath),
                        f'{logging}_{self.job_type}.txt')
                outlogging = open(logpath, 'w', 1)

                # TODO: may need to deal with piping in the future, for now
                # it is not needed
                p = Popen(cmd, text=True, bufsize=1,
                        stdin=subprocess.PIPE,
                        stdout=outlogging, stderr=outlogging)
                self.pid = p.pid
                stdout, stderr = p.communicate(input=stdin)
                # stdout and stderr are both written to outlogging
                # hence, assign them to be empty strings
                stdout, stderr = '', ''
                outlogging.close()
            else:
                p = Popen(cmd, text=True, bufsize=1,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.pid = p.pid
                stdout, stderr = p.communicate(input=stdin)
            self.returncode = p.returncode

            # successful run and write to log
            if self.returncode == 0:
                if lock:
                    try:
                        lock.acquire()
                        _LOG.debug(f'{self.job_type} completed, output: {outpath}')
                    finally:
                        lock.release()
                else:
                    _LOG.debug(f'{self.job_type} completed, output: {outpath}')
                return outpath
            else:
                error_msg = ' '.join([f'Error occurred running {self.job_type}.',
                    f'returncode: {self.returncode}'])
                if logging != None:
                    logpath = '\nLOGPATH: ' + os.path.join(
                            os.path.dirname(outpath),
                            f'{logging}_{self.job_type}.txt')
                else:
                    logpath = ''
                if lock:
                    try:
                        lock.acquire()
                        _LOG.error(error_msg + '\nSTDOUT: ' + stdout +
                                '\nSTDERR: ' + stderr + logpath)
                    finally:
                        lock.release()
                else:
                    _LOG.error(error_msg + '\nSTDOUT: ' + stdout +
                            '\nSTDERR: ' + stderr + logpath)
                exit(self.returncode)
        except Exception:
            log_exception(_LOG)

    # implemented in subclass
    # return: (cmd, outpath)
    @abstractmethod
    def get_invocation(self):
        raise NotImplementedError(
            'get_invocation() should be implemented by subclasses.')

'''
Generic job that runs the given command, represented as a list of strings
'''
class GenericJob(Job):
    def __init__(self, cmd=[], job_type='external'):
        Job.__init__(self)
        self.job_type = job_type
        self.cmd = cmd

    def get_invocation(self):
        return self.cmd, None

'''
A EPA-ng job that runs EPA-ng with given parameters
'''
class EPAngJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'epa-ng'

        self.path = ''
        self.info_path = ''
        self.tree_path = ''
        self.aln_path = ''
        self.qaln_path = ''
        self.outdir = ''
        #self.molecule = ''
        #self.model = ''
        self.num_cpus = 1

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'epa_result.jplace') 
        cmd = [self.path,
                '-m', self.info_path,
                '-t', self.tree_path, '-s', self.aln_path,
                '-q', self.qaln_path, '-w', self.outdir,
                '-T', str(self.num_cpus), '--redo']
        return cmd, self.outpath

'''
A taxtastic job that create a refpkg based on given parameters
'''
class TaxtasticJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'taxit'

        self.path = ''
        self.outdir = ''
        self.name = ''
        self.tree_path = ''
        self.aln_path = ''
        self.info_path = ''

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir)
        cmd = [self.path, 'create', '-P', self.outdir,
                '-l', self.name, '--aln-fasta', self.aln_path,
                '--tree-file', self.tree_path]
        # check which model file is provided
        if 'bestModel' in self.info_path:
            cmd.extend(['--model-file', self.info_path])
        else:
            cmd.extend(['--tree-stats', self.info_path])
        return cmd, self.outpath

'''
A pplacer job that uses taxtastic refpkg to place sequences
'''
class PplacerTaxtasticJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'pplacer'

        self.path = ''
        self.refpkg_dir = ''
        self.qaln_path = ''
        self.outdir = ''
        self.outpath = ''
        #self.model = 'GTR'
        self.num_cpus = 1

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        # outpath defined 
        cmd = [self.path,
                #'-m', self.model,
                '-c', self.refpkg_dir, '-o', self.outpath,
                '-j', str(self.num_cpus), self.qaln_path]
        return cmd, self.outpath
