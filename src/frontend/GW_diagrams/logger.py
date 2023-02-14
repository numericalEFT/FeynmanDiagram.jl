#!/usr/bin/env python
import logging
import signal
import os
import sys
import inspect
from color import *
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
# add parentdir into PYTHONPATH, where IO module can be found
#from memory_profiler import profile
workspace = parentdir

log = logging.getLogger()
log.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
fh = logging.FileHandler(os.path.join(workspace, 'diagram.log'))

ch.setLevel(logging.INFO)
fh.setLevel(logging.INFO)
formatter = logging.Formatter(fmt="[calc][%(asctime)s][%(levelname)s]:\n%(message)s",
                              datefmt='%y/%m/%d %H:%M:%S')
ch.setFormatter(formatter)
fh.setFormatter(formatter)
log.addHandler(ch)
log.addHandler(fh)


def Assert(condition, info):
    if not condition:
        log.error(red(info))
        raise AssertionError


def Abort(info):
    log.error(red(info))
    raise AssertionError


class DelayedInterrupt(object):
    def __enter__(self):
        self.signal_received = False
        self.old_handler_int = signal.getsignal(signal.SIGINT)
        self.old_handler_term = signal.getsignal(signal.SIGTERM)
        signal.signal(signal.SIGINT, self.handler)
        signal.signal(signal.SIGTERM, self.handler)

    def handler(self, signal, frame):
        self.signal_received = (signal, frame)
        log.info('SIG {0} received. Delaying Interrupt...'.format(signal))

    def __exit__(self, type, value, traceback):
        signal.signal(signal.SIGINT, self.old_handler_int)
        signal.signal(signal.SIGTERM, self.old_handler_term)
        if self.signal_received:
            log.info('Interrupt successfully dealyed!')
            self.old_handler_int(*self.signal_received)


def memory_usage():
    import subprocess
    out = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())],
                           stdout=subprocess.PIPE).communicate()[0].split(b'\n')
    vsz_index = out[0].split().index(b'RSS')
    mem = float(out[1].split()[vsz_index]) / 1024
    return mem
