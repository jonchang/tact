# fastMRCA functions

from __future__ import division

import functools
import logging
import math
import multiprocessing
import random
from time import time
logger = logging.getLogger(__name__)

from .lib import get_tip_labels 

global tree
global cores
global pool
global maxtax

def initialize(phy, max_singlethread_taxa=None, nproc=multiprocessing.cpu_count()):
    """
    Initialize the fastmrca singleton with a tree, tunes when the fastmrca
    parallel algorithm should kick in, and starts the multiprocessing pool.
    """
    global tree
    global cores
    global pool
    global maxtax
    tree = phy
    cores = 1
    pool = multiprocessing.Pool(processes=cores)
    maxtax = 1000000
    if not maxtax and cores == 1:
        maxtax = float('inf')
    if not maxtax:
        maxtax = autotune()

def autotune():
    """
    Figure how many tips a clade needs before multithreading
    beats singlethreaded math.
    """
    global tree
    global cores
    tn = tree.taxon_namespace
    ntax = max(cores * cores, 1024)
    while True:
        if ntax > len(tn):
            return len(tn)
        st = []
        mt = []
        for i in range(5):
            logger.debug("FastMRCA iteration={}, ntax={}".format(i, ntax))
            labels = [tx.label for tx in random.sample(tn, ntax)]
            start_time = time()
            tn.taxa_bitmask(labels=labels)
            st.append(time() - start_time)
            start_time = time()
            bitmask(labels)
            mt.append(time() - start_time)
        # Get median times
        st_s = sorted(st)[2]
        mt_s = sorted(mt)[2]
        if st_s - mt_s > 0.75:
            # Single-thread performs ~0.75s worse
            return ntax
        else:
            ntax = ntax * 4

def bitmask(labels):
    """
    Gets a bitmask for the taxa in `labels`, potentially in parallel.
    """
    global tree
    global pool
    global cores
    global maxtax
    tn = tree.taxon_namespace
    if len(labels) < maxtax:
        return tn.taxa_bitmask(labels=labels)
    f = functools.partial(fastmrca_getter, tn)
    full_bitmask = 0
    for res in pool.map(f, chunks(labels, int(math.ceil(len(labels) / cores))), chunksize=1):
        full_bitmask |= res
    return full_bitmask

def get(labels):
    """Pulls a MRCA node out for the taxa in `labels`."""
    global tree
    labels = set(labels)
    mrca = tree.mrca(leafset_bitmask=bitmask(labels))
    if not mrca:
        return None
    if mrca and labels.issuperset(get_tip_labels(mrca)):
        return mrca

def cleanup():
    """Destroy the fastmrca multiprocessing pool."""
    pool.terminate()

def chunks(l, n):
    """Yield successive `n`-sized chunks from `l`."""
    l = list(l)
    for i in range(0, len(l), n):
        yield l[i:i + n]

def fastmrca_getter(tn, x):
    """Helper function for submitting stuff to the pool."""
    taxa = tn.get_taxa(labels=x)
    bitmask = 0
    for taxon in taxa:
        bitmask |= tn.taxon_bitmask(taxon)
    return bitmask
