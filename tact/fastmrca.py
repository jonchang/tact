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
global pool
global maxtax

def initialize(phy, max_singlethread_taxa=None, nproc=multiprocessing.cpu_count()):
    """
    Initialize the fastmrca singleton with a tree, tunes when the fastmrca
    parallel algorithm should kick in, and starts the multiprocessing pool.
    """
    global tree
    global pool
    global maxtax
    tree = phy
    maxtax = 1000000

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
    for res in map(f, labels):
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

def fastmrca_getter(tn, x):
    """Helper function for submitting stuff to the pool."""
    taxa = tn.get_taxa(labels=x)
    bitmask = 0
    for taxon in taxa:
        bitmask |= tn.taxon_bitmask(taxon)
    return bitmask
