# fastMRCA functions

from __future__ import division

from .lib import get_tip_labels

import functools
import logging

logger = logging.getLogger(__name__)

global tree


def initialize(phy):
    """
    Initialize the fastmrca singleton with a tree.
    """
    global tree
    tree = phy


def bitmask(labels):
    """
    Gets a bitmask for the taxa in `labels`, potentially in parallel.
    """
    global tree
    tn = tree.taxon_namespace
    return tn.taxa_bitmask(labels=labels)


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
    """Helper function for submitting stuff."""
    taxa = tn.get_taxa(labels=x)
    bitmask = 0
    for taxon in taxa:
        bitmask |= tn.taxon_bitmask(taxon)
    return bitmask
