# -*- coding: utf-8 -*-

"""Module that provides a generic monophyly singleton data object."""

import re
from collections import defaultdict

from .tree_util import get_tip_labels
from .tree_util import get_monophyletic_node

_valid_monophyly = {}


def init(tree):
    if len(_valid_monophyly):
        return _valid_monophyly

    _genera_map = defaultdict(list)
    tips = get_tip_labels(tree)
    for tip in tips:
        genus, _ = re.split("[_ ]+", tip, maxsplit=1)
        _genera_map[genus].append(tip)

    for genus, species in _genera_map.items():
        node = get_monophyletic_node(tree, species)
        if node:
            _valid_monophyly[genus] = species
