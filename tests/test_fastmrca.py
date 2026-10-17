"""Parity tests for the label index in `tact.fastmrca` with DendroPy."""

import os
import random

import dendropy
import numpy as np
import pytest
from click.testing import CliRunner

from tact import cli_add_taxa, fastmrca
from tact.cli_taxonomy import build_taxonomic_tree
from tact.tree_util import get_monophyletic_node, get_tip_labels, update_tree_view

DATASETS = ["weirdness", "intrusion", "short_branch", "stem", "stem2"]


def dendropy_bitmask(labels):
    """DendroPy reference implementation."""
    return fastmrca.tree.taxon_namespace.taxa_bitmask(labels=labels)


def load_pair(datadir, stem):
    """Load a taxonomy and backbone the way `tact_add_taxa` does, sharing one taxon namespace."""
    taxonomy = dendropy.Tree.get(
        path=os.path.join(datadir, stem + ".taxonomy.tre"), schema="newick", rooting="default-rooted"
    )
    backbone = dendropy.Tree.get(
        path=os.path.join(datadir, stem + ".backbone.tre"),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=taxonomy.taxon_namespace,
    )
    update_tree_view(backbone)
    return taxonomy, backbone


def label_sets(taxonomy, rng):
    """Every rank to be checked plus random subsets."""
    for node in taxonomy.preorder_internal_node_iter():
        yield get_tip_labels(node)
    all_labels = sorted(get_tip_labels(taxonomy))
    for _ in range(50):
        yield rng.sample(all_labels, rng.randint(1, len(all_labels)))


@pytest.mark.parametrize("stem", DATASETS)
def test_bitmask_matches_dendropy(datadir, stem):
    taxonomy, backbone = load_pair(datadir, stem)
    fastmrca.initialize(backbone)
    for labels in label_sets(taxonomy, random.Random(stem)):
        assert fastmrca.bitmask(labels) == dendropy_bitmask(labels)
        # any iterable of labels is accepted
        assert fastmrca.bitmask(iter(list(labels))) == dendropy_bitmask(list(labels))


@pytest.mark.parametrize("stem", DATASETS)
def test_get_matches_reference(datadir, stem):
    taxonomy, backbone = load_pair(datadir, stem)
    fastmrca.initialize(backbone)
    backbone_tips = get_tip_labels(backbone)
    checked = 0
    for labels in label_sets(taxonomy, random.Random(stem)):
        extant = backbone_tips.intersection(labels)
        if not extant:
            continue
        assert fastmrca.get(extant) is get_monophyletic_node(backbone, extant)
        checked += 1
    assert checked > 0


def test_bitmask_ignores_unknown_labels(datadir):
    _taxonomy, backbone = load_pair(datadir, "weirdness")
    fastmrca.initialize(backbone)
    known = sorted(get_tip_labels(backbone))[:3]
    with_unknown = [*known, "Not a taxon"]
    assert fastmrca.bitmask(with_unknown) == dendropy_bitmask(with_unknown) == fastmrca.bitmask(known)
    assert fastmrca.bitmask(["Not a taxon"]) == dendropy_bitmask(["Not a taxon"]) == 0
    assert fastmrca.bitmask([]) == 0


@pytest.mark.parametrize("is_case_sensitive", [False, True])
def test_bitmask_respects_namespace_case_sensitivity(is_case_sensitive):
    tn = dendropy.TaxonNamespace(is_case_sensitive=is_case_sensitive)
    tree = dendropy.Tree.get(
        data="((Aa_bb:1,Cc_dd:1):1,Ee_ff:2);",
        schema="newick",
        taxon_namespace=tn,
        case_sensitive_taxon_labels=is_case_sensitive,
    )
    update_tree_view(tree)
    fastmrca.initialize(tree)
    for labels in (["Aa bb"], ["aa bb"], ["AA BB", "Ee ff"], ["cc DD", "Cc dd"]):
        assert fastmrca.bitmask(labels) == dendropy_bitmask(labels)
    assert (fastmrca.bitmask(["aa bb"]) == 0) is is_case_sensitive


def test_bitmask_combines_taxa_sharing_a_label():
    tree = dendropy.Tree.get(data="((A:1,B:1):1,C:2);", schema="newick")
    tn = tree.taxon_namespace
    duplicate = tn.new_taxon("A")
    update_tree_view(tree)
    fastmrca.initialize(tree)
    assert fastmrca.bitmask(["A"]) == dendropy_bitmask(["A"])
    assert fastmrca.bitmask(["A"]) & tn.taxon_bitmask(duplicate)


def test_bitmask_follows_namespace_changes():
    tree = dendropy.Tree.get(data="((A:1,B:1):1,C:2);", schema="newick")
    tn = tree.taxon_namespace
    update_tree_view(tree)
    fastmrca.initialize(tree)
    assert fastmrca.bitmask(["D"]) == 0
    new_taxon = tn.new_taxon("D")
    assert fastmrca.bitmask(["D"]) == dendropy_bitmask(["D"]) == tn.taxon_bitmask(new_taxon)
    tn.remove_taxon(new_taxon)
    assert fastmrca.bitmask(["D"]) == dendropy_bitmask(["D"]) == 0


def test_requires_initialization(monkeypatch):
    monkeypatch.setattr(fastmrca, "tree", None)
    with pytest.raises(RuntimeError):
        fastmrca.bitmask(["A"])
    with pytest.raises(RuntimeError):
        fastmrca.get(["A"])


def run_add_taxa(backbone, taxonomy, output, seed):
    random.seed(seed)
    # SciPy's annealing fallback draws from NumPy's global state
    np.random.seed(seed)  # noqa: NPY002
    # rates are kept in a module-level dict that would otherwise leak between runs
    cli_add_taxa.mrca_rates.clear()
    handlers = list(cli_add_taxa.logger.handlers)
    try:
        result = CliRunner().invoke(
            cli_add_taxa.main,
            ["--backbone", str(backbone), "--taxonomy", str(taxonomy), "--output", str(output)],
            catch_exceptions=False,
        )
    finally:
        for handler in cli_add_taxa.logger.handlers[:]:
            if handler not in handlers:
                cli_add_taxa.logger.removeHandler(handler)
                handler.close()
    assert result.exit_code == 0, result.output
    return {ext: open(f"{output}.{ext}", "rb").read() for ext in ("newick.tre", "nexus.tre", "rates.csv")}


def assert_output_unchanged(monkeypatch, tmp_path, backbone, taxonomy, seed):
    """`tact_add_taxa` must write the same files with the label index as with DendroPy's lookup."""
    indexed = run_add_taxa(backbone, taxonomy, tmp_path / "indexed", seed)
    with monkeypatch.context() as patch:
        patch.setattr(fastmrca, "bitmask", dendropy_bitmask)
        reference = run_add_taxa(backbone, taxonomy, tmp_path / "reference", seed)
    assert indexed == reference
    # a different seed gives a different tree
    assert run_add_taxa(backbone, taxonomy, tmp_path / "other", seed + 1)["newick.tre"] != indexed["newick.tre"]


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("stem", DATASETS)
def test_add_taxa_output_unchanged(monkeypatch, datadir, tmp_path, stem, seed):
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    assert_output_unchanged(monkeypatch, tmp_path, backbone, taxonomy, seed)


def test_add_taxa_output_unchanged_carangaria(monkeypatch, datadir, tmp_path):
    examples = os.path.join(os.path.dirname(os.path.dirname(datadir)), "examples")
    taxonomy = tmp_path / "Carangaria.taxonomy.tre"
    build_taxonomic_tree(os.path.join(examples, "Carangaria.csv")).write_to_path(str(taxonomy), schema="newick")
    assert_output_unchanged(monkeypatch, tmp_path, os.path.join(examples, "Carangaria.tre"), taxonomy, 0)
