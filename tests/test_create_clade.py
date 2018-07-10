from __future__ import division

from hypothesis import given, assume
import hypothesis.strategies as st

from dendropy import TaxonNamespace

from tact.cli_add_taxa import create_clade, count_locked, edge_iter

@given(st.lists(
    st.tuples(
        st.from_regex(r"\A[A-Z][a-z]+ [a-z]+\Z"),
        st.floats(min_value=1e-3, max_value=1000, allow_infinity=False, allow_nan=False)
        ),
    min_size=1,
    unique_by=lambda x: x[0]
    ).map(dict))
def test_create_clade(data):
    spp = data.keys()
    ages = data.values()
    assume(len(set(ages)) == len(ages))
    tn = TaxonNamespace(spp, label="taxa")
    clade = create_clade(tn, spp, ages)
    xx = [x.label == "locked" for x in edge_iter(clade.seed_node)]
    cnt = sum(xx)
    tot = len(list(edge_iter(clade.seed_node)))
    assert(tot == cnt + 1)
