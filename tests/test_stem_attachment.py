import pytest

from dendropy import Tree

pytest_plugins = ["pytester"]

@pytest.fixture
def run(testdir):
    def do_run(*args):
        args = ["tact_add_taxa"] + list(args)
        return testdir._run(*args)
    return do_run

def test_pyconv_latin1_to_utf8(tmpdir, run):
    backbone = tmpdir.join("backbone.tre")
    with backbone.open("w") as f:
        f.write("((a1:0.2199330096,a2:0.2199330096):0.4707234795,b1:0.6906564891);")
    taxonomy = tmpdir.join("taxonomy.tre")
    with taxonomy.open("w") as f:
        f.write("((a1,a2)cladeA,(b1,b2,b3,b4,b5,b6,b7,b8,b9)cladeB)everything;")
    output = tmpdir.join("res.newick.tre")
    result = run("--taxonomy", taxonomy, "--backbone", backbone, "--output", "res")
    assert result.ret == 0
    with output.open("rb") as f:
        tacted = Tree.get(path=output, schema="newick")
    labs = set(["b" + str(ii) for ii in range(1, 10)])
    mrca = tacted.mrca(taxon_labels=labs)
    newlabs = set([x.label for x in mrca.leaf_nodes()])
    assert labs == newlabs
