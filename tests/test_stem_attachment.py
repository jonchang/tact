import pytest

from dendropy import Tree

@pytest.mark.parametrize('execution_number', range(10))
def test_stem_attachment(script_runner, execution_number):
    backbone = "backbone.tre"
    with open(backbone, "w") as f:
        f.write("((a1:0.2199330096,a2:0.2199330096):0.4707234795,b1:0.6906564891);")
    taxonomy = "taxonomy.tre"
    with open(taxonomy, "w") as f:
        f.write("((a1,a2,a3,a4,a5,a6,a7,a8,a9)cladeA,(b1,b2,b3,b4,b5,b6,b7,b8,b9)cladeB)everything;")
    output = "res.newick.tre"
    result = script_runner.run("tact_add_taxa", "--taxonomy", taxonomy, "--backbone", backbone, "--output", "res", "-vvv")
    assert result.returncode == 0
    with open(output, "rb") as f:
        tacted = Tree.get(path=output, schema="newick")

    # is clade A monophyletic?
    labA = set(["a" + str(ii) for ii in range(1, 10)])
    mrcaA = tacted.mrca(taxon_labels=labA)
    newA = set([x.taxon.label for x in mrcaA.leaf_nodes()])
    assert labA == newA

    # is clade B monophyletic?
    labB = set(["b" + str(ii) for ii in range(1, 10)])
    mrcaB = tacted.mrca(taxon_labels=labB)
    newB = set([x.taxon.label for x in mrcaB.leaf_nodes()])
    assert labB == newB
