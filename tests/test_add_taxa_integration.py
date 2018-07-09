import pytest
import sys

from dendropy import Tree

@pytest.mark.parametrize('execution_number', range(10))
@pytest.mark.script_launch_mode('subprocess')
def test_stem_attachment(script_runner, execution_number):
    backbone = "backbone.tre"
    with open(backbone, "w") as f:
        f.write("((a1:0.2,a2:0.2):0.4,b1:0.6);")
    taxonomy = "taxonomy.tre"
    with open(taxonomy, "w") as f:
        f.write("((a1,a2,a3,a4,a5,a6,a7,a8,a9)cladeA,(b1,b2,b3,b4,b5,b6,b7,b8,b9)cladeB)everything;")
    output = "res.newick.tre"
    result = script_runner.run("tact_add_taxa", "--taxonomy", taxonomy, "--backbone", backbone, "--output", "res", "-vvv", "--cores=1")
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

@pytest.mark.parametrize('execution_number', range(10))
@pytest.mark.script_launch_mode('subprocess')
def test_clade_intrusion(script_runner, execution_number):
    backbone = "backbone.tre"
    with open(backbone, "w") as f:
        f.write("(Citharinus_congicus:119.251137,(Xenocharax_spilurus:89.531635,(Nannaethiops_unitaeniatus:57.887979,(Hemistichodus_vaillanti:40.503843,((Nannocharax_ansorgii:22.822841,(Nannocharax_elongatus:16.095353,Nannocharax_intermedius:16.095353):6.727488):8.678133,((Phago_boulengeri:8.133028,Microstomatichthyoborus_bashforddeani:8.133028):20.486394,((Paradistichodus_dimidiatus:20.792505,(Distichodus_sexfasciatus:12.26958,Distichodus_fasciolatus:12.26958):8.522925):1.762066,(Distichodus_hypostomatus:19.305748,(Distichodus_notospilus:11.30667,(Distichodus_decemmaculatus:7.311092,Distichodus_affinis:7.311092):3.995579):7.999078):3.248822):6.064851):2.881552):9.002869):17.384136):31.643657):29.719502);")
    taxonomy = "taxonomy.tre"
    with open(taxonomy, "w") as f:
        f.write("(((Citharidium_ansorgii)Citharidium,(Citharinops_distichodoides)Citharinops,(Citharinus_citharus,Citharinus_congicus,Citharinus_eburneensis,Citharinus_gibbosus,Citharinus_latus,Citharinus_macrolepis)Citharinus)Citharinidae,((Belonophago_hutsebouti,Belonophago_tinanti)Belonophago,(Congocharax_olbrechtsi,Congocharax_spilotaenia)Congocharax,(Distichodus_affinis,Distichodus_altus,Distichodus_antonii,Distichodus_atroventralis,Distichodus_brevipinnis,Distichodus_decemmaculatus,Distichodus_engycephalus,Distichodus_fasciolatus,Distichodus_hypostomatus,Distichodus_kolleri,Distichodus_langi,Distichodus_lusosso,Distichodus_maculatus,Distichodus_mossambicus,Distichodus_nefasch,Distichodus_noboli,Distichodus_notospilus,Distichodus_petersii,Distichodus_rostratus,Distichodus_rufigiensis,Distichodus_schenga,Distichodus_sexfasciatus,Distichodus_teugelsi)Distichodus,Dundocharax_bidentatus,(Eugnathichthys_eetveldii,Eugnathichthys_macroterolepis,Eugnathichthys_virgatus)Eugnathichthys,(Hemigrammocharax_angolensis,Hemigrammocharax_lineostriatus,Hemigrammocharax_machadoi,Hemigrammocharax_minutus,Hemigrammocharax_monardi,Hemigrammocharax_multifasciatus,Hemigrammocharax_ocellicauda,Hemigrammocharax_rubensteini,Hemigrammocharax_uniocellatus,Hemigrammocharax_wittei)Hemigrammocharax,(Hemistichodus_lootensi,Hemistichodus_mesmaekersi,Hemistichodus_vaillanti)Hemistichodus,(Ichthyborus_besse,Ichthyborus_congolensis,Ichthyborus_monodi,Ichthyborus_ornatus,Ichthyborus_quadrilineatus)Ichthyborus,Mesoborus_crocodilus,(Microstomatichthyoborus_bashforddeani,Microstomatichthyoborus_katangae)Microstomatichthyoborus,(Nannaethiops_bleheri,Nannaethiops_unitaeniatus)Nannaethiops,(Nannocharax_altus,Nannocharax_ansorgii,Nannocharax_brevis,Nannocharax_elongatus,Nannocharax_fasciatus,Nannocharax_fasciolaris,Nannocharax_gracilis,Nannocharax_hollyi,Nannocharax_intermedius,Nannocharax_latifasciatus,Nannocharax_lineomaculatus,Nannocharax_luapulae,Nannocharax_macropterus,Nannocharax_maculicauda,Nannocharax_micros,Nannocharax_niloticus,Nannocharax_occidentalis,Nannocharax_ogoensis,Nannocharax_parvus,Nannocharax_procatopus,Nannocharax_pteron,Nannocharax_reidi,Nannocharax_rubrolabiatus,Nannocharax_schoutedeni,Nannocharax_signifer,Nannocharax_taenia,Nannocharax_usongo,Nannocharax_zebra)Nannocharax,(Neolebias_ansorgii,Neolebias_axelrodi,Neolebias_gossei,Neolebias_gracilis,Neolebias_kerguennae,Neolebias_lozii,Neolebias_philippei,Neolebias_powelli,Neolebias_trewavasae,Neolebias_trilineatus,Neolebias_unifasciatus)Neolebias,Paradistichodus_dimidiatus,Paraphago_rostratus,(Phago_boulengeri,Phago_intermedius,Phago_loricatus)Phago,Xenocharax_spilurus)Distichodontidae)Characiformes;")
    output = "res.newick.tre"
    result = script_runner.run("tact_add_taxa", "--taxonomy", taxonomy, "--backbone", backbone, "--output", "res", "-vvv", "--cores=1")
    assert result.returncode == 0
    with open(output, "rb") as f:
        tacted = Tree.get(path=output, schema="newick")

    ss = tacted.as_ascii_plot()
    sys.stderr.write(ss)

    clade = set([x.replace("_", " ") for x in "Citharinus_citharus Citharinus_congicus Citharinus_eburneensis Citharinus_gibbosus Citharinus_latus Citharinus_macrolepis".split()])
    mrcaA = tacted.mrca(taxon_labels=clade)
    newA = set([x.taxon.label for x in mrcaA.leaf_nodes()])
    assert clade == newA


@pytest.mark.parametrize('execution_number', range(10))
@pytest.mark.script_launch_mode('subprocess')
def test_short_branch_problems(script_runner, execution_number):
    backbone = "backbone.tre"
    with open(backbone, "w") as f:
        f.write("(((((Coccotropsis_gymnoderma:24.102991,Synanceia_verrucosa:24.102991):2.782215,Aetapcus_maculatus:26.885206):1.989681,(Minous_trachycephalus:17.454997,Inimicus_didactylus:17.454997):11.41989):1.964048,(Ocosia_zaspilota:28.128453,Paracentropogon_rubripinnis:28.128453):2.710482):1.919125,Aploactisoma_milesii:32.75806);")
    taxonomy = "taxonomy.tre"
    with open(taxonomy, "w") as f:
        f.write("((((Acanthosphex_leurynnis)Acanthosphex,(Adventor_elongatus)Adventor,(Aploactis_aspera)Aploactis,(Aploactisoma_milesii)Aploactisoma,(Bathyaploactis_curtisensis,Bathyaploactis_ornatissima)Bathyaploactis,(Cocotropus_altipinnis,Cocotropus_astakhovi,Cocotropus_dermacanthus,Cocotropus_echinatus,Cocotropus_eksae,Cocotropus_izuensis,Cocotropus_keramaensis,Cocotropus_larvatus,Cocotropus_masudai,Cocotropus_microps,Cocotropus_monacanthus,Cocotropus_possi,Cocotropus_richeri,Cocotropus_roseomaculatus,Cocotropus_roseus,Cocotropus_steinitzi)Cocotropus,(Erisphex_aniarus,Erisphex_philippinus,Erisphex_pottii,Erisphex_simplex)Erisphex,(Kanekonia_florida,Kanekonia_pelta,Kanekonia_queenslandica)Kanekonia,(Matsubarichthys_inusitatus)Matsubarichthys,(Neoaploactis_tridorsalis)Neoaploactis,(Paraploactis_hongkongiensis,Paraploactis_intonsa,Paraploactis_kagoshimensis,Paraploactis_obbesi,Paraploactis_pulvinus,Paraploactis_taprobanensis,Paraploactis_trachyderma)Paraploactis,(Peristrominous_dolosus)Peristrominous,(Prosoproctus_pataecus)Prosoproctus,(Pseudopataecus_carnatobarbatus,Pseudopataecus_taenianotus)Pseudopataecus,(Ptarmus_gallus,Ptarmus_jubatus)Ptarmus,(Sthenopus_mollis)Sthenopus,(Xenaploactis_anopta,Xenaploactis_asperrima,Xenaploactis_cautes)Xenaploactis)Aploactinidae,((Choridactylus_lineatus,Choridactylus_multibarbus,Choridactylus_natalensis,Choridactylus_striatus)Choridactylus,(Dampierosa_daruma)Dampierosa,(Erosa_erosa)Erosa,(Inimicus_brachyrhynchus,Inimicus_caledonicus,Inimicus_cuvieri,Inimicus_didactylus,Inimicus_filamentosus,Inimicus_gruzovi,Inimicus_japonicus,Inimicus_joubini,Inimicus_sinensis,Inimicus_smirnovi)Inimicus,(Leptosynanceia_asteroblepa)Leptosynanceia,(Minous_andriashevi,Minous_coccineus,Minous_dempsterae,Minous_inermis,Minous_longimanus,Minous_monodactylus,Minous_pictus,Minous_pusillus,Minous_quincarinatus,Minous_trachycephalus,Minous_usachevi,Minous_versicolor)Minous,(Pseudosynanceia_melanostigma)Pseudosynanceia,(Synanceia_alula,Synanceia_horrida,Synanceia_nana,Synanceia_platyrhyncha,Synanceia_verrucosa)Synanceia,(Trachicephalus_uranoscopus)Trachicephalus)Synanceiidae,((Aetapcus_maculatus)Aetapcus,(Neopataecus_waterhousii)Neopataecus,(Pataecus_fronto)Pataecus)Pataecidae,(((Ablabys_binotatus,Ablabys_macracanthus,Ablabys_taenianotus)Ablabys,(Centropogon_australis,Centropogon_latifrons,Centropogon_marmoratus)Centropogon,(Coccotropsis_gymnoderma)Coccotropsis,(Cottapistus_cottoides)Cottapistus,(Glyptauchen_panduratus)Glyptauchen,(Gymnapistes_marmoratus)Gymnapistes,(Liocranium_pleurostigma,Liocranium_praepositum)Liocranium,(Neocentropogon_aeglefinus,Neocentropogon_affinis,Neocentropogon_japonicus,Neocentropogon_mesedai,Neocentropogon_profundus,Neocentropogon_trimaculatus)Neocentropogon,(Neovespicula_depressifrons)Neovespicula,(Notesthes_robusta)Notesthes,(Ocosia_apia,Ocosia_fasciata,Ocosia_possi,Ocosia_ramaraoi,Ocosia_spinosa,Ocosia_vespa,Ocosia_zaspilota)Ocosia,(Paracentropogon_longispinis,Paracentropogon_rubripinnis,Paracentropogon_vespa,Paracentropogon_zonatus)Paracentropogon,(Pseudovespicula_dracaena)Pseudovespicula,(Richardsonichthys_leucogaster)Richardsonichthys,(Snyderina_guentheri,Snyderina_yamanokami)Snyderina,(Tetraroge_barbata,Tetraroge_niger)Tetraroge,(Vespicula_bottae,Vespicula_cypho,Vespicula_trachinoides,Vespicula_zollingeri)Vespicula)Tetrarogidae))Scopan);")
    output = "res.newick.tre"
    result = script_runner.run("tact_add_taxa", "--taxonomy", taxonomy, "--backbone", backbone, "--output", "res", "-vvv", "--cores=1")
    assert result.returncode == 0
    with open(output, "rb") as f:
        tacted = Tree.get(path=output, schema="newick")

    n_short = 0
    for leaf in tacted.leaf_node_iter():
        if leaf.edge.length < 0.5:
            n_short += 1
    assert n_short < 10
