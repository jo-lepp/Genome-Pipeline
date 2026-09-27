process DREP {
    tag "${root}:${label}"
    label 'binning'
    conda "bioconda::drep"
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(genomes, stageAs: 'genomes_in/*'), val(ext), val(pub_subpath)

    output:
    tuple val(root), val(label), path("drep/dereplicated_genomes"), emit: derep_dir

    script:
    """
    dRep dereplicate drep -g genomes_in/*.${ext}
    echo "dereplication done (${label})"
    """
}

workflow DEREPLICATION {
    take:
    metabat_bins   // tuple(root, [fa files])          - from BINNING.out.metabat_bins
    semibin_bins   // tuple(root, model, [fa.gz files]) - from BINNING.out.semibin_bins

    main:
    ch_drep_metabat = metabat_bins.map { root, fa ->
        tuple(root, 'metabat', fa, 'fa', "${root}_metabat/drep")
    }
    ch_drep_semibin = semibin_bins.map { root, model, fagz ->
        tuple(root, model, fagz, 'fa.gz', "${root}_semibin/${model}_model/drep")
    }

    ch_drep_in = ch_drep_metabat.mix(ch_drep_semibin)
    DREP(ch_drep_in)

    emit:
    derep_dir = DREP.out.derep_dir   // tuple(root, label, dereplicated_genomes dir)
}