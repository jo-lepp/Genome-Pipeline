process GVCLASS_SETUP {
    label 'small'
    storeDir "${params.db_dir}/gvclass_install"
    conda "conda-forge::pixi conda-forge::git"

    output:
    path "gvclass", emit: install_dir

    script:
    """
    git clone ${params.gvclass_repo_url} gvclass
    cd gvclass && pixi install --frozen && pixi run setup-db
    """
}

/*
 * FNA_PREP now takes dRep's single output directly - no more separate
 * metabat/remag/semibin-unzipped inputs to merge, since there's only one
 * bin source per sample now. Oversized-bin subsampling logic is unchanged.
 */
process FNA_PREP {
    tag "$root"
    label 'gvclass'
    conda "bioconda::genometools-genometools"

    input:
    tuple val(root), path(derep_dir)

    output:
    tuple val(root), path("fna_prepped"), emit: fna_dir
    path "subsampled_bins.tsv", emit: subsampled, optional: true

    script:
    """
    mkdir -p fna_prepped
    : > subsampled_bins.tsv
    for bin in "${derep_dir}"/*.fa; do
        size_mb=\$(du -m "\$bin" | cut -f1)
        if [ "\$size_mb" -gt ${params.gvclass_max_bin_size_mb} ]; then
            gt splitfasta -targetsize ${params.gvclass_max_bin_size_mb}000000 "\$bin"
            mv "\$bin".* fna_prepped/ 2>/dev/null || true
            n_parts=\$(ls "\$bin".* 2>/dev/null | wc -l)
            echo -e "${root}\\t\$(basename "\$bin")\\t\${size_mb}MB\\t\${n_parts} parts" >> subsampled_bins.tsv
        else
            cp "\$bin" fna_prepped/
        fi
    done
    """
}

process GVCLASS {
    tag "$root"
    label 'gvclass'
    conda "conda-forge::pixi"
    publishDir { "${params.outdir}/${root}_gvclass" }, mode: params.publish_mode

    input:
    tuple val(root), path(fna_dir)
    path install_dir

    output:
    tuple val(root), path("gvclass_summary.tsv"), emit: summary
    path "gvclass_out", emit: gvclass_dir

    script:
    """
    cd "${install_dir}"
    pixi run ./gvclass -o "\$(pwd)/../gvclass_out" -t ${params.threads_gvclass} "\$(pwd)/../${fna_dir}"
    cd ..
    cp gvclass_out/gvclass_summary.tsv .
    """
}

process COLLECT_SUBSAMPLE_REPORT {
    label 'small'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path reports

    output:
    path "gvclass_subsampled_bins_summary.tsv"

    script:
    """
    echo -e "sample\\tbin\\toriginal_size\\tparts" > gvclass_subsampled_bins_summary.tsv
    cat ${reports.join(' ')} >> gvclass_subsampled_bins_summary.tsv 2>/dev/null || true
    """
}

workflow GVCLASS_TAXONOMY {
    take:
    derep_dir   // tuple(root, dir)

    main:
    GVCLASS_SETUP()
    FNA_PREP(derep_dir)
    GVCLASS(FNA_PREP.out.fna_dir, GVCLASS_SETUP.out.install_dir.first())
    COLLECT_SUBSAMPLE_REPORT(FNA_PREP.out.subsampled.collect())

    emit:
    gvclass_dir = GVCLASS.out.gvclass_dir
    summary     = GVCLASS.out.summary   // feeds REPORT directly
}
