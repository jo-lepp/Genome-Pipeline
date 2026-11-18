
process GVCLASS_SETUP {
    label 'small'
    conda "conda-forge::pixi conda-forge::git"

    output:
    val(install_dir), emit: install_dir

    script:
    install_dir = file("${params.db_dir}/gvclass_install/gvclass").toAbsolutePath().toString()

    """
    INSTALL_DIR="${install_dir}"

    # .setup_complete is created after pixi installation and database setup completes
    # if it already exists, ssuextract will not be reinstalled. If it doesn't
    # exist, creating the folder and cloning the repository will proceed.
    if [ ! -f "\$INSTALL_DIR/.setup_complete" ]; then
        mkdir -p "\$(dirname "\$INSTALL_DIR")"

        if [ ! -d "\$INSTALL_DIR/.git" ]; then
            git clone ${params.gvclass_repo_url} "\$INSTALL_DIR"
        fi

        cd "\$INSTALL_DIR"

        pixi install --frozen
        pixi run setup-db

        touch .setup_complete
    fi
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
    conda "${projectDir}/envs/genometools.yml"

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
    conda "${projectDir}/envs/pixi_tools.yml"
    publishDir { "${params.outdir}/${root}_gvclass" }, mode: params.publish_mode

    input:
    tuple val(root), path(fna_dir)
    val install_dir

    output:
    tuple val(root), path("gvclass_summary.tsv"), emit: summary
    path "gvclass_out", emit: gvclass_dir

     script:
    """
    export GVCLASS_DB="${params.db_dir}/gvclass"
    task_dir="\$(pwd)"
    cd "${install_dir}"
    pixi run ./gvclass -o "\${task_dir}/gvclass_out" -t ${params.threads_gvclass} "\${task_dir}/${fna_dir}"
    cd "\${task_dir}"
    cp gvclass_out/gvclass_summary.tsv .
    """
}

process COLLECT_SUBSAMPLE_REPORT {
    label 'small'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path reports, stageAs: 'subsample_report??.tsv'

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
    ch_install_dir = Channel.value(file("${params.software_dir}/gvclass", checkIfExists: true).toAbsolutePath().toString())
    FNA_PREP(derep_dir)
    GVCLASS(FNA_PREP.out.fna_dir, ch_install_dir)
    COLLECT_SUBSAMPLE_REPORT(FNA_PREP.out.subsampled.collect())

    emit:
    gvclass_dir = GVCLASS.out.gvclass_dir
    summary     = GVCLASS.out.summary   // feeds REPORT directly
}
