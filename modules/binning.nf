
// ---- jgi_summarize_bam_contig_depths (ships with the metabat2 package) ----
process DEPTH {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(root), path(bam), path(bai)

    output:
    tuple val(root), path("${root}_depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth "${root}_depth.txt" "${bam}"
    """
}

// ---- metabat2 ----
process METABAT2 {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"
    publishDir "${params.outdir}/${root}_metabat", mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly), path(depth_path)

    output:
    tuple val(root), path("bins/*.fa"), emit: bins
    path "bins", emit: bins_dir
    path "*.txt", emit: metabat_stats

    script:
    """
    mkdir -p bins
    metabat2 -i "${mylo_assembly}" -a "${depth_path}" -o bins/metabat

    mv bins/*.txt . 2>/dev/null || true
    echo "metabat2 binning done"
    """
}

// ---- SemiBin2, run once per selected model ----
process SEMIBIN2 {
    tag "${root}:${model}"
    label 'binning'
    conda "bioconda::semibin"
    publishDir "${params.outdir}/${root}_semibin/${model}_model", mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly), path(bam), path(bai), val(model), val(mode_flag)

    output:
    tuple val(root), val(model), path("output_bins/*.fa.gz"), emit: bins
    tuple val(root), val(model), path("output_bins"), emit: bins_dir

    script:
    """
    SemiBin2 single_easy_bin --sequencing-type long_read ${mode_flag} -i "${mylo_assembly}" -b "${bam}" -o .
    echo "Semibin2 binning done (${model})"
    """
}

// ---- REMAG ----
process REMAG {
    tag "$root"
    label 'binning'
    conda "bioconda::remag"
    publishDir "${params.outdir}/${root}_remag", mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly), path(bam), path(bai)

    output:
    tuple val(root), path("bins/*.fa"), emit: bins
    tuple val(root), path("bins"), emit: bins_dir

    script:
    """
    remag "${mylo_assembly}" -c "${bam}" -o . --save-filtered-contigs
    echo "remag binning done"
    """
}

/*
 * Subworkflow: takes the assembly + sorted BAM for one sample and runs
 * whichever binners params.run_metabat2 / run_semibin2 / run_remag ask
 * for. When SemiBin2 is on, params.semibin_models controls *which* of
 * the three models actually run.
 */
workflow BINNING {
    take:
    ch_assembly_bam   // tuple(root, mylo_assembly, bam, bai)

    main:
    ch_metabat_bins = Channel.empty()
    ch_semibin_bins = Channel.empty()
    ch_remag_bins   = Channel.empty()
    ch_remag_bins_dir = Channel.empty()

    if (params.run_metabat2) {
        DEPTH(ch_assembly_bam.map { root, assembly, bam, bai -> tuple(root, bam, bai) })

        ch_metabat_in = ch_assembly_bam
            .map { root, assembly, bam, bai -> tuple(root, assembly) }
            .join(DEPTH.out.depth)

        METABAT2(ch_metabat_in)
        ch_metabat_bins = METABAT2.out.bins
    }

    if (params.run_semibin2) {
        // the full set of models this pipeline knows about and their flags
        def all_models = [
            soil:   '--environment soil',
            self:   '--self-supervised',
            global: '--environment global'
        ]

        // params.semibin_models is a comma-separated string, e.g. "soil,global"
        def selected = params.semibin_models.split(',')*.trim()
        def unknown  = selected - all_models.keySet().toList()
        if (unknown) {
            error("Unknown semibin_models entry: ${unknown.join(', ')}. Valid options: ${all_models.keySet().join(', ')}")
        }

        ch_semibin_models = Channel.fromList(selected.collect { model -> [model, all_models[model]] })

        ch_semibin_in = ch_assembly_bam.combine(ch_semibin_models)
        SEMIBIN2(ch_semibin_in)
        ch_semibin_bins = SEMIBIN2.out.bins
    }

        if (params.run_remag) {
        REMAG(ch_assembly_bam)
        ch_remag_bins     = REMAG.out.bins
        ch_remag_bins_dir = REMAG.out.bins_dir     // <-- add this line
    }

    emit:
    metabat_bins     = ch_metabat_bins
    semibin_bins     = ch_semibin_bins
    remag_bins       = ch_remag_bins
    remag_bins_dir   = ch_remag_bins_dir   // <-- add this line
}
}