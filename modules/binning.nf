/*
 * Mirrors the "#binning (metabat2, semibin, remag)" block from the .def
 * file. Each binner now gets its own conda env (built automatically by
 * Nextflow) instead of one shared "binners" env, and each one can be
 * turned on/off independently via params. SemiBin2's three models can
 * also be selected individually, instead of always running all three.
 */

// ---- jgi_summarize_bam_contig_depths (ships with the metabat2 package) ----
process DEPTH {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"
    publishDir { "${params.outdir}" }, mode: params.publish_mode

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
    publishDir { "${params.outdir}/${root}_metabat" }, mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly), path(depth_path)

    output:
    tuple val(root), path("bins/*.fa"), emit: bins
    tuple val(root), path("bins"), emit: bins_dir
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
    publishDir { "${params.outdir}/${root}_semibin/${model}_model" }, mode: params.publish_mode

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
    publishDir { "${params.outdir}/${root}_remag" }, mode: params.publish_mode

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
 * whichever binners run_metabat2 / run_semibin2 / run_remag ask for
 * (passed in as plain arguments, not read from params - see the
 * boolean-coercion note in main.nf for why). When SemiBin2 is on,
 * params.semibin_models controls *which* of the three models run.
 *
 * Also builds a unified (root, label, dir) stream of every raw bin set
 * that actually ran, for DAS_Tool to consume independently of dRep.
 */
workflow BINNING {
    take:
    ch_assembly_bam   // tuple(root, mylo_assembly, bam, bai)
    run_metabat2
    run_semibin2
    run_remag

    main:
    ch_metabat_bins     = Channel.empty()
    ch_metabat_bins_dir = Channel.empty()
    ch_semibin_bins     = Channel.empty()
    ch_semibin_bins_dir = Channel.empty()
    ch_remag_bins       = Channel.empty()
    ch_remag_bins_dir   = Channel.empty()

    if (run_metabat2) {
        DEPTH(ch_assembly_bam.map { root, assembly, bam, bai -> tuple(root, bam, bai) })

        ch_metabat_in = ch_assembly_bam
            .map { root, assembly, bam, bai -> tuple(root, assembly) }
            .join(DEPTH.out.depth)

        METABAT2(ch_metabat_in)
        ch_metabat_bins     = METABAT2.out.bins
        ch_metabat_bins_dir = METABAT2.out.bins_dir
    }

    if (run_semibin2) {
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
        ch_semibin_bins     = SEMIBIN2.out.bins
        ch_semibin_bins_dir = SEMIBIN2.out.bins_dir
    }

    if (run_remag) {
        REMAG(ch_assembly_bam)
        ch_remag_bins     = REMAG.out.bins
        ch_remag_bins_dir = REMAG.out.bins_dir
    }

    // ---- unified (root, label, dir) stream for DAS_Tool, built only from
    // the already-safe ch_* variables above - never touches PROCESS.out
    // directly, so a disabled binner just contributes an empty channel
    // via .mix() instead of crashing on an uninvoked process ----
    ch_metabat_labeled = ch_metabat_bins_dir.map { root, dir -> tuple(root, 'metabat', dir) }
    ch_semibin_labeled = ch_semibin_bins_dir.map { root, model, dir -> tuple(root, model, dir) }
    ch_remag_labeled   = ch_remag_bins_dir.map   { root, dir -> tuple(root, 'remag', dir) }
    ch_all_bin_sets    = ch_metabat_labeled.mix(ch_semibin_labeled).mix(ch_remag_labeled)

    emit:
    metabat_bins   = ch_metabat_bins
    semibin_bins   = ch_semibin_bins
    remag_bins     = ch_remag_bins
    remag_bins_dir = ch_remag_bins_dir
    all_bin_sets   = ch_all_bin_sets
}
