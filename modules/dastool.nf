workflow BINNING {
    take:
    ch_assembly_bam
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
        def all_models = [
            soil:   '--environment soil',
            self:   '--self-supervised',
            global: '--environment global'
        ]

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

    // ---- unified (root, label, dir) stream for DAS_Tool, built from the
    // already-safe ch_* variables above (never touches PROCESS.out
    // directly) - so a disabled binner just contributes an empty channel
    // via .mix(), instead of crashing on an uninvoked process ----
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
