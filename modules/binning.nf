/*
 * Exactly ONE binner runs per sample now (params.binner), not all three in
 * parallel. REMAG is a special case: it only recovers eukaryotic bins by
 * design, so its rejected (bacterial-looking) contigs - written to
 * *_non_eukaryotic.fasta via --save-filtered-contigs - get run through a
 * SECOND, user-chosen binner (params.remag_reject_binner). Both sets of
 * bins are then merged into one list for dRep; since they represent
 * different organisms (eukaryotic vs bacterial), they won't share ANI
 * overlap, so dRep just keeps both untouched rather than wastefully
 * trying to dereplicate them against each other.
 */

process DEPTH {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"

    input:
    tuple val(root), path(bam), path(bai)

    output:
    tuple val(root), path("${root}_depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth "${root}_depth.txt" "${bam}"
    """
}

process METABAT2 {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"
    publishDir { "${params.outdir}/${root}_metabat" }, mode: params.publish_mode

    input:
    tuple val(root), path(assembly), path(depth_path)

    output:
    tuple val(root), path("bins/*.fa"), emit: bins

    script:
    """
    mkdir -p bins
    metabat2 -i "${assembly}" -a "${depth_path}" -o bins/metabat
    """
}

process SEMIBIN2 {
    tag "${root}:${model}"
    label 'binning'
    conda "bioconda::semibin"
    publishDir { "${params.outdir}/${root}_semibin/${model}_model" }, mode: params.publish_mode

    input:
    tuple val(root), path(assembly), path(bam), path(bai), val(model), val(mode_flag)

    output:
    tuple val(root), path("output_bins/*.fa.gz"), emit: bins

    script:
    """
    SemiBin2 single_easy_bin --sequencing-type long_read ${mode_flag} -i "${assembly}" -b "${bam}" -o .
    """
}

process REMAG {
    tag "$root"
    label 'binning'
    conda "bioconda::remag"
    publishDir { "${params.outdir}/${root}_remag" }, mode: params.publish_mode

    input:
    tuple val(root), path(assembly), path(bam), path(bai)

    output:
    tuple val(root), path("bins/*.fa"),                emit: bins
    tuple val(root), path("*_non_eukaryotic.fasta"),    emit: rejects, optional: true

    script:
    """
    remag "${assembly}" -c "${bam}" -o . --save-filtered-contigs
    """
}

/*
 * Depth + binning for the reject-rebinning path. METABAT2/SEMIBIN2 are
 * reused via aliasing (include ... as ...) so the reject path shares the
 * exact same, already-tested binning logic - the only difference is the
 * input assembly is REMAG's non_eukaryotic.fasta instead of the full
 * myloasm assembly, and depth is recalculated against just those contigs
 * using the SAME original BAM (jgi_summarize_bam_contig_depths only
 * reports depth for contigs actually present in the given fasta, so
 * there's no need to re-align reads against the reduced contig set).
 */
include { DEPTH as DEPTH_REJECT; METABAT2 as METABAT2_REJECT; SEMIBIN2 as SEMIBIN2_REJECT } from './binning.nf'

workflow BINNING {
    take:
    ch_assembly_bam       // tuple(root, assembly, bam, bai)
    binner                // 'metabat' | 'semibin' | 'remag'
    remag_reject_binner   // 'metabat' | 'semibin' - only used when binner == 'remag'

    main:
    def all_models = [
        soil:   '--environment soil',
        self:   '--self-supervised',
        global: '--environment global'
    ]
    def selected_models = params.semibin_models.split(',')*.trim()
    def unknown = selected_models - all_models.keySet().toList()
    if (unknown) {
        error("Unknown semibin_models entry: ${unknown.join(', ')}. Valid options: ${all_models.keySet().join(', ')}")
    }
    ch_model_list = Channel.fromList(selected_models.collect { m -> [m, all_models[m]] })

    ch_bins = Channel.empty()

    if (binner == 'metabat') {
        DEPTH(ch_assembly_bam.map { root, assembly, bam, bai -> tuple(root, bam, bai) })
        ch_in = ch_assembly_bam.map { root, assembly, bam, bai -> tuple(root, assembly) }.join(DEPTH.out.depth)
        METABAT2(ch_in)
        ch_bins = METABAT2.out.bins

    } else if (binner == 'semibin') {
        ch_in = ch_assembly_bam.combine(ch_model_list)
        SEMIBIN2(ch_in)
        // merge every selected model's bins into one flat list per sample
        ch_bins = SEMIBIN2.out.bins.map { root, model, bins -> tuple(root, bins) }
            .groupTuple(by: 0)
            .map { root, binlists -> tuple(root, binlists.flatten()) }

    } else if (binner == 'remag') {
        REMAG(ch_assembly_bam)

        // pair each sample's rejected-contigs fasta back with its own BAM
        ch_reject_assembly_bam = REMAG.out.rejects
            .join(ch_assembly_bam.map { root, assembly, bam, bai -> tuple(root, bam, bai) })
            .map { root, rejects, bam, bai -> tuple(root, rejects, bam, bai) }

        if (remag_reject_binner == 'metabat') {
            DEPTH_REJECT(ch_reject_assembly_bam.map { root, rejects, bam, bai -> tuple(root, bam, bai) })
            ch_reject_in = ch_reject_assembly_bam.map { root, rejects, bam, bai -> tuple(root, rejects) }.join(DEPTH_REJECT.out.depth)
            METABAT2_REJECT(ch_reject_in)
            ch_reject_bins = METABAT2_REJECT.out.bins
        } else {
            ch_reject_in = ch_reject_assembly_bam.combine(ch_model_list)
            SEMIBIN2_REJECT(ch_reject_in)
            ch_reject_bins = SEMIBIN2_REJECT.out.bins.map { root, model, bins -> tuple(root, bins) }
                .groupTuple(by: 0)
                .map { root, binlists -> tuple(root, binlists.flatten()) }
        }

        // combine REMAG's own (eukaryotic) bins with the secondary binner's
        // (bacterial) bins into one flat list per sample
        ch_bins = REMAG.out.bins
            .join(ch_reject_bins)
            .map { root, euk_bins, bac_bins -> tuple(root, ([euk_bins, bac_bins].flatten())) }
    }

    emit:
    all_bins = ch_bins   // tuple(root, [bin fasta paths])
}
