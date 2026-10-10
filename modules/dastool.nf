/*
 * DAS_Tool - integrates bins from multiple binning methods (metabat2,
 * SemiBin2 models, REMAG) and selects a non-redundant, optimized set
 * using its own single-copy-gene scoring - a cross-method consensus
 * step, distinct from dRep (which only dereplicates WITHIN one method's
 * own output). Consumes BINNING's raw bin sets directly, independent
 * of DEREPLICATION/dRep.
 *
 * Reference: Sieber et al. 2018, Nature Microbiology.
 * https://github.com/cmks/DAS_Tool
 */
process DAS_TOOL {
    tag "$root"
    label 'binning'
    conda "bioconda::das_tool"
    publishDir { "${params.outdir}/${root}_dastool" }, mode: params.publish_mode

    input:
    tuple val(root), val(labels), path(bin_dirs, stageAs: 'binset_*'), path(assembly)

    output:
    tuple val(root), path("dastool_DASTool_summary.tsv"),        emit: summary
    tuple val(root), path("dastool_DASTool_contigs2bin.tsv"),    emit: contigs2bin
    path "dastool_allBins.eval",                                 emit: bin_evals, optional: true
    path "DASTool_bins",                                         emit: bins_dir,  optional: true

    script:
    def label_str = labels.join(',')
    """
    mkdir -p contig2bin_tables

    # DAS_Tool's own helper script needs plain (uncompressed) fasta per
    # bin set - metabat2/remag bins are already .fa, but SemiBin2's are
    # .fa.gz, so decompress a working copy of any bin set that needs it
    # before converting to a contigs2bin table.
    i=0
    tables=""
    for d in ${bin_dirs.join(' ')}; do
        label=\$(echo "${label_str}" | cut -d',' -f\$((i+1)))

        if ls "\$d"/*.fa.gz >/dev/null 2>&1; then
            mkdir -p "unzipped_\$i"
            for f in "\$d"/*.fa.gz; do
                gunzip -c "\$f" > "unzipped_\$i/\$(basename "\${f%.gz}")"
            done
            src_dir="unzipped_\$i"
        else
            src_dir="\$d"
        fi

        Fasta_to_Contigs2Bin.sh -i "\$src_dir" -e fa > "contig2bin_tables/\${label}.tsv"
        tables="\${tables}\${tables:+,}contig2bin_tables/\${label}.tsv"
        i=\$((i+1))
    done

    DAS_Tool -i "\$tables" \\
        -l "${label_str}" \\
        -c "${assembly}" \\
        -o dastool \\
        --search_engine diamond \\
        --threads ${params.threads_dastool} \\
        --write_bins \\
        --write_bin_evals
    """
}

/*
 * Subworkflow: groups whichever bin sets BINNING actually produced by
 * sample, pairs each sample with its assembly, and runs DAS_Tool once
 * per sample across all of that sample's bin sets at once.
 */
workflow DASTOOL_INTEGRATION {
    take:
    all_bin_sets   // tuple(root, label, bin_dir) - from BINNING.out.all_bin_sets
    assembly       // tuple(root, assembly)       - e.g. ch_assembly

    main:
    ch_grouped = all_bin_sets
        .groupTuple(by: 0)   // -> tuple(root, [label, label, ...], [dir, dir, ...])
        .join(assembly)      // -> tuple(root, labels, dirs, assembly)

    DAS_TOOL(ch_grouped)

    emit:
    summary     = DAS_TOOL.out.summary
    contigs2bin = DAS_TOOL.out.contigs2bin
    bins_dir    = DAS_TOOL.out.bins_dir
}
