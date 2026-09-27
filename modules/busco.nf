
// ---- gunzip dereplicated SemiBin2 bins so BUSCO can read them ----
// Only needed for self/soil/global - metabat and remag bins are
// already plain .fa and skip straight to BUSCO.
process UNZIP_BINS {
    tag "${root}:${label}"
    label 'small'
    publishDir { "${params.outdir}/${root}_semibin/${label}_model" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(derep_dir)

    output:
    tuple val(root), val(label), path("output_bins_unzipped"), emit: unzipped_dir

    script:
    """
    mkdir -p output_bins_unzipped
    for f in "${derep_dir}"/*.fa.gz; do
        gunzip -c "\$f" > "output_bins_unzipped/\$(basename "\${f%.gz}")"
    done
    """
}

// ---- BUSCO, run once per bin set (metabat/remag/self/soil/global) ----
process BUSCO {
    tag "${root}:${label}"
    label 'qc'
    conda "bioconda::busco"
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(bins_dir), val(pub_subpath)

    output:
    tuple val(root), val(label), path("busco"), emit: busco_dir

    script:
    """
    busco -i "${bins_dir}" --out_path . -o busco -m genome -f -l ${params.busco_lineage}
    echo "busco done (${label})"
    """
}

/*
 * Subworkflow: takes DEREPLICATION's output (metabat + whichever
 * SemiBin2 models ran) and BINNING's REMAG directory, gunzips the
 * SemiBin2 bins, and runs BUSCO across everything that's actually present.
 */
workflow BUSCO_QC {
    take:
    derep_dir       // tuple(root, label, dereplicated_genomes dir) - from DEREPLICATION.out.derep_dir
    remag_bins_dir  // tuple(root, bins dir)                        - from BINNING.out.remag_bins_dir

    main:
    ch_derep_metabat = derep_dir.filter { root, label, dir -> label == 'metabat' }
    ch_derep_semibin = derep_dir.filter { root, label, dir -> label in ['self', 'soil', 'global'] }

    UNZIP_BINS(ch_derep_semibin)

    ch_busco_metabat = ch_derep_metabat.map { root, label, dir -> tuple(root, label, dir, "${root}_metabat") }
    ch_busco_remag   = remag_bins_dir.map    { root, dir       -> tuple(root, 'remag', dir, "${root}_remag") }
    ch_busco_semibin = UNZIP_BINS.out.unzipped_dir.map { root, label, dir -> tuple(root, label, dir, "${root}_semibin/${label}_model") }

    ch_busco_in = ch_busco_metabat.mix(ch_busco_remag).mix(ch_busco_semibin)
    BUSCO(ch_busco_in)

    emit:
    busco_dir = BUSCO.out.busco_dir
    unzipped_dir = UNZIP_BINS.out.unzipped_dir
}