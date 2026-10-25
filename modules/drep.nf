/*
 * NOTE: dRep is never run on REMAG bins specifically as their own category
 * anymore - with the single-binner architecture, BINNING.out.all_bins
 * already contains whichever bins actually resulted (one binner's output,
 * or REMAG's euk bins + the secondary binner's bacterial bins merged).
 * DREP runs once per sample on that combined set.
 *
 * IMPORTANT: dRep's default filtering step internally shells out to
 * checkM v1 (a separate, older tool from CheckM2, despite the similar
 * name) to score completeness/contamination before clustering. checkM1
 * needs its own reference database registered via `checkm data setRoot`
 * before it will run at all - without it, dRep fails at Step 1 with
 * "!!! checkM failed !!!" and silently produces zero dereplicated
 * genomes, even though dRep itself exits with status 0. This is skipped
 * entirely when drep_ignore_genome_quality is true.
 */

process CHECKM1_DB {
    label 'binning'
    conda "bioconda::drep conda-forge::python=3.10 conda-forge::setuptools=75.8.2"
    storeDir "${params.db_dir}/checkm1_db"

    output:
    path "checkm_data", emit: db_dir

    script:
    """
    wget -O checkm_data.tar.gz "${params.checkm1_data_url}"
    mkdir -p checkm_data
    tar -xzf checkm_data.tar.gz -C checkm_data
    rm checkm_data.tar.gz
    """
}

process DREP {
    tag "$root"
    label 'binning'
    conda "bioconda::drep conda-forge::python=3.10 conda-forge::setuptools=75.8.2"
    publishDir { "${params.outdir}/${root}_drep" }, mode: params.publish_mode

    input:
    tuple val(root), path(bins, stageAs: 'bin_??.fa')
    path checkm1_db

    output:
    tuple val(root), path("dereplicated_genomes"), emit: derep_dir

    script:
    def quality_flag = params.drep_ignore_genome_quality ? '--ignoreGenomeQuality' : ''
    def checkm_setup = params.drep_ignore_genome_quality ? '' : "checkm data setRoot \"\$(pwd)/${checkm1_db}\""
    """
    ${checkm_setup}
    # bins may arrive gzipped (SemiBin2) or plain (metabat2/REMAG) - normalize
    # to plain fasta before dRep ever sees them, so every downstream consumer
    # (BUSCO, CheckM2, GTDB-Tk, GVClass) can assume uniform, uncompressed input
    mkdir -p plain_bins
    for f in bin_*.fa; do
        if file "\$f" | grep -q gzip; then
            gunzip -c "\$f" > "plain_bins/\$(basename "\$f")"
        else
            cp "\$f" "plain_bins/"
        fi
    done

    dRep dereplicate drep_out -g plain_bins/*.fa --S_algorithm ANImf ${quality_flag}
    cp -r drep_out/dereplicated_genomes .
    """
}

workflow DEREPLICATION {
    take:
    all_bins                    // tuple(root, [bin fasta paths]) - from BINNING.out.all_bins
    drep_ignore_genome_quality  // boolean

    main:
    ch_checkm1_db = Channel.value(file("NO_FILE"))
    if (!drep_ignore_genome_quality) {
        CHECKM1_DB()
        ch_checkm1_db = CHECKM1_DB.out.db_dir
    }

    DREP(all_bins, ch_checkm1_db.first())

    emit:
    derep_dir = DREP.out.derep_dir
}
