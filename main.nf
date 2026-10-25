#!/usr/bin/env nextflow
/*
 *      DEFUALT STEPS:
 *   1. MYLOASM            - long-read assembly
 *   2. MAPPING            - minimap2 + samtools sort/index
 *   3. BINNING            - metabat2 / semibin2 (configurable models) / remag ; DAS_Tool
 *   4. DEREPLICATION      - dRep (metabat + semibin only, NOT remag)
 *   5. CHECKM2_QC         - bacterial completeness/contamination (NOT remag)
 *   6. BUSCO_QC           - eukaryotic QC (ALL bin sets, including remag)
 *   7. GTDBTK_CLASSIFY    - bacterial taxonomy (NOT remag)
 *   8. GVCLASS_TAXONOMY   - eukaryotic taxonomy (ALL bin sets, including remag)
 *
 *      Optional Changes to Steps:
 *   1. can toggle the semibin models, remag, and metabat2 binning
 *   2. can toggle gvclass and gtdbtk, option to replace with ssuExtract/cmsearch
 */

nextflow.enable.dsl = 2

include { MYLOASM }          from './modules/myloasm.nf'
include { MAPPING }          from './modules/mapping.nf'
include { BINNING }          from './modules/binning.nf'
include { DEREPLICATION }    from './modules/drep.nf'
include { CHECKM2_QC }       from './modules/checkm2.nf'
include { BUSCO_QC }         from './modules/busco.nf'
include { GTDBTK_CLASSIFY }  from './modules/gtdbtk.nf'
include { GVCLASS_TAXONOMY } from './modules/gvclass.nf'
include { CMSEARCH_EUK }     from './modules/cmsearch.nf'
include { SSU_EXTRACT_ALL }  from './modules/cmsearch.nf'
include { REPORT } from './modules/report.nf'

workflow {
    def drep_ignore_genome_quality = (params.drep_ignore_genome_quality instanceof Boolean)
        ? params.drep_ignore_genome_quality
        : params.drep_ignore_genome_quality.toString().trim().toLowerCase() == 'true'
    def run_gtdbtk   = (params.run_gtdbtk instanceof Boolean) ? params.run_gtdbtk : params.run_gtdbtk.toString().trim().toLowerCase() == 'true'
    def fast_mode       = (params.fast_mode instanceof Boolean) ? params.fast_mode : params.fast_mode.toString().trim().toLowerCase() == 'true'
    def run_gvclass     = (params.run_gvclass instanceof Boolean) ? params.run_gvclass : params.run_gvclass.toString().trim().toLowerCase() == 'true'
    def run_ssu_screen  = (params.run_ssu_screen instanceof Boolean) ? params.run_ssu_screen : params.run_ssu_screen.toString().trim().toLowerCase() == 'true'
    def fast_mode_tool  = params.fast_mode_tool

    if (!(params.binner in ['metabat', 'semibin', 'remag'])) {
        error("params.binner must be 'metabat', 'semibin', or 'remag'. Got: '${params.binner}'")
    }
    if (params.binner == 'remag' && !(params.remag_reject_binner in ['metabat', 'semibin'])) {
        error("params.remag_reject_binner must be 'metabat' or 'semibin' when params.binner is 'remag'. Got: '${params.remag_reject_binner}'")
    }

    if (!(params.read_type in ['hifi', 'ont'])) {
        error("params.read_type must be 'hifi' or 'ont', got: '${params.read_type}'")
    }
    if (!(params.read_type in ['hifi', 'ont'])) {
        exit 1, "params.read_type must be 'hifi' or 'ont', got: '${params.read_type}'"
    }
    if (!fast_mode && !(fast_mode_tool in ['ssu_extract', 'cmsearch'])) {
        error("params.fast_mode_tool must be 'ssu_extract' or 'cmsearch'. Got: '${fast_mode_tool}'")
    }

    if (fast_mode) {
        fast_mode_tool = 'cmsearch'
        run_gvclass    = false
        run_ssu_screen = true
    }

    // ================= 0. BUILD PER-SAMPLE CHANNELS =================
    if (!(params.read_type in ['hifi', 'ont'])) {
        error("params.read_type must be 'hifi' or 'ont', got: '${params.read_type}'")
    }

    if (params.samplesheet?.trim()) {
        // ---- multi-sample mode ----
        def sheet_file = file(params.samplesheet)
        if (!sheet_file.exists()) error("samplesheet not found: ${params.samplesheet}")

        // Parsed eagerly as a plain list (not a lazy channel), so every
        // validation error below surfaces immediately and halts the run
        // before any process launches. A channel .subscribe{} callback
        // fires asynchronously and can't reliably guarantee that.
        def rows = sheet_file.splitCsv(header: true)

        def parsed = rows.collect { row ->
            if (!row.fastq?.trim()) error("samplesheet row missing 'fastq': ${row}")

            def fastq_file = file(row.fastq)
            if (!fastq_file.exists()) error("fastq not found: ${row.fastq}")

            def sample_name = fastq_file.simpleName

            def assembly_file = row.assembly?.trim() ? file(row.assembly) : null
            if (assembly_file && !assembly_file.exists()) {
                error("assembly not found for ${sample_name}: ${row.assembly}")
            }

            [sample_name, fastq_file, assembly_file]
        }

        def dupes = parsed.collect { it[0] }.countBy { it }.findAll { k, v -> v > 1 }.keySet()
        if (dupes) {
            error("Duplicate sample name(s) derived from fastq filename: ${dupes.join(', ')}. Rename the input files to disambiguate.")
        }

        ch_samples = Channel.fromList(parsed)

        ch_input = ch_samples.map { sample, fastq, assembly -> tuple(sample, fastq) }

        ch_needs_assembly = ch_samples
            .filter { sample, fastq, assembly -> assembly == null }
            .map    { sample, fastq, assembly -> tuple(sample, fastq) }

        ch_has_assembly = ch_samples
            .filter { sample, fastq, assembly -> assembly != null }
            .map    { sample, fastq, assembly -> tuple(sample, assembly) }

    } else {
        // ---- single-sample mode ----
        def input_file = file(params.input_fgz)
        def root_name  = params.nametag?.trim() ? params.nametag.trim() : input_file.simpleName

        ch_input = Channel.of( tuple(root_name, input_file) )

        if (params.assembly?.trim()) {
            def assembly_file = file(params.assembly)
            if (!assembly_file.exists()) error("params.assembly was set but file not found: ${params.assembly}")
            ch_needs_assembly = Channel.empty()
            ch_has_assembly   = Channel.of( tuple(root_name, assembly_file) )
        } else {
            ch_needs_assembly = ch_input
            ch_has_assembly   = Channel.empty()
        }
    }

    // ================= 1. ASSEMBLY (only for samples that need it) =================
    MYLOASM(ch_needs_assembly)
    ch_assembly = MYLOASM.out.assembly.mix(ch_has_assembly)

    // ================= 2. MAPPING =================
    ch_map_in = ch_assembly.join(ch_input)
    MAPPING(ch_map_in)
    ch_assembly_bam = ch_assembly.join(MAPPING.out.bam)

    // ================= 3. BINNING (metabat2 / semibin2 / remag) =================
    BINNING(ch_assembly_bam, params.binner, params.remag_reject_binner)

    // ================= 4. DEREPLICATION (metabat + semibin, remag) =================
    DEREPLICATION(BINNING.out.all_bins, drep_ignore_genome_quality)

    // ================= 5. CHECKM2 (metabat + semibin, NOT remag) =================
    CHECKM2_QC(DEREPLICATION.out.derep_dir)

    // ================= 6. BUSCO (all bin sets, including remag) =================
    BUSCO_QC(DEREPLICATION.out.derep_dir)

    // ================= 7. GTDB-TK (metabat + semibin, NOT remag) =================
    ch_gtdbtk_summary = Channel.empty()
    if (run_gtdbtk) {
        GTDBTK_CLASSIFY(DEREPLICATION.out.derep_dir)
        ch_gtdbtk_summary = GTDBTK_CLASSIFY.out.summary
    }

    // ================= 8. GVCLASS (all bin sets, including remag) =================
    ch_gvclass_summary = Channel.empty()
    if (run_gvclass) {
        GVCLASS_TAXONOMY(DEREPLICATION.out.derep_dir)
        ch_gvclass_summary = GVCLASS_TAXONOMY.out.summary
    }

    // ================= fast mode: cmsearch instead of full taxonomy =================
    ch_ssu_summary = Channel.empty()
    if (run_ssu_screen) {
        if (params.fast_mode_tool == 'ssu_extract') {
            SSU_EXTRACT_ALL(ch_assembly)
            ch_ssu_summary = SSU_EXTRACT_ALL.out.summary
        } else {
            CMSEARCH_EUK(ch_assembly)
            ch_ssu_summary = CMSEARCH_EUK.out.summary
        }
    }
    if (run_ssu_screen) {
        if (fast_mode_tool == 'ssu_extract') {
            SSU_EXTRACT_ALL(ch_assembly)
        } else {
            CMSEARCH_EUK(ch_assembly)
        }
    }

    // ---- compiled QC + taxonomy spreadsheet ----

    def no_file = file("NO_FILE")

    ch_gtdbtk_for_report = run_gtdbtk
        ? ch_gtdbtk_summary
        : CHECKM2_QC.out.quality_report.map { root, quality_report -> tuple(root, no_file) }
    
    ch_gvclass_for_report = run_gvclass
        ? ch_gvclass_summary
        : CHECKM2_QC.out.quality_report.map { root, quality_report -> tuple(root, no_file) }
    
    ch_ssu_for_report = run_ssu_screen
        ? ch_ssu_summary
        : CHECKM2_QC.out.quality_report.map { root, quality_report -> tuple(root, no_file) }
    
    ch_report_in = CHECKM2_QC.out.quality_report
        .join(DEREPLICATION.out.derep_dir)
        .join(ch_gtdbtk_for_report)
        .join(ch_gvclass_for_report)
        .join(ch_ssu_for_report)
    
    REPORT(ch_report_in)
    REPORT(
        CHECKM2_QC.out.quality_report,
        DEREPLICATION.out.derep_dir,
        ch_gtdbtk_summary,
        ch_gvclass_summary,
        ch_ssu_summary
    )

    workflow.onComplete = {
        try {
            log.info """
            Pipeline finished at : ${workflow.complete}
            Total run time        : ${workflow.duration}
            Success               : ${workflow.success}
            Exit status           : ${workflow.exitStatus}
            """.stripIndent()
        } catch (Exception e) {
            log.warn "Could not print completion summary: ${e.message}"
    }
}
}

