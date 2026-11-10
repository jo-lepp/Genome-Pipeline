process REPORT {
    tag "$root"
    label 'small'
    conda "conda-forge::pandas conda-forge::openpyxl"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(root), path(checkm2_report), path(derep_dir), path(busco_summary), path(gtdbtk_summary), path(gvclass_summary), path(ssu_summary)

    output:
    path "${root}_bin_report.xlsx", emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob, os

    def is_placeholder(p):
        return os.path.basename(str(p)) == "NO_FILE"

    df = pd.read_csv("${checkm2_report}", sep="\\t")
    df = df.rename(columns={"Name": "Bin"})

    def safe_read(path, label):
        if is_placeholder(path):
            print(f"[REPORT] {label} not run for this sample - leaving those columns as N/A")
            return None
        matches = glob.glob(str(path))
        if not matches:
            print(f"[REPORT] no {label} file found at {path} - leaving those columns as N/A")
            return None
        return pd.read_csv(matches[0] if os.path.isdir(path) else path, sep="\\t")

    gtdb = safe_read("${gtdbtk_summary}", "GTDB-Tk summary")
    if gtdb is not None:
        gtdb_cols = gtdb[["user_genome", "classification", "closest_genome_reference",
                           "closest_genome_ani", "closest_genome_af", "red_value"]]
        gtdb_cols = gtdb_cols.rename(columns={"user_genome": "Bin", "classification": "GTDB_classification"})
        df = df.merge(gtdb_cols, on="Bin", how="left")
    else:
        for c in ["GTDB_classification", "closest_genome_reference", "closest_genome_ani", "closest_genome_af", "red_value"]:
            df[c] = "N/A"

    busco = safe_read("${busco_summary}", "BUSCO summary")
    if busco is not None:
        busco.columns = [c.strip() for c in busco.columns]
    
        # Column names are matched flexibly rather than hardcoded exactly, since
        # this file's literal header wasn't independently confirmed the way the
        # other three tools' were - verify against your real file (see note below)
        # and tell me if any of these end up matching the wrong column.
        def find_col(candidates):
            for cand in candidates:
                for c in busco.columns:
                    if c.lower().replace(" ", "_") == cand:
                        return c
            return None
    
        input_col      = find_col(["input_file", "input"])
        complete_col   = find_col(["complete"])
        single_col     = find_col(["single", "complete_single", "complete_sc"])
        dup_col        = find_col(["duplicated", "complete_duplicated"])
        frag_col       = find_col(["fragmented"])
        missing_col    = find_col(["missing"])
        nmarkers_col   = find_col(["n_markers", "total_buscos"])
    
        busco["Bin"] = busco[input_col].apply(lambda f: os.path.splitext(os.path.basename(str(f)))[0])
    
        rename_map = {"Bin": "Bin"}
        keep_cols = ["Bin"]
        for col, label in [
            (complete_col,  "BUSCO_Complete"),
            (single_col,    "BUSCO_Complete_Single"),
            (dup_col,       "BUSCO_Complete_Duplicated"),
            (frag_col,      "BUSCO_Fragmented"),
            (missing_col,   "BUSCO_Missing"),
            (nmarkers_col,  "BUSCO_n_markers"),
        ]:
            if col is not None:
                rename_map[col] = label
                keep_cols.append(col)
    
        busco_cols = busco[keep_cols].rename(columns=rename_map)
        df = df.merge(busco_cols, on="Bin", how="left")
    else:
        for c in ["BUSCO_Complete", "BUSCO_Complete_Single", "BUSCO_Complete_Duplicated",
                  "BUSCO_Fragmented", "BUSCO_Missing", "BUSCO_n_markers"]:
            df[c] = "N/A"

    gv = safe_read("${gvclass_summary}", "GVClass summary")
    if gv is not None:
        gv_cols = gv[["query", "taxonomy_majority"]].rename(
            columns={"query": "Bin", "taxonomy_majority": "GVclass_classification"})
        df = df.merge(gv_cols, on="Bin", how="left")
    else:
        df["GVclass_classification"] = "N/A"
    
    df["ssuextract_tax"] = "N/A"
    df["ssuextract_tax_reference"] = "N/A"
    df["SSU_seq"] = "N/A"
    ssu = safe_read("${ssu_summary}", "SSU screening summary")
    if ssu is not None and len(ssu) > 0:
        bin_fastas = glob.glob(os.path.join("${derep_dir}", "*.fa")) + glob.glob(os.path.join("${derep_dir}", "*.fasta"))
        for _, row in ssu.iterrows():
            contig_id = str(row.get("contig_name", row.iloc[0]))
            for bf in bin_fastas:
                with open(bf) as fh:
                    headers = [l[1:].split()[0] for l in fh if l.startswith(">")]
                if contig_id in headers:
                    bin_name = os.path.splitext(os.path.basename(bf))[0]
                    df.loc[df["Bin"] == bin_name, "ssuextract_tax"] = row.get("taxonomy", "N/A")
                    df.loc[df["Bin"] == bin_name, "ssuextract_tax_reference"] = row.get("reference_taxonomy", "N/A")
                    df.loc[df["Bin"] == bin_name, "SSU_seq"] = row.get("query_sequence", "N/A")
                    break

    df.insert(0, "Sample", "${root}")
    df.to_excel("${root}_bin_report.xlsx", index=False)
    """
}
