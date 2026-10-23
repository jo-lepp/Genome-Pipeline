process REPORT {
    tag "$root"
    label 'small'
    conda "conda-forge::pandas conda-forge::openpyxl"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(root), path(checkm2_report)
    tuple val(root2), path(derep_dir)
    tuple val(root3), path(gtdbtk_summary), optional: true
    tuple val(root4), path(gvclass_summary), optional: true
    tuple val(root5), path(ssu_summary), optional: true

    output:
    path "${root}_bin_report.xlsx", emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob, os

    df = pd.read_csv("${checkm2_report}", sep="\\t")
    df = df.rename(columns={"Name": "Bin"})

    def safe_read(path_glob, label):
        matches = glob.glob(path_glob)
        if not matches:
            print(f"[REPORT] no {label} file found ({path_glob}) - leaving those columns as N/A")
            return None
        return pd.read_csv(matches[0], sep="\\t")

    gtdb = safe_read("${gtdbtk_summary ?: ''}", "GTDB-Tk summary") if "${gtdbtk_summary ?: ''}" else None
    if gtdb is not None:
        gtdb["closest_genome_reference"] = gtdb["fastani_reference"].where(
            gtdb["fastani_reference"].ne("N/A"), gtdb["closest_placement_reference"])
        gtdb["closest_genome_ani"] = gtdb["fastani_ani"].where(
            gtdb["fastani_ani"].ne("N/A"), gtdb["closest_placement_ani"])
        gtdb["closest_genome_af"] = gtdb["fastani_af"].where(
            gtdb["fastani_af"].ne("N/A"), gtdb["closest_placement_af"])
        gtdb_cols = gtdb[["user_genome", "classification", "closest_genome_reference",
                           "closest_genome_ani", "closest_genome_af", "red_value"]]
        gtdb_cols = gtdb_cols.rename(columns={"user_genome": "Bin", "classification": "GTDB_classification"})
        df = df.merge(gtdb_cols, on="Bin", how="left")
    else:
        for c in ["GTDB_classification", "closest_genome_reference", "closest_genome_ani", "closest_genome_af", "red_value"]:
            df[c] = "N/A"

    gv = safe_read("${gvclass_summary ?: ''}", "GVClass summary") if "${gvclass_summary ?: ''}" else None
    if gv is not None:
        tax_col = next((c for c in gv.columns if "taxonomy" in c.lower() or "classification" in c.lower()), None)
        name_col = gv.columns[0]
        gv_cols = gv[[name_col, tax_col]].rename(columns={name_col: "Bin", tax_col: "GVclass_classification"})
        df = df.merge(gv_cols, on="Bin", how="left")
    else:
        df["GVclass_classification"] = "N/A"

    # SSUextract hits aren't bin-specific - attach a hit to whichever bin
    # actually contains that contig, by checking each bin fasta's headers
    df["ssuextract_tax"] = "N/A"
    df["SSU_seq"] = "N/A"
    if "${ssu_summary ?: ''}":
        ssu = safe_read("${ssu_summary}", "SSUextract summary")
        if ssu is not None and len(ssu) > 0:
            bin_fastas = glob.glob(os.path.join("${derep_dir}", "*.fa")) + glob.glob(os.path.join("${derep_dir}", "*.fasta"))
            for _, row in ssu.iterrows():
                contig_id = str(row.get("contig_id", row.iloc[0]))
                for bf in bin_fastas:
                    with open(bf) as fh:
                        headers = [l[1:].split()[0] for l in fh if l.startswith(">")]
                    if contig_id in headers:
                        bin_name = os.path.splitext(os.path.basename(bf))[0]
                        df.loc[df["Bin"] == bin_name, "ssuextract_tax"] = row.get("taxonomy", "N/A")
                        df.loc[df["Bin"] == bin_name, "SSU_seq"] = row.get("sequence", "N/A")
                        break

    df.insert(0, "Sample", "${root}")
    df.to_excel("${root}_bin_report.xlsx", index=False)
    """
}
