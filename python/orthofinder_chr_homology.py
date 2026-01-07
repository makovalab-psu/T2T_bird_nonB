import pandas as pd
from collections import Counter, defaultdict

# -------------------------
# INPUT FILES
# -------------------------
ORTHOGROUPS = "orthofinder/OrthoFinder/Results_Jan05/Orthogroups/Orthogroups.tsv"

gene2chr_files = {
    "bTaeGut7v0.4_MT_rDNA.noDots.proteins": "orthofinder/bTaeGut7v0.4_MT_rDNA_gene2chr.tsv",
    "chicken.v23.noDots.proteins": "orthofinder/chicken.v23_gene2chr.tsv",
    "bCalAnn1_v1.p.noDots.proteins": "orthofinder/bCalAnn1_v1.p_gene2chr.tsv",
    "CAU-Wild1.1.noDots.proteins": "orthofinder/CAU-Wild1.1_gene2chr.tsv",
    "OTswu.noDots.proteins": "orthofinder/OTswu_gene2chr.tsv",
    "bStrUra1.noDots.proteins": "orthofinder/bStrUra1_gene2chr.tsv",
}

# -------------------------
# LOAD GENE → CHR MAPS
# -------------------------
gene2chr = {}
for sp, fn in gene2chr_files.items():
    df = pd.read_csv(fn, sep="\t", names=["gene", "chr"])
    gene2chr[sp] = dict(zip(df["gene"], df["chr"]))

# -------------------------
# LOAD ORTHOGROUPS
# -------------------------
og = pd.read_csv(ORTHOGROUPS, sep="\t")

# -------------------------
# COUNT SHARED ORTHOGROUPS
# -------------------------
# homology_counts[sp1][chr1][sp2][chr2] += 1
homology_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(Counter)))

for _, row in og.iterrows():
    chr_by_species = {}

    for sp in gene2chr:
        if pd.isna(row[sp]):
            continue
        genes = row[sp].split(", ")
        chrs = [gene2chr[sp].get(g) for g in genes if g in gene2chr[sp]]
        if chrs:
            # Use most common chromosome if multiple genes
            chr_by_species[sp] = Counter(chrs).most_common(1)[0][0]

    for sp1, chr1 in chr_by_species.items():
        for sp2, chr2 in chr_by_species.items():
            if sp1 != sp2:
                homology_counts[sp1][chr1][sp2][chr2] += 1

# -------------------------
# OUTPUT HOMOLOGY TABLES
# -------------------------
for sp1 in homology_counts:
    for sp2 in homology_counts[sp1][next(iter(homology_counts[sp1]))]:
        rows = []
        for chr1 in homology_counts[sp1]:
            for chr2, count in homology_counts[sp1][chr1][sp2].items():
                rows.append((chr1, chr2, count))

        df = pd.DataFrame(rows, columns=[f"{sp1}_chr", f"{sp2}_chr", "shared_orthogroups"])
        df.sort_values("shared_orthogroups", ascending=False, inplace=True)

        out = f"orthofinder/{sp1}_vs_{sp2}_chr_homology.tsv"
        df.to_csv(out, sep="\t", index=False)

