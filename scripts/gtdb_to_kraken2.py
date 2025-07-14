import pandas as pd
from pathlib import Path

# === 配置路径 ===
summary_file = "all.summary.tsv"
fasta_dir = Path("../07.nr_95/dereplicated_genomes")
output_fasta_dir = Path("kraken2_MAG_db/library/MAGs")
taxonomy_dir = Path("kraken2_MAG_db/taxonomy")

# === 创建输出目录 ===
output_fasta_dir.mkdir(parents=True, exist_ok=True)
taxonomy_dir.mkdir(parents=True, exist_ok=True)

# === 初始化 taxid 表 ===
taxid_counter = 900000000
taxon_to_id = {"d__root": 1}
nodes = []
names = []

# === 读取 summary 文件 ===
df = pd.read_csv(summary_file, sep="\t")
df[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = df['classification'].str.split(';', expand=True)
df = df.apply(lambda col: col.str.strip() if col.dtype == 'object' else col)

# === 处理每条记录 ===
for _, row in df.iterrows():
    genome = row["user_genome"]
    lineage = [row["Domain"], row["Phylum"], row["Class"], row["Order"], row["Family"], row["Genus"], row["Species"]]
    lineage = [x for x in lineage if isinstance(x, str) and x != ""]

    parent_id = 1
    for clade in lineage:
        if clade not in taxon_to_id:
            taxid_counter += 1
            taxon_to_id[clade] = taxid_counter
            rank = clade.split("__")[0].lower()
            nodes.append(f"{taxid_counter}\t|\t{parent_id}\t|\t{rank}\t|")
            names.append(f"{taxid_counter}\t|\t{clade}\t|\t\t|\tscientific name\t|")
        parent_id = taxon_to_id[clade]

    final_taxid = taxon_to_id[lineage[-1]]

    # === 查找 fasta 文件（.fa）===
    input_fa = fasta_dir / f"{genome}.fa"
    output_fa = output_fasta_dir / f"{genome}_tax.fa"

    if input_fa.exists():
        with open(input_fa) as fin, open(output_fa, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    fout.write(f">kraken:taxid|{final_taxid} {genome}\n")
                else:
                    fout.write(line)
    else:
        print(f"[WARN] Fasta not found: {input_fa}")

# === 写出 taxonomy 文件 ===
with open(taxonomy_dir / "nodes.dmp", "w") as f:
    f.write("1\t|\t1\t|\tno rank\t|\n")
    for line in nodes:
        f.write(line + "\n")

with open(taxonomy_dir / "names.dmp", "w") as f:
    f.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
    for line in names:
        f.write(line + "\n")

print("[✔] Kraken2 input fasta and taxonomy files generated successfully.")
