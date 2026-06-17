# SEGC workflow

This document describes the reproducible workflow used to construct the Standardized Extremophilic Genomic Catalogue (SEGC) and to analyse biosynthetic gene clusters (BGCs) from paired-end metagenomic reads.

The workflow starts from a single input directory:

```text
00.rawdata/
├── SampleA_1.fq.gz
├── SampleA_2.fq.gz
├── SampleB_1.fq.gz
└── SampleB_2.fq.gz
```

Read files must be named as `<sample_id>_1.fq.gz` and `<sample_id>_2.fq.gz`. All output directories are created by the commands below.

## Software

The versions below match the manuscript where applicable.

```text
Trim Galore v0.5.0
MEGAHIT v1.1.3
Seqtk v1.5.0
Bowtie2 v2.3.5 or later
SAMtools v1.9.0 or later
BEDTools
MetaBAT2 v2.12.1
MaxBin2 v2.2.6
CONCOCT v1.0.0
DAS Tool v1.1.7
CheckM v1.2.3
CheckM2 v1.0.2
Barrnap v0.9
tRNAscan-SE v2.0.12
dRep v3.4.5
GTDB-Tk v2.4.0 with GTDB Release R226
BMGE v1.12
FastTree v2.1.10
antiSMASH v7.0
BiG-SCAPE v2.0 / BiG-SLiCE v2.0
skani v0.3.0
fastANI v1.33
Prodigal v2.6.3
PfamScan
CLEAN v1.0.1
GNU parallel
```

Before running, set project-level variables:

```bash
set -euo pipefail

THREADS=${THREADS:-60}
JOBS=${JOBS:-6}

CHECKM2_DB=${CHECKM2_DB:-/path/to/CheckM2_database/uniref100.KO.1.dmnd}
GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:-/path/to/gtdbtk/release226}
PFAM_DB=${PFAM_DB:-/path/to/Pfam-A.hmm}
PFAMSCAN_DB=${PFAMSCAN_DB:-/path/to/pfam_scan_db}

export THREADS JOBS CHECKM2_DB GTDBTK_DATA_PATH PFAM_DB PFAMSCAN_DB
```

## 1. Prepare sample list

```bash
mkdir -p metadata

find 00.rawdata -name "*_1.fq.gz" -type f \
  | sed 's#^00.rawdata/##; s#_1.fq.gz$##' \
  | sort > metadata/sample_ids.txt

while read -r id; do
  test -s "00.rawdata/${id}_1.fq.gz"
  test -s "00.rawdata/${id}_2.fq.gz"
done < metadata/sample_ids.txt
```

## 2. Quality control

Raw reads are quality-filtered with Trim Galore to remove adapters and low-quality bases.

```bash
mkdir -p 01.clean_reads

parallel -j "${JOBS}" \
  'trim_galore --paired --gzip -q 20 -o 01.clean_reads \
    00.rawdata/{1}_1.fq.gz 00.rawdata/{1}_2.fq.gz' \
  :::: metadata/sample_ids.txt
```

## 3. Per-sample assembly

Each sample is assembled independently with MEGAHIT.

```bash
mkdir -p 02.assembly_megahit

parallel -j "${JOBS}" \
  'megahit \
    -1 01.clean_reads/{1}_1_val_1.fq.gz \
    -2 01.clean_reads/{1}_2_val_2.fq.gz \
    -t "${THREADS}" \
    -o 02.assembly_megahit/{1}' \
  :::: metadata/sample_ids.txt
```

## 4. Filter contigs

Only contigs longer than 1,500 bp are retained for binning and downstream analyses.

```bash
mkdir -p 03.contigs_1500

parallel -j "${JOBS}" \
  'seqtk seq -L 1500 02.assembly_megahit/{1}/final.contigs.fa \
    > 03.contigs_1500/{1}.contigs.1500.fa' \
  :::: metadata/sample_ids.txt
```

## 5. Map clean reads to contigs

Clean reads are mapped back to their own assembly once. The sorted BAM, MetaBAT2 depth table, MaxBin2 coverage table and CONCOCT coverage table generated here are reused by the three binning algorithms.

```bash
mkdir -p \
  04.mapping/index \
  04.mapping/sam \
  04.mapping/bam \
  04.mapping/depth \
  04.mapping/maxbin_coverage \
  04.mapping/stat

while read -r id; do
  bowtie2-build \
    -f "03.contigs_1500/${id}.contigs.1500.fa" \
    --threads "${THREADS}" \
    "04.mapping/index/${id}.contigs.1500"

  bowtie2 \
    -x "04.mapping/index/${id}.contigs.1500" \
    -1 "01.clean_reads/${id}_1_val_1.fq.gz" \
    -2 "01.clean_reads/${id}_2_val_2.fq.gz" \
    -p "${THREADS}" \
    -S "04.mapping/sam/${id}.sam" \
    2> "04.mapping/stat/${id}.bowtie2.stat"

  samtools view \
    -@ "${THREADS}" \
    -b \
    -S "04.mapping/sam/${id}.sam" \
    -o "04.mapping/bam/${id}.bam"

  samtools sort \
    -@ "${THREADS}" \
    -l 9 \
    -O BAM \
    "04.mapping/bam/${id}.bam" \
    -o "04.mapping/bam/${id}.sorted.bam"

  samtools index \
    "04.mapping/bam/${id}.sorted.bam" \
    -@ "${THREADS}"

  jgi_summarize_bam_contig_depths \
    --outputDepth "04.mapping/depth/${id}.depth.txt" \
    "04.mapping/bam/${id}.sorted.bam"

  genomeCoverageBed \
    -ibam "04.mapping/bam/${id}.sorted.bam" \
    > "04.mapping/maxbin_coverage/${id}.histogram.tab"

  python ./calculate-contig-coverage.py \
    "04.mapping/maxbin_coverage/${id}.histogram.tab"

  rm -f "04.mapping/sam/${id}.sam" "04.mapping/bam/${id}.bam"
done < metadata/sample_ids.txt
```

## 6. Genome binning

Three complementary binning algorithms are used: MetaBAT2, MaxBin2 and CONCOCT.

### 6.1 MetaBAT2

```bash
mkdir -p 05.binning/metabat2

while read -r id; do
  mkdir -p "05.binning/metabat2/${id}"

  metabat2 \
    -m 1500 \
    -t "${THREADS}" \
    -i "03.contigs_1500/${id}.contigs.1500.fa" \
    -a "04.mapping/depth/${id}.depth.txt" \
    -o "05.binning/metabat2/${id}/${id}.metabat2" \
    -v
done < metadata/sample_ids.txt
```

### 6.2 MaxBin2

```bash
mkdir -p 05.binning/maxbin2

while read -r id; do
  mkdir -p "05.binning/maxbin2/${id}"

  run_MaxBin.pl \
    -contig "03.contigs_1500/${id}.contigs.1500.fa" \
    -abund "04.mapping/maxbin_coverage/${id}.histogram.tab.coverage.tab" \
    -max_iteration 50 \
    -out "05.binning/maxbin2/${id}/${id}.maxbin2" \
    -thread "${THREADS}"
done < metadata/sample_ids.txt
```

### 6.3 CONCOCT

```bash
mkdir -p 05.binning/concoct 05.binning/concoct_work

while read -r id; do
  mkdir -p "05.binning/concoct/${id}" "05.binning/concoct_work/${id}"

  cut_up_fasta.py \
    "03.contigs_1500/${id}.contigs.1500.fa" \
    -c 10000 \
    -o 0 \
    --merge_last \
    -b "05.binning/concoct_work/${id}/${id}.contigs_10K.bed" \
    > "05.binning/concoct_work/${id}/${id}.contigs_10K.fa"

  concoct_coverage_table.py \
    "05.binning/concoct_work/${id}/${id}.contigs_10K.bed" \
    "04.mapping/bam/${id}.sorted.bam" \
    > "05.binning/concoct_work/${id}/${id}.coverage_table.tsv"

  concoct \
    --composition_file "05.binning/concoct_work/${id}/${id}.contigs_10K.fa" \
    --coverage_file "05.binning/concoct_work/${id}/${id}.coverage_table.tsv" \
    -b "05.binning/concoct_work/${id}/${id}.concoct_output" \
    --threads "${THREADS}"

  merge_cutup_clustering.py \
    "05.binning/concoct_work/${id}/${id}.concoct_output_clustering_gt1000.csv" \
    > "05.binning/concoct_work/${id}/${id}.clustering_merged.csv"

  extract_fasta_bins.py \
    "03.contigs_1500/${id}.contigs.1500.fa" \
    "05.binning/concoct_work/${id}/${id}.clustering_merged.csv" \
    --output_path "05.binning/concoct/${id}"
done < metadata/sample_ids.txt
```

## 7. Integrate bins with DAS Tool

```bash
mkdir -p 06.das_tool/scaffolds2bin 06.das_tool/results 07.MAGs/raw_bins

parallel -j "${JOBS}" \
  'Fasta_to_Contig2Bin.sh -i 05.binning/maxbin2/{1} -e fasta \
    > 06.das_tool/scaffolds2bin/{1}.maxbin2.tsv
   Fasta_to_Contig2Bin.sh -i 05.binning/metabat2/{1} -e fa \
    > 06.das_tool/scaffolds2bin/{1}.metabat2.tsv
   Fasta_to_Contig2Bin.sh -i 05.binning/concoct/{1} -e fa \
    > 06.das_tool/scaffolds2bin/{1}.concoct.tsv
   DAS_Tool \
    -i 06.das_tool/scaffolds2bin/{1}.maxbin2.tsv,06.das_tool/scaffolds2bin/{1}.metabat2.tsv,06.das_tool/scaffolds2bin/{1}.concoct.tsv \
    -l maxbin2,metabat2,concoct \
    -c 03.contigs_1500/{1}.contigs.1500.fa \
    -o 06.das_tool/results/{1}.DASTool \
    --threads "${THREADS}" \
    --write_bins \
    --score_threshold 0' \
  :::: metadata/sample_ids.txt

find 06.das_tool/results -path "*_DASTool_bins/*.fa" -type f \
  | while read -r bin; do
      sample=$(basename "$(dirname "$bin")" _DASTool_bins)
      cp "$bin" "07.MAGs/raw_bins/${sample}__$(basename "$bin")"
    done
```

## 8. MAG quality assessment

MAG quality is assessed by both CheckM1 and CheckM2. Medium-quality MAGs are retained only if both tools report completeness >=50% and contamination <=10%.

```bash
mkdir -p 08.quality/checkm1 08.quality/checkm2 09.MAGs/medium_quality

checkm lineage_wf \
  -x fa \
  -t "${THREADS}" \
  --tmpdir 08.quality/checkm1/tmp \
  07.MAGs/raw_bins \
  08.quality/checkm1

checkm qa \
  08.quality/checkm1/lineage.ms \
  08.quality/checkm1 \
  --tab_table \
  -o 2 \
  -f 08.quality/checkm1/checkm1_qa.tsv

checkm2 predict \
  --threads "${THREADS}" \
  -x fa \
  -i 07.MAGs/raw_bins \
  -o 08.quality/checkm2 \
  --database_path "${CHECKM2_DB}"

python3 - <<'PY'
from pathlib import Path
import csv
import shutil

raw_dir = Path("07.MAGs/raw_bins")
out_dir = Path("09.MAGs/medium_quality")
out_dir.mkdir(parents=True, exist_ok=True)

checkm1 = {}
with open("08.quality/checkm1/checkm1_qa.tsv", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        name = row.get("Bin Id") or row.get("Bin")
        if not name:
            continue
        checkm1[name] = (
            float(row["Completeness"]),
            float(row["Contamination"]),
        )

checkm2 = {}
with open("08.quality/checkm2/quality_report.tsv", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        name = row.get("Name") or row.get("Bin Id")
        if not name:
            continue
        checkm2[name] = (
            float(row["Completeness"]),
            float(row["Contamination"]),
        )

kept = []
for fasta in raw_dir.glob("*.fa"):
    name = fasta.stem
    if name not in checkm1 or name not in checkm2:
        continue
    c1, x1 = checkm1[name]
    c2, x2 = checkm2[name]
    if c1 >= 50 and x1 <= 10 and c2 >= 50 and x2 <= 10:
        shutil.copy2(fasta, out_dir / fasta.name)
        kept.append(name)

Path("09.MAGs/medium_quality_MAGs.txt").write_text("\n".join(sorted(kept)) + "\n")
print(f"Retained {len(kept)} medium-quality MAGs")
PY
```

High-quality MAGs can be identified with the MIMAG criteria: completeness >90%, contamination <5%, at least 18 tRNAs, and 5S/16S/23S rRNAs.

```bash
mkdir -p 08.quality/barrnap 08.quality/trnascan

find 09.MAGs/medium_quality -name "*.fa" -type f \
  | sed 's#^09.MAGs/medium_quality/##; s#.fa$##' \
  | sort > metadata/medium_quality_MAG_ids.txt

parallel -j "${JOBS}" \
  'barrnap --kingdom bac --reject 0.01 --evalue 1e-3 \
    09.MAGs/medium_quality/{1}.fa \
    > 08.quality/barrnap/{1}.bac.gff || true
   barrnap --kingdom arc --reject 0.01 --evalue 1e-3 \
    09.MAGs/medium_quality/{1}.fa \
    > 08.quality/barrnap/{1}.arc.gff || true
   tRNAscan-SE -B -o 08.quality/trnascan/{1}.bac.tsv \
    09.MAGs/medium_quality/{1}.fa || true
   tRNAscan-SE -A -o 08.quality/trnascan/{1}.arc.tsv \
    09.MAGs/medium_quality/{1}.fa || true' \
  :::: metadata/medium_quality_MAG_ids.txt
```

## 9. Species-level genome bins and taxonomy

Medium-quality MAGs are dereplicated into species-level genome bins (SGBs) at 95% ANI.

```bash
mkdir -p 10.SGBs

dRep dereplicate \
  10.SGBs/dRep95 \
  -g 09.MAGs/medium_quality/*.fa \
  -p "${THREADS}" \
  --ignoreGenomeQuality \
  -pa 0.9 \
  -sa 0.95 \
  -nc 0.3

mkdir -p 10.SGBs/representatives
cp 10.SGBs/dRep95/dereplicated_genomes/*.fa 10.SGBs/representatives/

gtdbtk classify_wf \
  --genome_dir 10.SGBs/representatives \
  --out_dir 11.taxonomy_gtdbtk \
  --cpus "${THREADS}" \
  --pplacer_cpus "${THREADS}" \
  --skip_ani_screen \
  --extension fa
```

Novel SGBs are defined as representative genomes with ANI <95% and alignment fraction <30% relative to GTDB reference genomes. Habitat-specific SGBs are SGBs whose member MAGs all originate from the same habitat.

## 10. Phylogenetic trees

GTDB-Tk marker alignments are trimmed with BMGE and trees are inferred with FastTree.

```bash
mkdir -p 12.phylogeny

for aln in 11.taxonomy_gtdbtk/align/gtdbtk.bac120.user_msa.fasta.gz \
           11.taxonomy_gtdbtk/align/gtdbtk.ar53.user_msa.fasta.gz; do
  if [ -s "$aln" ]; then
    gzip -dc "$aln" > "12.phylogeny/$(basename "$aln" .gz)"
  fi
done

if [ -s 12.phylogeny/gtdbtk.bac120.user_msa.fasta ]; then
  bmge -i 12.phylogeny/gtdbtk.bac120.user_msa.fasta \
    -t AA -g 0.5 -h 1 -b 1 -w 1 \
    -of 12.phylogeny/gtdbtk.bac120.trimmed.fasta
  FastTree 12.phylogeny/gtdbtk.bac120.trimmed.fasta \
    > 12.phylogeny/bac120.SEGC.tree
fi

if [ -s 12.phylogeny/gtdbtk.ar53.user_msa.fasta ]; then
  bmge -i 12.phylogeny/gtdbtk.ar53.user_msa.fasta \
    -t AA -g 0.5 -h 1 -b 1 -w 1 \
    -of 12.phylogeny/gtdbtk.ar53.trimmed.fasta
  FastTree 12.phylogeny/gtdbtk.ar53.trimmed.fasta \
    > 12.phylogeny/ar53.SEGC.tree
fi
```

## 11. BGC discovery with antiSMASH

antiSMASH is run on all retained medium-quality MAGs.

```bash
mkdir -p 13.BGC/antismash 13.BGC/gbk

parallel -j "${JOBS}" \
  'antismash \
    09.MAGs/medium_quality/{1}.fa \
    --taxon bacteria \
    --output-dir 13.BGC/antismash/{1} \
    --genefinding-tool prodigal \
    --cb-knownclusters \
    --cc-mibig \
    --fullhmmer \
    -c 1' \
  :::: metadata/medium_quality_MAG_ids.txt

find 13.BGC/antismash -name "*.region*.gbk" -type f \
  | while read -r gbk; do
      mag=$(basename "$(dirname "$gbk")")
      cp "$gbk" "13.BGC/gbk/${mag}__$(basename "$gbk")"
    done
```

## 12. GCF clustering with BiG-SCAPE

BGCs are grouped into gene cluster families (GCFs). The manuscript used c = 0.3 and c = 0.7 for GCF-level analyses.

```bash
mkdir -p 14.GCF

for cutoff in 0.3 0.7; do
  bigscape cluster \
    -i 13.BGC/gbk \
    -o "14.GCF/bigscape_c${cutoff}" \
    -p "${PFAM_DB}" \
    -c "${THREADS}" \
    --gcf-cutoffs "${cutoff}" \
    --mix \
    --alignment-mode local \
    --extend-strategy greedy \
    --classify category \
    --include-singletons
done
```

GCFs containing BGCs from one habitat are defined as habitat-specific GCFs. GCFs containing BGCs from more than one habitat are defined as multi-habitat GCFs.

## 13. Comparison with external genome catalogues

This step compares SEGC SGB representatives with external genome catalogues using skani for candidate search and fastANI for final ANI estimation.

Prepare one genome list per external catalogue:

```text
external_catalogues/
├── aquatic.genome_list.txt
├── crop_root.genome_list.txt
├── human.genome_list.txt
├── ocean.genome_list.txt
├── soil.genome_list.txt
└── EEMC.genome_list.txt
```

Each list should contain absolute or project-relative FASTA paths, one genome per line.

```bash
mkdir -p 15.catalogue_compare/skani_refs 15.catalogue_compare/skani_hits 15.catalogue_compare/fastani

find 10.SGBs/representatives -name "*.fa" -type f \
  | sort > 15.catalogue_compare/SEGC_SGB_representatives.txt

for list in external_catalogues/*.genome_list.txt; do
  name=$(basename "$list" .genome_list.txt)
  skani sketch -l "$list" -o "15.catalogue_compare/skani_refs/${name}" -t "${THREADS}" --slow
  skani search \
    -d "15.catalogue_compare/skani_refs/${name}" \
    -l 15.catalogue_compare/SEGC_SGB_representatives.txt \
    -o "15.catalogue_compare/skani_hits/${name}.skani.tsv" \
    -t "${THREADS}"
done
```

For final species assignment, recompute ANI for candidate matches with fastANI using `--minFraction 0.3 --fragLen 1500`. SGBs sharing >=95% ANI with any genome in a reference catalogue are considered present in that catalogue.

## 14. Terpene pathway reconstruction

Terpene BGC proteins are annotated with CLEAN at medium confidence (confidence >=0.2). Predicted EC numbers are mapped to curated terpene marker genes from KEGG and MetaCyc to identify terpene classes and specific pathways such as retinal biosynthesis.

```bash
mkdir -p 17.terpene/proteins 17.terpene/clean

# Example: run CLEAN on a FASTA file of terpene BGC proteins.
CLEAN_infer_fasta.py \
  --fasta_data 17.terpene/proteins/SEGC.terpene.proteins.fasta \
  --out_dir 17.terpene/clean
```

## 15. Metatranscriptomic validation

For matched metatranscriptomes, raw RNA reads are processed with Trim Galore. Clean reads are mapped to beta-carotene dioxygenase and bacteriorhodopsin reference genes using Bowtie2, and expression is quantified as FPKM normalized by total raw read counts.

Input layout:

```text
18.metatranscriptome/
├── 00.rawdata/
│   ├── RNA_sampleA_1.fastq.gz
│   └── RNA_sampleA_2.fastq.gz
└── references/
    ├── beta_carotene_dioxygenase.fa
    └── bacteriorhodopsin.fa
```

Commands:

```bash
cd 18.metatranscriptome

mkdir -p 01.references 02.clean_reads 03.bowtie 04.stat 05.counts
cat references/beta_carotene_dioxygenase.fa references/bacteriorhodopsin.fa \
  > 01.references/retinal_phototrophy_genes.fa

find 00.rawdata -name "*_1.fastq.gz" -type f \
  | sed 's#^00.rawdata/##; s#_1.fastq.gz$##' \
  | sort > sample_ids.txt

parallel -j "${JOBS}" \
  'trim_galore --paired --gzip -q 20 -o 02.clean_reads \
    00.rawdata/{1}_1.fastq.gz 00.rawdata/{1}_2.fastq.gz' \
  :::: sample_ids.txt

bowtie2-build -f 01.references/retinal_phototrophy_genes.fa 01.references/retinal_phototrophy_genes

while read -r id; do
  raw_lines=$(zcat "02.clean_reads/${id}_1_val_1.fq.gz" | wc -l)
  total_reads=$((raw_lines / 4))

  bowtie2 \
    -x 01.references/retinal_phototrophy_genes \
    -1 "02.clean_reads/${id}_1_val_1.fq.gz" \
    -2 "02.clean_reads/${id}_2_val_2.fq.gz" \
    -S "03.bowtie/${id}.sam" \
    -p "${THREADS}" \
    --very-sensitive-local \
    --no-mixed \
    --no-discordant \
    2> "04.stat/${id}.bowtie2.stat"

  samtools view -@ "${THREADS}" -bS "03.bowtie/${id}.sam" \
    | samtools sort -@ "${THREADS}" -o "03.bowtie/${id}.sorted.bam"
  samtools index "03.bowtie/${id}.sorted.bam"
  samtools idxstats "03.bowtie/${id}.sorted.bam" > "05.counts/${id}.idxstats.tsv"

  awk -v total="${total_reads}" 'BEGIN{OFS="\t"; print "gene","length_bp","count","FPKM"}
    $1 != "*" {
      fpkm = ($3 * 1e9) / (total * $2)
      print $1, $2, $3, fpkm
    }' "05.counts/${id}.idxstats.tsv" > "05.counts/${id}.FPKM.tsv"
done < sample_ids.txt

cd ..
```

## Output summary

```text
01.clean_reads/                 quality-filtered reads
02.assembly_megahit/            per-sample assemblies
03.contigs_1500/                contigs longer than 1,500 bp
04.mapping/                     read mapping, BAM files and depth files
05.binning/                     MetaBAT2, MaxBin2 and CONCOCT bins
06.das_tool/                    integrated bins
09.MAGs/medium_quality/         MAGs passing CheckM1 and CheckM2 thresholds
10.SGBs/representatives/        representative species-level genome bins
11.taxonomy_gtdbtk/             GTDB-Tk taxonomy
12.phylogeny/                   bacterial and archaeal marker-gene trees
13.BGC/                         antiSMASH BGC predictions
14.GCF/                         BiG-SCAPE GCF clustering
15.catalogue_compare/           comparison with external genome catalogues
16.terpene/                     terpene pathway reconstruction
17.metatranscriptome/           RNA validation of retinal-based phototrophy
```

## Citation and data

SEGC was constructed from 1,462 metagenomic samples spanning seven extreme habitats: acid mine, cryosphere, deep sea, hot spring, hydrothermal plume, saline-alkaline and subsurface. The manuscript reports 54,661 medium-quality MAGs, 21,805 SGBs and 162,855 BGCs.

Raw metagenomic sequencing data for the 115 in-house saline-alkaline soil samples are available under NCBI accession PRJNA1285087. MAG FASTA files and BGC GenBank files are available from Zenodo: https://doi.org/10.5281/zenodo.15788452.
