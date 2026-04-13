# Standardized extremophilic genomic catalogueue (SEGC)


## 1. Quality control

```sh
## Only paired-end raw sequencing data were retained. The raw reads were quality-filtered using Trim Galore to obtain high-quality clean reads
for i in $(cat sample_id); do
trim_galore -o 01.clean.data --gzip --paired 00.raw_data/${i}_1.fastq.gz  00.raw_data/${i}_2.fastq.gz
done

```

## 2. *de novo* assembly

```SH
##MEGAHIT was used to assemble the clean data, as it is particularly suitable for large-sample analyses
for i in $(cat sample_id); do
megahit -1 01.clean.data/${i}_1_val_1.fq.gz -2 01.clean.data/${i}_2_val_2.fq.gz -t 80 -o 02.assembly_megahit/${i}_assembly
done
```

## 3. Contig length filtering

```sh
#First, the assembled contigs were extracted, and then length filtering was performed using Seqtk with a minimum length threshold of 1,500 bp
##
for i in $(cat sample_id); do 
cp 02.assembly_megahit/${i}_assembly/final.contigs.fa 03.seqtk_results/${i}.fa
done
##
for i in $(cat sample_id); do
seqtk seq -L 1500 03.seqtk_results/${i}.fa > 03.seqtk_results/${i}_l.fa
done
```

## 4. Binning analysis

```sh
## Metagenomic binning was performed using three different methods: MetaBAT2, MaxBin, and CONCOCT.
for i in $(cat sample_id); do 
mkdir 05.binning/metabat2_bin/${i}
mkdir 05.binning/MaxBin_bin/${i}
mkdir 05.binning/concoct_bin/${i}
##metabat2
bowtie2-build -f 03.seqtk_results/${i}_l.fa --threads 10 05.binning/MaxBin_bin/${i}/${i}_final
bowtie2 -1 01.clean.data/${i}_1_val_1.fq.gz -2 01.clean.data/${i}_2_val_2.fq.gz -p 10 -x 05.binning/MaxBin_bin/${i}/${i}_final -S 05.binning/MaxBin_bin/${i}/${i}_final.sam 2>0_SAM/${i}.bowtie2.stat
samtools view -@ 10 -b -S 05.binning/MaxBin_bin/${i}/${i}_final.sam -o 05.binning/MaxBin_bin/${i}/${i}_final.bam
samtools sort -@ 10 -l 9 -O BAM 05.binning/MaxBin_bin/${i}/${i}_final.bam -o 05.binning/MaxBin_bin/${i}/${i}_final.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth 05.binning/MaxBin_bin/${i}/${i}_final.depth.txt 05.binning/MaxBin_bin/${i}/${i}_final.sorted.bam
metabat2 -m 1500 -t 10 -i 03.seqtk_results/${i}_l.fa -a 05.binning/MaxBin_bin/${i}/${i}_final.depth.txt -o 05.binning/metabat2_bin/${i}/${i}_metabat2 -v
##MaxBin
genomeCoverageBed -ibam 05.binning/MaxBin_bin/${i}/${i}_final.sorted.bam > 05.binning/MaxBin_bin/${i}/${i}.histogram.tab
python ./calculate-contig-coverage.py 05.binning/MaxBin_bin/${i}/${i}.histogram.tab
run_MaxBin.pl -contig 03.seqtk_results/${i}_l.fa -abund 05.binning/MaxBin_bin/${i}/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 05.binning/MaxBin_bin/${i}/${i}_MaxBin -thread 10
##concoct
cut_up_fasta.py 03.seqtk_results/${i}_l.fa -c 10000 -o 0 --merge_last -b 05.binning/MaxBin_bin/${i}/${i}.contigs_10K.bed > 05.binning/MaxBin_bin/${i}/${i}.contigs_10K.fa
samtools index 05.binning/MaxBin_bin/${i}/${i}_final.sorted.bam -@ 10
concoct_coverage_table.py 05.binning/MaxBin_bin/${i}/${i}.contigs_10K.bed 05.binning/MaxBin_bin/${i}/${i}_final.sorted.bam > 05.binning/MaxBin_bin/${i}/${i}.coverage_table.tsv
concoct --composition_file 05.binning/MaxBin_bin/${i}/${i}.contigs_10K.fa --coverage_file 05.binning/MaxBin_bin/${i}/${i}.coverage_table.tsv -b 05.binning/MaxBin_bin/${i}/${i}.concoct_output --threads 10
merge_cutup_clustering.py 05.binning/MaxBin_bin/${i}/${i}.concoct_output_clustering_gt1000.csv > 05.binning/MaxBin_bin/${i}/${i}.clustering_merged.csv
extract_fasta_bins.py 03.seqtk_results/${i}_l.fa 05.binning/MaxBin_bin/${i}/${i}.clustering_merged.csv --output_path 05.binning/concoct_bin/${i}
done

```

## 5. Refinement

```sh
## DAS_Tool was used to integrate and de-replicate the binning results obtained from the different tools
for i in $(cat sample_id); do
Fasta_to_Contig2Bin.sh -i 04.binning/MaxBin_bin/${i} -e fasta > 06.das_tool/${i}_maxbin.scaffolds2bin.tsv
Fasta_to_Contig2Bin.sh -i 04.binning/metabat2_bin/${i} -e fa > 06.das_tool/${i}_metabat2.scaffolds2bin.tsv
Fasta_to_Contig2Bin.sh -i 04.binning/concoct_bin/${i} -e fa > 06.das_tool/${i}_concoct.scaffolds2bin.tsv
done

for i in $(cat sample_id); do
DAS_Tool -i 06.das_tool/${i}_maxbin.scaffolds2bin.tsv,06.das_tool/${i}_metabat2.scaffolds2bin.tsv,06.das_tool/${i}_concoct.scaffolds2bin.tsv -l maxbin,metabat,concoct -c 03.seqtk_results/${i}_l.fa -o 06.das_tool_results/${i} --threads 40 --write_bins --score_threshold 0
done

```

## 6. MAG Quality Assessment

```sh
# MAG quality was assessed using both CheckM and CheckM2. Initially, MAGs were filtered based on quality estimates using CheckM2. Those that met the thresholds were subsequently validated with CheckM1. Only MAGs that passed the quality criteria of both CheckM2 and CheckM1 were retained for downstream analyses

##checkM2
checkm2 predict --threads 40 -x fa -i 07.MAG_all -o 08.checkM2
#checkM1
concat () {
checkm lineage_wf -x fa 09.MAG_checkM2 10.checkm1/${1} -t 20 --tmpdir bin_checkm.tmp --pplacer_threads 20
}
export -f concat
cat tmp | parallel -j 4 concat {}
```

## 7. Species-level Dereplication with dRep

```sh
##Species-level dereplication was conducted using dRep with a 95% average nucleotide identity (ANI) threshold
dRep dereplicate 11.nr_95/ -g genome_list.txt -p 90 --ignoreGenomeQuality -pa 0.90 -sa 0.95 --S_algorithm fastANI -nc 0.3 --multiround_primary_clustering --primary_chunksize 10000 --genomeInfo genomeInfo.csv
```

## 8. Taxonomic Assignment Using GTDB-Tk

```sh
#The taxonomic affiliation of the dereplicated MAGs was determined using GTDB-Tk (Genome Taxonomy Database Toolkit), which assigns standardized taxonomy based on the GTDB reference database.
gtdbtk classify_wf --genome_dir 11.nr_95/dereplicated_genomes/ --out_dir 12.taxa_gtdb --cpus 85 --pplacer_cpus 85 --skip_ani_screen --extension fa
```

## 9.tRNA and rRNA Identification

```sh
###rRNA
for i in $(cat bac.MAG.id); do
barrnap --kingdom bac --threads 50 --outseq 13.RNA/01.rRNA/${i}_rRNA.fasta --quiet 04.MAG/${i}.fa --reject 0.01 –e-value 1e-3  > 13.RNA/01.rRNA/${i}_rRNA.gff3
done
for i in $(cat arc.MAG.id); do
barrnap --kingdom arc --threads 50 --outseq 13.RNA/01.rRNA/${i}_rRNA.fasta --quiet 04.MAG/${i}.fa --reject 0.01 –e-value 1e-3  > 13.RNA/01.rRNA/${i}_rRNA.gff3
done
###tRNA
concat () {
tRNAscan-SE -A -o 13.RNA/02.tRNA/${1}_tRNA.out -f 13.RNA/02.tRNA/${1}_tRNA.ss -m 13.RNA/02.tRNA/${1}_tRNA.stats 04.MAG/${1}.fa 
}
export -f concat
cat MAG_id | parallel --tmpdir ~/autodl-tmp/123 -j 30 concat {}

```

## 10.Prediction of Microbial Optimal Growth Temperatures

```sh
##To infer the thermal preferences of the recovered MAGs, the optimal growth temperature (OGT) for each genome was predicted using MetaThermo.
concat () {
python ./meta-thermo-main/metathermo.py -f 04.MAG/${1}.fa -t fna
mv MPT_output.csv 14.metathermo.out/${1}.thermo.csv
}
export -f concat
cat MAG_id | parallel -j 50 concat {}
```

## 11.Phylogenetic Analysis

```sh
#Phylogenetic trees were constructed based on concatenated protein sequences generated by GTDB-Tk. The phylogenies included archaeal and bacterial MAGs from the SEMGC dataset, as well as reference genomes from the GTDB. Separate phylogenetic trees were generated for bacteria and archaea by combining SEMGC MAGs with corresponding GTDB reference species.
##cope file
cp 12.taxa_gtdb/align/gtdbtk.bac120.user_msa.fasta.gz 15.tree/
cp 12.taxa_gtdb/align/gtdbtk.ar53.user_msa.fasta.gz 15.tree/
cp 12.taxa_gtdb/align/gtdbtk.ar53.msa.fasta.gz 15.tree/
cp 12.taxa_gtdb/align/gtdbtk.bac120.msa.fasta.gz 15.tree/
##Alignment trimming was performed using BMGE
bmge -i gtdbtk.ar53.user_msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of gtdbtk.ar53_trimmed.fasta
bmge -i gtdbtk.bac120.user_msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of gtdbtk.bac120_trimmed.fasta
bmge -i gtdbtk.ar53.msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of gtdbtk.ar53_msa_trimmed.fasta
bmge -i gtdbtk.bac120.msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of gtdbtk.bac120_msa_trimmed.fasta
##Phylogenetic trees were constructed using FastTree
FastTree gtdbtk.ar53_trimmed.fasta > ar53_SEMGC_gtdb.tree
FastTree gtdbtk.bac120_trimmed.fasta > bac120_SEMGC_gtdb.tree
FastTree gtdbtk.ar53_msa_trimmed.fasta > ar53_msa_gtdb.tree
FastTree gtdbtk.bac120_msa_trimmed.fasta > bac120_msa_gtdb.tree
```

## 12.Pangenome Analysis

```sh
#Gene families were analyzed through a pangenome approach and subsequent clustering. Coding sequences were first predicted for each genome using Prodigal (meta mode), generating protein, nucleotide, and GFF files for downstream analyses. Homologous gene families were then clustered using OrthoFinder for each species group. To examine gene family evolution, CAFE5 was used to estimate gene family expansion and contraction events based on gene count matrices and ultrametric species trees. For single-nucleotide variant (SNV) detection, whole-genome alignments were performed using nucmer with an all-vs-all strategy. The resulting SNV positions were summarized with custom scripts to generate reference-based SNP matrices. To assess selection pressure at the gene level, SNV sites were annotated relative to gene models, and pN/pS ratios were calculated for each gene using a custom pipeline. For each genome, annotated SNVs were mapped onto coding regions using the GFF annotations and reference FASTA files. Subsequently, gene IDs and corresponding protein sequences were appended to the pN/pS tables. Genes with elevated pN/pS ratios were extracted and functionally annotated using eggNOG-mapper.
##Gene prediction
concat () {
prodigal -i ../04.MAG/${1}.fa -a 01.protein/${1}.proteins.faa -d 02.gene/${1}.gene.fasta -p meta -f gff > 03.gff/${1}.gff
}
export -f concat
cat ../genome.id | parallel -j 86 concat {}
##OrthoFinder: homologous gene family clustering
concat () {
orthofinder -t 8 -f 04.species.protein/${1}/ -o 05.orthofinder/${1} -a 8
}
export -f concat
cat species.3_id | parallel -j 10 concat {}
##CAFE5: gene family expansion and contraction analysis
for i in $(cat species.3_id); do
cafe5 -i 06.add_desc/${i}.orthogroups_counts_with_desc.tsv -t 07.tree/${i}.SpeciesTree_ultrametric.txt -o 08.cafe5/${i}.out -c 0 -p -k 2
done
##Whole-genome alignment using nucmer
concat () {
./all_vs_all_nucmer.py -g ../03.genome/${1}/ -o 09.mummer/${1}_output
}
export -f concat
cat species.3_id | parallel -j 80 concat {}
##Summarize SNPs
for i in $(cat species.3_id); do
python generate_snp_matrix_by_ref.py -i 09.mummer/${i}_output/snps -o 08.snp_tables_by_ref/
done
###Calculate pN/pS
for i in $(cat ../genome.id); do
python annotate_snv_by_position.py --snv 08.snp_tables_by_ref/${i}.snp.matrix.tsv --gff 03.gff/${i}.gff --fasta ../04.MAG/${i}.fa -o 09.pspn/${i}
done
###Add gene IDs and protein sequences
for i in $(cat ../genome.id); do
python ../add_geneinfo_to_pnps.py -i 08.pnps/${i}_gene_level_pnps.tsv -f ../04.protein/${i}.proteins.faa -o 09.final.results/${i}.output.tsv
done
###Data extraction
python parse_pnps.py -i 09.final.results -o 10.analysis
##Functional annotation
emapper.py -i 10.analysis/pN2_or_inf_proteins.faa -o 11.eggnog/
```

## 13.CRISPR-Cas system identification

```sh
#CRISPR-Cas systems were identified using CRISPRCasTyper
concat () {
cctyper 04.MAG/${1}.fa 15.cas/${1}_cctyper_output
}
export -f concat
cat MAG_id | parallel -j 30 concat {}
```

## 14.Secondary Metabolite Potential Analysis

```sh
#Potential secondary metabolite biosynthetic gene clusters (BGCs) were predicted using antiSMASH.
concat () {
antismash 04.MAG/${i}.fa --taxon bacteria --output-dir 16.BGC/01.antismash/${1} --genefinding-tool prodigal --cb-knownclusters -c 1 --cc-mibig  --fullhmmer
}
export -f concat
cat MAG_id | parallel -j 30 concat {}
##Completeness
bigslice -i 02.gbk 03.bigslice.SEMGC.out
##BiG-SLiCE distance analysis
#human
tar -xzvf BGC_All.tar.gz #ABC-HuMi: the Atlas of Biosynthetic Gene Clusters in the Human Microbiome
bigslice -i input_folder_template_human human --num_threads 60
bigslice --query 04.complete.BGCs --n_ranks 1 human -t 40
#ocean
tar -xzvf antismash-bgcs-genomes-unfiltered.tar.gz
bigslice -i input_folder_template_ocean ocean --num_threads 60
bigslice --query 04.complete.BGCs --n_ranks 1 ocean -t 40
#fermentation food
bigslice -i input_folder_template_food fermentation --num_threads 60
bigslice --query 04.complete.BGCs --n_ranks 1 fermentation -t 40

```

## 15.AMP

```sh
#1.smORF prediction
concat () {
macrel get-smorfs -f 04.MAG/${1}.fa -o 17.amp/01.smorf/${1}
}
export -f concat
cat MAG_id | parallel -j 90 concat {}
#2.AMP prediction using Macrel
find 01.smorf -type f -name "macrel.out.smorfs.faa" | xargs cat > all.macrel.out.smorfs.faa
macrel peptides -f all.macrel.out.smorfs.faa -o 02.macrel
#3.AMP-SEMiner-main
python AMP-SEMiner-main/End_to_end_Tok_CLS/pre_Token_Byfa.py --input ../02.macrel/macrel_output.fasta --output AMPSEMiner.out_pred.tsv --model_name model_weights/Tok_CLS/epoch15 --batch_size 4 --max_len 300
#4.ampscannerv2
python amp_scanner_predict.py ../../02.macrel/macrel_output.fasta TrainedModels/021820_FULL_MODEL.h5
#5.apex
python ./apex-main/apex_predict.py -i ../02.macrel/macrel_output.fasta -o apex.out_pred.csv
#6.PepNet
python prott5_embedder.py --input ../../../02.macrel/macrel_output.fasta --output macrel_output.h5 --model Rostlab/prot_t5_xl_uniref50
python predict.py -type AMP -output_path ./ -test_fasta ../../../02.macrel/macrel_output.fasta -feature_file path to macrel_output.h5

```

