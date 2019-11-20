# Coat color change transcriptomics
Pipeline and scripts necessary for the analysis of RNA-Sequencing data from the mountain hare (Lepus timidus) fall moult

This pipeline assumes that raw reads were already cleaned and filtered with tools such as Cutadapt and Trimmomatic.

All python scripts use python 2.7

- [De novo transcriptome assembly](https://github.com/evochange/coat_color_change_transcriptomics/#de-novo-transcriptome-assembly)
- [Mapping and relative abundance estimation](https://github.com/evochange/coat_color_change_transcriptomics/#mapping-and-relative-abundance-estimation)

## De novo transcriptome assembly
Reads 1 and reads 2 for all individuals were concatenated in two files ```timidus_R1.fq.gz```and ```timidus_R2.fq.gz```

#### Generating a raw assembly
We run Trinity 2.6.6 to generate the assembly

```
/bin/trinityrnaseq-Trinity-v2.6.6/Trinity -seqType fq --max_memory 100G --CPU 8 --left timidus_R1.fq.gz --right timidus_R2.fq.gz --full_cleanup --SS_lib_type RF --output APT_transcriptome_trinity
```

#### Assessing assembly quality
We used Transrate 1.0.3 to assess raw assembly quality and generate a filtered assembly

```
transrate --threads 8 --assembly ~/timidus/raw_transcriptome_APT/APT_transcriptome_trinity.Trinity.fasta --left ~/timidus/filtered_reads/timidus_R1.fq --right ~/timidus/filtered_reads/timidus_R2.fq --output ~/timidus/transrate_good_APT/
```

#### Transcriptome annotation
We annotated our filtered assembly to Rabbit (OryCun2.0) and Mouse (GRCm38) obtained from ENSEMBLE 92 protein references, using the built-in Transrate reciprocal blast algorythm.

```
transrate --assembly ~/timidus/transrate_good_APT/good.APT_transcriptome_trinity.Trinity.fasta --reference ~/timidus/references/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa --threads 4 --output ~/timidus/transrate_good_APT/transrate_blast_ORY
```
```
transrate --assembly ~/timidus/transrate_good_APT/good.APT_transcriptome_trinity.Trinity.fasta --reference ~/timidus/references/Mus_musculus.GRCm38.pep.all.fa --threads 4 --output ~/timidus/transrate_good_APT/transrate_blast_MUS
```

#### Filtering the transcriptome by annotation
From the "good" transcriptome generated with Transrate, we also removed Trinity genes that annotate to multiple ENSEMBL genes. This happens because Trinity will assemble contigs that should correspond to isoforms. In some cases, these contigs will annotate to ENSEMBL proteins that originate from correspond to different genes and, in that case, we removed them from the analysis. The following are steps we took to find Trinity genes with multiple annotations.

In the same foulder with the blast results run ```do_annot_table.py```

```
do_annot_table.py good.APT_transcriptome_trinity.Trinity_into_Oryctolagus_cuniculus.OryCun2.0.pep.all.1.blast good.APT_transcriptome_trinity.Trinity_into_Mus_musculus.GRCm38.pep.all.1.blast annotation_ORY_MUS.txt annotation_ORY-or-MUS.txt
```

Retrieve a table with ENSEMBL gene and transcript annotations for each ENSEMBL protein from the original fasta references

```
grep '^>' Mus_musculus.GRCm38.pep.all.fa | cut -f1,4,5,8 -d' ' > Mus_musculus.GRCm38.pep.protein_information.txt
grep '^>' Oryctolagus_cuniculus.OryCun2.0.pep.all.fa | cut -f1,4,5,8 -d' ' > Oryctolagus_cuniculus.OryCun2.0.pep.protein_information.txt
```
Then run the scripts in the same foulder where the last three steps were run. ```/give_multiannotated_genes_step1.py``` will create an intermediate file with the ENSEMBL gene rabbit and mus annotation of each Trinity gene. This file will include multiple annotated Trinity genes.  ```/give_multiannotated_genes_step2.py``` will separate this file into single and multiple annotated genes.
```
./give_multiannotated_genes_step1.py
./give_multiannotated_genes_step2.py
```
The output of ```give_multiannotated_genes_step2.py``` are three files:

- single annotated genes (needed to filter RSEM's output): single_annotated_genes.txt
- multipleannotated genes, either to Ory or Mus: multi_annotated_genes.txt
- the genes where one isoform is not annotated: NA_annotated_genes.txt 

## Mapping and relative abundance estimation

Create bowtie2 indices first, using the RSEM wrapper:

```
~/my_programs/RSEM-1.3.0/extract-transcript-to-gene-map-from-trinity good.APT_transcriptome_trinity.Trinity.fasta APT_transcript-to-gene-map
```
```
~/my_programs/RSEM-1.3.0/rsem-prepare-reference --bowtie2 --bowtie2-path ~/my_programs/bowtie2-2.3.4.1-linux-x86_64/ --transcript-to-gene-map APT_transcript-to-gene-map good.APT_transcriptome_trinity.Trinity.fasta APT
```

Map reads using bowtie2:

```
for f in $(ls *1.fastq.gz); do ~/my_programs/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -q --phred33 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --dpad 0 --gbar 9999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant --nofw -p 1 -k 200 --fr -x ../transrate_good_APT/APT -1 $f -2 ${f/1.fastq.gz/2.fastq.gz} -S ${f/R1.fastq.gz/APT.sam}; done
```

Use RSEM to calculate expression:
```
for f in $(ls *.sam); do ~/my_programs/RSEM-1.3.0/rsem-calculate-expression --calc-ci --seed-length 23 --paired-end --sam $f ~/timidus/transrate_good_APT/APT_transcriptome_trinity.Trinity/APT ~/timidus/rsem_APT/${f/.sam/_rsem_results}; done
```

Finally, filter multiple annotated genes from RSEM's output with ```filter_RSEM_results.sh``` and ```single_annotated_genes.txt```generated previously.

```
cut -f1 single_annotated_genes.txt | grep -v 'CONTIG' > sag_to_remove
```
The script ```filter_RSEM_results.py``` will accept any list of Trinity codes. Here, I use a list with single annotated genes.

```
for i in $(ls *.genes.results); do filter_RSEM_results.py $i sag_to_remove; done
```

