# Coat color change transcriptomics
This is the pipeline and scripts used to analyze RNA-Sequencing data from skin samples of mountain hare (Lepus timidus) individuals undergoing autumn molt, by Mafalda Sousa Ferreira (mafalda_sferreira (at) hotmail.com). The dataset in this analysis comprises three skin samples representing three molt stages (early molt or "brown", middle molt or "intermediate" and late molt or "white") taken from four individuals (total of 12 skin samples). We had HiSeq 2500 and 2000 paired-end data for each skin sample.

All python scripts use python 2.7.

This pipeline assumes that raw reads were already cleaned and filtered with tools such as Cutadapt and Trimmomatic. 

- [De novo transcriptome assembly](https://github.com/evochange/coat_color_change_transcriptomics/#de-novo-transcriptome-assembly)
- [Mapping and relative abundance estimation](https://github.com/evochange/coat_color_change_transcriptomics/#mapping-and-relative-abundance-estimation)
- [Gene Ontology analysis with a costume annotation in Ontologizer](https://github.com/evochange/coat_color_change_transcriptomics/blob/master/README.md#gene-ontology-analysis-with-a-costume-annotation-in-ontologizer)

## De novo transcriptome assembly
Reads 1 and reads 2 for all individuals were concatenated in two files ```timidus_R1.fq.gz```and ```timidus_R2.fq.gz```.

#### Generating a raw assembly
We run Trinity 2.6.6 to generate the assembly.

```
/bin/trinityrnaseq-Trinity-v2.6.6/Trinity -seqType fq --max_memory 100G --CPU 8 --left timidus_R1.fq.gz --right timidus_R2.fq.gz --full_cleanup --SS_lib_type RF --output APT_transcriptome_trinity
```

#### Assessing assembly quality
We used Transrate 1.0.3 to assess raw assembly quality and generate a "good" (filtered) assembly.

```
transrate --threads 8 --assembly ~/timidus/raw_transcriptome_APT/APT_transcriptome_trinity.Trinity.fasta --left ~/timidus/filtered_reads/timidus_R1.fq --right ~/timidus/filtered_reads/timidus_R2.fq --output ~/timidus/transrate_good_APT/
```

#### Transcriptome annotation
We annotated our filtered assembly to the ENSEMBL 92 rabbit (OryCun2.0) and mouse (GRCm38) protein references, using the built-in Transrate reciprocal blast algorythm.

```
transrate --assembly ~/timidus/transrate_good_APT/good.APT_transcriptome_trinity.Trinity.fasta --reference ~/timidus/references/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa --threads 4 --output ~/timidus/transrate_good_APT/transrate_blast_ORY
```
```
transrate --assembly ~/timidus/transrate_good_APT/good.APT_transcriptome_trinity.Trinity.fasta --reference ~/timidus/references/Mus_musculus.GRCm38.pep.all.fa --threads 4 --output ~/timidus/transrate_good_APT/transrate_blast_MUS
```

#### Obtaining gene annotations and filtering multiple annotated genes.
Our gene expression analysis was done at the gene level, using TRINITY_DN10004_c0_g1 codes, instead of the isoform information. The following steps allowed us to derive an ENSEMBL **gene** annotation for each Trinity **gene** and also filter our assembly for Trinity **genes** with multiple annotations. This happens because Trinity will assemble contigs that should correspond to isoforms (e.g. TRINITY_DN10004_c0_g1_i1). Since we blasted our contigs to a protein reference, these contigs may annotate to ENSEMBL proteins that correspond to different ENSEMBL genes. If so, we opted to remove these contigs from the gene expression analyis (see [Mapping and relative abundance estimation](https://github.com/evochange/coat_color_change_transcriptomics/#mapping-and-relative-abundance-estimation)). We took the following steps to obtain ENSEMBL gene annotations and detect Trinity genes with multiple annotations.

In the same folder with the blast results run ```do_annot_table.py```. This script will output two files (you can name them something else):
- ```annotation_ORY_MUS.txt``` (output 1) will have three columns with the rabbit and mouse **protein** annotation for each Trinity **isoform** (contig);
- ```annotation_ORY-or-MUS.txt``` (output 2) will have two columns with the rabbit or mouse (if rabbit doesn't exist) **protein** annotation for each Trinity **isoform**.

```
do_annot_table.py good.APT_transcriptome_trinity.Trinity_into_Oryctolagus_cuniculus.OryCun2.0.pep.all.1.blast good.APT_transcriptome_trinity.Trinity_into_Mus_musculus.GRCm38.pep.all.1.blast annotation_ORY_MUS.txt annotation_ORY-or-MUS.txt
```

Now we have a protein annotation for each Trinity "isoform", but we want a **gene** annotation for each Trinity **gene**. Trinity **isoforms** can be converted to **gene** codes very simply by looking at the code of the contig. For instance, TRINITY_DN10004_c0_g1_i1 isoform corresponds to TRINITY_DN10004_c0_g1 gene. However, we need the reference to obtain ENSEMBL gene codes for ENSEMBL protein codes. Retrieve a table with ENSEMBL gene and transcript annotations for each ENSEMBL protein from the original fasta references.

```
grep '^>' Mus_musculus.GRCm38.pep.all.fa | cut -f1,4,5,8 -d' ' > Mus_musculus.GRCm38.pep.protein_information.txt
grep '^>' Oryctolagus_cuniculus.OryCun2.0.pep.all.fa | cut -f1,4,5,8 -d' ' > Oryctolagus_cuniculus.OryCun2.0.pep.protein_information.txt
```

These tables look like this:
```
>ENSOCUP00000010526.2 gene:ENSOCUG00000012231.3 transcript:ENSOCUT00000012228.2
>ENSOCUP00000021714.1 gene:ENSOCUG00000021536.1 transcript:ENSOCUT00000025500.1 gene_symbol:BOK
>ENSOCUP00000016504.1 gene:ENSOCUG00000021116.2 transcript:ENSOCUT00000027952.1
>ENSOCUP00000009605.2 gene:ENSOCUG00000011164.2 transcript:ENSOCUT00000011165.2 gene_symbol:THAP4
>ENSOCUP00000001915.2 gene:ENSOCUG00000002218.3 transcript:ENSOCUT00000002222.2
>ENSOCUP00000027137.1 gene:ENSOCUG00000002218.3 transcript:ENSOCUT00000034088.1
```

Then, run the following scripts in the same folder where you have ```annotation_ORY_MUS.txt```,```Mus_musculus.GRCm38.pep.protein_information.txt```and ```Oryctolagus_cuniculus.OryCun2.0.pep.protein_information.txt```. 

```give_multiannotated_genes_step1.py``` (you will have to change the script to read other inputs) will use these three files to create an intermediate file ```new_annotation_DONOTUSETHISFILE.txt```with the ENSEMBL **gene** rabbit and mouse annotation of each Trinity **gene**. The script convertes the Trinity isoform code to gene code by removing the isoform information from the contig name, and thus will contain duplicate Trinity gene codes that may have different ENSEMBL gene codes. For example:

```
TRINITY_DN10025_c0_g1	ENSOCUG00000029166.1	NA
TRINITY_DN10025_c0_g1	ENSOCUG00000029166.1	NA
TRINITY_DN10004_c0_g1 ENSOCUG00000009225.3	ENSMUSG00000020542.18
TRINITY_DN10004_c0_g1 ENSOCUG00000003826.3	ENSMUSG00000031027.15
```

The next script ```give_multiannotated_genes_step2.py``` will create a consensus file with uniq Trinity gene codes per line and a consensus annotation, which will be written to one file. If it finds that the same Trinity gene is annotated to more than one ENSEMBL gene code, it will write it to a separate file. The output of ```give_multiannotated_genes_step2.py``` will be:

- single_annotated_genes.txt: single annotated genes (needed to filter RSEM's output).
- multi_annotated_genes.txt: multipleannotated genes, either to rabbit or mouse.
- NA_annotated_genes.txt: the genes where one isoform is not annotated.

You will need to change the scripts to read your input files. Call both scripts like so:

```
give_multiannotated_genes_step1.py
give_multiannotated_genes_step2.py
```


## Mapping and relative abundance estimation

Create bowtie2 indices first, using the RSEM wrapper.

```
~/my_programs/RSEM-1.3.0/extract-transcript-to-gene-map-from-trinity good.APT_transcriptome_trinity.Trinity.fasta APT_transcript-to-gene-map
```
```
~/my_programs/RSEM-1.3.0/rsem-prepare-reference --bowtie2 --bowtie2-path ~/my_programs/bowtie2-2.3.4.1-linux-x86_64/ --transcript-to-gene-map APT_transcript-to-gene-map good.APT_transcriptome_trinity.Trinity.fasta APT
```

Map reads using bowtie2.

```
for f in $(ls *1.fastq.gz); do ~/my_programs/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -q --phred33 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --dpad 0 --gbar 9999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant --nofw -p 1 -k 200 --fr -x ../transrate_good_APT/APT -1 $f -2 ${f/1.fastq.gz/2.fastq.gz} -S ${f/R1.fastq.gz/APT.sam}; done
```

Use RSEM to calculate expression.
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

## Gene Ontology analysis with a costume annotation in Ontologizer
Here, we have already performed the gene expression analysis between molting stages using edgeR. We now have a set of differentially expressed genes, whose expression is high at early or late molt. We would like to understand if there are specific and different biological functions that these "early" or "late" expressed genes are performing. To do so, we performed a Gene Ontology enrichment analysis for each set of genes. To take advantage of our mixed annotation where we have rabbit and mouse ENSEMBL codes we used [Ontologizer](http://ontologizer.de/) because it allows using a costume Gene Association File (GAF) file where we have ENSEMBL rabbit and mouse gene annotations. There are perhaps better ways to perform this analysis now, and we would advise interested users to check the [Ontologizer website](http://ontologizer.de/), or [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) that now supports non-model organisms as well. We used the following pipeline to generate a costume GAF file.

#### Prepare the necessary inputs
From ENSEMBL [Biomart](https://www.ensembl.org/biomart/martview/ab1eb25c5409ef27a44f88d46e7448f9), retrieve the corresponding GO term accession for each ENSEMBL gene in the rabbit and mouse genome. The output should look like this:
```
Gene stable ID  GO term accession
ENSOCUG00000012231      
ENSOCUG00000021536      GO:0042981
ENSOCUG00000021536      GO:0005102
ENSOCUG00000021536      GO:0031625
ENSOCUG00000021536      GO:0042803
ENSOCUG00000021536      GO:0001836
ENSOCUG00000021536      GO:0006915
...
```

Store them in two separate files, like ```ENSEMBL92_Mmusculus_GRCm38_GO_annotation.txt``` and ```ENSEMBL92_OryCun2_GO_annotation.txt```.

Also, store a list with the ENSEMBL rabbit and mouse annotations in the transcriptome in two separate files (in our case ```ltimidus_APT_annot_ory.txt``` and ```ltimidus_APT_annot_mus.txt```). The files should look like this:

```
ENSOCUG00000004892
ENSOCUG00000015449
ENSOCUG00000026197
ENSOCUG00000024333
...
```

Finally, retrieve a Gene Ontology File (obo) from the Gene Ontology [website](http://geneontology.org/). In our case we named it ```go-basic_6Set2018.obo```.

#### Generate the GAF file
Use the scripts ```make_input_4_gaf.py```, ```parse_obo_file.py```and ```do_Gaf_file_new_version.py```. ```make_input_4_gaf.py``` and ```do_Gaf_file_new_version.py``` were published previously for a similar work on snowshoe hares that you can find [here](https://github.com/MafaldaSFerreira/Snowshoe-hare-transcriptome).

Use them like so:
```
python make_input_4_Gaf.py ENSEMBL92_OryCun2_GO_annotation.txt ENSEMBL92_Mmusculus_GRCm38_GO_annotation.txt ltimidus_APT_annot_ory.txt ltimidus_APT_annot_mus.txt ltimidus_GO_annot_ENSEMBL92.txt
```
```
parse_obo_file.py go-basic_6Set2018.obo go-basic_6Set2018.tab
```
```
python do_GAF_file_new_version.py ltimidus_GO_annot_ENSEMBL92.txt go-basic_6Set2018.tab > ltimidus_transcriptome_6Sep18.gaf
```

```ltimidus_transcriptome_6Sep18.gaf``` can now be used as input in Ontologizer, along with ```go-basic_6Set2018.obo``` and your gene sets made of the ENSEMBL annotations of your differentially expressed genes. 


