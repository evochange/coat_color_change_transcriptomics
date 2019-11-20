#!/usr/bin/python

##------------------------------------------------------------------------------------------------------------
# This script creates two dictionaries from an input created with give_multiannotated_gene_step1.py. For each
# Trinity gene, it will create a dictionary with ALL the EnSEMBL ORY or ENSEMBL MUS (or HOMO). Then, 
# it will find multiple annotated genes, by looping through the keys of each entry in these dictionaries
#
# Assumption 1: Genes that have one isoform not annotated will be removed, even if the remaining isoforms have 
# concordant annotations to the same gene. I am conservative and do not assume that the annotation of 
# non-annotated isoforms is equal to the remainder of the isoforms.
# Assumption 2: We always prefer the ORY annotation, so the script reads that annotation first, checks
# for NAs, and if there are NAs, it will look at the HUMAN annotation for a consistent annotation across 
# isoforms. If there are no NAs in the ory annotation, the script does not look at the Human annotation at 
# all. It looks for multipleannotation in the ORY genes.
# Assumption 3: In the end, we only accept genes with NA in ory but not NA in Homo, and that are not 
# multipleannotated to neither.
#------------------------------------------------------------------------------------------------------------
# input1: new_annotation_DONOTUSETHISFILE.txt that has to be in the folder where the script is run
# output1: single annotated genes (needed to filter RSEM's output): single_annotated_genes.txt
# output2: multipleannotated genes, either to Ory or Mus: multi_annotated_genes.txt
# output3: the genes where one isoform is not annotated: NA_annotated_genes.txt
#------------------------------------------------------------------------------------------------------------
# USAGE: give_multiannotated_genes_step2.py
#------------------------------------------------------------------------------------------------------------
# BY: MafaldaSF 15/05/2018: Modified to include this header.
#------------------------------------------------------------------------------------------------------------


def checkequal(iterator):
	try:
		iterator=iter(iterator)
		first=next(iterator)
		return all(first == rest for rest in iterator)
	except StopIteration:
		return len(set(iterator))<=1
	
#let's read the annotation again... (I know I'm lazy...)
readagainannotation=open("new_annotation_DONOTUSETHISFILE.txt","r")

#create a dictionary where to write the annotation
orydict={}
homodict={}


for line in readagainannotation:
	output2=''
	line=line.strip().split('\t')
	ncontig=line[0]
	orygene=line[1]
	homogene=line[2]
	
	#is the contig already in the ory dictionary?
	if ncontig in orydict:
		#if it's already there, let's add the new annotation to a list:
		#get the previous values
		oryvalue=orydict.get(ncontig)
		#append the new value to the old values
		oryvalue.append(orygene)
		#add the new key+values to the dictionary:
		orydict.setdefault(ncontig,oryvalue)
	else:
		#else create an empty list
		oryvalue=[]
		#add the annotation to the list
		oryvalue.append(orygene)
		#add both to the dictionary:
		orydict.setdefault(ncontig,oryvalue)
	
	if ncontig in homodict:
		homovalue=homodict.get(ncontig)
		homovalue.append(homogene)
		homodict.setdefault(ncontig,homovalue)
	else:
		homovalue=[]
		homovalue.append(homogene)
		homodict.setdefault(ncontig,homovalue)
		
	orykeys=orydict.keys()
	homokeys=homodict.keys()

#open the outputs:
multioutput=open("multi_annotated_genes.txt","w")
singleoutput=open("single_annotated_genes.txt","w")
NAoutput=open("NA_annotated_genes.txt","w") #this is only to check if the python is doing the right thing

multioutput.write('CONTIG'+'\t'+'Annotation'+'\n')
singleoutput.write("CONTIG"+'\t'+'Annotation'+'\n')
NAoutput.write("CONTIG"+'\n')

#now let's iterate the keys:
for key in orykeys:
	values=orydict.get(key)
	#1) is any of the isoforms not annotated?
	if any(a=='NA' for a in values):
		#if some of the isoform or all are NA, let's check the homo annotation:
		#clear the variable from the ory keys
		values=[] 
		values=homodict.get(key)
		#1a) is any or all of the values NA?
		if any(b=='NA' for b in values):
			NAoutput.write(key+'\n')
		#1b) are the values all equal between themselves?
		else:
			destiny=str(checkequal(values))
			if destiny=='True':
				#you only need one of the values, so index the first (values[0])
				singleoutput.write(key+'\t'+values[0]+'\n')
			else:
				for c in values: #write one line for each c
					multioutput.write(key+'\t'+c+'\n')
	#2) if all are annotated	
	else:
		#2a) do they all have the same annotation?
		destiny=str(checkequal(values))
		if destiny=='True':
			singleoutput.write(key+'\t'+values[0]+'\n')
		else:
			for d in values: #write one line for each d
				multioutput.write(key+'\t'+d+'\n')
		
					
multioutput.close()
singleoutput.close()
NAoutput.close()