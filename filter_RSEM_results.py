#!/usr/bin/python

#------------------------------------------------------------------------------------------------------------
# use the file sag_to_remove (or any file with a list of genes!) to
# filter the RSEM gene.results
# the output will be name_rsem_file.filtered
#------------------------------------------------------------------------------------------------------------
# USAGE: filter_RSEM_results.py rsem_file list_of_genes_to_keep
#------------------------------------------------------------------------------------------------------------
# BY: MafaldaSFerreira 20/04/2016
#------------------------------------------------------------------------------------------------------------

import sys

input1=sys.argv[1]
input2=sys.argv[2]

def makedict(input1):
	
	infile=open(input1,'r')
	dict={}
	
	for line in infile:
		line=str(line.strip())
		dict.setdefault(line,line)
	
	return dict

#make the filter with makedict()
filter=makedict(input2)

#open the rsem file
rsem_input=open(input1,'r')

#create the output
output_name=str(input1+'.filtered')
output=open(output_name,'w')

#read lines from the rsem file
lines=rsem_input.readlines()
#write the header to the output
output.write(lines[0])
#retrieve number of lines
number=len(lines)

#skip the header line
for line in lines[1:number]:
	newline=line.strip().split('\t')
	gene=str(newline[1])

	if gene in filter:
		output.write(line)
		
output.close()
	
	
		
	
	
