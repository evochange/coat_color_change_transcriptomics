#!/usr/bin/python

#------------------------------------------------------------------------------------------------------------
# This script parses a obo file and makes a table that has the GO category for each 
# GO:accession number. This is the input for do_Gaf_file_new_version.py 
#------------------------------------------------------------------------------------------------------------
# USAGE: parse_obo_file.py go-basic_6Set2018.obo go-basic_6Set2018.tab
#------------------------------------------------------------------------------------------------------------
# BY: MafaldaSFerreira 15/05/2018: Modified to include this header.
#------------------------------------------------------------------------------------------------------------

#make a dictionary with the correspondances between the GO IDs and GO categories
import re
import sys

input_file=sys.argv[1]
output_file=sys.argv[2]

def extractGOcategory(obo):
	obo=open(obo,'r')
	
	dict2={}
	
	lines=obo.readlines()
	
	number=len(lines)
	
	regexp=re.compile("obsolete")
	
	for i in range(0,number):
		if lines[i].startswith('[Term]'):
			
			GO=str(lines[i+1].strip('id:').strip())
			name=str(lines[i+2].strip('name: ').strip())
			category=str(lines[i+3].strip('namespace:').strip())
			
			if regexp.search(name) is not None: #if the term is obs
				dict2.setdefault(GO,"obs")
			elif category=="molecular_function":
				dict2.setdefault(GO,"M")
			elif category=="biological_process":
				dict2.setdefault(GO,"B")
			elif category =="cellular_component":
				dict2.setdefault(GO,"C")
		
	return dict2

mydict=extractGOcategory(input_file)
output=open(output_file,"w")

for GO in mydict:
	value=mydict[GO]
	output.write(''.join(GO+'\t'+value+'\n'))
	
output.close()