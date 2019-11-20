#A script to produce a correspondence file between the ENSEMBL codes of my transcriptome and their respective GO terms
#Because I have both MUS and ORYC ENSEMBL annotations, this script will work with both
# AUTHOR: Mafalda S. Ferreira, 2014 @ Porto
# ---------------------------------------------------------------------------------------------------------------------------
# USAGE: python make_input_4_gaf.py corr_btw_ory_ensembl&go corr_btw_mus_ensembl&go annot_ory2ensembl annot_mus2ensembl output
# ---------------------------------------------------------------------------------------------------------------------------
# FILES: 
# corr_btw_ory_ensembl&go = file obtained from biomart that has the correspondence between ALL ORY ensembl gene codes and GO terms
# corr_btw_mus_ensembl&go = file obtained from biomart that has the correspondence between ALL MUS ensembl gene codes and GO terms
# annot_ory2ensembl = a list with the ORY ENSEMBL codes annotated in my transcriptome  ATTENTION: Only one column with ENSEMBL gene codes
# annot_mus2ensembl = a list with the MUS ENSEMBL codes annotated in my transcriptome. ATTENTION: Only one column with ENSEMBL gene codes
# ----------------------------------------------------------------------------------------------------------------------------
import sys

ref_ory=sys.argv[1]
ref_mus=sys.argv[2]
infile_ory=sys.argv[3]
infile_mus=sys.argv[4]
output=sys.argv[5]

print 'Hello! Check if your files are well formated before usage.'

def makeorydict(ory):
	oryref=open(ory,'r')
	
	orydict={}
	
	for line in oryref:
		# print line
		line=line.strip('\n')
		fields=line.split('\t')
		#print fields
		key=fields[0]
		value=fields[1]+','
		#print key,value
		if key not in orydict:
			orydict.setdefault(key,value)
		elif key in orydict:
			orydict[key]=orydict[key]+value
	
	return orydict
	
def makemusdict(mus):
	
	musref=open(mus,'r')
	musdict={}

	for line in musref:
		# print line
		line=line.strip('\n')
		fields=line.split('\t')
		#print fields
		key=fields[0]
		value=fields[1]+','
		#print key,value
		if key not in musdict:
			musdict.setdefault(key,value)
		elif key in musdict:
			musdict[key]=musdict[key]+value
	
	return musdict
	
dictory=makeorydict(ref_ory)
dictmus=makemusdict(ref_mus)

annotmus=open(infile_mus,'r')
annotory=open(infile_ory,'r')

outfile=open(output,'w')

for line in annotory:
	myfields=line.strip()
	#print myfields
	if myfields in dictory:
		new_line_ory=''.join(myfields.strip()+' '+dictory[myfields])
		outfile.write(new_line_ory+'\n')

for line in annotmus:
	myfields=line.strip()
	#print myfields
	if myfields in dictmus:
		new_line_mus=''.join(myfields.strip()+' '+dictmus[myfields])
		#print new_line_mus
		outfile.write(new_line_mus+'\n')
	
outfile.close()
