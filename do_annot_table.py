#!/usr/bin/python

##------------------------------------------------------------------------------------------------------------
# This scripts takes the blast outputs from Transrate with annotations of every Trinity gene to 
# the Rabbit or Mouse reference.
#------------------------------------------------------------------------------------------------------------
# input1: blast file em formato 6 com anotacao a ORY
# input2: blast file em formato 6 com anotacao a MUS
# main output1 (three columns): contig; rabbit annotation; mouse annotation
# main output2 (two columns): contig; rabbit annotation or mouse annotation when rabbit does not exist
#------------------------------------------------------------------------------------------------------------
# USAGE: do_annotation_table.py input1 input2 output1 output2
#------------------------------------------------------------------------------------------------------------
# BY: MafaldaSF 19/04/2016 16:47am
# BY: MafaldasSF 14/05/2018; changed the script to work with the mouse annotation. "HOMO" or
# "homo" refer to a previous script where I used human annotation
#------------------------------------------------------------------------------------------------------------


import os
import sys

input1=sys.argv[1]
input2=sys.argv[2]
output1=sys.argv[3]
output2=sys.argv[4]

#Definir como fazer os dicionarios
def makeorydict(oryblast):

#abrir o dicionario:
	orydict={}

	#retirar os best blast hits das tabelas
	os.system("sort -k1,1 -k12,12gr -k11,11g -k3,3gr %s | sort -u -k1,1 --merge > bestHits_ory" % (oryblast))
	
	#ler os outputs do comando:
	BH_ory=open('bestHits_ory','r')

	for line in BH_ory: #para cada linha em BH_ory
		fields=line.strip().split('\t') #separar as linhas
		contig=''.join(fields[0]) #contig e o primeiro campo
		protein=''.join(fields[1]) #protein e o segundo campo
		
		orydict.setdefault(contig,protein) #escrever output para o dicionario
	
	return orydict
	

def makehomodict(homoblast):
	homodict={}
	
	os.system("sort -k1,1 -k12,12gr -k11,11g -k3,3gr %s | sort -u -k1,1 --merge > bestHits_homo" % (homoblast))

	BH_homo=open('bestHits_homo','r')

	for line in BH_homo:
		fields=line.strip().split('\t')
		contig=''.join(fields[0])
		protein=''.join(fields[1])
	
		homodict.setdefault(contig,protein)
		
	return homodict
	
#chamar as funcoes definidas em cima para criar os dicionarios:
orydict=makeorydict(input1)
homodict=makehomodict(input2)

#criar um ficheiro annotated_contigs.txt com todos os contigs anotados nos tres blasts:
os.system("cut -f1 bestHits_ory > contigs_2_ORY.txt")
os.system("sort contigs_2_ORY.txt > contigs_2_ORY_2.txt")
os.system("cut -f1 bestHits_homo > contigs_2_HOMO.txt")
os.system("sort contigs_2_HOMO.txt > contigs_2_HOMO_2.txt")
os.system("cat contigs_2_ORY.txt contigs_2_HOMO.txt > contigs_2_all.txt")
os.system("sort contigs_2_all.txt | uniq > annotated_contigs.txt")

#verificar que os ficheiros foram bem criados:
#se foram bem criados, o numero de contigs comuns entre annotated_contigs.txt e os ficheiros
#individuais com os contigs anotados em cada blast deve ser igual!
print "Are these the same numbers?"
print "How many contigs does contigs_2_ORY.txt have?:"
os.system("wc -l contigs_2_ORY.txt")
print "How many commons between contigs_2_ORY.txt and annotated_contigs.txt?"
os.system("comm -12 contigs_2_ORY_2.txt annotated_contigs.txt | wc -l")

print "How many contigs does contigs_2_HOMO.txt have?:"
os.system("wc -l contigs_2_HOMO.txt")
print "How many commons between contigs_2_HOMO.txt and annotated_contigs.txt?"
os.system("comm -12 contigs_2_HOMO_2.txt annotated_contigs.txt | wc -l")

#chamar annotated_contigs.txt
annotated=open("annotated_contigs.txt","r")

#abrir os outputs e escrever a primeira linha
table=open(output1,"w")
table2=open(output2,"w")
table.write("CONTIG"+"\t"+"ORY"+"\t"+"MUS"+"\n")
table2.write("CONTIG"+"\t"+"ANNOTATION"+"\n")


for contig in annotated:
	output=''
	output_2=''
	contig=str(contig.strip())
	
	#se o contig existe no dicionario ORY a=anotacao ory
	if contig in orydict:
		a=''.join(orydict.get(contig))
	#se nao existe, escrever um NA no output na coluna do ORY
	else:
		a=''.join('NA')
	#se o contig existe no dicionario HOMO a=anotacao homo
	if contig in homodict:
		c=''.join(homodict.get(contig))
	#se nao existe, escrever um NA no output na coluna do HOMO
	else:
		c=''.join('NA')
	
	#criar a linha de output:
	output=''.join(contig+'\t'+a+'\t'+c)
	table.write(output+'\n')
	
	#se nao existe uma anotacao ORY
	if a=='NA':
		#se nao existe uma anotacao HOMO
		if c=='NA':
			#nao deve acontecer porque so estamos a considerar contigs com pelo menos uma anotacao
			print "this can't happen!!"
			print c
			#fiquemos entao com a anotacao HOMO
		else:
			output_2=''.join(contig+'\t'+c)	
	#ficamos com a anotacao a ORY:
	else:
		output_2=''.join(contig+'\t'+a)
	#escrever a linha de output para a segunda tabela:
	table2.write(output_2+'\n')

#TCHANAM! 	
table.close()
table2.close()
		
		
	
	
	
