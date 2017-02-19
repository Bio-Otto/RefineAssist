#!/usr/bin/python

# Author: Sari Sabban
# Email:  sari.sabban@gmail.com
# URL:    https://github.com/sarisabban
#
# Created By:   	Sari Sabban
# Created Date: 	13 January 2017
#
# Modified By:   	Sari Sabban
# Modified Date: 	3 February 2017

import re
import itertools
from Bio.PDB import *

print 'File must be in the same directory as this script.'
filename=raw_input('Enter File Name>')

#This script calculates and prints out the SASA of each amino acid using DSSP from within biopython. 
y=None
w=None
sasalist=list()
p=PDBParser()
structure=p.get_structure('X',filename)
model=structure[0]
dssp=DSSP(model,filename,acc_array='Wilke')
for x in dssp:
	if x[1]=='A':
		y=129*(x[3])
		w='ALA'
	elif x[1]=='V':
		y=174*(x[3])
		w='VAL'
	elif x[1]=='I':
		y=197*(x[3])
		w='ILE'
	elif x[1]=='L':
		y=201*(x[3])
		w='LEU'
	elif x[1]=='M':
		y=224*(x[3])
		w='MET'
	elif x[1]=='P':
		y=159*(x[3])
		w='PRO'
	elif x[1]=='Y':
		y=263*(x[3])
		w='TYR'
	elif x[1]=='F':
		y=240*(x[3])
		w='PHE'
	elif x[1]=='W':
		y=285*(x[3])
		w='TRP'
	elif x[1]=='R':
		y=274*(x[3])
		w='ARG'
	elif x[1]=='N':
		y=195*(x[3])
		w='ASN'
	elif x[1]=='C':
		y=167*(x[3])
		w='CYC'
	elif x[1]=='Q':
		y=225*(x[3])
		w='GLN'
	elif x[1]=='E':
		y=223*(x[3])
		w='GLU'
	elif x[1]=='G':
		y=104*(x[3])
		w='GLY'
	elif x[1]=='H':
		y=224*(x[3])
		w='HIS'
	elif x[1]=='K':
		y=236*(x[3])
		w='LYS'
	elif x[1]=='S':
		y=155*(x[3])
		w='SER'
	elif x[1]=='T':
		y=172*(x[3])
		w='THR'
	elif x[1]=='D':
		y=193*(x[3])
		w='ASP'
	z=str(x[0])+' '+str(w)+' '+str(y)
	sasalist.append(z)
sasainfo='\n'.join(sasalist)

#This script prints out the secondary structure each amino acid belongs to using DSSP from within biopython. 
sslist=list()
p = PDBParser()
structure = p.get_structure('X', filename)
model = structure[0]
dssp = DSSP(model, filename)
for x in dssp:
	if x[2]=='G' or x[2]=='H' or x[2]=='I':
		y='H'
	elif x[2]=='B' or x[2]=='E':
		y='S'
	else:
		y='L'
	sslist.append(y)
aa_ss=''.join(sslist)

sasainf=sasainfo.splitlines()
somelist=list()
for x,z in zip(sasainf,aa_ss):
	y=x.strip().split()
	somelist.append(y[0]+' '+z+' '+y[1]+' '+y[2])	
data='\n'.join(somelist)

#Make a list of tuples, each containing an amino acid's secondary structure as key, and its corresponding SASA result as value. The amino acid sequence order is maintained.
lis=list()
value=None
key=None
for line in data.splitlines():
	line=line.strip()
	find_value=re.findall('\d+\.\d+$',line)
	for x in find_value:
		value=float(x)
	find_key=re.findall('\s[HSL]\s',line)
	for y in find_key:
		key=y.strip()
	lis.append((key,value))

#Label each amino acid depending on its SASA position according to the parameters highlighted in the paper by (Koga et.al., 2012 - PMID: 23135467). The parameters are as follows:
#Surface:
#	Helix or Sheet:	SASA=>60
#	Loop:		SASA=>40
#
#Boundry:
#	Helix or Sheet:	15<SASA<60
#	Loop:		25<SASA<40
#
#Core:
#	Helix or Sheet:	SASA=<15
#	Loop:		SASA=<25
core=list()
boundery=list()
surface=list()
count=0
SASA=list()
for x,y in lis:
	count=count+1

	if y <=25 and x=='L':
		core.append(count)
		SASA.append('C')
	elif 25<y<40 and x=='L':
		boundery.append(count)
		SASA.append('B')
	elif y>40 and x=='L':
		surface.append(count)
		SASA.append('S')

	elif y <=15 and x=='H':
		core.append(count)
		SASA.append('C')
	elif 15<y<60 and x=='H':
		boundery.append(count)
		SASA.append('B')
	elif y>60 and x=='H':
		surface.append(count)
		SASA.append('S')

	elif y <=15 and x=='S':
		core.append(count)
		SASA.append('C')
	elif 15<y<60 and x=='S':
		boundery.append(count)
		SASA.append('B')
	elif y>60 and x=='S':
		surface.append(count)
		SASA.append('S')

#Convert the amino acid short names in the variable sasa into an amino acid sequence.
aminos1=list()
for line in data.splitlines():
	line=line.strip()
	AA=re.findall('\S[A-Z]\S',line)
	aminos1.append(AA)

aminos2=list()
for x in aminos1:
	for y in x:
		if y=='ALA':
			aminos2.append('A')
		if y=='VAL':
			aminos2.append('V')
		if y=='ILE':
			aminos2.append('I')
		if y=='LEU':
			aminos2.append('L')
		if y=='MET':
			aminos2.append('M')
		if y=='PRO':
			aminos2.append('P')
		if y=='TYR':
			aminos2.append('Y')
		if y=='PHE':
			aminos2.append('F')
		if y=='TRP':
			aminos2.append('W')
		if y=='GLY':
			aminos2.append('G')
		if y=='CYS':
			aminos2.append('C')
		if y=='GLN':
			aminos2.append('Q')
		if y=='ASN':
			aminos2.append('N')
		if y=='THR':
			aminos2.append('T')
		if y=='SER':
			aminos2.append('S')
		if y=='ARG':
			aminos2.append('R')
		if y=='HIS':
			aminos2.append('H')
		if y=='LYS':
			aminos2.append('K')
		if y=='ASP':
			aminos2.append('D')
		if y=='GLU':
			aminos2.append('E')

#Print the custom PYMOL commands to select the different protein layers according to each amino acid's SASA.
print 'Surface Amino Acids Command:'
print 'select Surf, resi','+'.join(str(z) for z in surface),'\n'

print 'Boundery Amino Acids Command:'
print 'select Bound, resi','+'.join(str(z) for z in boundery),'\n'

print 'Core Amino Acids Command:'
print 'select Core, resi','+'.join(str(z) for z in core),'\n'

#Calculate which amino acids are in the wrong layer. The rules this script are from the following Rosetta LayerDesign Protocol (source: goo.gl/NsQubf) and they are as follows:
#Surface
#    Loop:		PGNQSTDERKH
#    Helix:		QEKH
#    Strand:		QTY
#
#Boundary
#    Loop:		AVILFYWGNQSTPDEHR
#    Helix:		AVILWQEKFM
#    Strand:		AVILFYWQTM
#
#Core:
#    Loop:		AVILPFWM
#    Helix:		AVILFW
#    Strand:		AVILFWM
d=list()
def pairwise(iterable):
	a, b = itertools.tee(iterable)
	next(b, None)
	return zip(a, b)

mutate=list()
for a,(b1,b2),c in zip(SASA,(pairwise(lis)),aminos2):
	if a=='S' and b2[0]=='L' and (c=='P' or c=='G' or c=='N' or c=='Q' or c=='S' or c=='T' or c=='D' or c=='E' or c=='R' or c=='K' or c=='H'):
		mutate.append(' ')
	elif a=='B' and b2[0]=='L' and (c=='A' or c=='V' or c=='I' or c=='L' or c=='F' or c=='Y' or c=='W' or c=='G' or c=='N' or c=='Q' or c=='S' or c=='T' or c=='P' or c=='D' or c=='E' or c=='H' or c=='R'):
		mutate.append(' ')
	elif a=='C' and b2[0]=='L' and (c=='A' or c=='V' or c=='I' or c=='L' or c=='P' or c=='F' or c=='W' or c=='M'):
		mutate.append(' ')
	elif a=='S' and b2[0]=='H' and (c=='Q' or c=='E' or c=='K' or c=='H'):
		mutate.append(' ')
	elif a=='B' and b2[0]=='H' and (c=='A' or c=='V' or c=='I' or c=='L' or c=='W' or c=='Q' or c=='E' or c=='K' or c=='F' or c=='M'):
		mutate.append(' ')
	elif a=='C' and b2[0]=='H' and (c=='A' or c=='V' or c=='I' or c=='L' or c=='F' or c=='W'):
		mutate.append(' ')
	elif a=='S' and b2[0]=='S' and (c=='Q' or c=='T' or c=='Y'):
		mutate.append(' ')
	elif a=='B' and b2[0]=='S' and (c=='A' or c=='V' or c=='I' or c=='L' or c=='F' or c=='Y' or c=='W' or c=='Q' or c=='T' or c=='M'):
		mutate.append(' ')
	elif a=='C' and b2[0]=='S' and (c=='A' or c=='V' or c=='I' or c=='L' or c=='F' or c=='W' or c=='M'):
		mutate.append(' ')
	else:
		mutate.append('*')

#Print alignment of SASA sequence, secondary structure sequence, amino acids sequence, and indicate with * which amino acids are in the wrong layer and should be mutated.
print '---------------\n'
print 'SASA:\t',''.join(str(x) for x in SASA)
print 'SS:\t',''.join(str(x) for x,y in lis)
print 'AA:\t',''.join(str(x) for x in aminos2)
print 'Mutate:\t',''.join(str(x) for x in mutate)
