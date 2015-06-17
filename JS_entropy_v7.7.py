#! /usr/bin/python
import numpy as np
import sys
import os
import random
from string import split

#	This version describes totCDR3 = aa+vj
#	INCLUDES SAMPLED CONTROL FOR JS DIVERGENCE:  File1 downsampled to the size of File2
#	CREATES NEW DIRECTORY FOR OUTPUT FILES BASED ON INPUT FILE NAMES
#	GENERATES 3 ABUNDANCE OUTPUT FILES:  *_cdr3.out, *_VJ.out, *_tot.out, all tab-delimited, which contain reads for each identifier in each population and blanks for absence
#	GENERATES 1 STATSTICS OUTPUT FILE:  *_H_CL_JS.out, tab-delimited, formatted as follows:
#		Row1: Header
#		Rows2-3:  Col0 = input file name, entropy for each identifier type (Col1-3), clonality for each identifier type (Col4-6), max entropy for each identifier type (Col7-9), total number of population members for each identifier type (Col10-12)
#		Row4:  Divergence header
#		Row5:  Col0 = Joined sample name, JSD for each identifier type (Col1-3), number of total members for each identifier type (Col4-6), JSD of File1 with control population for each identifier type (Col7-9), variance in JSD of File1 with control population for each identifier type across downsampling iterations (Col10-12) - number of iterations currently 10
#		Row6:  Control header
#		Row7:  Col0 = control population name, entropy (Col1-3), clonality (Col4-6), max entropy (Col7-9), total number of population members (Col10-12) in control population
#		Row8:  Variance in the values of Row7 over downsampling iterations

#	NOTE: the naming conventions of the rows is conducive to extracting lines from batches of output directories
#	NOTE: Output includes Jensen-Shannon Divergence, NOT Jensen-Shannon Divergence Metric!!!!

def wholenum(x):
	return x >= 1

def div_limited(dictname):	# dictname is a dictionary whose values are tuples of a p and q, zeroes will be removed for some functions
	a = 0
	c = 0
	for key in dictname.keys():
		if c == 0:
			if type(key) == type('') and key.startswith('TR'):
				combo1 = key.split('.')
				cassette = str(combo1[len(combo1)-1])
#				print cassette
				end = cassette.find('J')
				chain = cassette[end-1]
				if chain == 'A' or chain == 'B':
					cassettetype = chain
				else:
					cassettetype = 'P'
			else:
				cassettetype = 'P'
				c = 1
		else:
			break
	list1 = []
	list2 = []
	for key in dictname.keys():
		list1.append(dictname[key][0])
		list2.append(dictname[key][1])
	pvalues1 = np.array(list1).astype(float)
	qvalues1 = np.array(list2).astype(float)

	px1=pvalues1/pvalues1.sum()
	qx1=qvalues1/qvalues1.sum()
       
	a = len(pvalues1)	#this is the total number of members in the population, whether or not they are zeroes in one population

	# for entropy calculations, get rid of zeroes
      
	pvalues2=np.array(filter(wholenum,list1)).astype(float)
	px2=pvalues2/pvalues2.sum()
	px2=px2[px2.nonzero()]
	Hp2=sum(-px2*np.log2(px2))
	Hpmax2 = -np.log2(1/float(len(pvalues2)))

	qvalues2=np.array(filter(wholenum,list2)).astype(float)
	qx2=qvalues2/qvalues2.sum()
	qx2=qx2[qx2.nonzero()]
	Hq2=sum(-qx2*np.log2(qx2))
	Hqmax2 = -np.log2(1/float(len(qvalues2)))

	mx1=(px1+qx1)/2					#	The M distribution is generated including those zeroes 
	Hm1=sum(-mx1*np.log2(mx1))
	JSpq = Hm1-0.5*(Hp2+Hq2)

	p_ans = Hp2,Hpmax2,len(pvalues2),sum(pvalues2)
	q_ans = Hq2,Hqmax2,len(qvalues2),sum(qvalues2)
	pq_JS = JSpq,a

	return p_ans,q_ans,pq_JS
	
#	infile1 and infile2 in typical output file format:  header row, Col0 = vjcombo, Col1 = reads, Col2 = V, Col3 = J, Col4 = aa, Col5 = nt

infile1 = sys.argv[1]	#First input file, usually PBMC
infile2 = sys.argv[2]	#Second input file, usually BRAIN
instrip1 = infile1.rstrip('_.productive.tsv')	#	ADJUST TO SUIT INPUT FILE TYPE - instrip NAMES LEAVE MINIMUM UNIQUE IDENTIFYING INFORMATION FOR ECONOMY OF FILE NAMING
instrip2 = infile2.rstrip('_.productive.tsv')

dirname = instrip1 +'_'+ instrip2 + '_v7.7'
cmd1 = 'mkdir %(dirname)s' % vars()
os.system(cmd1)


#	FILE 1 (PBMC) INPUT:  Makes libraries for identifiers, counts unique occurrences of identifiers, counts total rows

dCDR3_1 = {}
dVJ_1 = {}
dTOT_1 = {}

ct1 = 0
input1 = open(infile1,'r')

k1 = 0
i1 = 0
j1 = 0
m1 = 0
for line in input1.readlines():
	if k1 != 0:
		llist = split(line)
		reads1 = int(llist[1])
		ct1 = ct1 + reads1
		combo1 = llist[0]
		aa1 = llist[4]
		try:
			dCDR3_1[aa1] = dCDR3_1[aa1] + reads1
			i1 = i1
		except KeyError:
			dCDR3_1[aa1] = reads1
			i1 = i1+1
		try:
			dVJ_1[combo1] = dVJ_1[combo1] + reads1
			j1 = j1
		except KeyError:
			dVJ_1[combo1] = reads1
			j1 = j1+1
		totid1 = aa1,combo1
		try:
			dTOT_1[totid1] = dTOT_1[totid1] + reads1
			m1 = m1
		except KeyError:
			dTOT_1[totid1] = reads1
			m1 = m1+1
	k1 = k1+1
input1.close()
print infile1,m1,'totCDR3',ct1,'reads'

#	FILE 2 (BRAIN) INPUT:  Makes libraries for identifiers, counts unique occurrences of identifiers, counts total rows

dCDR3_2 = {}
dVJ_2 = {}
dTOT_2 = {}

ct2 = 0
input2 = open(infile2,'r')
k2 = 0
i2 = 0
j2 = 0
m2 = 0
for line in input2.readlines():
	if k2 != 0:
		llist = split(line)
		reads2 = int(llist[1])
		ct2 = ct2 + reads2
		combo2 = llist[0]
		aa2 = llist[4]
		try:
			dCDR3_2[aa2] = dCDR3_2[aa2] + reads2
			i2 = i2
		except KeyError:
			dCDR3_2[aa2] = reads2
			i2 = i2+1
		try:
			dVJ_2[combo2] = dVJ_2[combo2] + reads2
			j2 = j2
		except KeyError:
			dVJ_2[combo2] = reads2
			j2 = j2+1
		totid2 = aa2,combo2
		try:
			dTOT_2[totid2] = dTOT_2[totid2] + reads2
			m2 = m2
		except KeyError:
			dTOT_2[totid2] = reads2
			m2 = m2+1
	k2 = k2+1
input2.close()
print infile2,m2,'totCDR3',ct2,'reads'

#	MAKE LIBRARIES OF MERGED POPULATIONS

dCDR3_share = {}
dCDR3_share.update(dCDR3_1)
dCDR3_share.update(dCDR3_2)
dVJ_share = {}
dVJ_share.update(dVJ_1)
dVJ_share.update(dVJ_2)
dTOT_share = {}
dTOT_share.update(dTOT_1)
dTOT_share.update(dTOT_2)

#	MERGED OUTPUT FILES

outfile1 = dirname + '/' + str(instrip1) + '_' + str(instrip2) + '_cdr3.out'

output1 = open(outfile1,'w')
h1 = 'aaCDR3'
h2 = instrip1
h3 = instrip2
output1.write('%(h1)s\t%(h2)s\t%(h3)s  \n' % vars())

for aa in dCDR3_share.keys():
	pt1 = aa
	try:
		pt2 = int(dCDR3_1[aa])
	except KeyError:
		pt2 = ''
	try:
		pt3 = int(dCDR3_2[aa])
	except KeyError:
		pt3 = ''
	output1.write('%(pt1)s\t%(pt2)s\t%(pt3)s  \n' % vars())
	if pt2 == '':
		pt2 = 0
	if pt3 == '':
		pt3 = 0
	dCDR3_share[aa] = pt2,pt3
output1.close()

outfile2 = dirname + '/' + str(instrip1) + '_' + str(instrip2) + '_VJ.out'

output2 = open(outfile2,'w')
h1 = 'VJcombo'
h2 = instrip1
h3 = instrip2
output2.write('%(h1)s\t%(h2)s\t%(h3)s  \n' % vars())

for combo in dVJ_share.keys():
	pt1 = combo
	try:
		pt2 = int(dVJ_1[combo])
	except KeyError:
		pt2 = ''
	try:
		pt3 = int(dVJ_2[combo])
	except KeyError:
		pt3 = ''
	output2.write('%(pt1)s\t%(pt2)s\t%(pt3)s  \n' % vars())
	if pt2 == '':
		pt2 = 0
	if pt3 == '':
		pt3 = 0
	dVJ_share[combo] = pt2,pt3
output2.close()

outfile3 = dirname + '/' + str(instrip1) + '_' + str(instrip2) + '_tot.out'

output3 = open(outfile3,'w')
h1 = 'aaCDR3'
h2 = 'VJcombo'
h3 = instrip1
h4 = instrip2
output3.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s  \n' % vars())

for totid in dTOT_share.keys():
	pt1 = totid[0]
	pt2 = totid[1]
	try:
		pt4 = int(dTOT_1[totid])
	except KeyError:
		pt4 = ''
	try:
		pt5 = int(dTOT_2[totid])
	except KeyError:
		pt5 = ''
	output3.write('%(pt1)s\t%(pt2)s\t%(pt4)s\t%(pt5)s  \n' % vars())
	if pt4 == '':
		pt4 = 0
	if pt5 == '':
		pt5 = 0
	dTOT_share[totid] = pt4,pt5
output3.close()

print instrip1,instrip2,'merged'

#	PERFORM ENTROPY & DIVERGENCE CALCULATIONS ON THE MERGED POPULATION DICTIONARIES

aa_ans = div_limited(dCDR3_share)
vj_ans = div_limited(dVJ_share)
tot_ans = div_limited(dTOT_share)

#	DIVERGENCE CONTROL: randomly selects reads from the File1 totCDR3 population until the control population contains as many totCDR3 as the File2 population

tot_list = []

for totid in dTOT_1.keys():		#	RE-CREATE THE FULL POPULATION BASED ON READ COUNTS IN dTOT_1
	reads = int(dTOT_1[totid])
	for e in range(0,reads):
		tot_list.append(totid)
		
ctrlH_aa = []
ctrlHmax_aa =[]
ctrlCL_aa =[]
ctrln_aa = []
ctrlJS_aa = []

ctrlH_vj = []
ctrlHmax_vj =[]
ctrlCL_vj =[]
ctrln_vj = []
ctrlJS_vj = []

ctrlH_tot = []
ctrlHmax_tot =[]
ctrlCL_tot =[]
ctrln_tot = []
ctrlJS_tot = []

for h in range(0,10):					#	CHANGE THE NUMBER OF ITERATIONS FOR GENERATING A DOWNSAMPLED CONTROL POPULATION HERE
	dTOT_ctrl = {}
	dCDR3_ctrl = {}
	dVJ_ctrl = {}
	g = 0
	for w in range(0,len(tot_list)):
		if g < m2:				#	CRITERIA FOR SAMPLING CAN BE CHANGED HERE
			totid = random.choice(tot_list)	#	CHANGE THE SAMPLING METHOD HERE. INCLUDING A RANDOM.SHUFFLE() STEP DID NOT CHANGE PROPERTIES OF THE RESULTING POPs.
			try:
				dTOT_ctrl[totid] = dTOT_ctrl[totid] + 1		#	NOTE: tot_list IS THE FULL SYNTHETIC POPULATION, NOT THE HISTOGRAM
				g = g
			except KeyError:
				dTOT_ctrl[totid] = 1
				g = g+1
	#print str(h) + ': ' + str(g)			#	DE-COMMENT TO PRINT AN UPDATE ON THE NUMBER OF totCDR3 FOR EACH ITERATION
	f = 0
	e = 0
	for totid in dTOT_ctrl.keys():			#	NOW CHARACTERIZE THE SAMPLES POPULATION
		aa = totid[0]
		combo = totid[1]
		reads = int(dTOT_ctrl[totid])
		try:
			dCDR3_ctrl[aa] = dCDR3_ctrl[aa] + reads
			f = f
		except KeyError:
			dCDR3_ctrl[aa] = reads
			f = f+1
		try:
			dVJ_ctrl[combo] = dVJ_ctrl[combo] + reads
			e = e
		except KeyError:
			dVJ_ctrl[combo] = reads
			e = e+1
	for aa in dCDR3_1.keys():
		reads = dCDR3_1[aa]
		try:
			ctrlreads = dCDR3_ctrl[aa]
			dCDR3_ctrl[aa] = reads,ctrlreads
		except KeyError:
			ctrlreads = 0
			dCDR3_ctrl[aa] = reads,ctrlreads
	for combo in dVJ_1.keys():
		reads = dVJ_1[combo]
		try:
			ctrlreads = dVJ_ctrl[combo]
			dVJ_ctrl[combo] = reads,ctrlreads
		except KeyError:
			ctrlreads = 0
			dVJ_ctrl[combo] = reads,ctrlreads
	for totid in dTOT_1.keys():
		reads = dTOT_1[totid]
		try:
			ctrlreads = dTOT_ctrl[totid]
			dTOT_ctrl[totid] = reads,ctrlreads
		except KeyError:
			ctrlreads = 0
			dTOT_ctrl[totid] = reads,ctrlreads
	
	#	PERFORM THE ENTROPY/DIVERGENCE CALCULATIONS ON File1 vs CONTROL
		
	aa_ctrl = div_limited(dCDR3_ctrl)
	vj_ctrl = div_limited(dVJ_ctrl)
	tot_ctrl = div_limited(dTOT_ctrl)

	#	EXTRACT THE OUTPUT AND APPEND FOR EACH ITERATION
	
	aactrl_diverge = aa_ctrl[2]
	ctrlJS_aa.append(float(aactrl_diverge[0]))
	
	aa_ctrlinfo = aa_ctrl[1]
	ctrlH_aa.append(float(aa_ctrlinfo[0]))
	ctrlHmax_aa.append(float(aa_ctrlinfo[1]))
	ctrlCL_aa.append(1-float(aa_ctrlinfo[0])/float(aa_ctrlinfo[1]))
	ctrln_aa.append(float(aa_ctrlinfo[2]))
	
	vjctrl_diverge = vj_ctrl[2]
	ctrlJS_vj.append(float(vjctrl_diverge[0]))
	
	vj_ctrlinfo = vj_ctrl[1]
	ctrlH_vj.append(vj_ctrlinfo[0])
	ctrlHmax_vj.append(float(vj_ctrlinfo[1]))
	ctrlCL_vj.append(1-float(vj_ctrlinfo[0])/float(vj_ctrlinfo[1]))
	ctrln_vj.append(vj_ctrlinfo[2])
	
	totctrl_diverge = tot_ctrl[2]
	ctrlJS_tot.append(totctrl_diverge[0])
	
	tot_ctrlinfo = tot_ctrl[1]
	ctrlH_tot.append(tot_ctrlinfo[0])
	ctrlHmax_tot.append(float(tot_ctrlinfo[1]))
	ctrlCL_tot.append(1-float(tot_ctrlinfo[0])/float(tot_ctrlinfo[1]))
	ctrln_tot.append(tot_ctrlinfo[2])
	

#	OUTPUT FILE WITH ALL THE STATISTICS

outfile4 = dirname + '/' + str(instrip1) + '_' + str(instrip2) + '_H_CL_JS.out'

output4 = open(outfile4,'w')
h1 = 'FileName'
h2 = 'Hcdr3'
h3 = 'Hvj'
h4 = 'Htot'
h5 = 'CLcdr3'
h6 = 'CLHvj'
h7 = 'CLtot'
h8 = 'Hcdr3_max'
h9 = 'Hvj_max'
h10 = 'Htot_max'
h11 = 'Num_CDR3'
h12 = 'Num_VJ'
h13 = 'Num_totCDR3'
output4.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s\t%(h5)s\t%(h6)s\t%(h7)s\t%(h8)s\t%(h9)s\t%(h10)s\t%(h11)s\t%(h12)s\t%(h13)s   \n' % vars())

for f in range(0,2):
	aa_info = aa_ans[f]
	vj_info = vj_ans[f]
	tot_info = tot_ans[f]
	if f == 0:
		pt1 = instrip1
	elif f == 1:
		pt1 = instrip2
	#Hs
	pt2 = aa_info[0]
	pt3 = vj_info[0]
	pt4 = tot_info[0]
	#Hmax
	pt8 = aa_info[1]
	pt9 = vj_info[1]
	pt10 = tot_info[1]
	#CL
	pt5 = 1-pt2/pt8
	pt6 = 1-pt3/pt9
	pt7 = 1-pt4/pt10
	#num
	pt11 = aa_info[2]
	pt12 = vj_info[2]
	pt13 = tot_info[2]
	output4.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\t%(pt5)s\t%(pt6)s\t%(pt7)s\t%(pt8)s\t%(pt9)s\t%(pt10)s\t%(pt11)s\t%(pt12)s\t%(pt13)s   \n' % vars())

h14 = 'JS_Divergence'
h15 = 'JScdr3'
h16 = 'JSvj'
h17 = 'JStotCDR3'
h18 = 'cdr3_COMBINEDnum'
h19 = 'vj_COMBINEDnum'
h20 = 'tot_COMBINEDnum'
h21 = 'JScdr3_ctrl'
h22 = 'JSvj_ctrl'
h23 = 'JStot_ctrl'
h24 = 'JScdr3_ctrlVAR'
h25 = 'JSvj_ctrlVAR'
h26 = 'JStot_ctrlVAR'

output4.write('%(h14)s\t%(h15)s\t%(h16)s\t%(h17)s\t%(h18)s\t%(h19)s\t%(h20)s\t%(h21)s\t%(h22)s\t%(h23)s\t%(h24)s\t%(h25)s\t%(h26)s   \n' % vars())

pt14 = 'JS_'+str(instrip1)+'_'+str(instrip2)

aa_diverge = aa_ans[2]
vj_diverge = vj_ans[2]
tot_diverge = tot_ans[2]
#	JS
pt15 = aa_diverge[0]
pt16 = vj_diverge[0]
pt17 = tot_diverge[0]
#	num
pt18 = aa_diverge[1]
pt19 = vj_diverge[1]
pt20 = tot_diverge[1]
#	JS_ctrl_means
pt21 = np.mean(ctrlJS_aa)
pt22 = np.mean(ctrlJS_vj)
pt23 = np.mean(ctrlJS_tot)
#	JS_ctrl_vars
pt24 = np.var(ctrlJS_aa)
pt25 = np.var(ctrlJS_vj)
pt26 = np.var(ctrlJS_tot)

output4.write('%(pt14)s\t%(pt15)s\t%(pt16)s\t%(pt17)s\t%(pt18)s\t%(pt19)s\t%(pt20)s\t%(pt21)s\t%(pt22)s\t%(pt23)s\t%(pt24)s\t%(pt25)s\t%(pt26)s   \n' % vars())

h27 = 'CTRL_Name'
h28 = 'Hcdr3_ctrl'
h29 = 'Hvj_ctrl'
h30 = 'Htot_ctrl'
h31 = 'CLcdr3_ctrl'
h32 = 'CLvj_ctrl'
h33 = 'CLtot_ctrl'
h34 = 'Hcdr3_max_ctrl'
h35 = 'Hvj_max_ctrl'
h36 = 'Htot_max_ctrl'
h37 = 'Num_CDR3_ctrl'
h38 = 'Num_VJ_ctrl'
h39 = 'Num_totCDR3_ctrl'

output4.write('%(h27)s\t%(h28)s\t%(h29)s\t%(h30)s\t%(h31)s\t%(h32)s\t%(h33)s\t%(h34)s\t%(h35)s\t%(h36)s\t%(h37)s\t%(h38)s\t%(h39)s   \n' % vars())

pt27 = 'CTRL_'+str(instrip1)+'_'+str(instrip2)

#	H_sample
pt28 = np.mean(ctrlH_aa)
pt29 = np.mean(ctrlH_vj)
pt30 = np.mean(ctrlH_tot)
#	CL_sample
pt31 = np.mean(ctrlCL_aa)
pt32 = np.mean(ctrlCL_vj)
pt33 = np.mean(ctrlCL_tot)
#	Hmax_sample
pt34 = np.mean(ctrlHmax_aa)
pt35 = np.mean(ctrlHmax_vj)
pt36 = np.mean(ctrlHmax_tot)
#	num_sample
pt37 = np.mean(ctrln_aa)
pt38 = np.mean(ctrln_vj)
pt39 = np.mean(ctrln_tot)

output4.write('%(pt27)s\t%(pt28)s\t%(pt29)s\t%(pt30)s\t%(pt31)s\t%(pt32)s\t%(pt33)s\t%(pt34)s\t%(pt35)s\t%(pt36)s\t%(pt37)i\t%(pt38)i\t%(pt39)i   \n' % vars())

pt40 = 'CTRL_VAR_'+str(instrip1)+'_'+str(instrip2)

#	H_sample
pt41 = np.var(ctrlH_aa)
pt42 = np.var(ctrlH_vj)
pt43 = np.var(ctrlH_tot)
#	CL_sample
pt44 = np.var(ctrlCL_aa)
pt45 = np.var(ctrlCL_vj)
pt46 = np.var(ctrlCL_tot)
#	Hmax_sample
pt47 = np.var(ctrlHmax_aa)
pt48 = np.var(ctrlHmax_vj)
pt49 = np.var(ctrlHmax_tot)
#	num_sample
pt50 = np.var(ctrln_aa)
pt51 = np.var(ctrln_vj)
pt52 = np.var(ctrln_tot)

output4.write('%(pt40)s\t%(pt41)s\t%(pt42)s\t%(pt43)s\t%(pt44)s\t%(pt45)s\t%(pt46)s\t%(pt47)s\t%(pt48)s\t%(pt49)s\t%(pt50)s\t%(pt51)s\t%(pt52)s   \n' % vars())

output4.close()

batmansays = ['BOOM!','KAPOW!','SHAZZAM!','WHAMMY!','VRONK!','SPLAT!','BANG!','WHAP!','ZOWIE!','SPLAT!','BAM!(just kidding)','CLANK!(just kidding)','AWK!(just kidding)']

print instrip1,instrip2,'...'+str(random.choice(batmansays))
