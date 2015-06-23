#!/usr/bin/env python

import sys
from operator import itemgetter

file=sys.argv[1]

# Cassette Action
D={'TRBV1':'del',
'TRBV10-1':'TRBV10-1',
'TRBV10-2':'TRBV10-2',
'TRBV10-3':'TRBV10-3',
'TRBV11-1':'TRBV11-1',
'TRBV11-2':'TRBV11-2',
'TRBV11-3':'TRBV11-3',
'TRBV12-2':'TRBV12-2',
'TRBV12-3':'TRBV12-3/4',
'TRBV12-4':'TRBV12-3/4',
'TRBV12-5':'TRBV12-5',
'TRBV13':'TRBV13',
'TRBV14':'TRBV14',
'TRBV15':'TRBV15',
'TRBV16':'TRBV16',
'TRBV17':'del',
'TRBV18':'TRBV18',
'TRBV19':'TRBV19',
'TRBV2':'TRBV2',
'TRBV20-1':'TRBV20-1',
'TRBV23-1':'del',
'TRBV24-1':'TRBV24-1',
'TRBV25-1':'TRBV25-1',
'TRBV26':'TRBV26',
'TRBV27':'TRBV27',
'TRBV28':'TRBV28',
'TRBV29-1':'TRBV29-1',
'TRBV3-1':'TRBV3-1',
'TRBV3-2':'TRBV3-1',
'TRBV4-1':'TRBV4-1',
'TRBV4-2':'TRBV4-2',
'TRBV4-3':'TRBV4-3',
'TRBV5-1':'TRBV5-1',
'TRBV5-2':'TRBV5-2',
'TRBV5-3':'TRBV5-5',
'TRBV5-4':'TRBV5-4',
'TRBV5-5':'TRBV5-5',
'TRBV5-6':'TRBV5-6',
'TRBV5-7':'TRBV5-6',
'TRBV5-8':'TRBV5-8',
'TRBV6-1':'TRBV6-1',
'TRBV6-2':'TRBV6-2/3',
'TRBV6-3':'TRBV6-2/3',
'TRBV6-4':'TRBV6-4',
'TRBV6-5':'TRBV6-5',
'TRBV6-6':'TRBV6-6',
'TRBV6-7':'TRBV6-7',
'TRBV6-8':'TRBV6-8',
'TRBV6-9':'TRBV6-9',
'TRBV7-1':'TRBV7-1',
'TRBV7-2':'TRBV7-2',
'TRBV7-3':'TRBV7-3',
'TRBV7-4':'TRBV7-4',
'TRBV7-5':'TRBV7-5',
'TRBV7-6':'TRBV7-6',
'TRBV7-7':'TRBV7-7',
'TRBV7-8':'TRBV7-8',
'TRBV7-9':'TRBV7-9',
'TRBV8-2':'TRBV8-2',
'TRBV9':'TRBV9',
'TRBJ1-1':'TRBJ1-1',
'TRBJ1-2':'TRBJ1-2',
'TRBJ1-3':'TRBJ1-3',
'TRBJ1-4':'TRBJ1-4',
'TRBJ1-5':'TRBJ1-5',
'TRBJ1-6':'TRBJ1-6',
'TRBJ2-1':'TRBJ2-1',
'TRBJ2-2':'TRBJ2-2',
'TRBJ2-3':'TRBJ2-3',
'TRBJ2-4':'TRBJ2-4',
'TRBJ2-5':'TRBJ2-5',
'TRBJ2-6':'TRBJ2-6',
'TRBJ2-7':'TRBJ2-7',
'Unknown':'del',
'TRAV10':'TRAV10',
'TRAV1-1':'TRAV1-1',
'TRAV1-2':'TRAV1-2',
'TRAV12-1':'TRAV12-1',
'TRAV12-2':'TRAV12-2',
'TRAV12-3':'TRAV12-3',
'TRAV13-1':'TRAV13-1',
'TRAV13-2':'TRAV13-2',
'TRAV14DV4':'TRAV14DV4',
'TRAV16':'TRAV16',
'TRAV17':'TRAV17',
'TRAV18':'TRAV18',
'TRAV19':'TRAV19',
'TRAV2':'TRAV2',
'TRAV20':'TRAV20',
'TRAV21':'TRAV21',
'TRAV22':'TRAV22',
'TRAV23DV6':'TRAV23DV6',
'TRAV24':'TRAV24',
'TRAV25':'TRAV25',
'TRAV26-1':'TRAV26-1',
'TRAV26-2':'TRAV26-2',
'TRAV27':'TRAV27',
'TRAV29DV5':'TRAV29DV5',
'TRAV3':'TRAV3',
'TRAV30':'TRAV30',
'TRAV32':'del',
'TRAV34':'TRAV34',
'TRAV35':'TRAV35',
'TRAV36DV7':'TRAV36DV7',
'TRAV38-1':'TRAV38-1',
'TRAV38-2DV8':'TRAV38-2DV8',
'TRAV39':'TRAV39',
'TRAV4':'TRAV4',
'TRAV40':'TRAV40',
'TRAV41':'TRAV41',
'TRAV5':'TRAV5',
'TRAV6':'TRAV6',
'TRAV7':'TRAV7',
'TRAV8-1':'TRAV8-1',
'TRAV8-2':'TRAV8-2',
'TRAV8-3':'TRAV8-3',
'TRAV8-4':'TRAV8-4',
'TRAV8-5':'TRAV8-3',
'TRAV8-6':'TRAV8-6',
'TRAV8-7':'TRAV8-3',
'TRAV9-1':'TRAV9-1',
'TRAV9-2':'TRAV9-2',
'TRAV28':'del',
'TRAJ1':'del',
'TRAJ10':'TRAJ10',
'TRAJ11':'TRAJ11',
'TRAJ12':'TRAJ12',
'TRAJ13':'TRAJ13',
'TRAJ14':'TRAJ14',
'TRAJ15':'TRAJ15',
'TRAJ16':'TRAJ16',
'TRAJ17':'TRAJ17',
'TRAJ18':'TRAJ18',
'TRAJ2':'del',
'TRAJ20':'TRAJ20',
'TRAJ21':'TRAJ21',
'TRAJ22':'TRAJ22',
'TRAJ23':'TRAJ23',
'TRAJ24':'TRAJ24',
'TRAJ25':'TRAJ25',
'TRAJ26':'TRAJ26',
'TRAJ27':'TRAJ27',
'TRAJ28':'TRAJ28',
'TRAJ29':'TRAJ29',
'TRAJ3':'TRAJ3',
'TRAJ30':'TRAJ30',
'TRAJ31':'TRAJ31',
'TRAJ32':'TRAJ32',
'TRAJ33':'TRAJ33',
'TRAJ34':'TRAJ34',
'TRAJ35':'TRAJ35',
'TRAJ36':'TRAJ36',
'TRAJ37':'TRAJ37',
'TRAJ38':'TRAJ38',
'TRAJ39':'TRAJ39',
'TRAJ4':'TRAJ4',
'TRAJ40':'TRAJ40',
'TRAJ41':'TRAJ41',
'TRAJ42':'TRAJ42',
'TRAJ43':'TRAJ43',
'TRAJ44':'TRAJ44',
'TRAJ45':'TRAJ45',
'TRAJ46':'TRAJ46',
'TRAJ47':'TRAJ47',
'TRAJ48':'TRAJ48',
'TRAJ49':'TRAJ49',
'TRAJ5':'TRAJ5',
'TRAJ50':'TRAJ50',
'TRAJ52':'TRAJ52',
'TRAJ53':'TRAJ53',
'TRAJ54':'TRAJ54',
'TRAJ56':'TRAJ56',
'TRAJ57':'TRAJ57',
'TRAJ58':'TRAJ58',
'TRAJ6':'TRAJ6',
'TRAJ61':'TRAJ61',
'TRAJ7':'TRAJ7',
'TRAJ8':'TRAJ8',
'TRAJ9':'TRAJ9'}

fin=open(file,'r')
head=fin.readline().strip()
data=fin.readlines()

data=[i.strip().split("\t") for i in data]

# replace by new designations
reassigned=[]
for line in data:
	reassigned.append(['@'.join([D[line[2]],D[line[3]],line[4],line[5]]),int(line[1])])


# remove lines containing 'del'
reassigned=[line for line in reassigned if "del" not in line[0]]

# merge identical lines and add up counts
merged={}
for line in reassigned:
	merged[line[0]]=merged.get(line[0],0)+line[1]  # add counts

# write output
S=sorted(merged.items(),key=itemgetter(1),reverse=True)
print head
for line in S:
	count=line[1]
	id=line[0].split("@")
	print('%s.%s\t%d\t%s\t%s\t%s\t%s') %(id[0],id[1],count,id[0],id[1],id[2],id[3])



	

