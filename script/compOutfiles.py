#! /bin/env python

from tnt.lib import *

outputtype = 2

import sys
variations = {}
samples = []
for pair in sys.argv[1:]:
    sample_id = pair.split(':')[0]
    f = open(pair.split(':')[1],'r')
    samples.append(sample_id)

    for var_id in VARIATIONINFO.keys():
	if var_id not in variations.keys():
	    variations[var_id] = {}
	try:		variations[var_id][sample_id] = {}
	except KeyError:variations[var_id] = { sample_id:{} }

    for line in f:
	line = line.rstrip()
	line = line.split('\t')
	variations[line[0]]['alt'] = line[1]
	variations[line[0]][sample_id]['reads'] = {}
	temp = 1
	base = None
	perc = None
	count = None
	for part in line[2:]:
	    if temp == 4: temp = 1; variations[line[0]][sample_id]['reads'][base] = [count, perc]
	    if temp == 1: base = part
	    if temp == 2: perc = part
	    if temp == 3: count = part
	    temp += 1
	variations[line[0]][sample_id]['reads'][base] = [count, perc]

if outputtype == 1:

    sys.stdout.write('Variation id\tAllele\t'+'\t'.join(['\t'+sample+'\t' for sample in samples])+'\n')
    sys.stdout.write('\t\t'+'\t'.join(['\tFrequency\tCount' for sample in samples])+'\n')

    for var_id, variation in variations.iteritems():

	total_reads = {}
	for sample in samples: total_reads[sample] = 0

	sys.stdout.write(var_id)

	for alt in variation['alt'].split('/'):
	    sys.stdout.write('\t'+alt)
	    for sample in samples:
		sys.stdout.write('\t')
		try :
		    sys.stdout.write('\t'+variations[var_id][sample]['reads'][alt][1]+'\t')
		    sys.stdout.write(variations[var_id][sample]['reads'][alt][0])
		    total_reads[sample] += int(variations[var_id][sample]['reads'][alt][0])
		except KeyError:
		    sys.stdout.write('\t0 %\t0')
	    sys.stdout.write('\n')
	sys.stdout.write('\t\t'+'\t'.join(['\ttot\t'+str(total_reads[sample]) for sample in samples])+'\n')
	sys.stdout.write('\n')

elif outputtype==2:

    for var_id, variation in variations.iteritems():

	total_reads = {}
	for sample in samples: total_reads[sample] = 0
	sys.stdout.write(var_id+'\tA\t\tG\t\tT\t\tC\t\t-\t\tTotal')

	for sample in samples:
	    sys.stdout.write('\n'+sample)
	    for alt in ['A','G','T','C','-']:
		try:
		    sys.stdout.write('\t'+variations[var_id][sample]['reads'][alt][1]+'\t')
		    sys.stdout.write(variations[var_id][sample]['reads'][alt][0])
		    total_reads[sample] += int(variations[var_id][sample]['reads'][alt][0])
		except KeyError:
		    sys.stdout.write('\t0 %\t0')
	    sys.stdout.write('\t'+str(total_reads[sample]))
	sys.stdout.write('\n\n')
