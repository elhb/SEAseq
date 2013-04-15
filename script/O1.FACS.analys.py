#!/usr/bin/env python
import os
import sys
import re
from time import time
import argparse
import multiprocessing

WELLS = {	'A1':	'TCTCTGTG',
		'A2':	'TGTACGTG',
		'A3':	'ATCGTCTG',
		'A4':	'TAGCTCTG',
		'A5':	'AGTATCTG',
		'A6':	'TCGAGCTG',
		'A7':	'TCATACTG',
		'A8':	'TACGACTG',
		'A9':	'ACTCACTG',
		'A10':	'AGAGTATG',
		'A11':	'AGCTGATG',
		'A12':	'TATCGATG',
		'B1':	'ATGCGATG',
		'B2':	'ACGTCATG',
		'B3':	'TCATGTCG',
		'B4':	'TAGCGTCG',
		'B5':	'TCTACTCG',
		'B6':	'ATGACTCG',
		'B7':	'ATCTATCG',
		'B8':	'ACAGATCG',
		'B9':	'ATACTGCG',
		'B10':	'TATATGCG',
		'B11':	'TGCTCGCG',
		'B12':	'ATCGCGCG',
		'C1':	'TAGTAGCG',
		'C2':	'AGATAGCG',
		'C3':	'TGTGAGCG',
		'C4':	'TCACAGCG',
		'C5':	'ACTGTACG',
		'C6':	'TGCGTACG',
		'C7':	'TCGCTACG',
		'C8':	'TACTGACG',
		'C9':	'AGACGACG',
		'C10':	'TGTAGACG',
		'C11':	'ACGAGACG',
		'C12':	'ATATCACG',
		'D1':	'TCAGCACG',
		'D2':	'TAGACACG',
		'D3':	'AGCACACG',
		'D4':	'ATGTGTAG',
		'D5':	'ACTCGTAG',
		'D6':	'TGCAGTAG',
		'D7':	'TGATCTAG',
		'D8':	'TACGCTAG',
		'D9':	'TCGTATAG',
		'D10':	'AGACATAG',
		'D11':	'AGCGTGAG',
		'D12':	'ATGATGAG',
		'E1':	'ACATCGAG',
		'E2':	'TCTGCGAG',
		'E3':	'ATAGAGAG',
		'E4':	'TATCAGAG',
		'E5':	'ACGCAGAG',
		'E6':	'ACAGTCAG',
		'E7':	'TCTATCAG',
		'E8':	'TAGTGCAG',
		'E9':	'TGACGCAG',
		'E10':	'ATCAGCAG',
		'E11':	'TGCTACAG',
		'E12':	'AGTGACAG',
		'F1':	'ACTGTGTC',
		'F2':	'TACATGTC',
		'F3':	'ATGACGTC',
		'F4':	'AGCGAGTC',
		'F5':	'TCGCAGTC',
		'F6':	'ATACAGTC',
		'F7':	'TGCGTCTC',
		'F8':	'TCACTCTC',
		'F9':	'ATCTGCTC',
		'F10':	'TGTAGCTC',
		'F11':	'ACGTACTC',
		'F12':	'TCTGACTC',
		'G1':	'ACGCTATC',
		'G2':	'ATCATATC',
		'G3':	'TCGTGATC',
		'G4':	'TGACGATC',
		'G5':	'TGCTCATC',
		'G6':	'TATGCATC',
		'G7':	'ACAGCATC',
		'G8':	'AGTACATC',
		'G9':	'AGTGCTGC',
		'G10':	'TGCGATGC',
		'G11':	'ATGCATGC',
		'G12':	'TCACATGC',
		'H1':	'AGAGTCGC',
		'H2':	'ACTATCGC',
		'H3':	'TAGATCGC',
		'H4':	'TCATGCGC',
		'H5':	'TACTACGC',
		'H6':	'ATATACGC',
		'H7':	'TGTCACGC',
		'H8':	'AGTCTAGC',
		'H9':	'ATGTGAGC',
		'H10':	'TAGCGAGC',
		'H11':	'ACACGAGC',
		'H12':	'TCTAGAGC'
		}
WELLSBYINDEX = {}
for well, index in WELLS.iteritems():
    WELLSBYINDEX[index] = well

def main():
    # commandline argument parsing
    argparser = argparse.ArgumentParser(	description='Doing N20 o1 facs analysis.', formatter_class=argparse.RawTextHelpFormatter,)
    argparser.add_argument('-r1',	dest='reads1',	metavar='FILE',	type=file,	required=True, help='Input "fastq"-file read1.')
    argparser.add_argument('-r2',	dest='reads2',	metavar='FILE',	type=file,	required=True, help='Input "fastq"-file read2.')
    argparser.add_argument('-hm',	dest='handleMM',metavar='N',	type=int,	required=False,default=0, help='specify number of missmatches allowed in handle, ie DLA or POSTR sequences (hamming distance).')
    argparser.add_argument('-wm',	dest='wbcmm',metavar='N',	type=int,	required=False,default=0, help='specify number of missmatches allowed in well barcode (hamming distance).')
    argparser.add_argument('-bm',	dest='bbcmm',metavar='N',	type=int,	required=False,default=0, help='specify number of missmatches allowed in bead barcode (hamming distance).')
    argparser.add_argument('-be',	dest='bbced',metavar='N',	type=int,	required=False,default=0, help='specify edit distance allowed in clustering of bead barcodes (levenshtein distance).')
    argparser.add_argument('-mr',	dest='minread',metavar='N',	type=int,	required=False,default=0, help='Minimum number of reads to output barcode (default 0).')
    argparser.add_argument('-mp',	dest='minperc',metavar='N',	type=int,	required=False,default=0, help='Minimum percentage of reads in well to output barcode (default 0).')
    argparser.add_argument('-pdf',	dest='pdf',metavar='N',		type=str,	required=True,default='multipage.pdf', help='Name of output pdffile.')
    indata = argparser.parse_args(sys.argv[1:])

    #	AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT NNNNNNNNNNNNNNNNNNNN CTTGATCCTCTCTGAGCGCAGCGGGCG ---GV--- CAAGATAACGCGTGCTGGTT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --IH-- ATCTCGTATGCCGTCTTCTGCTTG
    #	TTACTATGCCGCTGGTGGCTCTAGATGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGA NNNNNNNNNNNNNNNNNNNN GAACTAGGAGAGACTCGCGTCGCCCGC ---CB--- GTTCTATTGCGCACGACCAA TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTG --ID-- TAGAGCATACGGCAGAAGACGAAC

    #define sequences
    C = sequence('The C sequence','CTAAGTCCATCCGCACTCCT','####################')
    TjH = sequence('The TjH sequence','CCATGTCATACACCGCCTTCAGAGC','####################')
    DLA = sequence('The DLA primer sequence with spacer','CTTGATCCTCTCTGAGCGCAGCGGGCG','CTTGATCCTCTCTGAGCGCAGCGGGCG')
    POSTRH = sequence('The postr primer handle','CAAGATAACGCGTGCTGGTT','CAAGATAACGCGTGCTGGTT')
    MOLBCLEN = 20
    WELLBCLEN = 8
    handle = DLA
    indata.DLA = DLA
    indata.POSTRH = POSTRH
    indata.MOLBCLEN = MOLBCLEN
    
    uniqN20s = {}
    bothreadsmatch = {}
    handle_matches = 0

    #WorkerPool = multiprocessing.Pool(multiprocessing.cpu_count(),maxtasksperchild=100)
    #results = WorkerPool.imap_unordered(magicFunction,clusterGen(indata.reads1, indata.reads2,indata),chunksize=100)
    #results = WorkerPool.imap(magicFunction,clusterGen(indata.reads1, indata.reads2,indata),chunksize=100)
    
    #for each read pair start
    wells = {}
    for well in WELLS: wells[well] = Well(well)
    wells['OTHER'] = Well('OTHER')

    sys.stderr.write('cmd:'+' '.join(sys.argv)+'\n')
    sys.stderr.write('PART1: identification of barcodes and cluster to well classification\n')
    rc = 0
    progress = Progress(bufcount(indata.reads1.name)/4)
    with progress:
	for [cluster,X] in clusterGen(indata.reads1, indata.reads2,indata):   
	    cluster.findDLA(DLA,indata)
	    cluster.findpostr(POSTRH,indata)
	    cluster.getSeqs(MOLBCLEN,indata)
	   
        #for cluster in results:
	    progress.update()
	    
	    try: wells[str(cluster.well)].add(cluster)
	    except KeyError:wells['OTHER'].add(cluster)
    
	    if cluster.fwd.dla:
		# add to total counter of matches
		handle_matches += 1
		# add to counter of uniq N20
		try:	     uniqN20s[cluster.fwd.subseq(max(cluster.fwd.dla_start-MOLBCLEN,0),cluster.fwd.dla_start).seq] += 1
		except KeyError: uniqN20s[cluster.fwd.subseq(max(cluster.fwd.dla_start-MOLBCLEN,0),cluster.fwd.dla_start).seq] = 1
    
		if cluster.fwd.subseq(max(cluster.fwd.dla_start-MOLBCLEN,0),cluster.fwd.dla_start).seq == cluster.rev.subseq(cluster.rev.dla_end,cluster.rev.dla_end+MOLBCLEN).revcomp().seq:
		    # if paired end match add to the "bothreadsmatch"-counter
		    try: bothreadsmatch[cluster.fwd.subseq(0,20).seq] += 1
		    except KeyError: bothreadsmatch[cluster.fwd.subseq(0,20).seq] = 1

		    #try: wells[str(cluster.well)].add(cluster)
		    #except KeyError:wells['OTHER'].add(cluster)

    sys.stderr.write('PART2: within each well;\n 1.try to classify PE missmatching barcodes to other known barcodes in well\n2.reduce the number of barcodes in well by grouping those with low editdistance\n')
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(indata.pdf)
    wellids = []
    for i in range(1,13):
	for i2 in ['A','B','C','D','E','F','G','H']: wellids.append(i2+str(i))
    for wellid in wellids:
	well = wells[wellid]
	if well.id == 'OTHER': continue
	well.double2single(indata)
	well.reducebarcodes(indata,pp)
    pp.close()

    occurance1 = 0
    N = 0
    N20max = ['',0]
    for N20 in uniqN20s:
	if uniqN20s[N20] == 1: occurance1 +=1
	if re.match('N',N20): N += 1
	if uniqN20s[N20] > N20max[1]: N20max = [N20,uniqN20s[N20]]

    occuranceboth = 0
    Nboth = 0
    N20maxboth = ['',0]
    for N20 in bothreadsmatch:
	if bothreadsmatch[N20] == 1: occuranceboth +=1
	if re.match('N',N20): Nboth += 1
	if bothreadsmatch[N20] > N20maxboth[1]: N20maxboth = [N20,bothreadsmatch[N20]]
    

    print str(len(uniqN20s))+ ' unique N20 sequences in '+str(handle_matches) +' reads where DLA is found (out of '+str(cluster.number)+' reads in total, ie DLA is found in '+str(round(handle_matches/float(cluster.number)*100,2))+'%).'
    print 'Out of theese '+str(occurance1)+' N20s are observed only once and '+str(N)+' include a N base.'
    print 'The most frequently observed N20, '+N20max[0]+', has been observed '+str(N20max[1])+' times.'
    print 'In average number of observations per uniqe N20 are '+str(round(sum(uniqN20s.values())/float(len(uniqN20s.values())),2))+' times and the median number of observations is '+str(sorted(uniqN20s.values())[len(uniqN20s.values())/2])+'.'
    print 'In '+str(round(sum(bothreadsmatch.values())/float(handle_matches)*100,2))+'% of the clusters ('+str(sum(bothreadsmatch.values()))+'st) does read2 have the same N20, if only considering theese clusters there are:'
    print '\t' +str(len(bothreadsmatch))+ ' unique N20 sequences in '+str(sum(bothreadsmatch.values())) +' reads'
    print '\tOut of theese '+str(occuranceboth)+' N20s are observed only once and '+str(Nboth)+' include a N base'
    print '\tThe most frequently observed N20, '+N20maxboth[0]+', has been observed '+str(N20maxboth[1])+' times'
    print '\tIn average every N20 has been observed '+str(round(sum(bothreadsmatch.values())/float(len(bothreadsmatch.values())),2))+' times and the median number of observations is '+str(sorted(bothreadsmatch.values())[len(bothreadsmatch.values())/2])
    print str(len(uniqN20s))+ '\t'+str(handle_matches) +'\t'+str(cluster.number)+'\t'+str(occurance1)+'\t'+str(N)+'\t'+str(sum(bothreadsmatch.values()))+'\t'+str(len(bothreadsmatch))+'\t'+str(occuranceboth)+'\t'+str(Nboth)+'\n'

def magicFunction(X):
    [cluster, indata] = X
    cluster.findDLA(indata.DLA,indata)
    cluster.findpostr(indata.POSTRH,indata)
    cluster.getSeqs(indata.MOLBCLEN,indata)
    return cluster
	
class Progress():
    def __init__(self,total, verb='full'):
	import time
	self.total = total
	self.current = 0
	self.type = verb
	self.logfile = sys.stderr
	self.ltime = time.time()
	self.lcurrent = self.current
	self.lpercentage = 0
    
    def __enter__(self):
	if self.type == 'minimal':
	    self.logfile.write('0%                 50%                 100%\n')
    
    def update(self):
	import time
	self.current += 1
	self.percentage = int(round(100*float(self.current)/self.total))
	if self.percentage % 5 == 0 and self.percentage != self.lpercentage:
		self.stf=int(round((self.total-self.current)/((self.current-self.lcurrent)/(time.time()-self.ltime))))
		if self.type == 'full': self.logfile.write(
			'#Progress => '+str(self.percentage)+'%, '+
			str( round((self.current-self.lcurrent)/(time.time()-self.ltime),2) )+' reads/second, '+
			time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+
			', left:'+str(self.stf/60)+'min '+str(self.stf%60)+'s'+
			'\n'
			)
		if self.type == 'minimal': self.logfile.write('..')
		self.ltime = time.time()
		self.lcurrent = self.current
		self.lpercentage = self.percentage

    def __exit__(self, *args):
	self.logfile.write('\n')

class Well():
    def __init__(self, wellid):
	self.id = wellid
	self.barcodes = {}
	self.singlen20clusters = []
	self.doublen20clusters = []
    
    def add(self,cluster):
	if type(cluster.n20) == type('this is a string'):
	    self.singlen20clusters.append(cluster)
	    try: self.barcodes[cluster.n20] += 1
	    except KeyError: self.barcodes[cluster.n20] = 1
	elif type(cluster.n20) == type(['this','is','a','list']): self.doublen20clusters.append(cluster)
	else: raise ValueError
	
    def double2single(self,indata):
	#sys.stderr.write(self.id+' double2single ...')
	temp =[]
	for cluster in self.doublen20clusters:
		mindist=[25,'barcode_in_cluster','barcode_in_well']
		for bc in self.barcodes:
			try: dist1 = hamming_distance(cluster.n20[0],bc)
			except AssertionError: dist1= 25
			try: dist2 = hamming_distance(cluster.n20[1],bc)
			except AssertionError: dist2 = 25
			if dist1 <= indata.bbcmm or dist2 <= indata.bbcmm:
				if dist1 == dist2 and dist1 != 0: pass#print 'barcodes equally bad'
				elif dist1 > dist2:
					if dist2 < mindist[0]: mindist = [dist2, cluster.n20[1], bc]
				elif dist1 < dist2:
					if dist1 < mindist[0]: mindist = [dist1, cluster.n20[0], bc]
		if mindist[0] <= indata.bbcmm:
			cluster.n20 = mindist[2]
			if mindist[0] != 0: cluster.n20mm = mindist[1]
			self.singlen20clusters.append(cluster)
			self.barcodes[mindist[2]] += 1
		else: temp.append(cluster)
	self.doublen20clusters = temp
	#sys.stderr.write(' done.\n')

    @property
    def clustercount(self):
	return len(self.singlen20clusters)#+len(self.doublen20clusters)

    def reducebarcodes(self,indata, pp):
	""" Find most common barcodes in well ( > 10% ??), then try to place other barcodes to this cluster
	"""
	print ''
	print 'Well',self.id, self.clustercount,'read pairs (with 1 n20/cluster)'
	maxdist = indata.bbcmm
	matchfunc = hamming_distance
	if indata.bbced: maxdist = indata.bbced; matchfunc = levenshtein
	
	print self.id, 'precluster'
	for bc, count in self.barcodes.iteritems():
	    percentage = round(100*float(count)/self.clustercount,2)
	    if percentage >= indata.minperc and count >= indata.minread:print '\t',count,'\t',bc,'\t',percentage,'%'
	maxrounds = 100; temp = 0
	go = True
	print 'Doing "clustering":'
	print 'round 0', len(self.barcodes),'barcodes'
	while go:
	    
	    lenbcs = len(self.barcodes)
	    
	    # calculate the percentage of the read population that "support" the current barcode and get the top supported barcodes
	    percentages = {}
	    for bc, count in self.barcodes.iteritems():
		percentage = round(100*float(count)/self.clustercount,2)
		try: percentages[percentage].append(bc)
		except KeyError:percentages[percentage]=[bc]
	    highten = percentages.keys()
	    highten.sort()
	    try: highten = highten[:10]
	    except IndexError: pass
	    
	    # For each barcode in the top ten supported barcodes try to merge all other barcodes to that one
	    merged = []
	    for percentage in highten:
		for template_bc in percentages[percentage]:
		    if template_bc in merged: continue # if this bc already has been merged to another do not use it as template
		    for bc in self.barcodes.keys():
			if bc == template_bc: continue # do not merge to self
			dist = matchfunc(bc,template_bc)
			if dist <= maxdist:
			    self.barcodes[template_bc] += self.barcodes[bc]
			    merged.append(bc)
			    del self.barcodes[bc]
	    
	    temp += 1
	    if temp == maxrounds or lenbcs == len(self.barcodes): go = False
	    
	    print 'round', temp, len(self.barcodes), 'barcodes','(', len(self.barcodes)-lenbcs,'st)'
	
	print self.id, 'postcluster'
	for bc, count in self.barcodes.iteritems():
	    percentage = round(100*float(count)/self.clustercount,2)
	    if percentage >= indata.minperc and count >= indata.minread:print '\t',count,'\t',bc,'\t',percentage,'%'
	    
	print 'plotting'
	import numpy as np
	import matplotlib.pyplot as plt
	plt.figure()
	pos = np.arange(len(self.barcodes.keys()))
	width = 1.0     # gives histogram aspect to the bar diagram
	ax = plt.axes()
	ax.set_xticks(pos + (width / 2))
	ax.set_xticklabels(self.barcodes.keys(),rotation='horizontal')
	plt.bar(pos, self.barcodes.values(), width, color='r')
	#plt.show()
	#plt.savefig(pp,format='pdf',bbox_inches=0)
	plt.suptitle(self.id, fontsize=12)
	pp.savefig()
	plt.close()
	print 'done'


class beadbarcode():
    def __init__(self,seq=None):
	self.clusters = []
	self.seq = seq

class sequence():
    def __init__(self,header,seq,qual):
	self.header = header.rstrip()
	self.qual = qual.rstrip()
	self.seq = seq.rstrip()
	self.len = len(seq.rstrip())

    def subseq(self,start,end):
	return sequence(self.header,self.seq[start:end],self.qual[start:end])

    def revcomp(self):
	''' Takes a sequence and reversecomplements it'''
	complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
	revcompseq = "".join([complement.get(nt.upper(), '') for nt in self.seq[::-1]])
	return sequence(self.header,revcompseq,self.qual[::-1])

    def comp(self):
	''' Takes a sequence and complements it'''
	complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
	revcompseq = "".join([complement.get(nt.upper(), '') for nt in self.seq])
	return sequence(self.header,revcompseq,self.qual)
    
class cluster():
    def __init__(self,fwd,rev,number):
	self.fwd	= fwd
	self.rev	= rev
	self.number	= number
	
    def findDLA(self,DLA,indata):
	[dla,start,end] = self.findHandle(DLA,indata,self.fwd)
	self.fwd.dla = dla
	self.fwd.dla_start = start
	self.fwd.dla_end = end
	[dla,start,end] = self.findHandle(DLA.revcomp(),indata,self.rev)
	self.rev.dla = dla
	self.rev.dla_start = start
	self.rev.dla_end = end
	
    def findpostr(self,POSTRH,indata):
	[postr,start,end] = self.findHandle(POSTRH,indata,self.fwd)
	self.fwd.postr = postr
	self.fwd.postr_start = start
	self.fwd.postr_end = end
	[postr,start,end] = self.findHandle(POSTRH.revcomp(),indata,self.rev)
	self.rev.postr = postr
	self.rev.postr_start = start
	self.rev.postr_end = end
	
    def findHandle(self,handle,indata,read):
	handle_match = re.search(handle.seq,read.seq)
	if handle_match:
	    hstart = handle_match.start()
	    hend = handle_match.end()
	elif indata.handleMM:
	    MMs = [[10000,-1]]
	    for i in range(read.len):
		dist = hamming_distance(handle.seq,read.seq[i:i+handle.len])
		MMs.append([dist,i])
		if dist < indata.handleMM: handle_match = True
		else: handle_match = None
		if i == self.fwd.len-handle.len: break
	    minpair = [10000,0]
	    for pair in MMs:
		if pair[0] < minpair[0]: minpair = pair
	    if minpair[0] < indata.handleMM: handle_match = True
	    hstart = minpair[1]
	    hend = minpair[1]+handle.len
	else:
	    hstart =-1;hend=-1
	return [handle_match,hstart,hend]

    def getSeqs(self,MOLBCLEN,indata):
	wellbcfwd = self.fwd.subseq(self.fwd.dla_end , self.fwd.postr_start ).revcomp().seq
	wellbcrev = self.rev.subseq(self.rev.postr_end, self.rev.dla_start ).seq
	if wellbcfwd == wellbcrev:
	    try: self.well = WELLSBYINDEX[wellbcfwd];
	    except KeyError:
		self.well = wellbcfwd
		if len(wellbcfwd) == 8:
		    for bc in WELLSBYINDEX.keys():
			dist = hamming_distance(wellbcfwd,bc)
			if dist <= indata.wbcmm: self.well = WELLSBYINDEX[bc];break
	else:
	    self.well = [wellbcfwd,wellbcrev]
	    #if len(wellbcfwd) == 8 and len(wellbcrev) == 8: dist1 = hamming_distance(wellbcrev,wellbcfwd)
	    dists = []
	    if len(wellbcfwd) == 8:
		    for bc in WELLSBYINDEX.keys():
			dist = hamming_distance(wellbcfwd,bc)
			if dist <= indata.wbcmm:dists.append([dist,bc])
	    if len(wellbcrev) == 8:
		    for bc in WELLSBYINDEX.keys():
			dist = hamming_distance(wellbcrev,bc)
			if dist <= indata.wbcmm: dists.append([dist,bc])
	    mindist = [1000000000,'BC']
	    if dists:
		for dist, bc in dists:
		    if dist < mindist[0]: mindist = [dist,bc]
		if mindist[0] <= indata.wbcmm: self.well = WELLSBYINDEX[mindist[1]];
		#else:self.well = None
	    #else: self.well = None
	
	n20fwd = self.fwd.subseq(max(self.fwd.dla_start-MOLBCLEN,0),self.fwd.dla_start).seq 
	n20rev = self.rev.subseq(self.rev.dla_end,self.rev.dla_end+MOLBCLEN).revcomp().seq
	if n20fwd == n20rev: self.n20 = n20fwd
	else:self.n20 = [n20fwd,n20rev] # do a fuzzy check with others on the same barcode ie later stage

def clusterGen(reads1,reads2,indata):
    import gzip
    if reads1.name.split('.')[-1] in ['gz','gzip']: reads1 = gzip.open(reads1.name)
    if reads2.name.split('.')[-1] in ['gz','gzip']: reads2 = gzip.open(reads2.name)
    countreads=0
    break_true=False
    while  not break_true:
	countreads+=1
	try:
	    r1_1=reads1.next().rstrip()
	    r1_2=reads1.next().rstrip()
	    r1_3=reads1.next().rstrip()
	    r1_4=reads1.next().rstrip()

	    r2_1=reads2.next().rstrip()
	    r2_2=reads2.next().rstrip()
	    r2_3=reads2.next().rstrip()
	    r2_4=reads2.next().rstrip()
	except StopIteration:
	    break_true=True
	fwd	= sequence(r1_1,r1_2,r1_4)
	rev	= sequence(r2_1,r2_2,r2_4)
	yield [cluster(fwd,rev,countreads),indata]
    reads1.close()
    reads2.close()

def bufcount(filename):
    """ returns the number of lines in a file
    """
    f = open(filename)
    if f.name.split('.')[-1] in ['gz','gzip']:
		import gzip
		f = gzip.open(f.name)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
        f.close
    return lines

def hamming_distance(s1, s2):
    assert len(s1) == len(s2), 'Error: '+str(len(s1)) + ' != ' + str(len(s2))
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if not s1:
        return len(s2)
 
    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
	    #if c1 == 'N' or c2 == 'N': substitutions -= 1 #if N then no mismatch
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
 
    return previous_row[-1]

if __name__ == '__main__':
    main()