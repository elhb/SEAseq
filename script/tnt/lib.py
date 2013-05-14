import sys

def lib_main(): pass

class Progress():

	def __init__(self,total, verb='full', logfile=sys.stderr, unit='read'):
		import time
		self.total = total
		self.current = 0
		self.type = verb
		self.logfile = logfile
		self.ltime = time.time()
		self.lcurrent = self.current
		self.lpercentage = 0
		if verb == 'full': self.printint = 5
		elif verb == 'minimal':self.printint = 5
		self.unit = unit

	def __enter__(self):
		if self.type == 'minimal': self.logfile.write('0%                 50%                 100%\n')
		#                                              ....................................................................................


	def update(self):
		import time
		self.current += 1
		self.percentage = int(round(100*float(self.current)/self.total))
		if self.percentage % self.printint == 0 and self.percentage != self.lpercentage:
			self.stf=int(round((self.total-self.current)/((self.current-self.lcurrent)/(time.time()-self.ltime))))
			if self.type == 'full': self.logfile.write(
				'#Progress => '+str(self.percentage)+'%, '+
				str( round((self.current-self.lcurrent)/(time.time()-self.ltime),2) )+' '+self.unit+'s/second, '+
				time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+
				', left: '+str(self.stf/60/60)+'h '+str(self.stf/60%60)+'min '+str(self.stf%60)+'s'+
				'\n'
				)
			if self.type == 'minimal': self.logfile.write('..')
			self.ltime = time.time()
			self.lcurrent = self.current
			self.lpercentage = self.percentage

	def __exit__(self, *args):
		self.logfile.write('\n')

def hamming_distance(s1, s2):
	assert len(s1) == len(s2), 'Error: '+str(len(s1)) + ' != ' + str(len(s2))
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def bufcount(filename):
	""" returns the number of lines in a file
	"""
	import gzip
	if filename.split('.')[-1] in ['gz','gzip']: f = gzip.open(filename)
	else: f = open(filename)
	lines = 0
	buf_size = 1024 * 1024
	read_f = f.read # loop optimization
	
	buf = read_f(buf_size)
	while buf:
		lines += buf.count('\n')
		buf = read_f(buf_size)
		f.close
	return lines

class sequence():
	def __init__(self,header,seq,qual):
		self.header = header.rstrip()
		self.qual = qual.rstrip()
		self.seq = seq.rstrip()
		assert len(self.qual) == len(self.seq), 'Error: qual and seq has different lengths!\n'
		self.len = len(seq.rstrip())
	
	def subseq(self,start,end):
		return sequence(self.header,self.seq[start:end],self.qual[start:end])
	
	def revcomp(self):
		''' Takes a sequence and reversecomplements it'''
		complementary = self.comp()
		return sequence(complementary.header,complementary.seq[::-1],complementary.qual[::-1])
	
	def comp(self):
		''' Takes a sequence and complements it'''
		complement = {'A':'T','T':'A',
					  'C':'G','G':'C',
					  'N':'N',
					  'R':'Y','Y':'R',
					  'K':'M','M':'K',
					  'B':'V','V':'B',
					  'D':'H','H':'D',
					  }
		compseq = "".join([complement.get(nt.upper(), '') for nt in self.seq])
		return sequence(self.header,compseq,self.qual)

class read(sequence):
    "Represents one of several reads from a DNA fragment"

def getPairs(indata):
	""" yield one readpair at a time from indata
	"""
	import re
	
	# unpack infiles
	[file1,file2]=[indata.reads1,indata.reads2]

	# generate selftest read set
	if indata.selftest:
		raise NotImplementedError
		return
	
	#check if files are gzipped
	import gzip
	if file1.name.split('.')[-1] in ['gz','gzip']: file1 = gzip.open(file1.name)
	else: file1 = open(file1.name,'r')
	if file2.name.split('.')[-1] in ['gz','gzip']: file2 = gzip.open(file2.name)
	else: file2 = open(file2.name,'r')

	# choose random reads to analyze from fastq files
	if indata.n:
		import random
		numreads=bufcount(indata.reads1.name)/4
		if indata.stop: numreads = indata.stop
		indata.logfile.write('Choosing '+str(indata.n)+' random pairs to analyze ... ')
		readNumbersToPrint = {}
		while len(readNumbersToPrint) < indata.n: readNumbersToPrint[random.randint(indata.skip+1,numreads)] = True
		tempo = readNumbersToPrint.keys()
		tempo.sort()
		readNumbersToPrint = tempo 
		indata.logfile.write('done.\n')
	
	# set the counters to initial values
	counter = 0
	tmp=0
	header="NA"
	if indata.skip: skip =True

	# itarate through fastq files and return readpairs
	for r1line in file1:
		r2line = file2.readline()
		
		tmp+=1 # increment linecount
		
		# skip or stop if option is set on indata
		if indata.skip and tmp < (4*indata.skip) and skip: continue
		elif indata.skip and tmp == (4*indata.skip) and skip: skip=False; tmp =0;continue
		if indata.stop and counter == indata.stop: break

		# depending on line number (within entry) do ...	
		if tmp == 1: #header check match between files
			header=r1line
			if r1line.split(" ")[0] != r2line.split(" ")[0]: indata.logfile.write('Error mismatching headers!');raise ValueError #os.kill(MASTER,1);sys.exit(1);#REALLYNOTOPTIMAL
		elif tmp == 2: # sequence check that its DNA and save sequences till later
			if counter == 1:
				indata.logfile.write('Checking data type of read 1 in pair 1 ... ')
				match = re.match("^[AGTCN]+$",r1line.rstrip())
				if match: indata.logfile.write('this is DNA data.\n')
				else: indata.logfile.write(' this is not a DNA sequence ('+r1line.rstrip()+') could something be wrong with your indata file?.\n\n');raise ValueError#os.kill(MASTER);sys.exit();#REALLYNOTOPTIMAL
			r1seq = r1line.rstrip()
			r2seq = r2line.rstrip()
		elif tmp == 3: # "+"-line do some format check
				if counter in {1:True,67:True,438:True,9675:True,53678:True,864513:True,1337354:True,317955:True,1226844:True,20389:True,118261:True}:
					if r1line[0] != r2line[0] or r1line[0] != '+': indata.logfile.write('Error Format not fastq!');raise ValueError#os.kill(MASTER);sys.exit(1);#REALLYNOTOPTIMAL
		elif tmp == 4: # quality values and end of entry, reset counter and yeild readpair
				tmp=0 # reset line counter
				counter+=1 # increase entry counter
				r1qual = r1line.rstrip() #get qual strings
				r2qual = r2line.rstrip()
				
				#yield readpair
				if not indata.n: yield [readpair(header.rstrip(), read(header.rstrip(),r1seq,r1qual), read(header.rstrip(),r2seq,r2qual)),indata]
				else:
					if counter == readNumbersToPrint[0]:
						yield [readpair(header.rstrip(), read(header.rstrip(),r1seq,r1qual), read(header.rstrip(),r2seq,r2qual)),indata]
						readNumbersToPrint = readNumbersToPrint[1:]
						if len(readNumbersToPrint) == 0: break


class readpair():
    """ object representing an illumina cluster """
    
    def __init__(self,header,r1,r2):
		self.header = header
		self.r1 = r1 #first read
		self.r2 = r2 #second read
		self.threads = []

    def identifythreads(self,indata):
		import re
		# a thread should look like this:
		# gagctgctgcaccatattcctgaac GACCATCACTTAAATCAGGTCCTCC NNNNNNNNNNN AGAGTCAAGTTATTTAAAAAATCTGGCC gctctgaaggcggtgtatgacatgg
		# GAGCTGCTGCACCATATTCCTGAAC CAATCTCCCCTATTATTTCTATCCTATGCCACTCCTGCTCATATCTCTAGTG GCTCTGAAGGCGGTGTATGACATGGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGT
		# gagctgctgcaccatattcctgaac                                                      gctctgaaggcggtgtatgacatgg
		
		exthandle_seq = 'gagctgctgcaccatattcctgaac'.upper()
		tjhandle_seq  = 'gctctgaaggcggtgtatgacatgg'.upper()
		
		matchfunk = hamming_distance
	
		#READ1
		exthandle = None
		tjhandle =  None
		exthandle = re.search('^'+exthandle_seq,self.r1.seq)
		tjhandle =  re.search(tjhandle_seq,self.r1.seq)
		if exthandle and tjhandle:
			thread = tntthread(self.r1.seq[exthandle.end():tjhandle.start()])
			self.threads.append(thread)
		elif indata.handlemm: # do some kind of fuzzy matching to allow for missmatches
			if not exthandle: 
				mindist = [10000,-1]
				for i in range(len(self.r1.seq)):
					if i+len(exthandle_seq) <= len(self.r1.seq): dist = matchfunk(exthandle_seq,self.r1.seq[i:i+len(exthandle_seq)])
					else: dist = 1000
					if dist < mindist[0]: mindist =[dist,i]
				if mindist[0] < indata.handlemm: exthandle_end = i+len(exthandle_seq)
				else: exthandle_end = None
			else: exthandle_end = exthandle.end()
			if not tjhandle: 
				mindist = [10000,-1]
				for i in range(len(self.r1.seq)):
					if i+len(tjhandle_seq) <= len(self.r1.seq): dist = matchfunk(tjhandle_seq,self.r1.seq[i:i+len(tjhandle_seq)])
					else: dist = 1000
					if dist < mindist[0]: mindist =[dist,i]
				if mindist[0] < indata.handlemm: tjhandle_start = i
				else: tjhandle_start = None
			else: tjhandle_start = tjhandle.start()
			if tjhandle_start and exthandle_end:
				thread = tntthread(self.r1.seq[exthandle_end:tjhandle_start])
				self.threads.append(thread)
		
		#READ2
		exthandle = None
		tjhandle =  None
		exthandle = re.search(   revcomp(exthandle_seq),self.r2.seq)
		tjhandle =  re.search('^'+revcomp(tjhandle_seq),self.r2.seq)
		if exthandle and tjhandle:
			thread = tntthread(revcomp(self.r2.seq[tjhandle.end():exthandle.start()]))
			self.threads.append(thread)
		elif indata.handlemm: #### do some kind of fuzzy matching to allow for missmatches
			if not exthandle: 
				mindist = [10000,-1]
				for i in range(len(self.r2.seq)):
					if i+len(exthandle_seq) <= len(self.r2.seq): dist = matchfunk(exthandle_seq,revcomp(self.r2.seq)[i:i+len(exthandle_seq)])
					else: dist = 1000
					if dist < mindist[0]: mindist =[dist,i]
				if mindist[0] < indata.handlemm: exthandle_end = i+len(exthandle_seq)
				else: exthandle_end = None
			else: exthandle_end = exthandle.end()
			if not tjhandle: 
				mindist = [10000,-1]
				for i in range(len(self.r2.seq)):
					if i+len(tjhandle_seq) <= len(self.r2.seq): dist = matchfunk(tjhandle_seq,revcomp(self.r2.seq)[i:i+len(tjhandle_seq)])
					else: dist = 1000
					if dist < mindist[0]: mindist =[dist,i]
				if mindist[0] < indata.handlemm: tjhandle_start = i
				else: tjhandle_start = None
			else: tjhandle_start = tjhandle.start()
			if tjhandle_start and exthandle_end:
				thread = tntthread(revcomp(self.r2.seq)[exthandle_end:tjhandle_start])
				self.threads.append(thread)

class tntthread(sequence):
    #Represents the specific part of a tnt thread without the general handles"""

	def __init__(self,seq):
		self.seq = seq

	def identifyspecific(self,indata,extbyseq,ext_primers,tjbyseq,tj_primers):
		import re
		self.extprimer = None
		self.tjprimer = None
		try:
			for ID in extbyseq[self.seq[:19]]:
				match = re.search('^'+ext_primers[ID], self.seq)
				if match:
					self.extprimer= ID;
		except KeyError: pass
		try:
			for ID in tjbyseq[self.seq[-19:]]:
				match = re.search(tj_primers[ID]+'$', self.seq)
				if match:
					self.tjprimer= ID;
		except KeyError: pass
		
		matchfunk = levenshtein
		
		if not self.tjprimer and indata.primeredit:
			minedit = [300, 'NA']
			for primer_id, seq in tj_primers.iteritems():
				dist = matchfunk(seq,self.seq[-len(seq)::])
				if dist < minedit[0]: minedit = [dist, primer_id]
				if minedit[0] <= indata.primeredit: break
			if minedit[0] <= indata.primeredit: self.tjprimer = primer_id
		if not self.extprimer and indata.primeredit:
			minedit = [300, 'NA']
			for primer_id, seq in ext_primers.iteritems():
				dist = matchfunk(seq,self.seq[:len(seq)])
				if dist < minedit[0]: minedit = [dist, primer_id]
				if minedit[0] <= indata.primeredit: break
			if minedit[0] <= indata.primeredit: self.extprimer = primer_id
		
		#if not self.extprimer:
		#	mm=3
		#	vote=[]
		#	for i in range(19):
		#		#print i, self.seq[i], extprimerpositions[i][self.seq[i]]
		#		for ID in extprimerpositions[i][self.seq[i]]:
		#			vote.append(ID)
		#	if vote.count(max(vote)) > 19-mm: print 'update ext:',max(vote)
		return

def comp(str):
	return str.replace("A","X").replace("T","A").replace("X","T").replace("G","X").replace("C","G").replace("X","C")

def revcomp(str):
	return comp(str[::-1])

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
			if c1 == 'N' or c2 == 'N': substitutions -= 1 #if N then no mismatch
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]

class Variation():
	
	def __init__(self,var_id):
		self.id = var_id
		self.threads = []
		self.reference = sequence(self.id,VARIATIONINFO[self.id]['ref'],VARIATIONINFO[self.id]['ref'])
		self.snp = VARIATIONINFO[self.id]['snp']
		self.output = 'No alignment performed'

	
	def add(self, thread):
		self.threads.append(thread)
	
	@property
	def len(self):
		return len(self.threads)
	
	def alignthreads(self):
		
		from Bio import pairwise2
		from Bio.SubsMat import MatrixInfo as matlist
		import re
		
		matrix = matlist.genetic
		gap_open = -10
		gap_extend = -1.5

		var_pos_bases = {}
		
		for thread in self.threads:
			
			if 'N' in self.reference.seq: raise ValueError
			
			alns = pairwise2.align.localds(self.reference.seq.replace(self.snp,'N'), thread.seq, matrix, gap_open, gap_extend)
			top_aln = alns[0]
			#read.aligned_ref, read.aligned_seq, read.aln_score, begin, end = top_aln
			aligned_ref = top_aln[0]
			aligned_seq = top_aln[1]
			aln_score   = top_aln[2]
			
			var_aln_pos = re.search('N',aligned_ref)
			if not var_aln_pos: raise Error

			try: 			var_pos_bases[ aligned_seq[var_aln_pos.start()] ] += 1
			except KeyError:var_pos_bases[ aligned_seq[var_aln_pos.start()] ] = 1
		total = sum([count for base,count in var_pos_bases.iteritems()])

		self.output = self.id +'\t' + '/'.join(uipac(self.snp,'bases'))
		for base,count in var_pos_bases.iteritems():
			self.output += '\t'+ base + '\t' + str(round(100*float(count)/total,2)) + ' %\t' + str(count)
		self.output += '\n'

def uipac(bases, back='uipac'): #U	Uracil NOT SUPPORTED!!!
	if back == 'uipac':
		if 'N' in bases: return 'N'
		uniqbases={}
		for i in bases:
			uniqbases[i]=True
		bases = uniqbases.keys()
		if 'U' in bases: sys.stderr.write('WARNING in function "uipac(bases)": Uracil NOT SUPPORTED!')
		if len(bases)==1:
			if 'A' in bases: return 'A' #A	Adenine
			if 'C' in bases: return 'C' #C	Cytosine
			if 'G' in bases: return 'G' #G	Guanine
			if 'T' in bases: return 'T' #T	Thymine
			#U	Uracil NOT SUPPORTED!!!
		if len(bases)==2:
			if 'A' in bases and 'G' in bases: return 'R' #R	Purine (A or G)
			if 'C' in bases and 'T' in bases: return 'Y' #Y	Pyrimidine (C, T, or U)
			if 'A' in bases and 'C' in bases: return 'M' #M	C or A
			if 'T' in bases and 'G' in bases: return 'K' #K	T, U, or G
			if 'A' in bases and 'T' in bases: return 'W' #W	T, U, or A
			if 'C' in bases and 'G' in bases: return 'S' #S	C or G
		if len(bases)==3:
			if 'C' in bases and 'T' in bases and 'G' in bases: return 'B' #B	C, T, U, or G (not A)
			if 'A' in bases and 'T' in bases and 'G' in bases: return 'D' #D	A, T, U, or G (not C)
			if 'A' in bases and 'T' in bases and 'C' in bases: return 'H' #H	A, T, U, or C (not G)
			if 'A' in bases and 'C' in bases and 'G' in bases: return 'V' #V	A, C, or G (not T, not U)
		if len(bases)==4:
			return 'N' #N	Any base (A, C, G, T, or U)
	elif back == 'bases':
		if   bases == 'R': return ['A','G'] 
		elif bases == 'Y': return ['C','T']
		elif bases == 'M': return ['A','C']
		elif bases == 'K': return ['G','T']
		elif bases == 'W': return ['A','T']
		elif bases == 'S': return ['C','G']
		elif bases == 'B': return ['C','T','G']
		elif bases == 'D': return ['A','T','G']
		elif bases == 'V': return ['A','C','G']
		elif bases == 'H': return ['A','C','T']
		elif bases == 'N': return ['A','G','T','C']

class SEAseqpair(readpair):
	
	def getN15(self):
		if self.handle_start:self.n15 = self.r1.subseq(0,self.handle_start)
		else: self.n15 = None
		return 0

    
	def identify(self, handle, indata):
		[handle_start, handle_end] = self.matchHandle(handle, indata, self.r1)
		self.handle_start = handle_start
		self.handle_end   = handle_end
		return 0

	def identifyIllumina(self, indata):
		handle = sequence('illuminaUniversal','AGATCGGAAGAGC','AGATCGGAAGAGC')
		indata.handlemm = 2
		#handle = sequence('illumina','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
		#indata.handlemm = 5
		[handle_start, handle_end] = self.matchHandle(handle, indata, self.r1)
		if handle_start: self.r1.illuminaadapter = True
		else: self.r1.illuminaadapter = False
		
		#handle = sequence('illumina','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
		#indata.handlemm = 5
		[handle_start, handle_end] = self.matchHandle(handle, indata, self.r2)
		if handle_start: self.r2.illuminaadapter = True
		else: self.r2.illuminaadapter = False
		return 0

	def matchHandle(self, handle, indata, read, matchfunk=hamming_distance):
		
		import re
		#matchfunk = hamming_distance
	
		handle_start = None
		handle_end   = None
	
		perfect_match = re.search(handle.seq, read.seq)
		if perfect_match:
		    handle_start = perfect_match.start()
		    handle_end = perfect_match.end()
		
		elif indata.handlemm:
			mindist = [10000,-1]
			for i in range(len(read.seq)):
			    
				if i+len(handle.seq) <= len(read.seq):
					dist = matchfunk(handle.seq,read.seq[i:i+len(handle.seq)])
				else: dist = 1000
				
				if dist < mindist[0]: mindist =[dist,i]

			if mindist[0] < indata.handlemm:
				handle_start = i
				handle_end = i+len(handle.seq)
			else:
				handle_start = None
				handle_end = None

		return [handle_start, handle_end]

	def get_cid(self,indata):
		if self.n15 and self.n15.len == 15:
			try: self.cid = indata.cid_by_bc[self.n15.seq]
			except KeyError: self.cid = False
		else: self.cid = None

class SEAseqSummary():
	
	def __init__(self):
		self.barcodes = {}
		self.readcount = 0
		self.handlefound = 0
		self.pairs = {}
		
	def add(self, pair):
		self.readcount += 1
		if pair.handle_start: self.handlefound += 1
		if pair.n15 and pair.n15.len == 15:
			try:
				self.barcodes[pair.n15.seq] += 1
#				self.pairs[pair.n15.seq].append(pair)
			except KeyError:
				self.barcodes[pair.n15.seq] = 1
#				self.pairs[pair.n15.seq] = [pair]
		return
		
	def part1(self):
		perc_c = str(round(100*float(self.handlefound)/self.readcount,2))+'%'
		uniq_n15s = str(len(self.barcodes.keys()))
		return 'in '+perc_c+' of the reads can the chandle be found, there are '+uniq_n15s+' uniq n15s'

	def loadclusters(self,filename):
		f = open(filename,'r',)
		self.clusters = eval(f.read())
		f.close()
		#print len(self.clusters)
		return

	def reducebarcodes(self,indata):
		""" Find most common barcodes in well ( > 10% ??), then try to place other barcodes to this cluster
		"""
		
		indata.logfile.write( '\n' )
		
		indata.minperc = 0
		maxdist = indata.bcmm
		matchfunc = hamming_distance

		percentages={}
		for bc, count in self.barcodes.iteritems():
			percentage = round(100*float(count)/self.readcount,4)
			try: percentages[percentage].append(bc)
			except KeyError:percentages[percentage] = [bc]
		highest = []
		perc = percentages.keys()
		perc.sort(reverse=True)
		#print perc
		while len(highest) < indata.seed:
			try:
				for bc in percentages[perc[0]]: highest.append(bc)
			except IndexError: pass
			perc=perc[1:]
		tempfile = open(indata.outfolder+'/seed_bcs.tempfile','w')
		for bc in highest: tempfile.write('>'+bc+'\n'+bc+'\n')
		tempfile.close()

		tempfile = open(indata.outfolder+'/raw_bcs.tempfile','w')
		for bc in self.barcodes: tempfile.write('>'+bc+'\n'+bc+'\n')
		tempfile.close()
		del percentages

		import subprocess
		from cStringIO import StringIO
		import time
		tempo = time.time()
		indata.logfile.write('starting '+' '.join(['dnaclust','--similarity',str(1-(float(indata.bcmm)/15)),'--input-file',indata.outfolder+'/raw_bcs.tempfile','-t',str(indata.cpus),'--predetermined-cluster-centers',indata.outfolder+'/seed_bcs.tempfile'])+'\n')
		dnaclust =               subprocess.Popen(['dnaclust','--similarity',str(1-(float(indata.bcmm)/15)),'--input-file',indata.outfolder+'/raw_bcs.tempfile','-t',str(indata.cpus),'--predetermined-cluster-centers',indata.outfolder+'/seed_bcs.tempfile'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		dnaclust_out, errdata = dnaclust.communicate()
		if dnaclust.returncode != 0:
			print 'dnaclust view Error code', dnaclust.returncode, errdata
			sys.exit()
		dnaclust_out = StringIO(dnaclust_out)
		indata.logfile.write('dnaclust done after '+str(round(time.time()-tempo,2))+'s, parsing result ... ')
		del dnaclust

		clusters={}
		cc=0
		for line in dnaclust_out:
			cc+=1
			clusters[cc] ={'total':0,'barcodes':{},'highest':[0,'XXXXXXXXXXXXXX']}
			line = line.rstrip().split('\t')
			for bc in line:
				clusters[cc]['barcodes'][bc]=self.barcodes[bc]
				clusters[cc]['total']+=self.barcodes[bc]
				if self.barcodes[bc] > clusters[cc]['highest'][0]: clusters[cc]['highest']=[self.barcodes[bc],bc]
		indata.logfile.write('almost done ... ')
		del dnaclust_out

		counter = 0
		reads_in_clusters={}
		for cc in clusters:
			try: reads_in_clusters[clusters[cc]['total']]+=1
			except KeyError: reads_in_clusters[clusters[cc]['total']]=1
			if int(clusters[cc]['total']) > 1:
				counter+=1
				#print cc,clusters[cc]['total'],clusters[cc]['highest'][1],clusters[cc]['highest'][0]

		indata.logfile.write('ok done now I\'ll just print and plot some info ... then done ... for real!\n\n')
		indata.outfile.write(str( cc)+' clusters whereof '+str(counter)+' has more than one read\n\n')

		temp_x=reads_in_clusters.keys()
		temp_x.sort()
		y=[];x =[]
		for i in xrange(max(temp_x)+1):
			x.append(i)
			try: y.append(reads_in_clusters[i])
			except KeyError: y.append(0)
		x=x
		y=y
		import numpy as np
		import matplotlib.pyplot as plt
		plt.figure()
		plt.axis([0,5000,0,20])
		plt.xlabel('Total Number of Reads per Barcode Cluster')
		plt.ylabel('Number of Clusters')
		#pos = np.arange(len(x))
		#width = 1.0     # gives histogram aspect to the bar diagram
		#ax = plt.axes()
		#ax.set_xticks(pos + (width / 2))
		#ax.set_xticklabels(x,rotation='horizontal')
		#plt.bar(pos, y, width, color='r')
		##plt.show()
		##plt.savefig(pp,format='pdf',bbox_inches=0)
		plt.plot(x,y)
		plt.suptitle('Diststibution', fontsize=12)
		plt.savefig('pelle.pdf')
		plt.close()
		indata.logfile.write( 'done\n')

		self.clusters = clusters
		del clusters

		return

def classify_cluser(indata,infile="temporary.cluster.files/1.reads",database="reference/4amplicons/4amplicons.fasta"):

	#database="../reference/4amplicons/4amplicons.fasta"
	
	from Bio.Blast.Applications import NcbiblastnCommandline
	from Bio.Blast import NCBIXML
	from cStringIO import StringIO
	import time
	
	#setting up blast
	database="reference/4amplicons/4amplicons.fasta"
	if indata.blastsetting == 'strict':cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, out=infile+'.'+indata.blastid+'.blastout')
	elif indata.blastsetting == 'sloppy':cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn', out=infile+'.'+indata.blastid+'.blastout')
	
	#Blasting all reads in cluster 
	blast_handle = cline.__call__()

#	blast_handle = StringIO(blast_handle[0])
#	blast_handle.seek(0)

	f = open(infile+'.'+indata.blastid+'.blastout')
	records = NCBIXML.parse(f)

	results = {'total':0}
	records_counter = 0
	#checking blast results
	for blast_record in records:
		records_counter +=1

		#print blast_record.query+'\t',
#		if indata.printblast: f.write(blast_record.query+'\n')
		readnumber = int(blast_record.query.split('_r')[1].split('_')[0])
		if readnumber == 1: r1_header = blast_record.query
		elif readnumber == 2 and records_counter%2 == 0: r2_header = blast_record.query
		else: indata.logfile.write('Warning: readnumber is funky!\n')
		
		if blast_record.alignments:
			#		if len(hsp.query) < 40 : f.write('\tTo short will ba counted as "No Hit"\n')
			if readnumber == 1:
				r1_header = blast_record.query
				r1_subj_name = blast_record.alignments[0].title.split(' ')[1]
				if len(blast_record.alignments[0].hsps[0].query) <= 40: r1_subj_name = 'No Hits'
			elif readnumber == 2 and records_counter%2 == 0:
				r2_header = blast_record.query
				r2_subj_name = blast_record.alignments[0].title.split(' ')[1]
				if len(blast_record.alignments[0].hsps[0].query) <= 40: r2_subj_name = 'No Hits'
			else: indata.logfile.write('Warning: readnumber is funky!\n')
		else:
			if readnumber == 1: r1_subj_name = 'No Hits'
			elif readnumber == 2: r2_subj_name = 'No Hits'

		if records_counter%2 == 0:
			if r1_header.split('_')[0] == r2_header.split('_')[0]:
				[junk0,junk1,rawrn,cluster_id,n15] = r1_header.split('_') #@M00275:102:000000000-A33TB:1:1101:16119:1648 _ 1:N:0:7 _ 205 _ AAGAGTCAGACTGAA
				cluster_id = int(cluster_id)
				results['total']+=1
				
				try: results[cluster_id]['total']+=1
				except KeyError: results[cluster_id]= {'total':1}
				
				if r1_subj_name == r2_subj_name:
					hit = r1_subj_name
				else:
					hit = 'Pair Disagree'
				
				try: results[cluster_id][hit]+=1
				except KeyError:results[cluster_id][hit]=1
			
			else: indata.logfile.write('WARNING: read pair headers missmatch!\n')
	f.close()
	
	import os
	os.remove(infile)
	os.remove(infile+'.'+indata.blastid+'.blastout')
	
	return results


VARIATIONINFO = {
					'c1':{'ext':'GACCATCACTTAAATCAGGTCCTCC', 'tj':'AGAGTCAAGTTATTTAAAAAATCTGGCC', 'ref':'ATTCAACAGTGATGGGACCATCACTTAAATCAGGTCCTCCKCCTCAGAGTCAAGTTATTTAAAAAATCTGGCCGATTTTGA','snp':'K'},
					'c2':{'ext':'CTCTTCAGCTCCCAGAGTCACCA', 'tj':'CCACTTACATACTTGTCCCTCTAAG', 'ref':'CCAAGGAATGGCCACTGCTCTTCAGCTCCCAGAGTCACCARTGGTTCCACTTACATACTTGTCCCTCTAAGTCATTTAGCC','snp':'R'},
					'c3':{'ext':'TGGATTTTTACATATGAGCCTTCAATG', 'tj':'CCTGGTTACTATTATTAAAGAATTTCTC', 'ref':'AGAAAGGAAGGCATGGATTTTTACATATGAGCCTTCAATGRTGTTGCCTGGTTACTATTATTAAAGAATTTCTCAGCTTTT','snp':'R'},
					'c4':{'ext':'GGCCAACAAGATCAGTCTGTTCTCT', 'tj':'CTGCTTATCTGACAACCCTCAGATC', 'ref':'CCTGTGCTGGCAGAAGGCCAACAAGATCAGTCTGTTCTCTRGTAGACTGCTTATCTGACAACCCTCAGATCCCAGCAGCTG','snp':'R'},
					'c5':{'ext':'ACCCTAGATTGTATTGTAGGAGGCAT', 'tj':'ATGGATGGCTGCTGGAAACCCCT', 'ref':'ATGGATCATACTTTACCCTAGATTGTATTGTAGGAGGCATYGTGGATGGATGGCTGCTGGAAACCCCTTGCCATAGCCAGC','snp':'Y'},
					'c6':{'ext':'CTTTCATTCTGCTCAAGTCTTCGCC', 'tj':'CCGATCCCTATCTACTTTCTCTCC', 'ref':'TTTCCCCTAGATTTCCTTTCATTCTGCTCAAGTCTTCGCCKGTGTCCGATCCCTATCTACTTTCTCTCCTCTTGTAGGCAA','snp':'K'},
					'c7':{'ext':'CCGGCCTCTGGCGTTCTACTCA', 'tj':'ACCTCTTTGTCCTGCAGCAGTTTTC', 'ref':'TGGTTTGTGCGAGTCGTCCCGGCCTCTGGCGTTCTACTCAYGTGCTGACCTCTTTGTCCTGCAGCAGTTTTCCAGCTTGCA','snp':'Y'},
					'c8':{'ext':'TTGCACAAGGTAAGTTTATTCTAGCTT', 'tj':'ATTTTGTGGTGAAAGTGCCTAAATTTG', 'ref':'GAGCCATTTACAATTGCACAAGGTAAGTTTATTCTAGCTTYCTTTCATTTTGTGGTGAAAGTGCCTAAATTTGTGAAAATG','snp':'Y'},
					'c9':{'ext':'CCTGAATATCAGGTAGGAATGTTTGT', 'tj':'ATCGCGCTCCCTGAGGATACTCA', 'ref':'CGGGAGTGGGCTTCCCTGAATATCAGGTAGGAATGTTTGTYCCTCATCGCGCTCCCTGAGGATACTCATGCCTGTGGTGGG','snp':'Y'},
					'c10':{'ext':'ACAAAGGCAGCCTTTCATACCTTAG', 'tj':'CCAGCTTCCCAGAAGTTGTGCATA', 'ref':'GGTCTGGGTTGAGGTACAAAGGCAGCCTTTCATACCTTAGKATAAGTCCAGCTTCCCAGAAGTTGTGCATATTTCCACTGA','snp':'K'},
					'c11':{'ext':'CAAACCTTACAGCATCGATAAGTTTC', 'tj':'AAAGAAACCAATGTTTGCAAAGTAGATA', 'ref':'AAATTAAATAACCCCAAACCTTACAGCATCGATAAGTTTCYGCCTAAAGAAACCAATGTTTGCAAAGTAGATAGGAGATGG','snp':'Y'},
					'c12':{'ext':'GCGGTGGGCGCGCAGTGCGTTCTC', 'tj':'AGGGCAGCTGTTCCGCCTGCGAT', 'ref':'GTACGCGCGTGGCGTGGCGGTGGGCGCGCAGTGCGTTCTCKGTGTGGAGGGCAGCTGTTCCGCCTGCGATGATTTATACTC','snp':'K'},
					'c13':{'ext':'CCAGGGTAAGGGAGAGGGCCACA', 'tj':'CCTTGGAGTTTGAAAGACAAAGGGA', 'ref':'GTCTCCAAGAGGCCTTCCCAGGGTAAGGGAGAGGGCCACARGGGAACCTTGGAGTTTGAAAGACAAAGGGAACACATGACA','snp':'R'},
					'c14':{'ext':'ATATGTTCCCTTCCTGAATTTGTCTAA', 'tj':'CTGAGTTGGGGGGTGGTTGCTG', 'ref':'CCTGCCTTTCAGGATATGTTCCCTTCCTGAATTTGTCTAARAAGGCTGAGTTGGGGGGTGGTTGCTGATTTTTATAACATA','snp':'R'},
					'c15':{'ext':'AAAAGCCTAAGAACTCTAAGCATTGC', 'tj':'CGCCAAGTTTGAAGGAACTCGAATT', 'ref':'TCTGAGACCAAAGCAAAAGCCTAAGAACTCTAAGCATTGCRGAGACGCCAAGTTTGAAGGAACTCGAATTCCACACCTGGT','snp':'R'},
					'c16':{'ext':'GAACACTGCCCTCGCACGGCCCC', 'tj':'CCCCCACATTTGTGCCGACACTGTT', 'ref':'GACTGCACCCAGGAGGGGAACACTGCCCTCGCACGGCCCCRAAGGGCCCCCACATTTGTGCCGACACTGTTCTCCCTTCGG','snp':'R'},
					'c17':{'ext':'CAGCTAGAGATAGGGGAAGTGTATA', 'tj':'CCCAATGCTTCCAAAGAGGAAAGG', 'ref':'TTCTACAAGATCTTCCAGCTAGAGATAGGGGAAGTGTATARGAATCCCAATGCTTCCAAAGAGGAAAGGAAAAGGTGGCAG','snp':'R'},
					'c18':{'ext':'CTTCCACTCTCAAAGGGCTTCTGAT', 'tj':'ACATTTGAATCTAATGGATCAGTATCAT', 'ref':'TTGGAGATTTTGTCACTTCCACTCTCAAAGGGCTTCTGATKTGCTACATTTGAATCTAATGGATCAGTATCATTTGGTTCC','snp':'K'},
					'c19':{'ext':'TGGATTCTGAATGTGCTTAATTTAAAAG', 'tj':'ATTTTTACAAAGTGATCGAAAGTTTTATC', 'ref':'TTGTCTTTCCCATGGATTCTGAATGTGCTTAATTTAAAAGSCTTTGATTTTTACAAAGTGATCGAAAGTTTTATCAAAGCA','snp':'S'},
					'c20':{'ext':'CTAATTCTGTTTCATTTCTATAGGCGA', 'tj':'CATGATCAGGAAGGCCGGGTGAT', 'ref':'ATAGTTTTTTATGCTAATTCTGTTTCATTTCTATAGGCGAKGAGGAGCATGATCAGGAAGGCCGGGTGATTGTGGCTGAAT','snp':'K'},
					'c21':{'ext':'CTTCCCTCTCCTCAGCGCCTAAGAG', 'tj':'CCCAGTGCGGGTGAGGAGTCGCGA', 'ref':'AGCAAGGTCCCAGGTCTTCCCTCTCCTCAGCGCCTAAGAGRGAGGCCCAGTGCGGGTGAGGAGTCGCGAGGAAGAGGCGGA','snp':'R'},
					'c22':{'ext':'GTGAGAATCCCTGAGCCCAGGCGGT', 'tj':'CCAGGAGGTGTCTGGACTGGCTGG', 'ref':'GGACCAGGCTTGGTTGTGAGAATCCCTGAGCCCAGGCGGTWGATGCCAGGAGGTGTCTGGACTGGCTGGGCCATGCCTGGG','snp':'W'},
					'c23':{'ext':'GGCACCTTTCCTCTAATCCAGCAAA', 'tj':'CCCTGCACACCAGAGACAAGCAG', 'ref':'TTGGAAACAGACGGCGGCACCTTTCCTCTAATCCAGCAAARTGATTCCCTGCACACCAGAGACAAGCAGAGTAACAGGATC','snp':'R'},
					'c24':{'ext':'CATTTCAAGGCAGTTTTTAAAGAAGTC', 'tj':'CAGACGGCGGCACCTTTCCTCT', 'ref':'GATGAACTAATTTCATTTCAAGGCAGTTTTTAAAGAAGTCWTGGAAACAGACGGCGGCACCTTTCCTCTAATCCAGCAAAG','snp':'W'},
					'c25':{'ext':'CTAGGACTCCTGCTAAGCTCTCCT', 'tj':'ACGCTTTTGCTAAAAACAGCAGAACT', 'ref':'TGTATGGGTGAAAGGGCTAGGACTCCTGCTAAGCTCTCCTYTCTGGACGCTTTTGCTAAAAACAGCAGAACTTTCCTTAAT','snp':'Y'},
					'c26':{'ext':'GTTTCAAAGCGCCAGTCATTTGCTC', 'tj':'AAATCCAGGAAATGCAGAAGAGGAAT', 'ref':'CAGAATACATTCAAGGTTTCAAAGCGCCAGTCATTTGCTCHGTTTTCAAATCCAGGAAATGCAGAAGAGGAATGTGCAACA','snp':'H'},
					'c27':{'ext':'GGAGCTGCTGGTGCAGGGGCCACG', 'tj':'AGCAGCCTCTGGCATTCTGGGAGCT', 'ref':'AGGGGCCGCCGGTGTAGGAGCTGCTGGTGCAGGGGCCACGSGGGGAGCAGCCTCTGGCATTCTGGGAGCTTCATCTGGACC','snp':'S'},
					'c28':{'ext':'GTGCCTGGACATTTGCCTTGCTGGA', 'tj':'ACTGGGGATGTGGGAGGGAGCAGA', 'ref':'TTTAGGGCCAGGCAGGTGCCTGGACATTTGCCTTGCTGGAYGGGGACTGGGGATGTGGGAGGGAGCAGACAGGAGGAATCA','snp':'Y'},
					'c29':{'ext':'GTGTGAGAGAAAAGAATGGAATAAGC', 'tj':'CTGCCATGCTCAGAGAATCCTAGA', 'ref':'AATGCTGATCCCCTGTGTGAGAGAAAAGAATGGAATAAGCRGAAACTGCCATGCTCAGAGAATCCTAGAGATACTGAAGAT','snp':'R'},
					'c30':{'ext':'TTCCAATGTCCCCCGGTTGGAGTC', 'tj':'ATCTGCAATGCTCCAGAGGGCAAGA', 'ref':'GCTGAGCATCGGTAAATTCCAATGTCCCCCGGTTGGAGTCSGCTGATCTGCAATGCTCCAGAGGGCAAGACCACCACTCGG','snp':'S'},
					'c31':{'ext':'GGTGGAGCCACACGAAGCGGTG', 'tj':'AATTAGTGGTCGGATTTCCAAAGACA', 'ref':'ACAGGCACAGAACATCCAGGTGGAGCCACACGAAGCGGTGYTTGGCAATTAGTGGTCGGATTTCCAAAGACAGGAGTTTTG','snp':'Y'},
					'c32':{'ext':'AGGCAGGGCAGGTCCAGGCCCTGC', 'tj':'CTGAGTGCCGGGGACGTACAGTGGC', 'ref':'GGCCACACCTCAACACAGGCAGGGCAGGTCCAGGCCCTGCRTGTGCTGAGTGCCGGGGACGTACAGTGGCCCCCACAGCTT','snp':'R'},
					'c33':{'ext':'GAGGTATCAGAGGTAATAAATATTCTAT', 'tj':'CAATAAGGTTTCTCAAGGGGCTGG', 'ref':'GAAATGGGGGTTGAGGTATCAGAGGTAATAAATATTCTATRAGAGAGGTACAATAAGGTTTCTCAAGGGGCTGGGTCAGCT','snp':'R'},
					'c34':{'ext':'CTCCCCCACCCCGCTTCACTCTG', 'tj':'ACTGGAGCCCCACAGCCCCCAC', 'ref':'GCCCCTCCACCTCTACACTCCCCCACCCCGCTTCACTCTGYTCGTTCTGACTGGAGCCCCACAGCCCCCACTCCCTCCAGT','snp':'Y'},
					'c35':{'ext':'TCCCCGCCTGGGGAAGGAGAGGCC', 'tj':'AGCGGCGGCACGGCCGCCTTCG', 'ref':'CTCGCTCCAGGCCGCATCCCCGCCTGGGGAAGGAGAGGCCYCGGGGGCAGCGGCGGCACGGCCGCCTTCGCGTGAGGCCCT','snp':'Y'},
					'c36':{'ext':'GAATTCATCGGACATGTTACTGTTTTT', 'tj':'ATGTACCACCAACTTTACGTTTGCAT', 'ref':'GTGAAAGACTAATGAATTCATCGGACATGTTACTGTTTTTSCTCCCTGATGTACCACCAACTTTACGTTTGCATGAAGGTG','snp':'S'},
					'c37':{'ext':'TAAAGTTTGGTCCCTTTCGCTCCG', 'tj':'AAAATGCTCTAACGGCAGGAGGTC', 'ref':'AAGGTGAACTGCTGACTAAAGTTTGGTCCCTTTCGCTCCGKCTCCTTGCGAAAATGCTCTAACGGCAGGAGGTCACGCGAG','snp':'K'},
					'c38':{'ext':'GTCAGGATTGCTGGTGTTATATTCTG', 'tj':'CAAAACAAGGAGAGAAGGAGTTGGA', 'ref':'ATGACTAGCAACAAGTCAGGATTGCTGGTGTTATATTCTGWTGGGGAAACAAAACAAGGAGAGAAGGAGTTGGATGTGGTA','snp':'W'},
					'c39':{'ext':'GTTTGATTCACCTTTGAGCACTCTTA', 'tj':'CAGAAGAGGCGTAACTGGTTACAAA', 'ref':'GGAATGATATACAAGTTTGATTCACCTTTGAGCACTCTTABAGAGTTTACAGAAGAGGCGTAACTGGTTACAAAACAAAAG','snp':'B'},
					'c40':{'ext':'ATTCCTCTGGTATTTTCTAAAACAGAAA', 'tj':'CATCAACTTTGTGTATCTATTCAGAAAA', 'ref':'TAAATACATTAAATTCCTCTGGTATTTTCTAAAACAGAAARTTTAATACATCAACTTTGTGTATCTATTCAGAAAATAGTA','snp':'R'},
					'c41':{'ext':'CTAGCTGGCCTTATTCTACTGACTT', 'tj':'AAGGCTGCTCCAGCCAGCACTG', 'ref':'CCAATTCAAAAACTCCTAGCTGGCCTTATTCTACTGACTTSGTGCGTGGAAGGCTGCTCCAGCCAGCACTGGTCCTATGGA','snp':'S'},
					'c42':{'ext':'GTACTGGACGTTGATGCCACTGAA', 'tj':'AGTTTTCCACACAGACTTTCCCGAT', 'ref':'CAGAGGAATGCAGCTTGTACTGGACGTTGATGCCACTGAAYGCCTGGCAGTTTTCCACACAGACTTTCCCGATCTAGGGGG','snp':'Y'},
					'c43':{'ext':'AGAAAGTGAGGACAGACAGAAGGC', 'tj':'CACCATGTCACCCTCAGGCGTG', 'ref':'GGAGGGACAGGGCAGGAGAAAGTGAGGACAGACAGAAGGCRGTGAGGGCACCATGTCACCCTCAGGCGTGCCAGGGAAGCT','snp':'R'},
					'c44':{'ext':'CATCTGCAGTCCCTGGGGGCCA', 'tj':'AGCTGAGCAGTGCCTTCCAGAGC', 'ref':'TCCGCCCCCAGGGCCAGGCATCTGCAGTCCCTGGGGGCCAYGCTGCGTGAGCTGAGCAGTGCCTTCCAGAGCCCTGTGCTG','snp':'Y'},
					'c45':{'ext':'GGATGGCGAGGCAGCAGGGGCA', 'tj':'AGCTTGGGATGCCCTAGGAAGGG', 'ref':'TACAGCCAGCTTTGGAGCGGATGGCGAGGCAGCAGGGGCAKCCCCTGGGAGCTTGGGATGCCCTAGGAAGGGCGAGGACTC','snp':'K'},
					'c46':{'ext':'TGCGTCGGCGGCTGCCCTCCC', 'tj':'CCTCACACGCCAACCCTGCTCCT', 'ref':'GTGCTGGACTGTCACCGCATGCGTCGGCGGCTGCCCTCCCRGAGGTAAGGCCTCACACGCCAACCCTGCTCCTTATCCTGT','snp':'R'},
					'c47':{'ext':'GATTCTCACGACTCTACCATGTCCA', 'tj':'CAGAGATCGTTCCTATACATTTCTGT', 'ref':'CTCTGGAAAACTCAGGATTCTCACGACTCTACCATGTCCARTGGAGTTCAGAGATCGTTCCTATACATTTCTGTTCATCTT','snp':'R'},
					'c48':{'ext':'GGAGGGCTGGGAGTCCGGAATG', 'tj':'AGTTGTGGGGCCCAGACTCCTAG', 'ref':'TTAAGATCCTGGCAGAGAGGAGGGCTGGGAGTCCGGAATGYTTGCGCTGAGTTGTGGGGCCCAGACTCCTAGTTCTTAGGG','snp':'Y'},
					'c49':{'ext':'CCAGTCCTGTCATGTCAGGGTTG', 'tj':'CTCAGGGGAAAATATGAGGGGCAC', 'ref':'TGCAAGTCTGGGATGGACCAGTCCTGTCATGTCAGGGTTGRGAATGAAGGCTCAGGGGAAAATATGAGGGGCACTGGAGCC','snp':'R'},
					'c50':{'ext':'AGCGGATGGTGGATTTCGCTGGC', 'tj':'CAAGGTGTGCATGCCTGACCCGT', 'ref':'CTGTGCCGCCATCACCCAGCGGATGGTGGATTTCGCTGGCRTGAAGGACAAGGTGTGCATGCCTGACCCGTTGTCAGACCT','snp':'R'},
					'c51':{'ext':'TTCCTGGCTCCATTTTCTCCACCAG', 'tj':'AATACATGAGCCAGACCCGCAGCTT', 'ref':'GATCTTAGCAATGACTTCCTGGCTCCATTTTCTCCACCAGYCCGAATACATGAGCCAGACCCGCAGCTTTCTCTTCCACAG','snp':'Y'},
					'c52':{'ext':'CCTGCCCAGGCTGAGTGCGGA', 'tj':'AGAGCTCCACCGAGGGGTGGG', 'ref':'GGCAGAGGGGCCTGTGTCTCCTGCCCAGGCTGAGTGCGGAYGGCTCGGTCTCAGAGCTCCACCGAGGGGTGGGTAGGTAAC','snp':'Y'},
					'c53':{'ext':'GCGAGGATGTCGTCCACCTCTG', 'tj':'ACCTTGGACACCGGGACAAAGGT', 'ref':'GTGGTCGGGCGGACTGCCGCGAGGATGTCGTCCACCTCTGYCTGCCCGCTCACCTTGGACACCGGGACAAAGGTGACCGCT','snp':'Y'},
					'c54':{'ext':'TGATGCTTTCGAAGTTTCAGTTGAAC', 'tj':'AAAAAGTTTCTTTAAATGTAAGAGCAGG', 'ref':'GTCTCTGCAAATAATGATGCTTTCGAAGTTTCAGTTGAACYGTCCCTCGCGAAAAAGTTTCTTTAAATGTAAGAGCAGGTC','snp':'Y'},
					'c55':{'ext':'TCGGACCGCATGGGTCGGACAGGT', 'tj':'ATCTAAGTGGGACAAGTTGACCCAGG', 'ref':'CCCAGGTCGGTGCGGTTCGGACCGCATGGGTCGGACAGGTYGCCATCTAAGTGGGACAAGTTGACCCAGGAACCGGGGCTG','snp':'Y'},
					'c56':{'ext':'GCAAGCGGTGAGTTTTCAGATGGGCA', 'tj':'AGAAGTGCATGTCCAGTGCCTGGC', 'ref':'CCGCTGGCCTGCGAGCAAGCGGTGAGTTTTCAGATGGGCAYGGGAGAAGTGCATGTCCAGTGCCTGGCTGCACCCGGGCCC','snp':'Y'},
					'c57':{'ext':'TTGATATTCTTCATTAGCTTGCCTGAT', 'tj':'CCATAGTAGAAACACATCAGTACTGA', 'ref':'AGTTCGCCAGTATTTGATATTCTTCATTAGCTTGCCTGATRAGAATAATGGGCCATAGTAGAAACACATCAGTACTGAAGC','snp':'R'},
					'c58':{'ext':'TGTCTGCTCCAAATATAGCAATGAAG', 'tj':'ACCCTCAACCTCTGCTGTCAAAGT', 'ref':'ACCAGTGCCAAGTGTGTCTGCTCCAAATATAGCAATGAAGYCCTCCTTCCCACCCTCAACCTCTGCTGTCAAAGTCAACCT','snp':'Y'},
					'c59':{'ext':'CTGGGGGGTGGCTCAGGGGAGGGTG', 'tj':'CAGGCGGTTCAAGCCGTTGGCTGGA', 'ref':'GCGAAGCGGACCACGCTGGGGGGTGGCTCAGGGGAGGGTGWGGGCAGGCGGTTCAAGCCGTTGGCTGGAGACACCTATTTA','snp':'W'},
					'c60':{'ext':'GTTCTGGCAGGCATTTGGCATCAGC', 'tj':'ACTCGTCCACATCCTCGGTACAGTA', 'ref':'GTGGCAGGTCCCGCCGTTCTGGCAGGCATTTGGCATCAGCYGGCACTCGTCCACATCCTCGGTACAGTACTGACCTGCAGG','snp':'Y'},
					'c61':{'ext':'GAGGGCTAAAAAGGTCCTGTAAGAAG', 'tj':'ATAGCACCCGTGGTGCATGGTATG', 'ref':'TTGGAAGGTCTTTTGAGGGCTAAAAAGGTCCTGTAAGAAGYGGCATAGCACCCGTGGTGCATGGTATGTGGGTGACCCTGG','snp':'Y'},
					'c62':{'ext':'CCTTCCCGGTCAGCTACTCCTCTTCC', 'tj':'CCCGCCGGCCCTCGCTGGACTCCA', 'ref':'GAGGCCCGCCCAGGCCTTCCCGGTCAGCTACTCCTCTTCCRGTGCCCGCCGGCCCTCGCTGGACTCCATGGAGAACCAGGT','snp':'R'},
					'c63':{'ext':'GTCGTTTCCAAGAGAGATCCTTTCTT', 'tj':'ATACACATCTGTACTGGGAGCTTGG', 'ref':'GAACTTTTTAACATGTCGTTTCCAAGAGAGATCCTTTCTTYGCCATACACATCTGTACTGGGAGCTTGGTCCCTCCTCTGG','snp':'Y'},
					'c64':{'ext':'TGTTTTCACTGTCTTGCTTCTGGTAA', 'tj':'ACCAGGTGTGGAATTCGAGTTCCT', 'ref':'CTTGGCCTCACTCTTGTTTTCACTGTCTTGCTTCTGGTAASGCCTTTTCTTCACCAGGTGTGGAATTCGAGTTCCTTCAAA','snp':'S'},
					'c65':{'ext':'GCACAGAAAGACCTGTGTGCTGC', 'tj':'ACTCACATTCAGTACATTTGGGTTCC', 'ref':'GTCTTGCAGGTAGGAATGCACAGAAAGACCTGTGTGCTGCKGGTTGGTGTGGACTCACATTCAGTACATTTGGGTTCCTGG','snp':'K'},
					'c66':{'ext':'CAGGCCATACTCTCCTTTACCATACTA', 'tj':'ACTGGACTCAGGCTGGAGGCAGA', 'ref':'CTTCCCTCTCCCTCAGGCCATACTCTCCTTTACCATACTAYTGGACTGGACTCAGGCTGGAGGCAGATAGACACAGTATAG','snp':'Y'},
					'c67':{'ext':'GCTGGACAGGAAGGGAGAATTCTGA', 'tj':'AACATGCAGCGAAGTATCATGTGAG', 'ref':'CAGACGAAGAAAAGTGCTGGACAGGAAGGGAGAATTCTGABGCCAACATGCAGCGAAGTATCATGTGAGCACCCTCCCTTG','snp':'B'},
					'c68':{'ext':'GGAGCTGCAGAGGTGTGGGCCCCTG', 'tj':'CCCAGAAGTCCAGCCACTGGGCTCC', 'ref':'ACTGGACAGTAACATGGAGCTGCAGAGGTGTGGGCCCCTGRGGACCCAGAAGTCCAGCCACTGGGCTCCTTCAGGACACAG','snp':'R'},
					'c69':{'ext':'ATATTTCAAGAGCTCCCATGTTCAGT', 'tj':'CCTGGAGAAGGAAAGCTCTTTTTGT', 'ref':'ACTTTTAATACTACATATTTCAAGAGCTCCCATGTTCAGTRTTCCTGGAGAAGGAAAGCTCTTTTTGTTTCTAAGTATGAG','snp':'R'},
					'c70':{'ext':'TTTTAATAATCGACTTTTTAAATGTGATCA', 'tj':'CTATAATATTGTACAGTTATTATAGGGC', 'ref':'AACCCTTGTTTTTTAATAATCGACTTTTTAAATGTGATCARGACTATAATATTGTACAGTTATTATAGGGCTTTTGGGGAA','snp':'R'},
					'c71':{'ext':'TAGGTAAAAATTATTAAGTGAAATTATTCAT', 'tj':'CTTTCAGTTTCATCTTTCTCCTGGG', 'ref':'AGTGGAGATTAGGTAAAAATTATTAAGTGAAATTATTCATRTGCTTTCAGTTTCATCTTTCTCCTGGGAAAGCCCCAAAAC','snp':'R'},
					'c72':{'ext':'GGCGGGTCTGGGTGCGGCCTGCCGCA', 'tj':'CTCGTGGTTCGGAGGCCCACGTGGC', 'ref':'TCCCGTGTGGGAGGGGCGGGTCTGGGTGCGGCCTGCCGCAWGACTCGTGGTTCGGAGGCCCACGTGGCCGGGGCGGGGACT','snp':'W'},
					'c73':{'ext':'TTTTAGTTTGGTGATAGAACAGCTCTT', 'tj':'CATGAGTACGTATCTTTTCTTTTAAAAG', 'ref':'TTTGTATTATAGCTTTTAGTTTGGTGATAGAACAGCTCTTRTTCATGAGTACGTATCTTTTCTTTTAAAAGAAAAAAGTAA','snp':'R'},
					'c74':{'ext':'GGCCAGGACCAGTTGGGCAACAAAAT', 'tj':'CTGGGCACTTGCTGCCAGTACTGGG', 'ref':'GTGAGTGAAGACCAGGCCAGGACCAGTTGGGCAACAAAATRTACTGGGCACTTGCTGCCAGTACTGGGAGCTGCCAGCACG','snp':'R'},
					'c75':{'ext':'ACCAGAGAGGCTCTGAGAAACCTC', 'tj':'AAACTTAGATCATCAGTCACCGAAGGT', 'ref':'ATCCCCGATTGAAAGAACCAGAGAGGCTCTGAGAAACCTCSGGAAACTTAGATCATCAGTCACCGAAGGTCCTACAGGGCC','snp':'S'},
					'f1':{'ext':'CATGGTCTGGAGTCTCGGAGTCC', 'tj':'TGGCCACGATGACCAGCAGGTTG', 'ref':'GTCACGAACACGTTGGTCATGGTCTGGAGTCTCGGAGTCCRGGCGATGGCCACGATGACCAGCAGGTTGCCTCCCACGGTG','snp':'R'},
					'f2':{'ext':'ACTCATTTTCCTGGTGTTTCAAATATTT', 'tj':'TTTTAGGTAGTAAAGGCTTTTGTCGTT', 'ref':'TACAAATGCATAACTCATTTTCCTGGTGTTTCAAATATTTRTTTTAGGTAGTAAAGGCTTTTGTCGTTCTAAATCCTGATT','snp':'R'},
					'f3':{'ext':'TTGAGTCGTGGTTTCCTGGTCATG', 'tj':'TAATAGCAGTAGAACAGCTCCCAGC', 'ref':'GCAGGACTCCGGGCCCTTGAGTCGTGGTTTCCTGGTCATGMCCGGGCAGAGCTAATAGCAGTAGAACAGCTCCCAGCAACA','snp':'M'},
					'f4':{'ext':'TGTCTAGGCCTTAGTTAATAATGAATG', 'tj':'TTCATATAGTTTATATCAGTGTAGGAGG', 'ref':'TTCATCACAGACCTCCTACACTGATATAAACTATATGAAGKCATTCATTATTAACTAAGGCCTAGACACAGGGAGAAAGCA','snp':'K'},
					'f5':{'ext':'GCGCCTTCTTGCTGGCACCCAAT', 'tj':'CCATGCGCCGGACCACGACGT', 'ref':'GCAACCCGGGAACGGCAGCGCCTTCTTGCTGGCACCCAATRGAAGCCATGCGCCGGACCACGACGTCACGCAGGAAAGGGA','snp':'R'},
					'f6':{'ext':'TGCATCTGAATGGGCAAGAAGGAG', 'tj':'CCTGACACAATCCACACCATCAGA', 'ref':'GGTGGCCCGGTACCAGTGCATCTGAATGGGCAAGAAGGAGRTAAGGCCTGACACAATCCACACCATCAGAATGATCACCCG','snp':'R'},
					'f7':{'ext':'CCATTGGGAGAGAAGGAAGAGTGG', 'tj':'CAGAAAAGACTTCTTGGGCTGGGT', 'ref':'GAGGGAACAACTAATGCCATTGGGAGAGAAGGAAGAGTGGRTGATCAGAAAAGACTTCTTGGGCTGGGTGTGGTGGCTCAC','snp':'R'},
					'f8':{'ext':'AAGGAGAAAGAGAGCCAGGACAAG', 'tj':'CTCTCTCCCTCCCTGAGCTGGA', 'ref':'CTTCTTCAGCACCTTCAAGGAGAAAGAGAGCCAGGACAAGWCTCTCTCCCTCCCTGAGCTGGAGCAACAGCAGGAACAGCA','snp':'W'},
					'f9':{'ext':'AGAAAGAGCTGGGACCCTGGGAA', 'tj':'TGGCCTCCAGGTAGTCTCAGGAG', 'ref':'TGAGCTCAGGGGCCTCTAGAAAGAGCTGGGACCCTGGGAASCCCTGGCCTCCAGGTAGTCTCAGGAGAGCTACTCGGGGTC','snp':'S'},
					'f10':{'ext':'CAACTAATCCAATGTGTGCCTTCACT', 'tj':'TTCAATTCATTGAGAAAACCTGCCCA', 'ref':'AAATGTTTCTGTGACAACTAATCCAATGTGTGCCTTCACTSAAGGAATTCAATTCATTGAGAAAACCTGCCCAGGTGCAAC','snp':'S'},
					'f11':{'ext':'GCGGATGGTGGATTTCGCTGGC', 'tj':'CAAGGTGTGCATGCCTGACCCG', 'ref':'CTGTGCCGCCATCACCCAGCGGATGGTGGATTTCGCTGGCRTGAAGGACAAGGTGTGCATGCCTGACCCGTTGTCAGACCT','snp':'R'},
					'f12':{'ext':'ACCCCAGGCCTTTTGCCAAACCA', 'tj':'TAGAAGGCCACAGAGGAAAGGGG', 'ref':'AGATGTAGTCGGTGTACACCCCAGGCCTTTTGCCAAACCARCACCTAGAAGGCCACAGAGGAAAGGGGAGGAAGGGTCATC','snp':'R'},
					'f13':{'ext':'CAACCTCAGCAGCACCTCCCAC', 'tj':'CCCGGAGGTGGACCCGGGCT', 'ref':'CCGGTGTGTGCTGCTCTCCAACCTCAGCAGCACCTCCCACRTCCCGGAGGTGGACCCGGGCTCCGCAGAACTGCAGAAGGT','snp':'R'},
					'f14':{'ext':'CTTTTAGGGATTGGGAAAAGTTTCAAAA', 'tj':'CTTCTTGCTGAAAGCGAGTTTTTACT', 'ref':'GAGGCCAAGTCCCTTTTAGGGATTGGGAAAAGTTTCAAAARTCTTCTTGCTGAAAGCGAGTTTTTACTTTCCTCGGGCCAC','snp':'R'},
					'f15':{'ext':'GAGGTACCGACCCGGACTTCCA', 'tj':'CATCCTCTGCCCCTTCTGCTTGA', 'ref':'CAGAAGCATAAAACTGCAGAGGTACCGACCCGGACTTCCAKTTCATCCTCTGCCCCTTCTGCTTGACCTCCATCTTCCGGG','snp':'K'},
					'f16':{'ext':'CACCAGCCAAGCTGCAGCCCC', 'tj':'TCTAGCAGGCTGGCTCCGCCC', 'ref':'AGAAAGGAGTCGAAGAAGCCACCAGCCAAGCTGCAGCCCCRAGCTCTAGCAGGCTGGCTCCGCCCGGAAGATGGAGGTCAA','snp':'R'},
					'f17':{'ext':'CCTGTTTGAGAAGTGGAAGCTCTAC', 'tj':'CCAGTGTCACCACAACCTGAGCC', 'ref':'TCAGGTGATGGACTTCCTGTTTGAGAAGTGGAAGCTCTACRGTGACCAGTGTCACCACAACCTGAGCCTGCTGCCCCCTCC','snp':'R'},
					'f18':{'ext':'AGACTGGGCCCTGCACCTCCC', 'tj':'TGCTAGCATTTGCAGGCCTACCC', 'ref':'CACTGGGGTCGAGATGGGCAGACTGGGCCCTGCACCTCCCRGGGCTGCTAGCATTTGCAGGCCTACCCGGGCAGTGCCCAG','snp':'R'},
					'f19':{'ext':'GCCACCGCCCAGGGCCCGG', 'tj':'TCATCGTTGTCCTCGGACACCGG', 'ref':'CTGCCCCCCGCCTCGGCCGTTGCCACCGCCCAGGGCCCGGRCGCCGCCTCATCGTTGTCCTCGGACACCGGGGACAATGGT','snp':'R'},
					'f20':{'ext':'CTTTTGCTGAGCCCGTTTGTGTTTT', 'tj':'CATCCACAGGGTCCAGGCTGCT', 'ref':'GGCAGTGAGAAGCATCTTTTGCTGAGCCCGTTTGTGTTTTRTCATCCACAGGGTCCAGGCTGCTGGTGCCACCCAGCCAGC','snp':'R'},
					'f21':{'ext':'CGTTATCTTCTCAAACTCTTTGGCC', 'tj':'CCATTCTGGAGGAGTTTGTTTTGTTC', 'ref':'GTGCCTCCAGCTACACGTTATCTTCTCAAACTCTTTGGCCRTCCATTCTGGAGGAGTTTGTTTTGTTCCTCTCTGCAGCAC','snp':'R'},
					'f22':{'ext':'GAGATGATCCCTTGCTGACTTACC', 'tj':'CCCACCAAAGTTCACCCTGAAGG', 'ref':'ATCAAGCGCCAGAATGGAGATGATCCCTTGCTGACTTACCRGTTCCCACCAAAGTTCACCCTGAAGGCTGGGCAGGTGGTG','snp':'R'},
					'f23':{'ext':'GAGAATAATTAATGGAGATACTATGAAAAA', 'tj':'TGTCACTTTACTTTGGAAGGTATTCC', 'ref':'CTGAATTTTGGAGAATAATTAATGGAGATACTATGAAAAASGAGAAAAATGTCACTTTACTTTGGAAGGTATTCCCAATTT','snp':'S'},
					'f24':{'ext':'AAGTGCTTCTCTAGCAAAAATTCTCC', 'tj':'CTTTCTCTAATAGCTCATGGGAGATA', 'ref':'ATAGTTCAGTCACCAAGTGCTTCTCTAGCAAAAATTCTCCRTTGAAGGATTCTTTCTCTAATAGCTCATGGGAGATAGAGG','snp':'R'},
					'f25':{'ext':'AATCACATCTGGTGGAGTAATTTTCC', 'tj':'CACCTCTAATGTCAGTTCAGCCCA', 'ref':'CTTATGTGTTTGAAAATCACATCTGGTGGAGTAATTTTCCRGTCACCTCTAATGTCAGTTCAGCCCATAAATATGGGTAAG','snp':'R'},
					'f26':{'ext':'CTTATGTGCAGACAACATTGAAGGAA', 'tj':'TTTGTTTCAACAGTAAATTCTTTAGTTTTT', 'ref':'GATAGAAACTGCTCCTTATGTGCAGACAACATTGAAGGAARGACATTTGTTTCAACAGTAAATTCTTTAGTTTTTCAACAA','snp':'R'},
					'f27':{'ext':'GCCAGAACATCCAAGGCCTCCC', 'tj':'TGCGATGAGGTGCACATAGCTGG', 'ref':'AGGAAGGCCACTGTCTTGGCCAGAACATCCAAGGCCTCCCRGCAGGTGCGATGAGGTGCACATAGCTGGACGGGCCGGTGC','snp':'R'},
					'f28':{'ext':'CGTTGTGTGACATCCCGACAGAAA', 'tj':'TTGTGAAATCGACAATACTGGCATTTG', 'ref':'GGGCCACCTTACCTACCGTTGTGTGACATCCCGACAGAAARGCACTTGTGAAATCGACAATACTGGCATTTGTTTCTGTTC','snp':'R'},
					'f29':{'ext':'AAAGCATGGTACAAAATGAACTATGGT', 'tj':'CAAGAAGGAGAGGTGGCTAGACC', 'ref':'AAAGCTCATTCTGAAAGCATGGTACAAAATGAACTATGGTRAGTCAAGAAGGAGAGGTGGCTAGACCCAGATGCAAGAGGA','snp':'R'},
					'f30':{'ext':'CTGACCAGTCTCTGGACATGAAGA', 'tj':'TGTTGGAAGTGGGAAGGCTCCCT', 'ref':'CCTGGGGCACCACCTCCTGACCAGTCTCTGGACATGAAGASGGAGGCCCTGTTGGAAGTGGGAAGGCTCCCTTCCTCCTCC','snp':'S'},
					'f31':{'ext':'TGCCAGGGATCAGTGAGGTCTCT', 'tj':'TTGGGGAGCTTGCCCAAGTGGG', 'ref':'GGTTGGAGTCAGGTCTGTGCCAGGGATCAGTGAGGTCTCTRAGACCCTTGGGGAGCTTGCCCAAGTGGGGGGTCCTCACTT','snp':'R'},
					'f32':{'ext':'GCTCCTTTGCCCCCATTAGCCAG', 'tj':'TCCTTCTCACCACCACCTCCGC', 'ref':'GCCAGCCCGCGCCTCTGGCTCCTTTGCCCCCATTAGCCAGRCGCCACCATCCTTCTCACCACCACCTCCGCTGGTCCCTCC','snp':'R'},
					'f33':{'ext':'GTTGCTTTGTTGTTGTCTCTGTATGA', 'tj':'TGTCTAGTGATGGGTGGATCAGTC', 'ref':'ATAAAGACTTAGATGTTGCTTTGTTGTTGTCTCTGTATGASAGCTGTCTAGTGATGGGTGGATCAGTCGCTGTTTGGTCAT','snp':'S'},
					'f34':{'ext':'AAAGGTGTCCTGGGAACAATCACC', 'tj':'TGTGGTAAAAACAGAAGGGCGGATG', 'ref':'CAGTGTTATTAGGTATAAAGGTGTCCTGGGAACAATCACCRCTGTGGTAAAAACAGAAGGGCGGATGAAACTCTACAGCGG','snp':'R'},
					'f35':{'ext':'ACACCGCGGTACTGGGCGCTG', 'tj':'TGTAGCGCGCACTGGCCCCTG', 'ref':'GGTCAGAATGGTGCCCATCACACCGCGGTACTGGGCGCTGRCTGTAGCGCGCACTGGCCCCTGACTTTCTCCTTGGATCTG','snp':'R'},
					'f36':{'ext':'GGTGTTTTCTCTGTTTATGGACAGAG', 'tj':'TGTGTGTGTAGTATGTGTTTTTTGCAT', 'ref':'ATACTGCTTTCTTCGGTGTTTTCTCTGTTTATGGACAGAGMAGAAACCAGTGTGTGTGTAGTATGTGTTTTTTGCATGGGG','snp':'M'},
					'f37':{'ext':'AAAAGTTATTGACCCCTAGAAAATGAC', 'tj':'CACATACTTCAATTAAAATGAGTAGTAC', 'ref':'TATAAATCAATACAAAAGTTATTGACCCCTAGAAAATGACRATCACATACTTCAATTAAAATGAGTAGTACAAAGAGACTT','snp':'R'},
					'f38':{'ext':'TGTGGGTGTATTTAACCAATAGCAGT', 'tj':'TTTCCAATCAGTTTGTCTTAAAGATTTG', 'ref':'AAATTATCTCAACATGTGGGTGTATTTAACCAATAGCAGTMGTTTCCAATCAGTTTGTCTTAAAGATTTGCCTTTTGTTTT','snp':'M'},
					'f39':{'ext':'TTAGAAGGCAAAGATTATCATGCTGG', 'tj':'CCCAATTATAAACAACCTACAAGAAATC', 'ref':'CCTTAACATCCCATTTAGAAGGCAAAGATTATCATGCTGGRTTAAAAAGACCCAATTATAAACAACCTACAAGAAATCAAT','snp':'R'},
					'f40':{'ext':'AAGATTTAAAATTAAGCATAGAAATATTGAA', 'tj':'CCCCATAATATTAATGGTTTTGTGTTAA', 'ref':'ATCCTAAAAAAGATTTAAAATTAAGCATAGAAATATTGAARAGAACCCCATAATATTAATGGTTTTGTGTTAAATCTGCAT','snp':'R'},
					'f41':{'ext':'TACTGTGATTGCTTTCAATCTATAAATC', 'tj':'CGGAAAATTGATGTCTCAACAATACTG', 'ref':'TACCGGGATTTTTACTGTGATTGCTTTCAATCTATAAATCWTGTTGCGGAAAATTGATGTCTCAACAATACTGAGTCCAAT','snp':'W'},
					'f42':{'ext':'CTTTGTAGGCTAGGGCAAGGAATTT', 'tj':'CTAAGCTCTTAAAAGGTTTCAAGCATG', 'ref':'ACCAGACCATAAAGGCTTTGTAGGCTAGGGCAAGGAATTTRGGGTTTATTCTAAGCTCTTAAAAGGTTTCAAGCATGCGTG','snp':'R'},
					'f43':{'ext':'CACTGAGCTTTAACACCAATAAGTGC', 'tj':'CCGTATTCATCCTCATTCTCATTTCT', 'ref':'TCAACTGCTGTGGACACTGAGCTTTAACACCAATAAGTGCRTGCCGTATTCATCCTCATTCTCATTTCTTTCTCAGTGCCT','snp':'R'},
					'f44':{'ext':'CTAGAGCAATCCGAATAAAACATAATCA', 'tj':'CCACTCCTGCTCATATCTCTAGTG', 'ref':'ATTAAAACACAGCTAGAGCAATCCGAATAAAACATAATCAKGCCACTCCTGCTCATATCTCTAGTGGCGTCTCATCTCTAA','snp':'K'},
					'f45':{'ext':'TTGCTCCTCCACATATATTAATTTTATG', 'tj':'TCCTTTATTGGTTACATTTGTCATTTCC', 'ref':'TCCTTTTCTTTTTTGCTCCTCCACATATATTAATTTTATGRGCTCCTTTATTGGTTACATTTGTCATTTCCACCACGATAT','snp':'R'},
					'f46':{'ext':'GCATAAAGTTATATAGAAAAGACCTACA', 'tj':'CGGAAGACAAGTACATCAACACGTA', 'ref':'CATTTAAGGGTGGCATAAAGTTATATAGAAAAGACCTACARGACGGAAGACAAGTACATCAACACGTAGTGAAATTAAGCA','snp':'R'},
					'f47':{'ext':'GATAGGATGCGTACAGGAATAAATCT', 'tj':'CCTAATTTTACGTGAGATATTTCTCAG', 'ref':'TTTTAATGATGTAAGATAGGATGCGTACAGGAATAAATCTRTAGCCTAATTTTACGTGAGATATTTCTCAGATATTGATCT','snp':'R'},
					'f48':{'ext':'ACATGTTTTTCTCTATACCTTTTGGTG', 'tj':'TCCTAAAAGCCTAAGTCATAAATGCAT', 'ref':'CTGAGCTCTTTATACATGTTTTTCTCTATACCTTTTGGTGRGAAAATCCTAAAAGCCTAAGTCATAAATGCATATAAACTA','snp':'R'},
					'f49':{'ext':'GGTAAATTAAAGAAATGTGAAGAATTTTTA', 'tj':'CAGTGGGTGTGAATGTTCTTATTTGT', 'ref':'TTCAAGTTCTGGTAAATTAAAGAAATGTGAAGAATTTTTAWAATACAGTGGGTGTGAATGTTCTTATTTGTAGGACTGATT','snp':'W'},
					'f50':{'ext':'CAATCCTGGATCTATGTAATGGATGT', 'tj':'TATTGATTGGATATCACTTTTTCACATC', 'ref':'AAACTGTTTTCTTTCAATCCTGGATCTATGTAATGGATGTRTATTGATTGGATATCACTTTTTCACATCTCAGATAACTAT','snp':'R'},
					'f51':{'ext':'AATTACCTTATGATAATTTATTATGGCAAT', 'tj':'CTTATTTTTCTTTGAAAAATGTTTTCAAG', 'ref':'CTGTAATTATAATTACCTTATGATAATTTATTATGGCAATRATGGGAATTTGCTTATTTTTCTTTGAAAAATGTTTTCAAG','snp':'R'},
					'f52':{'ext':'AGCACTTTATTCTCTGCAGGAGTATA', 'tj':'CAGCAAAACCCACCAGGCCTATTC', 'ref':'TTCACGAAGCACAAAGCACTTTATTCTCTGCAGGAGTATARTGGGATTTGACAGCAAAACCCACCAGGCCTATTCCCTGAC','snp':'R'},
					'f53':{'ext':'TAGTTTTGAGTATTGTAGATACCTTATG', 'tj':'CAAATCTTTTATCAGATAGATAATTGGC', 'ref':'ATTATTATTATTTAGTTTTGAGTATTGTAGATACCTTATGKATACAAATCTTTTATCAGATAGATAATTGGCAAATATGTT','snp':'K'},
					'f54':{'ext':'CATTTACTTGCTTGTCCCCTTCCAG', 'tj':'CATGTAGAGCAGTACCATAATTGTTAG', 'ref':'CCACCATCCACCATCCATTTACTTGCTTGTCCCCTTCCAGWGTACATGTAGAGCAGTACCATAATTGTTAGCCCATAGCCC','snp':'W'},
					'f55':{'ext':'ATACTCTTCCAATTGTTTTCTATCTTTTA', 'tj':'TGAATGAAAATTACCCACTTCTCATTAG', 'ref':'TTCTTTGATCTATACTCTTCCAATTGTTTTCTATCTTTTARCAATGAATGAAAATTACCCACTTCTCATTAGCCCTCTCAG','snp':'R'},
					'f56':{'ext':'GGATGCTTTAAGGAAGCATCTTGCC', 'tj':'CACAAAGCTAACAAGAGTGGAGCC', 'ref':'TAGGAGGAAAGGAAAGGATGCTTTAAGGAAGCATCTTGCCRGGAGACACAAAGCTAACAAGAGTGGAGCCTGCAGCTGGAG','snp':'R'},
					'f57':{'ext':'GAGGCTTGGCAGTCACCTGGGT', 'tj':'TACAAGGGCCTGAGCCAAAGTGG', 'ref':'ATCTCGGAGACCAGCTTAGAGGCTTGGCAGTCACCTGGGTRCAGGATACAAGGGCCTGAGCCAAAGTGGTGAGGGAGGGTG','snp':'R'},
					'f58':{'ext':'TTAACTGTGGTCATAGTCACTCTCAA', 'tj':'CCACCCTCTAATACCTAAGTGGAAA', 'ref':'AATTTATTACCAGATTAACTGTGGTCATAGTCACTCTCAARTAATGGCCACCCTCTAATACCTAAGTGGAAAAACAGTTTT','snp':'R'},
					'f59':{'ext':'TGCTGTATCAACATTTCAGCAGCTC', 'tj':'CTTTTAACTCTAAGAAAATGGCAGCTC', 'ref':'CAGAGGAAAGAATGATGCTGTATCAACATTTCAGCAGCTCRTCTTTTAACTCTAAGAAAATGGCAGCTCCTAAATGTTCAA','snp':'R'},
					'f60':{'ext':'GTTCCTGACCTCTGAATGAGAGGG', 'tj':'TGTGTAAGGCCAATGCCTGGGAG', 'ref':'TTCATTTTCCAGGGAAGTTCCTGACCTCTGAATGAGAGGGRCTGTGTAAGGCCAATGCCTGGGAGGAAGGCAAGGATGAGT','snp':'R'},
					'f61':{'ext':'CCAGTGGGTGGGAGAGAAAGGAC', 'tj':'CAAGTGGTAAAGGCCCTCACAGAG', 'ref':'TAGTGGGATTCAGGCTCCCAGTGGGTGGGAGAGAAAGGACRTAAGGAAGCAAGTGGTAAAGGCCCTCACAGAGTATCAGCA','snp':'R'},
					'f62':{'ext':'TTCTTACGGAGTTGTCTTAATAAAGAC', 'tj':'CTTCAGCGTCCAGGTCTCCCTC', 'ref':'GGCCCAGCAGGAATTCTTACGGAGTTGTCTTAATAAAGACWGAAAATGCTTCAGCGTCCAGGTCTCCCTCTCTTGGGGACC','snp':'W'},
					'f63':{'ext':'CCCACCTTGGGTGGTGTAATACAC', 'tj':'CACGGGAGCCCCGGGTAGTAAC', 'ref':'GATAGATCCCGAGTGCCCCACCTTGGGTGGTGTAATACACRATATCACGGGAGCCCCGGGTAGTAACCACGGAGGTGTCAG','snp':'R'},
					'f64':{'ext':'GGTCTTCGGGAACCCTAAAATCAC', 'tj':'TAATCCTATCTCTTGTCTCTGCTTCT', 'ref':'TGGGAGATATGTAGATGGTCTTCGGGAACCCTAAAATCACSAGTAATCCTATCTCTTGTCTCTGCTTCTCCATATGAGTCA','snp':'S'},
					'f65':{'ext':'GGTTTTAAATCCCTTACAAACATTAGATT', 'tj':'CCCTTTTAAAATATACTTTTTTACACATG', 'ref':'CATATTTATATGGTTTTAAATCCCTTACAAACATTAGATTKTATCCCTTTTAAAATATACTTTTTTACACATGGAAGCATT','snp':'K'},
					'f66':{'ext':'TATTTTTTAAACATCTCAAGTGCTGTCA', 'tj':'TGCAGCACTGGAAACCTCTGTAGT', 'ref':'AGAAGGGATCAGTATTTTTTAAACATCTCAAGTGCTGTCARCGTGCAGCACTGGAAACCTCTGTAGTAGAAGTCCCTGAGC','snp':'R'},
					'f67':{'ext':'AATCATATCGTCTGGCTAAAAGTTGTA', 'tj':'TAAAGAAAGCACTCTGAAGGCATGG', 'ref':'AGAAAAAGTAAGAAATCATATCGTCTGGCTAAAAGTTGTARGAGAAATAAAGAAAGCACTCTGAAGGCATGGCCGAACTAC','snp':'R'},
					'f68':{'ext':'CATGATCAGAAACTACCAAATCATGAT', 'tj':'TTTATGTGTGAATGTGTTAGAAAATTGTG', 'ref':'TCTCAGTACTACACATGATCAGAAACTACCAAATCATGATSACGTTTATGTGTGAATGTGTTAGAAAATTGTGATGTAGTG','snp':'S'},
					'f69':{'ext':'TTATTGCTTGAGAAATAAATAATATAAATAAA', 'tj':'CTTGATTTAGGAAATGAAGTTGTGCTA', 'ref':'TTTGGCATTTATTGCTTGAGAAATAAATAATATAAATAAARTAGCTTGATTTAGGAAATGAAGTTGTGCTAGCTGGAACAC','snp':'R'},
					'f70':{'ext':'CAATCTCCCCTATTATTTCTATCCTAT', 'tj':'CATTGGTCATAACCAATGAACCAATAT', 'ref':'TATACTTCACACCCAATCTCCCCTATTATTTCTATCCTATRTAAGTGAGGTGCATTGGTCATAACCAATGAACCAATATTG','snp':'R'},
					'f71':{'ext':'AAAAAAATTGCCATCTGGTTTCCAACA', 'tj':'TATTTTCTCTATCAGTACTTCAAAAATAAT', 'ref':'TTGCAATTATAAAAAAAAAATTGCCATCTGGTTTCCAACASTATTTTCTCTATCAGTACTTCAAAAATAATCTTTTTTTCT','snp':'S'},
					'f72':{'ext':'CCTCTCATGACCCAGGAAGGCC', 'tj':'CCCTCTTTGTGTTGTAGTCCATGC', 'ref':'AGATAACGCTGAAGCAGGCCTCTCATGACCCAGGAAGGCCRGGGTGGATCCCTCTTTGTGTTGTAGTCCATGCTATTAAAA','snp':'R'},
					'f73':{'ext':'CCACAATACTGTAATCACACAGATAC', 'tj':'CCCAGACAGCACCTACAAGCAGG', 'ref':'AATAATTCATATCACCACAATACTGTAATCACACAGATACRCCCAGACAGCACCTACAAGCAGGAGTTGGGATGAAGCCAC','snp':'R'},
					'f74':{'ext':'CGCAGAGCCTAATCACTACACCC', 'tj':'TCTCTGCTAGGGTTTCATGACTTCG', 'ref':'GAGCCTGCAGCTGGAGCCGCAGAGCCTAATCACTACACCCRCCCATCTCTGCTAGGGTTTCATGACTTCGTATCGGGGATT','snp':'R'},
					'f75':{'ext':'GGGAAGTGATTTACAGGCAGGAGTA', 'tj':'TGTGGAATACAGAATAATTTGACTTTAAT', 'ref':'CTGGTCATGAGCTAGGGGAAGTGATTTACAGGCAGGAGTARCATGTGGAATACAGAATAATTTGACTTTAATAAATTTCAA','snp':'R'}
				}

ext_primers = {}
tj_primers = {}
extbyseq = {}
tjbyseq = {}
for var_id, info in VARIATIONINFO.iteritems():

	ext_primers[var_id] = info['ext']
	tj_primers[var_id] = info['tj']

	try:			 extbyseq[info['ext'][:19]].append(var_id)
	except KeyError: extbyseq[info['ext'][:19]] = [var_id]

	try:			 tjbyseq[info['tj'][-19:]].append(var_id)
	except KeyError: tjbyseq[info['tj'][-19:]] = [var_id]

if __name__ == "__main__": lib_main()
