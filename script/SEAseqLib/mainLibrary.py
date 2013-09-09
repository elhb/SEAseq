import sys
import os
MASTER = os.getpid()
version = 'ALPHA 1.3'

######################### MAIN #########################

def lib_main():
	pass

######################### FUNCTIONS #########################

def initiate_file(filename, logfile, mode='w'):
    import os
    
    if type(logfile) == file and mode != 'r': logfile.write('Initiating '+filename+' ...\n')
    
    if mode == 'w' and os.path.exists(filename):
	tmp = filename
	filename = raw_input('WARNING: the file '+filename+' already excists. Enter an alternative filename (leave empty to overwrite):')
	if type(logfile) == file and mode != 'r': logfile.write('WARNING: the file '+filename+' already excists. Enter an alternative filename (leave empty to overwrite):')
	if not filename:
	    filename = tmp
	    if type(logfile) == file and mode != 'r': logfile.write('overwriting\n')
	#else:
	#    if type(indata.logfile) == file and mode != 'r': config.logfile.write('Creating file '+filename+'.\n')
    
    if mode =='ow': mode ='w'
    import re
    #if re.search('log',filename):	out = open(filename, mode,0)
    #else:				out = open(filename, mode,1)
    out = open(filename, mode,1)
    
    if type(logfile) == file and mode != 'r': logfile.write('File '+filename+' sucessfully initiated.\n')
    
    return out

def writelogheader(logfile):
    import sys
    import time
    import getpass
    username = getpass.getuser()
    logfile.write('----------------\n')
    logfile.write('Running program: '+' '.join(sys.argv)+'.\n')
    logfile.write('Version: '+version+'\n')
    logfile.write('time: '+time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+'\n')
    logfile.write('Master process id='+str(MASTER)+'\n')
    logfile.write('Started by user = '+username+'\n')

def gi2orgname(gi_number):
	from Bio import Entrez
	Entrez.email = "erik.borgstrom@scilifelab.se"
	handle = Entrez.efetch(db="nucleotide", id=gi_number, retmode="xml")
	records = Entrez.read(handle)
	assert len(records) == 1
	return records[0]['GBSeq_organism']

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
			if c1 == 'N' or c2 == 'N': substitutions -= 1 #if N then no mismatch
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]

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

def getPairs(config):
	""" yield one readpair at a time from indata
	"""
	import re
	import gzip	

	# choose random reads to analyze from fastq files
	if config.random:
		import random
		numreads=config.randomumreads
		if config.stop: numreads = config.stop
		config.logfile.write('Choosing '+str(config.random)+' random pairs to analyze ... ')
		readNumbersToPrint = {}
		while len(readNumbersToPrint) < config.random: readNumbersToPrint[random.randint(config.skip+1,numreads)] = True
		tempo = readNumbersToPrint.keys()
		tempo.sort()
		readNumbersToPrint = tempo
		config.logfile.write('done.\n')
	
	# set the counters to initial values
	counter = 0
	tmp=0
	header="NA"
	if config.skip: skip =True

	# unpack infiles
	for file1, file2 in zip(config.infiles['r1'],config.infiles['r2']):
	
		#check if files are gzipped
		if file1.split('.')[-1] in ['gz','gzip']: file1 = gzip.open(file1)
		else: file1 = open(file1,'r')
		if file2.split('.')[-1] in ['gz','gzip']: file2 = gzip.open(file2)
		else: file2 = open(file2,'r')

		# itarate through fastq files and return readpairs
		for r1line in file1:
			r2line = file2.readline() #get line from  read2 file
			
			tmp+=1 # increment linecount
			
			#random speedup
			if config.random and counter != readNumbersToPrint[0]:
				#print tmp
				if tmp == 4: tmp = 0;continue
				elif tmp != 1: continue
				#elif tmp == 1:
					#counter+=1 # increase entry counter
					#print '\t',counter,readNumbersToPrint[0]
					#continue
				
			
			# skip or stop if option is set on config
			if config.skip and tmp < (4*config.skip) and skip: continue
			elif config.skip and tmp == (4*config.skip) and skip: skip=False; tmp =0;continue
			if config.stop and counter == config.stop: break

			# depending on line number (within entry) do ...	
			if tmp == 1: #header check match between files
				counter+=1 # increase entry counter
				header=r1line
				if r1line.split(" ")[0] != r2line.split(" ")[0]: config.logfile.write('Error mismatching headers!');raise ValueError #os.kill(MASTER,1);sys.exit(1);#REALLYNOTOPTIMAL
			elif tmp == 2: # sequence check that its DNA and save sequences till later
				if counter == 1:
					config.logfile.write('Checking data type of read 1 in pair 1 ... ')
					match = re.match("^[AGTCN]+$",r1line.rstrip())
					if match: config.logfile.write('this is DNA data.\n')
					else: config.logfile.write(' this is not a DNA sequence ('+r1line.rstrip()+') could something be wrong with your fastq file?.\n\n');raise ValueError#os.kill(MASTER);sys.exit();#REALLYNOTOPTIMAL
				r1seq = r1line.rstrip()
				r2seq = r2line.rstrip()
			elif tmp == 3: # "+"-line do some format check
					if counter in {1:True,67:True,438:True,9675:True,53678:True,864513:True,1337354:True,317955:True,1226844:True,20389:True,118261:True}:
						if r1line[0] != r2line[0] or r1line[0] != '+': config.logfile.write('Error Format not fastq!');raise ValueError#os.kill(MASTER);sys.exit(1);#REALLYNOTOPTIMAL
			elif tmp == 4: # quality values and end of entry, reset counter and yeild readpair
					tmp=0 # reset line counter
					r1qual = r1line.rstrip() #get qual strings
					r2qual = r2line.rstrip()
					
					#yield readpair
					if not config.random: yield [readpair(header.rstrip(), read(header.rstrip(),r1seq,r1qual), read(header.rstrip(),r2seq,r2qual)),config]
					elif counter == readNumbersToPrint[0]:
						yield [readpair(header.rstrip(), read(header.rstrip(),r1seq,r1qual), read(header.rstrip(),r2seq,r2qual)),config]
						readNumbersToPrint = readNumbersToPrint[1:]
						if len(readNumbersToPrint) == 0: break

def comp(str):
	return str.replace("A","X").replace("T","A").replace("X","T").replace("G","X").replace("C","G").replace("X","C")

def revcomp(str):
	return comp(str[::-1])

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

def UIPAC2REGEXP(string):
    return string.replace('R','[AG]').replace('Y','[CT]').replace('S','[GC]').replace('W','[AT]').replace('K','[GT]').replace('M','[AC]').replace('B','[CGT]').replace('D','[AGT]').replace('H','[ACT]').replace('V','[ACG]').replace('N','.')

def classify_cluser(config,infile="temporary.cluster.files/1.reads",database="reference/4amplicons/4amplicons.fasta"):

	#database="../reference/4amplicons/4amplicons.fasta"
	
	from Bio.Blast.Applications import NcbiblastnCommandline
	from Bio.Blast import NCBIXML
	from cStringIO import StringIO
	import time
	
	#setting up blast
	database=config.blastdb
	if config.blastsetting == 'strict':cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, out=infile+'.'+config.blastid+'.blastout')
	elif config.blastsetting == 'sloppy':cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn', out=infile+'.'+config.blastid+'.blastout')
	
	#Blasting all reads in cluster 
	blast_handle = cline.__call__()

#	blast_handle = StringIO(blast_handle[0])
#	blast_handle.seek(0)

	f = open(infile+'.'+config.blastid+'.blastout')
	records = NCBIXML.parse(f)

	results = {'total':0}
	records_counter = 0
	#checking blast results
	if config.printblast:o=''
	for blast_record in records:
		records_counter +=1

		#print blast_record.query+'\t',
#		if config.printblast: f.write(blast_record.query+'\n')
		readnumber = int(blast_record.query.split('_r')[1].split('_')[0])
		if readnumber == 1: r1_header = blast_record.query
		elif readnumber == 2 and records_counter%2 == 0: r2_header = blast_record.query
		else: config.logfile.write('Warning: readnumber is funky!\n')
		
		if config.printblast:o+=blast_record.query+'\n'
		if blast_record.alignments:
			if config.printblast:
				alignment = blast_record.alignments[0]
				for hsp in alignment.hsps:
					o+='\t'+ '****Alignment****'+'\n'
					o+='\t'+ 'sequence: '+ alignment.title+'\n'
					o+='\t'+ 'length: '+ str(alignment.length)+'\n'
					o+='\t'+ 'e value: '+ str(hsp.expect)+'\n'
					o+='\tq: '+str(hsp.query_start)+'\t'+ hsp.query +'\n'
					o+='\tm:  \t'+ hsp.match +'\n'
					o+='\ts: '+str(hsp.sbjct_start)+'\t'+ hsp.sbjct +'\n'
					if len(hsp.query) < 40 : o+='\tTo short will ba counted as "No Hit"\n'

			#		if len(hsp.query) < 40 : f.write('\tTo short will ba counted as "No Hit"\n')
			if readnumber == 1:
				r1_header = blast_record.query
				r1_subj_name = blast_record.alignments[0].title.split(' ')[1]
				if len(blast_record.alignments[0].hsps[0].query) <= 40: r1_subj_name = 'No Hits'
			elif readnumber == 2 and records_counter%2 == 0:
				r2_header = blast_record.query
				r2_subj_name = blast_record.alignments[0].title.split(' ')[1]
				if len(blast_record.alignments[0].hsps[0].query) <= 40: r2_subj_name = 'No Hits'
			else: config.logfile.write('Warning: readnumber is funky!\n')
		else:
			if readnumber == 1: r1_subj_name = 'No Hits'
			elif readnumber == 2: r2_subj_name = 'No Hits'

		if records_counter%2 == 0:
			if r1_header.split('_')[0] == r2_header.split('_')[0]:
				[junk0,junk1,rawrn,cluster_id,n15] = r1_header.split('_') #@M00275:102:000000000-A33TB:1:1101:16119:1648 _ 1:N:0:7 _ 205 _ AAGAGTCAGACTGAA
				cluster_id = int(cluster_id)
				results['total']+=1
				
				try: results[cluster_id]['total']+=1
				except KeyError: results[cluster_id]= {'total':1,'output':''}
				
				if r1_subj_name == r2_subj_name:
					hit = r1_subj_name
					if config.printblast: o+='## Read Pair Agree!\n'
				else:
					hit = 'Pair Disagree'
					if config.printblast: o+='## Read Pair Disagree!\n'
				try: results[cluster_id][hit]+=1
				except KeyError:results[cluster_id][hit]=1
				
				if config.printblast:results[cluster_id]['output'] += o+'\n\n'; o=''
			
			else: config.logfile.write('WARNING: read pair headers missmatch!\n')
		if config.printblast:o+='\n'+'\n'
	if config.printblast:pass
		#f2= open(config.path+'/blastReport.txt','a')
		#f2.write(o+'\n\n');
		#f2.close()
	f.close()
	
	import os
	os.remove(infile)
	os.remove(infile+'.'+config.blastid+'.blastout')
	
	return results

######################### CLASSES #########################

class Progress():

	def __init__(self,total, verb='full', logfile=sys.stderr, unit='read' ,mem=False, printint=0):
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
		self.mem = mem
		if printint: self.printint = printint

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
				', left: '+str(self.stf/60/60)+'h '+str(self.stf/60%60)+'min '+str(self.stf%60)+'s')
			if self.mem:
				import resource
				self.logfile.write(', using '+str((resource.getrusage(resource.RUSAGE_SELF).ru_maxrss+resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)/1024)+' ('+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)+') MB.\n')
			else:	self.logfile.write('\n')
			if self.type == 'minimal': self.logfile.write('..')
			self.ltime = time.time()
			self.lcurrent = self.current
			self.lpercentage = self.percentage

	def __exit__(self, *args):
		self.logfile.write('\n')

class Configuration():
    
    def __init__ (self, path, cmd, stop=None, skip=None ,random=None ):

	# permanent
	self.path 		= path
	self.config		= self.path+'/'+'config'
	self.init_logfile	= self.path + '/' + 'init.log.txt'
	self.init_outfile	= self.path + '/' + 'init.out.txt'
	self.addfqs_logfile	= self.path + '/' + 'addfqs.log.txt'
	self.addfqs_outfile	= self.path + '/' + 'addfqs.out.txt'
	self.cluster_logfile	= self.path + '/' + 'cluster.log.txt'
	self.cluster_outfile	= self.path + '/' + 'cluster.out.txt'
	self.sortreads_logfile	= self.path + '/' + 'sort.log.txt'
	self.sortreads_outfile	= self.path + '/' + 'sort.out.txt'
	self.meta_logfile	= self.path + '/' + 'meta.log.txt'
	self.meta_outfile	= self.path + '/' + 'meta.out.txt'
	self.sbatch_logfile	= self.path + '/' + 'sbatch.log.txt'
	self.sbatch_outfile	= self.path + '/' + 'sbatch.out.txt'
	self.metagraph_logfile	= self.path + '/' + 'graph.log.txt'
	self.metagraph_outfile	= self.path + '/' + 'graph.out.txt'
	self.classifymeta_logfile= self.path + '/' + 'classify.log.txt'
	self.classifymeta_outfile= self.path + '/' + 'classify.out.txt'
	self.clusters_file	= self.path+'/barcode_clusters_dictionary'
	self.abspath		= None # not loaded or set

	# longlasting
	self.infiles		= {'r1':[],'r2':[]}
	self.readcounts		= []
	self.chandlemissmatch	= None
	self.barcodemissmatch	= None
	self.clustercount	= None
	self.numberofseeds	= None

	# for each run
	self.cmd		= cmd
	self.stop		= stop
	self.skip		= skip
	self.random		= random
	self.sortformat		= 'fq'
	self.read_count_per_barcode_cluster_cutoff = 1
	
	
	if cmd == 'init':
	    self.logfile = self.init_logfile
	    self.outfile = self.init_outfile
	elif cmd == 'addfqs':
	    self.logfile = self.addfqs_logfile
	    self.outfile = self.addfqs_outfile
	elif cmd == 'clusterbarcodes':
	    self.logfile = self.cluster_logfile
	    self.outfile = self.cluster_outfile
	elif cmd == 'sortreads':
	    self.logfile = self.sortreads_logfile
	    self.outfile = self.sortreads_outfile
	elif cmd == 'meta':
	    self.logfile = self.meta_logfile
	    self.outfile = self.meta_outfile
	elif cmd == 'sbatch':
	    self.logfile = self.sbatch_logfile
	    self.outfile = self.sbatch_outfile
	elif cmd == 'metagraph':
	    self.logfile = self.metagraph_logfile
	    self.outfile = self.metagraph_outfile
	elif cmd == 'classifymeta':
	    self.logfile = self.classifymeta_logfile
	    self.outfile = self.classifymeta_outfile

    def getreads2process(self, ):
	
	self.logfile.write('Getting readcounts ...\n')
	total = 0
	for i in range(len(self.readcounts)):
	    rc = self.readcounts[i]
	    self.logfile.write(self.infiles['r1'][i]+' -> '+str(rc) +' reads.\n')
	    total += rc
	self.logfile.write(str(total)+' read pairs in fastq files.\n');
	
	# calculate the number of reads to process
	self.reads2process = total
	if self.skip: 	self.reads2process -= self.skip
	if self.stop: 	self.reads2process = self.stop
	if self.random:	self.reads2process = self.random

    def set(self,varname, value ):
	if varname == 'chandlemissmatch':	self.chandlemissmatch	= value
	elif varname == 'barcodemissmatch':	self.barcodemissmatch	= value
	elif varname == 'clustercount':		self.clustercount	= value
	elif varname == 'numberofseeds':	self.numberofseeds	= value
	elif varname == 'sortformat':		self.sortformat		= value
	elif varname == 'read_count_per_barcode_cluster_cutoff':self.read_count_per_barcode_cluster_cutoff = value
	else: raise ValueError

    def load(self, ):
	self.config = initiate_file(self.config, self.logfile , mode='r')
	for line in self.config:
	    if line.rstrip() == '# Absolute path:':
		self.abspath = self.config.next().rstrip()
	    if line.rstrip() == '# Infiles dictionary:':
		self.infiles = eval(self.config.next())
	    if line.rstrip() == '# Read counts list:':
		self.readcounts = eval(self.config.next())
	    if line.rstrip() == '# Number of cluster seeds:':
		self.numberofseeds = eval(self.config.next().rstrip())
	    if line.rstrip() == '# Number of barcode clusters identified:':
		self.clustercount = eval(self.config.next().rstrip())
	self.config.close()
	self.config = self.config.name

    def openconnections(self, ):
	if self.cmd == 'init':
	    self.outfile = initiate_file(self.outfile, self.logfile)
	    self.logfile = initiate_file(self.logfile, self.logfile)
	else:
	    self.outfile = initiate_file(self.outfile, self.logfile, mode='a')
	    self.logfile = initiate_file(self.logfile, self.logfile, mode='a')

    def save(self, ):

	if not os.path.exists(self.config):
	    self.config = initiate_file(self.config, self.logfile)
	    self.abspath = os.path.abspath(self.path)
	else:
	    self.config = initiate_file(self.config, self.logfile, mode='ow')

	self.logfile.write('Writing settings to config file ...\n')
	self.config.write(
	    '# Absolute path:\n'+str(self.abspath)+'\n'+
	    '# Infiles dictionary:\n'+str(self.infiles)+'\n'+
	    '# Read counts list:\n'+str(self.readcounts)+'\n'+
	    '# Number of cluster seeds:\n'+str(self.numberofseeds)+'\n'+
	    '# Number of barcode clusters identified:\n'+str(self.clustercount)+'\n'
	    )
	self.config.close()
	self.config = self.config.name

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

class readpair():
	""" object representing an illumina cluster """

	def __init__(self,header,r1,r2):
		self.header = header
		self.r1 = r1 #first read
		self.r2 = r2 #second read
		self.matchingprimerpairs = []
		self.id = 0

	def getN15(self):
		if self.handle_start:self.n15 = self.r1.subseq(0,self.handle_start)
		else: self.n15 = None
		return 0

	def identify(self, handle, config):
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r1)
		self.handle_start = handle_start
		self.handle_end   = handle_end
		return 0

	def identifyIllumina(self, config):
		handle = sequence('illuminaUniversal','AGATCGGAAGAGC','AGATCGGAAGAGC')
		config.chandlemissmatch = 2
		#handle = sequence('illumina','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r1)
		if handle_start: self.r1.illuminaadapter = True
		else: self.r1.illuminaadapter = False
		
		#handle = sequence('illumina','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r2)
		if handle_start: self.r2.illuminaadapter = True
		else: self.r2.illuminaadapter = False
		
		if self.r1.illuminaadapter or self.r2.illuminaadapter:
			self.isillumina = True
		else:	self.isillumina = False
		return 0

	def matchHandle(self, handle, config, read, matchfunk=hamming_distance):
		
		import re
		#matchfunk = hamming_distance
	
		handle_start = None
		handle_end   = None
	
		perfect_match = re.search(handle.seq, read.seq)
		if perfect_match:
		    handle_start = perfect_match.start()
		    handle_end = perfect_match.end()
		
		elif config.chandlemissmatch:
			mindist = [10000,-1]
			for i in range(len(read.seq)):
			    
				if i+len(handle.seq) <= len(read.seq):
					dist = matchfunk(handle.seq,read.seq[i:i+len(handle.seq)])
				else: dist = 1000
				
				if dist < mindist[0]: mindist =[dist,i]

			if mindist[0] < config.chandlemissmatch:
				handle_start = i
				handle_end = i+len(handle.seq)
			else:
				handle_start = None
				handle_end = None

		return [handle_start, handle_end]

	def get_cid(self,config):
		if self.n15 and self.n15.len == 15:
			try: self.cid = config.cid_by_bc[self.n15.seq]
			except KeyError: self.cid = False
		else: self.cid = None

	def matchprimerpair(self, primerpair):
		forwardprimermatch = self.matchfwd(primerpair)
		reverseprimermatch = self.matchrev(primerpair)
		if forwardprimermatch and reverseprimermatch: self.matchingprimerpairs.append(primerpair)
	
	def matchfwd(self, primerpair):
		import re
		return re.match(primerpair.fwdReStr, self.r1.seq[self.handle_end:])
	
	def matchrev(self, primerpair):
		import re
		return re.match(primerpair.revReStr, self.r2.seq)
	
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

	def reducebarcodes(self,config):
		""" Find most common barcodes in well ( > 10% ??), then try to place other barcodes to this cluster
		"""
		
		config.logfile.write( '\n' )
		
		config.minperc = 0
		maxdist = config.barcodemissmatch
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
		while len(highest) < config.numberofseeds:
			try:
				for bc in percentages[perc[0]]: highest.append(bc)
			except IndexError: pass
			perc=perc[1:]
		tempfile = open(config.path+'/predetermined_cluster_centers.fa','w')
		for bc in highest: tempfile.write('>'+bc+' '+str(self.barcodes[bc])+' readpairs\n'+bc+'\n')
		tempfile.close()

		tempfile = open(config.path+'/raw_barcode_sequences.fa','w')
		for bc in self.barcodes: tempfile.write('>'+bc+' '+str(self.barcodes[bc])+' readpairs\n'+bc+'\n')
		tempfile.close()
		del percentages

		# alternatively we could use chdit, though homopolymers seems to be problematic
		# cd-hit-454 -i tasks/20130805.1_index10/raw_barcode_sequences.fa -o DELETEME -c 0.85 -g 1 -T 8 -gap -6 -gap-ext -2 -AS 2
		# cdhit-cluster-consensus DELETEME.clstr tasks/20130805.1_index10/raw_barcode_sequences.fa DELETEME.cons DELETEME.aln
		import subprocess
		from cStringIO import StringIO
		import time
		import multiprocessing
		tempo = time.time()
		config.logfile.write('starting '+' '.join(['dnaclust','--similarity',str(1-(float(config.barcodemissmatch)/15)),'--input-file',config.path+'/raw_barcode_sequences.fa','-t',str(multiprocessing.cpu_count()),'--predetermined-cluster-centers',config.path+'/predetermined_cluster_centers.fa'])+'\n')
		dnaclust =               subprocess.Popen(['dnaclust','--similarity',str(1-(float(config.barcodemissmatch)/15)),'--input-file',config.path+'/raw_barcode_sequences.fa','-t',str(multiprocessing.cpu_count()),'--predetermined-cluster-centers',config.path+'/predetermined_cluster_centers.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		dnaclust_out, errdata = dnaclust.communicate()
		if dnaclust.returncode != 0:
			print 'dnaclust view Error code', dnaclust.returncode, errdata
			sys.exit()
		dnaclust_out = StringIO(dnaclust_out)
		seconds = round(time.time()-tempo,2)
		config.logfile.write('dnaclust done after '+str(int(seconds/60/60))+'h '+str(int(seconds/60%60))+'min '+str(int(round(seconds%60)))+'s, parsing result ... ')
		del dnaclust

		clusters={}
		cc=0
		for line in dnaclust_out:
			cc+=1
			clusters[cc] ={'total':0,'barcodes':{},'highest':[0,'XXXXXXXXXXXXXX']}
			line = line.rstrip().split('\t')
			for bc in line:
				if bc not in clusters[cc]['barcodes']: clusters[cc]['total']+=self.barcodes[bc]
				clusters[cc]['barcodes'][bc]=self.barcodes[bc]
				if self.barcodes[bc] > clusters[cc]['highest'][0]: clusters[cc]['highest']=[self.barcodes[bc],bc]
		config.logfile.write('almost done ... ')
		del dnaclust_out

		counter = 0
		reads_in_clusters={}
		for cc in clusters:
			try: reads_in_clusters[clusters[cc]['total']]+=1
			except KeyError: reads_in_clusters[clusters[cc]['total']]=1
			if int(clusters[cc]['total']) > 1:
				counter+=1
				#print cc,clusters[cc]['total'],clusters[cc]['highest'][1],clusters[cc]['highest'][0]

		config.logfile.write('ok done now I\'ll just print and plot some info ... then done ... for real!\n\n')
		config.outfile.write(str( cc)+' clusters whereof '+str(counter)+' has more than one read\n\n')

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
		for scale in [[0,5000,0,20],[0,1000,0,20]]:
			plt.figure()
			plt.axis(scale)
			#plt.xlabel('Total Number of Reads per Barcode Cluster')
			#plt.ylabel('Number of Clusters')
			#pos = np.arange(len(x))
			#width = 1.0     # gives histogram aspect to the bar diagram
			#ax = plt.axes()
			#ax.set_xticks(pos + (width / 2))
			#ax.set_xticklabels(x,rotation='horizontal')
			#plt.bar(pos, y, width, color='r')
			##plt.show()
			##plt.savefig(pp,format='pdf',bbox_inches=0)
			plt.plot(x,y)
			plt.suptitle('Y is Total Number of Clusters with X Reads Pairs per Barcode Cluster.', fontsize=12)
			plt.savefig(config.path+'/read_pairs_per_barcode_cluster.x_scale_'+str(scale[0])+'-'+str(scale[1])+'.y_scale_'+str(scale[2])+'-'+str(scale[3])+'.pdf')
			plt.close()
		config.logfile.write( 'done\n')

		self.clusters = clusters
		del clusters

		return

class PrimerPair(object):
	
	def __init__(self, fwd, rev, name):
		self.fwd = fwd
		self.rev = rev
		self.name = name
		self.fwdReStr = UIPAC2REGEXP(self.fwd)
		self.revReStr = UIPAC2REGEXP(self.rev)
		self.referenceAlleles = []

class BarcodeCluster(object):

	def __init__(self, id_number, barcode_sequence=None):
		self.id = id_number
		self.barcodesequence = barcode_sequence
		self.readpairs = []
		self.consensuses = {}
		self.amplicons = {}
		self.definedamplicons = {}
		self.readpairbyheader = {}
		self.headerbyint = {}

	def addreadpair(self, pair):
		if pair.cid: cluster_id = pair.cid
		else: cluster_id = int(pair.header.split(':')[-1].split('_')[1])
		if cluster_id == self.id:
			if self.barcodesequence and pair.n15:
				if pair.n15 != self.barcodesequence:
					import sys
					sys.stderr.write('ERROR: Barcode sequence from read pair does not match the barcode in cluster.\n')
					raise ValueError
			self.readpairs.append(pair)
		else:
			import sys
			sys.stderr.write('ERROR: Cluster id from read pair does not match the barcode cluster id.\n')
			raise ValueError
	
	@property
	def readcount(self):
		return len(self.readpairs)

	@property
	def ampliconpairs(self):
		tmp_counter = 0
		for pair in self.readpairs:
			if pair.isillumina or pair.primererror: continue
			tmp_counter += 1
		return tmp_counter
	
	@property
	def adaptercount(self):
		tmp_counter = 0
		for pair in self.readpairs:
			if pair.isillumina: tmp_counter += 1
		return tmp_counter
	
	@property
	def primererrors(self):
		tmp_counter = 0
		for pair in self.readpairs:
			if pair.primererror: tmp_counter += 1
		return tmp_counter

	@property
	def ampliconcount(self):
		return len(self.amplicons)

	@property
	def definedamplicons(self):
		tmp_counter = 0
		for amplicon in self.amplicons.values:
			if amplicon.allelecount >= 1: tmp_counter += 1
		return tmp_counter

	def createtempfile(self, config, verb=True):

		tmpcounter = 0
		if verb: output = 'PAIRS:\n'
		f = open(config.path+'/sortedReads/temporary.'+str(self.id)+'.fa','w')
		
		for pair in self.readpairs:
		    
		    tmpcounter += 1
		    pair.id = tmpcounter
		    self.headerbyint[pair.id] = pair.header
		    self.readpairbyheader[pair.header] = pair
		    if verb: output += str(pair.id)+'\t'
		    pair.primererror = None
		    
		    #identify subparts of read
		    C_HANDLE = sequence('c handle',"CTAAGTCCATCCGCACTCCT","CTAAGTCCATCCGCACTCCT")
		    pair.identify(C_HANDLE, config)
		    pair.getN15()
		    pair.identifyIllumina(config)
		    
		    #check if read is adaptersequence
		    if pair.isillumina:
			if verb: output+='illumina adapter in read pair\n';
			continue
		    
		    # Match the primer pairs
		    for name, primerpair in config.primerpairs.iteritems(): pair.matchprimerpair(primerpair)
	
		    # If only one primer pair match
		    if   len(pair.matchingprimerpairs) == 1:
			pair.p1 = pair.matchingprimerpairs[0].name
			pair.p2 = pair.matchingprimerpairs[0].name
	
		    # if no or several pairs match match fwd and rev seperately
		    elif len(pair.matchingprimerpairs) != 1:
			pair.p1 = ''
			pair.p2 = ''
			for name, primerpair in config.primerpairs.iteritems():
			    if pair.matchfwd(primerpair):
				if pair.p1: pair.p1 += '/'
				pair.p1 += primerpair.name
			    if pair.matchrev(primerpair):
				if pair.p2: pair.p2 += '/'
				pair.p2 += primerpair.name
			if not pair.p1: pair.p1 = '???'
			if not pair.p2: pair.p2 = '???'
		    if verb: output += pair.p1+'\t'+pair.p2+'\t';
		    
		    # check that primers match
		    if pair.p1 == '???' or pair.p1 != pair.p2 or len(pair.matchingprimerpairs) != 1:
			if pair.p1 == '???':                        pair.primererror = 'primers-not-identifiable'
			elif pair.p1 != pair.p2:                    pair.primererror = 'fwd-rev-pair-missmatch'
			elif len(pair.matchingprimerpairs) != 1:    pair.primererror = 'more-than-one-pair-match'
			if verb:
			    output += pair.primererror+'\t'+pair.r1.seq +' '+ pair.r2.seq+'\n'
			continue
	
		    #look for unexpected primer sequences
		    import re
		    matched_primerpair = config.primerpairs[pair.p1]
		    fwd_in_r2 = None
		    rev_in_r1 = None
		    fwdcount = 0
		    revcount = 0
		    other_fwd_in_any = None
		    other_rev_in_any = None
		    for name, primerpair in config.primerpairs.iteritems():
			if primerpair.name == matched_primerpair.name:
			    fwd_in_r2  = re.search(     matched_primerpair.fwdReStr,    pair.r2.seq) #fwd in read 2
			    rev_in_r1  = re.search(     matched_primerpair.revReStr,	pair.r1.seq) #rev in read1
			    fwdcount = len(re.findall(  matched_primerpair.fwdReStr,    pair.r1.seq))
			    revcount = len(re.findall(  matched_primerpair.revReStr,    pair.r2.seq))
			else:
			    other_fwd_in_any = re.search( primerpair.fwdReStr,   	pair.r1.revcomp().seq + 'NNNNN' + pair.r1.seq + 'NNNNN' + pair.r2.seq + 'NNNNN' + pair.r2.revcomp().seq) # fwd in any read
			    other_rev_in_any = re.search( primerpair.revReStr,          pair.r1.revcomp().seq + 'NNNNN' + pair.r1.seq + 'NNNNN' + pair.r2.seq + 'NNNNN' + pair.r2.revcomp().seq) # rev in any read
		    if fwd_in_r2 or rev_in_r1 or other_fwd_in_any or other_rev_in_any or revcount != 1 or fwdcount != 1:
			    primererror+=1;
			    if verb:
				output+='PRIMER ODD COMBO\t';
				output += pair.r1.seq +' '+ pair.r2.seq+'\n'
			    continue
	
		    # Add sequences to output
		    if verb: output += pair.r1.seq +' '+ pair.r2.seq+'\n'
		    
		    # prepare for clustering by writing read pair to tempfile and saving temporary id number and mappinf header to pair-object
		    tem_seq = pair.r1.seq[pair.handle_end:][len( config.primerpairs[pair.p1].fwd )+1:]+'NNNNNNNNNN'+pair.r2.revcomp().seq[:-(len(    config.primerpairs[pair.p1].rev   )+1)]
		    f.write('>'+str(pair.id)+'\n'+ tem_seq +'\n')
		
		f.close()
		
		if verb:return output
	
	def clusterreadpairs(self, config, indata, verb=False):
		
		#check that there is data to work with
		if self.adaptercount+self.primererrors == self.readcount:
		    import os
		    os.remove( config.path+'/sortedReads/temporary.'+str(self.id)+'.fa' )
		    return 'All adapter and/or primer error.\n'#['ONLY JUNK',return_info]
	
		# Cluster Read pairs
		import subprocess
		from cStringIO import StringIO
		cdhit = subprocess.Popen(
			[	'cd-hit-454',
				'-i',config.path+'/sortedReads/temporary.'+str(self.id)+'.fa',
				'-o',config.path+'/sortedReads/cluster.'+str(self.id)+'.fa',
				'-g','1',
				'-c',str(indata.clustering_identity/100.0)
				],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE )
		cdhit_out, errdata = cdhit.communicate()
		if cdhit.returncode != 0:
			print 'cmd: '+' '.join( ['cd-hit-454','-i',config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa','-o',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa','-g','1','-c',str(indata.clustering_identity/100.0)])
			print 'cd-hit cluster='+str(cluster.id)+' view Error code', cdhit.returncode, errdata
			sys.exit()
	
		# Build consensus sequences for read pair clusters
		ccc = subprocess.Popen(
			[	'cdhit-cluster-consensus',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.fa.clstr',
				config.path+'/sortedReads/temporary.'+str(self.id)+'.fa',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.consensus',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.aligned'
				],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE )
		ccc_out, errdata = ccc.communicate()
		if ccc.returncode != 0:
			print 'cmd: '+' '.join( ['cdhit-cluster-consensus',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr',config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned'])
			print 'cluster='+str(cluster.id)+' cdhit-cluster-consensus view Error code', ccc.returncode, errdata
			print ccc_out
			sys.exit()
	
		# output info from temporary files
		if verb:
			output = ''
			for info in [
				[config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus.fasta','\nCD-HIT consensus sequences:\n'],
				[config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr','\nClustering details:\n'],
				[config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned','\nAlignment details:\n']
				]:
				filename , message = info
				f = open(filename)
				tmp = f.read()
				output += message
				output += tmp
				f.close()
			return output

	def removetempfiles(self, config):
		import os
		os.remove(config.path + '/sortedReads/temporary.' + str(self.id) + '.fa')
		os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.fa')
		os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.fa.clstr')
		os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.consensus.fasta')
		os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.aligned')

	def loadconsensuses(self, config):

	#	get identity from file
		f = open(config.path+'/sortedReads/cluster.'+str(self.id)+'.fa.clstr')
		data = f.read()
		f.close()
		for consensus_identities in data.split('>Cluster '):
			consensusid = consensus_identities.split('\n')[0]
			consensus = Consensus(consensusid)
			for line in consensus_identities.split('\n')[1:]:
			    line=line.rstrip()
			    if not line: continue
			    read_id     = int(line.split('>')[1].split('.')[0])
			    identity = line.split('/')[-1]
			    pair = self.readpairbyheader[self.headerbyint[read_id]]
			    consensus.addreadpair(pair,identity)
			self.consensuses[consensus.id] = consensus
			
	def consensusesToAmplicons(self, config):
		for consensusid, consensus in self.consensuses.iteritems():
			try: self.amplicons[consensus.type].addallele(consensus)
			except KeyError:
				self.amplicons[consensus.type] = Amplicon(consensus.type, consensus.primer)
				self.amplicons[consensus.type].addallele(consensus)

	def loadconsensusalignemnts(self, config):
		# get alignments from file
		f = open(config.path+'/sortedReads/cluster.'+str(self.id)+'.aligned')
		data = f.read()
		f.close()
		consensus = None
		for cluster_aln in data.split('===========================================\n'):
			consensusalignment = ''
			tmpseqs = {}
			for part in cluster_aln.split('\n\n'):
				for line in part.split('\n'):
				    line=line.rstrip()
				    if not line: continue
				    if line[0] == 'A':
					consensusid=line.split(' ')[-1].split(':')[0];
					consensus = self.consensuses[consensusid]
					consensus.alignmentStr = ''
					continue
				    read_id = line.split(':  +  ')[0]
				    seq     = line.split(':  +  ')[1]
				    if read_id == 'Consensus': consensus.alignmentStr += seq.split(' ')[0]
				    else:
					try: 		 tmpseqs[read_id] += seq
					except KeyError: tmpseqs[read_id]  = seq
			for readid, seq in tmpseqs.iteritems():
				consensus.readpairs[int(read_id)].alignmentStr = seq
	
	def loadconsensussequences(self, config):
		##load singelton consensus sequences
		f = open(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus.fasta')
		data = f.read()
		f.close()
		consensusid = None
		for consensus_identities in data.split('>')[1:]:
			if   consensus_identities[0] == 'c':
				consensusid = consensus_identities.split(' ')[0].split('_')[-1]
				readcount  = int(consensus_identities.split(' ')[2])
				assert self.consensuses[consensusid].readcount == readcount,	'ERRORSHMERROR 1'
			elif consensus_identities[0] == 's':
				consensusid = consensus_identities.split(' ')[0].split('_')[-1]
				read_id     = consensus_identities.split(' ')[1].split('\n')[0]
				assert self.consensuses[consensusid].readcount  == 1,		'ERRORSHMERROR 2'
				assert self.consensuses[consensusid].seedpairid == read_id,	'ERRORSHMERROR 3'
			consensus = self.consensuses[consensusid]
			tmpseq = ''
			for line in consensus_identities.split('\n')[1:]:
			    line=line.rstrip()
			    tmpseq += line
			consensus.sequence = sequence(consensus.id,tmpseq,tmpseq)

	def getDefinedAmplicons(self, ):
		for amplicon in self.amplicons.values():
			if amplicon.allelecount >= 1: self.definedamplicons[amplicon.type] = amplicon	
	

class Amplicon(object):

	def __init__(self, amplicon_type, primerpair):
		self.allels = []#hold the consensus sequences produced
		self.type = amplicon_type
		self.primer = primerpair
		self.allelecount = None

	@property
	def readcount(self):
		return sum([consensus.readcount for consensus in self.allels])
		
	def addallele(self,consensus ):
		if consensus.type == self.type: self.consensussses.append(consensus)
		else:
			import sys
			sys.stderr.write('ERROR: Amplicon type does not match the Consensus type.\n')
			raise ValueError

	def checkmono(self, indata):
		output = ''
		output += '\t'+self.type+' '+str(self.readcount)+' reads in total.\n'
		self.allelecount = 0
		self.monoclonal = None
		
		for consensus in self.allels:
			consensus.percentagesupport = 100*float(consensus.readcount)/float(self.readcount)
			output += '\t\tConsensus '+consensus.id+' supported by '+str(consensus.percentagesupport)+'% of readpop ('+str(consensus.readcount)+' reads)\t'+consensus.sequence+'\n'
			if  consensus.percentagesupport >= indata.minimum_support and consensus.readcount > indata.minimum_reads:
				self.allelecount += 1
		
		if self.allelecount == 1:
			self.monoclonal = True
			output += 'Monoclonal for '+self.type
		elif self.allelecount >1:
			self.monoclonal = False
			output += 'Polyclonal for '+self.type
		else:   output += 'No "good allels" found for '+self.type
		
		return output

class Consensus(object):

	def __init__(self, idnumber, amplicon_type=None, primerpair=None, sequence=None):
		self.readpairs = {}
		self.type = amplicon_type
		self.primer = primerpair
		self.id = idnumber
		self.sequence = sequence

	@property
	def readcount(self):
		return len(self.readpairs)

	def addreadpair(self, pair, identity):
		if identity[-1] == '*':
			identity = 'SEED'
			self.seedpairid = pair.id
		pair.consensusidentity = identity
		if self.type == None:
			self.type = pair.p1
			self.primer = pair.matchingprimerpairs[0]
		if pair.p1 == self.type:
			self.readpairs[pair.id] = pair
		else:
			import sys
			sys.stderr.write('ERROR: Consensus type does not match the readpair fwd primer.\n')
			sys.stderr.write('WARNING: mixed consensus clustering!\n')
			raise ValueError

	def alignmentoutput(self, config):
		#make output for alnignments
		output = ''
		output += 'Consensus number '+self.id+' from '+str(self.readcount)+' read pairs\t'+self.alignmentStr
		for readid, pair in self.readpairs.iteritems():
			output += 'Read pair id = '+str(pair.id)+'    \t'+pair.consensusidentity+'\t'+pair.p1+'\t'+pair.alignmentStr +'\n'		
		return output
	

if __name__ == "__main__": lib_main()
