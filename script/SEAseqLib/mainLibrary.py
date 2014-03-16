import sys
import os
MASTER = os.getpid()
version = 'ALPHA 1.5'

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
    if os.path.exists(filename) and os.path.islink(filename): os.unlink(filename)
    out = open(filename, mode,1)
    
    if type(logfile) == file and mode != 'r': logfile.write('File '+filename+' sucessfully initiated.\n')
    
    return out

def writelogheader(logfile):
    import sys
    import time
    import getpass
    import commands
    from socket import gethostname
    username = getpass.getuser()
    computer = gethostname()
    logfile.write('----------------\n')
    logfile.write('Running program: '+' '.join(sys.argv)+'.\n')
    logfile.write('Version: '+version+'\n')
    logfile.write('time: '+time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+'\n')
    logfile.write('Master process id='+str(MASTER)+'\n')
    logfile.write('Started by user = '+username+' on host = '+computer+'\n')
    if gethostname().split('.')[1] == 'uppmax':
        logfile.write('Program is run on uppmax, temporary files will be placed in '+commands.getoutput('echo $SNIC_TMP')+' .\n')

def gi2orgname(gi_number,database='/proj/b2011011/SEAseq/reference/NCBItaxonomy.db',lock=None):
	gi_number = int(gi_number)
	#print gi_number, 'local'
	import sqlite3
	# Connect to local database
	conn = sqlite3.connect(database)
	c = conn.cursor()
	t = (int(gi_number),)
	try:
		gi_number, name = c.execute('SELECT * FROM gi2name WHERE gi=?', t).fetchone()
		return name
	except TypeError:
		print gi_number, 'not in local db will fetch online and update db.'
		values      = []
		# do online
		from Bio import Entrez
		Entrez.email = "erik.borgstrom@scilifelab.se"
		handle = Entrez.efetch(db="nucleotide", id=str(gi_number), retmode="xml")
		records = Entrez.read(handle)
		assert len(records) == 1
		values.append((gi_number,records[0]['GBSeq_organism']))
		if lock: lock.acquire()
		try:
			gi_number, name = c.execute('SELECT * FROM gi2name WHERE gi=?', t).fetchone()
		except TypeError:
			c.executemany('INSERT INTO gi2name VALUES (?,?)', values)
			conn.commit();
		conn.close()
		if lock: lock.release()
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
		numreads=config.random
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
	for file1, file2 in zip(config.infilesDictionary['r1'],config.infilesDictionary['r2']):
	
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
			if config.stop and counter == config.stop+1: break

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
				if config.readLength:
					config.readLength = int(config.readLength)
					if config.readLength > len(r1seq) or config.readLength > len(r2seq):
						config.logfile.write('WARNING: readlength limit is '+str(config.readLength)+'bp, r1 is '+str(len(r1seq))+'bp and r2 is '+str(len(r2seq))+'bp.\n')
					else:
						if config.readLength < len(r1seq): r1seq = r1seq[:config.readLength]
						if config.readLength < len(r2seq): r2seq = r2seq[:config.readLength]
			elif tmp == 3: # "+"-line do some format check
					if counter in {1:True,67:True,438:True,9675:True,53678:True,864513:True,1337354:True,317955:True,1226844:True,20389:True,118261:True}:
						if r1line[0] != r2line[0] or r1line[0] != '+': config.logfile.write('Error Format not fastq!');raise ValueError#os.kill(MASTER);sys.exit(1);#REALLYNOTOPTIMAL
			elif tmp == 4: # quality values and end of entry, reset counter and yeild readpair
					tmp=0 # reset line counter
					r1qual = r1line.rstrip() #get qual strings
					r2qual = r2line.rstrip()

					if config.readLength:
						config.readLength = int(config.readLength)
						if config.readLength > len(r1qual) or config.readLength > len(r2qual):
							config.logfile.write('WARNING: readlength limit is '+str(config.readLength)+'bp, r1 is '+str(len(r1qual))+'bp and r2 is '+str(len(r2qual))+'bp.\n')
						else:
							if config.readLength < len(r1qual): r1qual = r1qual[:config.readLength]
							if config.readLength < len(r2qual): r2qual = r2qual[:config.readLength]
					
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
    return string.replace('R','[AG]').replace('Y','[CT]').replace('S','[GC]').replace('W','[AT]').replace('K','[GT]').replace('M','[AC]').replace('B','[CGT]').replace('D','[AGT]').replace('H','[ACT]').replace('V','[ACG]').replace('N','[AGTC]')

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

def clusterGenerator(config,indata):
    
    import cPickle
    import gzip
    
    import os
    import sys
    if   os.path.exists(config.path+'/clusters.pickle'):    filename = config.path+'/clusters.pickle'
    elif os.path.exists(config.path+'/clusters.pickle.gz'): filename = config.path+'/clusters.pickle.gz'
    else: config.logfile.write('Please run the "SEAseq meta" step first.\nNow exiting program.\n');sys.exit()
    
    if indata.tempFileFolder and not indata.debug and indata.tempFileFolder != config.path:
        config.logfile.write('Copying pickle file to temporary location for fast access:\n ');
        import os, shutil, glob
	pid = os.getpid()
        noPathFileName = filename.split('/')[-1]
        src = filename
        dst = indata.tempFileFolder+'/SEAseqtemp/'+str(pid)+'.'+noPathFileName
        if not os.path.exists(dst):
		config.logfile.write('\tcopying: '+src+' to '+dst+ '\n');
		shutil.copyfile(src, dst)
		config.logfile.write('\tDone.\n');
        filename = dst

    clusterundump = open(filename,mode='rb', buffering=1024*64)
    if clusterundump.name.split('.')[-1] in ['gz','gzip']:
        clusterundump.close()
        clusterundump = gzip.open(clusterundump.name)
    
    config_backup = config
    tmpCounter = 0
    while True:
	config = config_backup
        try:
            cluster = cPickle.load(clusterundump)
	    tmpCounter += 1
	    if indata.skip and tmpCounter <= indata.skip: config = 'skipped'
            yield [cluster,config,indata]
        except EOFError:
            config.logfile.write('All clusters read from file.\n')
            break

def getClassification(taxid=None,gi=None,database='/proj/b2011011/SEAseq/reference/NCBItaxonomy.db',lock=None):
    #
    # imports
    #
    import time
    import sys
    import sqlite3
    #
    # set initial values
    #
    start = time.time()
    classStr = []
    rankStr = []
    classification = {}
    #print 'database:';print database
    conn = sqlite3.connect(database)
    c = conn.cursor()
    try: # Try ready dist
	if taxid:
		t = (int(taxid),)
		taxid, classdist = c.execute('SELECT * FROM taxid2classdist WHERE taxid=?', t).fetchone()
		#print "got a prebuilt string"
		return eval(classdist)
	elif gi:
		t = (int(gi),)
		gi, classdist = c.execute('SELECT * FROM gi2classdist WHERE gi=?', t).fetchone()
		#print "got a prebuilt string"
		return eval(classdist)
	else: sys.stderr.write('Error: you need to supply either a gi number or taxid\n')
    except TypeError: # else build one
	#print 'no prebuilt found, bulding string:'
	name = None
	rank = None
	#
	# get taxid
	#
	if gi and not taxid:
	    #print 'fetching taxid ...'
	    t = (int(gi),)
	    try: gi, taxid = c.execute('SELECT * FROM gi2taxid WHERE gi=?', t).fetchone()
	    except TypeError:
		name = gi2orgname(gi,lock=lock)
		t = (name,)
		taxid, name = c.execute('SELECT * FROM taxid2name WHERE name=?', t).fetchone()
	    #print 'done'
	if not gi and not taxid: sys.stderr.write('Error: you need to supply either a gi number or taxid\n')
	#taxid = int(taxid)
	if taxid:classification['taxid'] = int(taxid)
	if gi: classification['gi'] = int(gi)
	originalTaxid = taxid
	#
	# While not rank root build on the classstring
	#
	itercounter = 0
	while name != 'root':
	    itercounter += 1
	    t = (int(taxid),)
	    taxid, name = c.execute('SELECT * FROM taxid2name WHERE taxid=?', t).fetchone()
	    taxid, rank = c.execute('SELECT * FROM taxid2rank WHERE taxid=?', t).fetchone()
	    oldtaxid, taxid = c.execute('SELECT * FROM taxid2parent WHERE taxid=?', t).fetchone()
	    classStr.append(name)
	    rankStr.append(rank)
	    classification[rank] = name
	    if itercounter > 20: break
	taxid = originalTaxid
	#
	#save the classtring to database
	#
	if lock: lock.acquire()
	if taxid:
		#print 'saving taxid to classstr db'
		try:
			t = (int(taxid),)
			taxid, classdist = c.execute('SELECT * FROM taxid2classdist WHERE taxid=?', t).fetchone()
			#print 'already in' , classdist
		except TypeError: #add
			values      = []
			values.append((int(taxid),str(classification)))
			c.executemany('INSERT INTO taxid2classdist VALUES (?,?)', values)
			conn.commit();
			#print 'saved'
	if gi:
		#print 'saving gi to classstr db'
		try:
			t = (int(gi),)
			gi, classdist = c.execute('SELECT * FROM gi2classdist WHERE gi=?', t).fetchone()
			#print 'already in'
		except TypeError: #add
			values      = []
			values.append((int(gi),str(classification)))
			c.executemany('INSERT INTO gi2classdist VALUES (?,?)', values)
			conn.commit();
			#print 'saved'
	if lock: lock.release()
    conn.close()
    end = round(time.time()-start,2)
    #print 'search took',int(end/60),'min',end%60,'seconds'
    return classification

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

	import sys
	import time
	import os
	if not os.path.isdir(path): os.mkdir(path)
	f = open(path+'/command.log.txt','a')
	f.write(time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+'\t'+' '.join(sys.argv)+'\n')
	f.close()

	# permanent
	self.path 		= path
	self.config		= self.path + '/'+'config'
	self.init_logfile	= self.path + '/' + 'init.log.txt'
	self.init_outfile	= self.path + '/' + 'init.out.txt'
	self.clone_logfile	= self.path + '/' + 'clone.log.txt'
	self.clone_outfile	= self.path + '/' + 'clone.out.txt'
	self.compare_logfile	= self.path + '/' + 'compare.log.txt'
	self.compare_outfile	= self.path + '/' + 'compare.out.txt'
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
	self.makegraphs_logfile	= self.path + '/' + 'graph.log.txt'
	self.makegraphs_outfile	= self.path + '/' + 'graph.out.txt'
	self.classifymeta_logfile= self.path + '/' + 'classify.log.txt'
	self.classifymeta_outfile= self.path + '/' + 'classify.out.txt'
	self.setVariables_logfile= self.path + '/' + 'setVariables.log.txt'
	self.setVariables_outfile= self.path + '/' + 'setVariables.out.txt'	
	self.clusters_file	= self.path+'/barcode_clusters_dictionary'
	self.absolutePath	= None # not loaded or set

	# longlasting
	self.infilesDictionary	= {'r1':[],'r2':[]}
	self.readCountsList	= []
	self.maxHandleMissMatch	= None
	self.maxBeadBarcodeMissMatch	= None
	self.numberOfBarcodeClustersIdentified	= None
	self.numberOfClusterSeeds	= None
	self.tempFilesFolder	= None
	self.skipReadCounting	= None
	self.handlePosition	= None
	self.sortFormat		= None
	self.trimmingRead1	= None
	self.trimmingRead2	= None
	self.minReadCountPerBead= None
	self.minReadCountPerConsensus= None
	self.minReadPopSupportConsensus=None
	self.minConsensusClusteringIdentity=None
	self.primerset		= None
	self.blastDb		= None
	self.gidatabase		= None
	self.minBlastIdentity	= None
	self.minBlastCoverage	= None
	self.shortTimeJobs	= None
	self.jobName		= None
	self.mostCommonToShow	= None
	self.subSpecies		= False
	self.skipPrevotella	= False
	self.readLength		= None
	self.allowedAllelLevelVariation = None
	self.rdpConfidence	= None

	# for each run
	self.cmd		= cmd
	self.stop		= stop
	self.skip		= skip
	self.random		= random
	self.sortFormat		= 'fq'
	self.read_count_per_barcode_cluster_cutoff = 1
	
	
	if cmd == 'init':
	    self.logfile = self.init_logfile
	    self.outfile = self.init_outfile
	elif cmd == 'setVariables':
	    self.logfile = self.setVariables_logfile
	    self.outfile = self.setVariables_outfile
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
	elif cmd == 'clone':
	    self.logfile = self.clone_logfile
	    self.outfile = self.clone_outfile
	elif cmd == 'compare':
	    self.logfile = self.compare_logfile
	    self.outfile = self.compare_outfile
	elif cmd == 'sbatch':
	    self.logfile = self.sbatch_logfile
	    self.outfile = self.sbatch_outfile
	elif cmd == 'makegraphs':
	    self.logfile = self.makegraphs_logfile
	    self.outfile = self.makegraphs_outfile
	elif cmd == 'classifymeta':
	    self.logfile = self.classifymeta_logfile
	    self.outfile = self.classifymeta_outfile

    def getreads2process(self, ):
	
	self.logfile.write('Getting readcounts ...\n')
	total = 0
	for i in range(len(self.readCountsList)):
	    rc = self.readCountsList[i]
	    self.logfile.write(self.infilesDictionary['r1'][i]+' -> '+str(rc) +' reads.\n')
	    total += rc
	self.logfile.write(str(total)+' read pairs in fastq files.\n');
	
	# calculate the number of reads to process
	self.reads2process = total
	if self.skip: 	self.reads2process -= self.skip
	if self.stop: 	self.reads2process = self.stop
	if self.random:	self.reads2process = self.random

    def set(self,varname, value, setDefaults=False ):
	if varname == 'chandlemissmatch':	self.maxHandleMissMatch	= value
	elif varname == 'barcodemissmatch':	self.maxBeadBarcodeMissMatch	= value
	elif varname == 'clustercount':		self.numberOfBarcodeClustersIdentified	= value
	elif varname == 'numberofseeds':	self.numberOfClusterSeeds	= value
	elif varname == 'sortformat':		self.sortFormat		= value
	elif varname == 'read_count_per_barcode_cluster_cutoff':self.read_count_per_barcode_cluster_cutoff = value
	else: raise ValueError

    def setMany(self,indata, setDefaults=False ):
	pass

    def load(self, ):
	self.config = initiate_file(self.config, self.logfile , mode='r')
	for line in self.config:
		if line[0] == '#': continue
		[varName, varValue, varComment] = line.rstrip().split('\t')
		try:	self.__dict__[varName] = eval(varValue)
		except: self.__dict__[varName] = varValue
		##general
		#if varName == 'absolutePath' :				self.absolutePath	= varValue.rstrip()
		#if varName == 'tempFilesFolder' :			self.tempFilesFolder	= varValue
		#if varName == 'infilesDictionary' :			self.infilesDictionary 	= eval(varValue)
		#if varName == 'readCountsList' :			self.readCountsList 	= eval(varValue)
		#if varName == 'skipReadCounting':			self.skipReadCounting	= eval(varValue)
		#if varName == 'jobName':				self.jobName		= varValue
		##bc clustering
		#if varName == 'numberOfClusterSeeds' :			self.numberOfClusterSeeds = eval(varValue.rstrip())
		#if varName == 'numberOfBarcodeClustersIdentified' :	self.numberOfBarcodeClustersIdentified 	= eval(varValue.rstrip())
		#if varName == 'maxBeadBarcodeMissMatch' :		self.maxBeadBarcodeMissMatch	= eval(varValue)
		#if varName == 'maxHandleMissMatch' :			self.maxHandleMissMatch	= eval(varValue)
		#if varName == 'handlePosition' :			self.handlePosition	= varValue
		##sorting
		#if varName == 'sortFormat' :				self.sortFormat		= varValue
		#if varName == 'trimmingRead1' :				self.trimmingRead1		= eval(varValue)
		#if varName == 'trimmingRead2' :				self.trimmingRead2		= eval(varValue)
	self.config.close()
	if self.primerset: self.primerset = open(self.primerset,'r')
	if type(self.jobName) == long: self.jobName = str(self.jobName)+'L'
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
	    self.absolutePath = os.path.abspath(self.path)
	else:
	    self.config = initiate_file(self.config, self.logfile, mode='ow')

	if type(self.primerset) == file: self.primerset = self.primerset.name

	self.logfile.write('Writing settings to config file ...\n')
	self.config.write(
		'#VariableName\t#Value\t#Comment\n'+

		'#\n# GENERAL VARIABLES:\n#\n'+
		'path'+					'\t'	+str(self.path)+		'\t'+	'# Relative analysis path'+	'\n'+
		'absolutePath'+				'\t'	+str(self.absolutePath)+	'\t'+	'# Absolute analysis path'+	'\n'+
		'tempFileFolder'+			'\t'	+str(self.tempFilesFolder)+	'\t'+	'# Folder used to store temporary files, host dependent, might change for each run'+	'\n'+
		'skipReadCounting'+			'\t'	+str(self.skipReadCounting)+	'\t'+	'# Bool skip counting reads in fastq files, debug option for speed'+	'\n'+
		'jobName'+				'\t'	+str(self.jobName)+		'\t'+	'# Name of job in sbatch files etc'+	'\n'+
		'infilesDictionary'+			'\t'	+str(self.infilesDictionary)+	'\t'+	'# Dictionary storing locations of pairs of input fastqfiles'+	'\n'+
		'readCountsList'+			'\t'	+str(self.readCountsList)+	'\t'+	'# List with the read pair count within each fastq file pair'+	'\n'+
		'readLength'+				'\t'	+str(self.readLength)+		'\t'+	'# Limits the readlengths to this value, if any read is longer it will be trimmed in the 3prime end, if read is shorter a warning will be displayed.'+	'\n'+

		'#\n# SETTINGS FOR BEAD BARCODES CLUSTERING AND IDENTIFICATION:\n#\n'+
		'numberOfClusterSeeds'+			'\t'	+str(self.numberOfClusterSeeds)+'\t'+	'# Number of sequences to use as seeds during clustering to identify bead barcodes '+	'\n'+
		'numberOfBarcodeClustersIdentified'+	'\t'	+str(self.numberOfBarcodeClustersIdentified)+	'\t'+	'# Number of Beads identified during clustering of barcode sequences'+	'\n'+
		'maxBeadBarcodeMissMatch'+		'\t'	+str(self.maxBeadBarcodeMissMatch)+	'\t'+	'# Number off missmatches allowed in clustering of bead barcodes'+	'\n'+
		'maxHandleMissMatch'+			'\t'	+str(self.maxHandleMissMatch)+	'\t'+	'# Number of missmatches allaowed when identifyng the handle'+	'\n'+
		'handlePosition'+			'\t'	+str(self.handlePosition)+	'\t'+	'# Position of handle, if set it overrides the sequence based identification'+	'\n'+

		'#\n# SETTINGS FOR SORTING FAST(A/Q) READS TO IDENTIFIED BEADS:\n#\n'+
		'sortFormat'+				'\t'	+str(self.sortFormat)+		'\t'+	'# File format to store sorted reads to'+	'\n'+
		'trimmingRead1'+			'\t'	+str(self.trimmingRead1)+	'\t'+	'# Number of bases to trim from read 1 during sorting (and in all downstream steps)'+	'\n'+
		'trimmingRead2'+			'\t'	+str(self.trimmingRead2)+	'\t'+	'# Number of bases to trim from read 2 during sorting (and in all downstream steps)'+	'\n'+

		'#\n# SETTINGS FOR CLUSTERING OF READS TO CONSENSUS SEQUENCES AND AMPLICONS:\n#\n'+
		'minReadCountPerConsensus'+		'\t'	+str(self.minReadCountPerConsensus)+	'\t'+	'# minimum number of reads needed for a consensus sequence to be considered as a variant of an amplicon'+	'\n'+
		'minReadPopSupportConsensus'+		'\t'	+str(self.minReadPopSupportConsensus)+	'\t'+	'# minimum %support of read pop needed for a consensus sequence to be considered as a variant of an amplicon'+	'\n'+
		'minConsensusClusteringIdentity'+	'\t'	+str(self.minConsensusClusteringIdentity)+'\t'+	'# minimum identity for two reads to cluster as one consensus sequence'+	'\n')
	if type(self.primerset) == file:
		self.config.write(
		'primerset'+				'\t'	+str(self.primerset.name)+		'\t'+	'# the primer set to be used when identifying amplicon supporting reads'+	'\n')
	else:
		self.config.write(
		'primerset'+				'\t'	+str(self.primerset)+		'\t'+	'# the primer set to be used when identifying amplicon supporting reads'+	'\n')
	self.config.write(
		'allowedAllelLevelVariation'+		'\t'	+str(self.allowedAllelLevelVariation)+'\t'+'# allowed "fraction variation" of allele readcount compared to most represented allele (readcount/mostrepresented >= 1-fraction).'+	'\n'+


		'#\n# SETTINGS FOR ALIGNMENT OF CONSENSUS SEQUENCES TO GENOMES:\n#\n'+
		'blastDb'+				'\t'	+str(self.blastDb)+		'\t'+	'# database to be used for the alignment of amplicon variants'+	'\n'+
		'gidatabase'+				'\t'	+str(self.gidatabase)+		'\t'+	'# dictionary of gi id to organism name mappings'+	'\n'+
		'minBlastIdentity'+			'\t'	+str(self.minBlastIdentity)+	'\t'+	'# minimum blast identity to consider the blast hit'+	'\n'+
		'minBlastCoverage'+			'\t'	+str(self.minBlastCoverage)+	'\t'+	'# minimum alignment length coverage to consider the blast hit'+	'\n'+
		'mostCommonToShow'+			'\t'	+str(self.mostCommonToShow)+	'\t'+	'# the number of most common genomes in hitlists for >1 defined all mono clusters'+	'\n'+
		'subSpecies'+				'\t'	+str(self.subSpecies)+		'\t'+	'# flag for showing the subspecies information or not'+	'\n'+
		'skipPrevotella'+			'\t'	+str(self.skipPrevotella)+	'\t'+	'# flag for skipping all Prevotella Hits or not'+	'\n'+
		'rdpConfidence'+			'\t'	+str(self.rdpConfidence)+	'\t'+	'# confidence treshold for rdp claissification'+	'\n'+
		#''+			'\t'	+str(self)+		'\t'+	'# '+	'\n'+
		'# A None value usually means that the variable is not yet set.\n'
		)
	#for varName in self.__dict__.keys(): print varName, str(self.__dict__[varName])
	self.config.close()
	self.config = self.config.name

    def loadPrimers(self, ):
	# set primerpairs
	#from SEAseqLib.mainLibrary import PrimerPair
	self.logfile.write('Loading primer pairs from '+self.primerset.name+'.\n')
	self.primerpairs = {}
	for line in self.primerset:
	    if line[0] != "#":
		line = line.rstrip().split('\t')
		try:self.primerpairs[line[0]] = PrimerPair(line[1],line[2],line[3])
		except IndexError: print line;print 'ERROR: Wrong number of columns in primerset file.'; sys.exit()

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
		if config.handlePosition:
			handle_start	= int(config.handlePosition.split(':')[0])
			handle_end	= int(config.handlePosition.split(':')[1])
		else:[handle_start, handle_end] = self.matchHandle(handle, config, self.r1)
		self.handle_start = handle_start
		self.handle_end   = handle_end
		return 0

	def identifyIllumina(self, config):
		#handle = sequence('illuminaUniversal','AGATCGGAAGAGC','AGATCGGAAGAGC')
		config.maxHandleMissMatch = 3
		#handle = sequence('illumina','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
		handle = sequence('illumina','AGATCGGAAGAGCACACGTCT','AGATCGGAAGAGCACACGTCT')
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r1)
		if handle_start: self.r1.illuminaadapter = True
		else: self.r1.illuminaadapter = False
		
		#handle = sequence('illumina','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
		handle = sequence('illumina','AGATCGGAAGAGCGTCGTGT','AGATCGGAAGAGCGTCGTGT')
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
		    self.missMatchesInTheHandle = 0
		
		elif config.maxHandleMissMatch:
			mindist = [10000,-1]
			for i in range(len(read.seq)):
			    
				if i+len(handle.seq) <= len(read.seq):
					dist = matchfunk(handle.seq,read.seq[i:i+len(handle.seq)])
				else: dist = 1000
				
				if dist < mindist[0]: mindist =[dist,i]

			if mindist[0] < config.maxHandleMissMatch:
				handle_start = i
				handle_end = i+len(handle.seq)
				self.missMatchesInTheHandle = mindist[0]
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
		return 0
	
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
		maxdist = config.maxBeadBarcodeMissMatch
		matchfunc = hamming_distance

		percentages={}
		for bc, count in self.barcodes.iteritems():
			percentage = round(100*float(count)/self.readcount,4)
			try: percentages[percentage].append(bc)
			except KeyError:percentages[percentage] = [bc]
		highest = []
		perc = percentages.keys()
		perc.sort(reverse=True)
		last_highest = 0
		#print perc
		while len(highest) < config.numberOfClusterSeeds:
			try:
				for bc in percentages[perc[0]]: highest.append(bc)
			except IndexError: pass
			perc=perc[1:]
			if len(highest) <= last_highest:
			    config.logfile.write('WARNING: could only find '+str(len(highest))+' uniqe barcodes, using all as cluster centers.\n')
			    break
			last_highest = len(highest)
		if len(highest) == 0:
			config.logfile.write('ERROR: could not find any uniqe barcodes, are you sure this dataset is correct, exciting.\n')
			import sys
			sys.exit(1)
		import os
		if os.path.exists(config.path+'/predetermined_cluster_centers.fa') and os.path.islink(config.path+'/predetermined_cluster_centers.fa'): os.unlink(config.path+'/predetermined_cluster_centers.fa')
		tempfile = open(config.path+'/predetermined_cluster_centers.fa','w')
		for bc in highest: tempfile.write('>'+bc+' '+str(self.barcodes[bc])+' readpairs\n'+bc+'\n')
		tempfile.close()

		if os.path.exists(config.path+'/raw_barcode_sequences.fa') and os.path.islink(config.path+'/raw_barcode_sequences.fa'): os.unlink(config.path+'/raw_barcode_sequences.fa')
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
		config.logfile.write('starting '+' '.join(['dnaclust','--similarity',str(1-(float(config.maxBeadBarcodeMissMatch)/15)),'--input-file',config.path+'/raw_barcode_sequences.fa','-t',str(multiprocessing.cpu_count()),'--predetermined-cluster-centers',config.path+'/predetermined_cluster_centers.fa'])+'\n')
		dnaclust =               subprocess.Popen(['dnaclust','--similarity',str(1-(float(config.maxBeadBarcodeMissMatch)/15)),'--input-file',config.path+'/raw_barcode_sequences.fa','-t',str(multiprocessing.cpu_count()),'--predetermined-cluster-centers',config.path+'/predetermined_cluster_centers.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

		#config.logfile.write('ok done now I\'ll just print and plot some info ... then done ... for real!\n\n')
		config.outfile.write(str( cc)+' clusters whereof '+str(counter)+' has more than one read\n\n')

		if os.path.exists(config.path+'/cluster.graphStats') and os.path.islink(config.path+'/cluster.graphStats'): os.unlink(config.path+'/cluster.graphStats')
		f = open(config.path+'/cluster.graphStats','w')
		f.write(str(reads_in_clusters))
		f.close()

		#temp_x=reads_in_clusters.keys()
		#temp_x.sort()
		#y=[];x =[]
		#for i in xrange(max(temp_x)+1):
		#	x.append(i)
		#	try: y.append(reads_in_clusters[i])
		#	except KeyError: y.append(0)
		#x=x
		#y=y
		#import numpy as np
		#import matplotlib.pyplot as plt
		#for scale in [[0,5000,0,20],[0,1000,0,20]]:
		#	plt.figure()
		#	plt.axis(scale)
		#	#plt.xlabel('Total Number of Reads per Barcode Cluster')
		#	#plt.ylabel('Number of Clusters')
		#	#pos = np.arange(len(x))
		#	#width = 1.0     # gives histogram aspect to the bar diagram
		#	#ax = plt.axes()
		#	#ax.set_xticks(pos + (width / 2))
		#	#ax.set_xticklabels(x,rotation='horizontal')
		#	#plt.bar(pos, y, width, color='r')
		#	##plt.show()
		#	##plt.savefig(pp,format='pdf',bbox_inches=0)
		#	plt.plot(x,y)
		#	plt.suptitle('Y is Total Number of Clusters with X Reads Pairs per Barcode Cluster.', fontsize=12)
		#	plt.savefig(config.path+'/read_pairs_per_barcode_cluster.x_scale_'+str(scale[0])+'-'+str(scale[1])+'.y_scale_'+str(scale[2])+'-'+str(scale[3])+'.pdf')
		#	plt.close()
		config.logfile.write( 'done\n')

		self.clusters = clusters
		del clusters

		return

	def reducebarcodes2(self,config):
		""" Find most common barcodes in well ( > 10% ??), then try to place other barcodes to this cluster
		"""
		
		config.logfile.write( '\n' )
		
		config.minperc = 0
		maxdist = config.maxBeadBarcodeMissMatch

		barcodesByCount = {}
		for bc, count in self.barcodes.iteritems():
		
			try: barcodesByCount[count].append(bc)
			except KeyError:barcodesByCount[count] = [bc]
		
		counts = barcodesByCount.keys()
		counts.sort()
		
		if os.path.exists(config.path+'/raw_barcode_sequences.fa') and os.path.islink(config.path+'/raw_barcode_sequences.fa'): os.unlink(config.path+'/raw_barcode_sequences.fa')
		tempfile = open(config.path+'/raw_barcode_sequences.fa','w')
		for count in counts:
			for bc in barcodesByCount[count]:
				for i in range(count):
					tempfile.write('>'+bc+' '+str(i)+'of'+str(self.barcodes[bc])+'\n'+bc+'\n')
		tempfile.close()

		# alternatively we could use chdit, though homopolymers seems to be problematic
		# cd-hit-454 -i tasks/20130805.1_index10/raw_barcode_sequences.fa -o DELETEME -c 0.85 -g 1 -T 8 -gap -6 -gap-ext -2 -AS 2
		# cdhit-cluster-consensus DELETEME.clstr tasks/20130805.1_index10/raw_barcode_sequences.fa DELETEME.cons DELETEME.aln
		
		# cd-hit-454 -i DELETEME.raw.bc.2.fa -o DELETEME -c 0.85 -M 0 -g 1 -n 5 -T 0 -gap -15 -gap-ext -2 -AS 2
		# cdhit-cluster-consensus DELETEME.clstr DELETEME.raw.bc.2.fa DELETEME.cons DELETEME.aln
		
		import subprocess
		from cStringIO import StringIO
		import time
		import multiprocessing
		tempo = time.time()
		command = ['dnaclust','--similarity',str(1-(float(config.maxBeadBarcodeMissMatch)/15)),'--input-file',config.path+'/raw_barcode_sequences.fa','-t',str(multiprocessing.cpu_count()),'--predetermined-cluster-centers',config.path+'/predetermined_cluster_centers.fa']		

		### COPY PASTE

		# Cluster Read pairs
		import subprocess
		from cStringIO import StringIO
		if indata.tempFileFolder:
			command = [	'cd-hit-454',
				'-i',indata.tempFileFolder+'/SEAseqtemp/temporary.'+str(self.id)+'.fa',
				'-o',indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.fa',
				'-g','1',#possibly change to 0 for more accurate though slower progress?? or other way around
				'-c',str(config.minConsensusClusteringIdentity/100.0)
				]
		else:	command = [	'cd-hit-454',
				'-i',config.path+'/sortedReads/temporary.'+str(self.id)+'.fa',
				'-o',config.path+'/sortedReads/cluster.'+str(self.id)+'.fa',
				'-g','1',#possibly change to 0 for more accurate though slower progress?? or other way around
				'-c',str(config.minConsensusClusteringIdentity/100.0)
				]

		cdhit = subprocess.Popen(
			command,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE )
		cdhit_out, errdata = cdhit.communicate()
		if cdhit.returncode != 0:
			print 'cmd: '+' '.join( command )
			print 'cd-hit cluster='+str(cluster.id)+' view Error code', cdhit.returncode, errdata
			sys.exit()
	
		# Build consensus sequences for read pair clusters
		if indata.tempFileFolder:
			command = ['cdhit-cluster-consensus',
				indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.fa.clstr',
				indata.tempFileFolder+'/SEAseqtemp/temporary.'+str(self.id)+'.fa',
				indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.consensus',
				indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.aligned'
				]
		else:	command =['cdhit-cluster-consensus',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.fa.clstr',
				config.path+'/sortedReads/temporary.'+str(self.id)+'.fa',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.consensus',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.aligned'
				] 
		ccc = subprocess.Popen(
			command,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE )
		ccc_out, errdata = ccc.communicate()
		if ccc.returncode != 0:
			print 'cmd: '+' '.join( command )
			print 'cluster='+str(self.id)+' cdhit-cluster-consensus view Error code', ccc.returncode, errdata
			print ccc_out
			sys.exit()

		### COPY PASTE

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

		config.outfile.write(str( cc)+' clusters whereof '+str(counter)+' has more than one read\n\n')

		if os.path.exists(config.path+'/cluster.graphStats') and os.path.islink(config.path+'/cluster.graphStats'): os.unlink(config.path+'/cluster.graphStats')
		f = open(config.path+'/cluster.graphStats','w')
		f.write(str(reads_in_clusters))
		f.close()

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
	def definedampliconcount(self):
		tmp_counter = 0
		for amplicon in self.amplicons.values():
			if amplicon.allelecount >= 1: tmp_counter += 1
		return tmp_counter

	def createtempfile(self, config, indata, verb=True):

		tmpcounter = 0
		if verb: output = 'PAIRS:\n'
		
		tofilestr = ''
		seqsForFasta = {}
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
			if verb: output+='---\t---\tillumina adapter in read pair\t'+pair.r1.seq +' '+ pair.r2.seq+'\n';
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
		    if pair.p1 == '???' or pair.p2 == '???' or pair.p1 != pair.p2 or len(pair.matchingprimerpairs) != 1:
			if pair.p1 == '???' and pair.p1 == pair.p2:	pair.primererror = 'primers-not-identifiable'
			if pair.p1 != pair.p2:				pair.primererror = 'fwd-rev-primer-missmatch'
			if len(pair.matchingprimerpairs) > 1:		pair.primererror = 'more-than-one-pair-match'
			if verb:
			    output += pair.primererror+'\t'+pair.r1.seq +' '+ pair.r2.seq+'\n'
			continue
	
		    #look for unexpected primer sequences
		    # this part could be improved and more verbose when searching for reverse complement primers
		    import re
		    matched_primerpair = config.primerpairs[pair.p1]
		    fwd_in_r2 = None
		    rev_in_r1 = None
		    fwdcount = 0
		    revcount = 0
		    other_fwd_in_any = None
		    other_rev_in_any = None
		    for name, primerpair in config.primerpairs.iteritems():
			if primerpair.name == matched_primerpair.name: # could imrpove here by searching for revcomp of fwd and rev though should be carefull with finding correct primers in small inserts, could be tagged as small insert amplicon?
				fwd_in_r2  = re.search(     matched_primerpair.fwdReStr,pair.r2.seq) #fwd in read 2
				rev_in_r1  = re.search(     matched_primerpair.revReStr,pair.r1.seq) #rev in read1
				fwdcount = len(re.findall(  matched_primerpair.fwdReStr,pair.r1.seq))
				revcount = len(re.findall(  matched_primerpair.revReStr,pair.r2.seq))
				if fwd_in_r2: 		pair.primererror = primerpair.name+'_'+'fwd_in_r2'+''.join(['-' for i in xrange(len('strange-primerpair-combo')-len(primerpair.name+'_'+'fwd_in_r2'))])
				if rev_in_r1: 		pair.primererror = primerpair.name+'_'+'rev_in_r1'+''.join(['-' for i in xrange(len('strange-primerpair-combo')-len(primerpair.name+'_'+'rev_in_r1'))])
				if revcount != 1: 	pair.primererror = primerpair.name+'_'+'revcount='+str(revcount)+''.join(['-' for i in xrange(len('strange-primerpair-combo')-len(primerpair.name+'_'+'revcount='+str(revcount)))])
				if fwdcount != 1: 	pair.primererror = primerpair.name+'_'+'fwdcount='+str(fwdcount)+''.join(['-' for i in xrange(len('strange-primerpair-combo')-len(primerpair.name+'_'+'fwdcount='+str(fwdcount)))])
			else: # could imrpove verbosity here by telling if we found "reverse complementery" and in which read
				other_fwd_in_any = re.search( primerpair.fwdReStr,   	pair.r1.revcomp().seq + 'NNNNN' + pair.r1.seq + 'NNNNN' + pair.r2.seq + 'NNNNN' + pair.r2.revcomp().seq) # fwd in any read
				other_rev_in_any = re.search( primerpair.revReStr,      pair.r1.revcomp().seq + 'NNNNN' + pair.r1.seq + 'NNNNN' + pair.r2.seq + 'NNNNN' + pair.r2.revcomp().seq) # rev in any read
				if other_fwd_in_any: 	pair.primererror = primerpair.name+'_fwd_in_any'+''.join(['-' for i in xrange(len('strange-primerpair-combo')-len(primerpair.name+'_fwd_in_any'))])
				if other_rev_in_any: 	pair.primererror = primerpair.name+'_rev_in_any'+''.join(['-' for i in xrange(len('strange-primerpair-combo')-len(primerpair.name+'_rev_in_any'))])
		    if pair.primererror:
			if verb:
				#pair.primererror = 'strange-primerpair-combo'
				output += pair.primererror+'\t'+pair.r1.seq +' '+ pair.r2.seq+'\n'
			continue
	
		    # Add sequences to output
		    if verb: output += '                        \t'+pair.r1.seq +' '+ pair.r2.seq+'\n'
		    
		    # prepare for clustering by writing read pair to tempfile and saving temporary id number and mappinf header to pair-object
		    tem_seq = pair.r1.seq[pair.handle_end:][len( config.primerpairs[pair.p1].fwd )+1:]+'NNNNNNNNNN'+pair.r2.revcomp().seq[:-(len(    config.primerpairs[pair.p1].rev   )+1)]
		    seqsForFasta[pair.id] = tem_seq
		    #tofilestr += '>'+str(pair.id)+'\n'+ tem_seq +'\n'
		
		#sort sequences by number of occurances of uniqe species
		import operator
		temporaryDict = {}
		for pairID, temporarySequence in seqsForFasta.iteritems():
			try: temporaryDict[temporarySequence].append(pairID)
			except KeyError: temporaryDict[temporarySequence] = [pairID]
		temporaryDict2 = {}
		for temporarySequence, pairIDList in temporaryDict.iteritems():
			try: temporaryDict2[len(pairIDList)] += pairIDList
			except KeyError: temporaryDict2[len(pairIDList)] = pairIDList
		for numberOfUniqeSequences, pairIDList in sorted(temporaryDict2.iteritems(), key=operator.itemgetter(0))[::-1]:
			for pairID in pairIDList: tofilestr += '>'+str(pairID)+'\n'+ seqsForFasta[pairID] +'\n'
		
		if indata.tempFileFolder:
			import os
			try: os.mkdir(indata.tempFileFolder+'/SEAseqtemp')
			except: pass
			f = open(indata.tempFileFolder+'/SEAseqtemp/temporary.'+str(self.id)+'.fa','w')
		else: 	f = open(config.path+'/sortedReads/temporary.'+str(self.id)+'.fa','w')
		f.write(tofilestr)
		f.close()
		
		if verb:return output
	
	def clusterreadpairs(self, config, indata, verb=False):
		
		#check that there is data to work with
		if self.adaptercount+self.primererrors == self.readcount:
		    import os
		    if indata.tempFileFolder: os.remove( indata.tempFileFolder+'/SEAseqtemp/temporary.'+str(self.id)+'.fa')
		    else: os.remove( config.path+'/sortedReads/temporary.'+str(self.id)+'.fa' )
		    return 'All adapter and/or primer error.\n'#['ONLY JUNK',return_info]
	
		# Cluster Read pairs
		import subprocess
		from cStringIO import StringIO
		if indata.tempFileFolder:
			command = [	'cd-hit-454',
				'-i',indata.tempFileFolder+'/SEAseqtemp/temporary.'+str(self.id)+'.fa',
				'-o',indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.fa',
				'-g','1',#possibly change to 0 for more accurate though slower progress?? or other way around
				'-c',str(config.minConsensusClusteringIdentity/100.0)
				]
		else:	command = [	'cd-hit-454',
				'-i',config.path+'/sortedReads/temporary.'+str(self.id)+'.fa',
				'-o',config.path+'/sortedReads/cluster.'+str(self.id)+'.fa',
				'-g','1',#possibly change to 0 for more accurate though slower progress?? or other way around
				'-c',str(config.minConsensusClusteringIdentity/100.0)
				]

		cdhit = subprocess.Popen(
			command,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE )
		cdhit_out, errdata = cdhit.communicate()
		if cdhit.returncode != 0:
			print 'cmd: '+' '.join( command )
			print 'cd-hit cluster='+str(cluster.id)+' view Error code', cdhit.returncode, errdata
			sys.exit()
	
		# Build consensus sequences for read pair clusters
		if indata.tempFileFolder:
			command = ['cdhit-cluster-consensus',
				indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.fa.clstr',
				indata.tempFileFolder+'/SEAseqtemp/temporary.'+str(self.id)+'.fa',
				indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.consensus',
				indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.aligned'
				]
		else:	command =['cdhit-cluster-consensus',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.fa.clstr',
				config.path+'/sortedReads/temporary.'+str(self.id)+'.fa',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.consensus',
				config.path+'/sortedReads/cluster.'+str(self.id)+'.aligned'
				] 
		ccc = subprocess.Popen(
			command,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE )
		ccc_out, errdata = ccc.communicate()
		if ccc.returncode != 0:
			print 'cmd: '+' '.join( command )
			print 'cluster='+str(self.id)+' cdhit-cluster-consensus view Error code', ccc.returncode, errdata
			print ccc_out
			sys.exit()
	
		# output info from temporary files
		if verb:
			output = ''
			if indata.tempFileFolder:
				datacombos = [
				[indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.consensus.fasta','\nCD-HIT consensus sequences:\n'],
				[indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.fa.clstr','\nClustering details:\n'],
				[indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.aligned','\nAlignment details:\n']
				]
			else:	datacombos = [
				[config.path+'/sortedReads/cluster.'+str(self.id)+'.consensus.fasta','\nCD-HIT consensus sequences:\n'],
				[config.path+'/sortedReads/cluster.'+str(self.id)+'.fa.clstr','\nClustering details:\n'],
				[config.path+'/sortedReads/cluster.'+str(self.id)+'.aligned','\nAlignment details:\n']
				]
			for info in datacombos:
				filename , message = info
				f = open(filename)
				tmp = f.read()
				output += message
				output += tmp
				f.close()
			return output
		return ''

	def removetempfiles(self, config, indata):
		import os
		if indata.tempFileFolder:
			os.remove(indata.tempFileFolder+'/SEAseqtemp/temporary.' + str(self.id) + '.fa')
			os.remove(indata.tempFileFolder+'/SEAseqtemp/cluster.'   + str(self.id) + '.fa')
			os.remove(indata.tempFileFolder+'/SEAseqtemp/cluster.'   + str(self.id) + '.fa.clstr')
			os.remove(indata.tempFileFolder+'/SEAseqtemp/cluster.'   + str(self.id) + '.consensus.fasta')
			os.remove(indata.tempFileFolder+'/SEAseqtemp/cluster.'   + str(self.id) + '.aligned')

		else: 
			os.remove(config.path + '/sortedReads/temporary.' + str(self.id) + '.fa')
			os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.fa')
			os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.fa.clstr')
			os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.consensus.fasta')
			os.remove(config.path + '/sortedReads/cluster.'   + str(self.id) + '.aligned')

	def loadconsensuses(self, config, indata):

	#	get identity from file
		if indata.tempFileFolder: f = open(indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.fa.clstr')
		else: f = open(config.path+'/sortedReads/cluster.'+str(self.id)+'.fa.clstr')
		data = f.read()
		f.close()
		for consensus_identities in data.split('>Cluster ')[1:]:
			consensusid = consensus_identities.split('\n')[0]
			if not consensusid: print 'WARNING!:',data.split('>Cluster ')
			consensus = Consensus(consensusid)
			for line in consensus_identities.split('\n')[1:]:
			    line=line.rstrip()
			    if not line: continue
			    read_id     = int(line.split('>')[1].split('.')[0])
			    identity = line.split('/')[-1]
			    pair = self.readpairbyheader[self.headerbyint[read_id]]
			    consensus.addreadpair(pair,identity)
			self.consensuses[consensus.id] = consensus
			if consensus.id == '': raise ValueError
			
	def consensusesToAmplicons(self, config):
		for consensusid, consensus in self.consensuses.iteritems():
			try: self.amplicons[consensus.type].addallele(consensus)
			except KeyError:
				self.amplicons[consensus.type] = Amplicon(consensus.type, consensus.primer)
				self.amplicons[consensus.type].addallele(consensus)

	def loadconsensusalignemnts(self, config, indata):
		# get alignments from file
		if indata.tempFileFolder: f = open(indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.aligned')
		else: f = open(config.path+'/sortedReads/cluster.'+str(self.id)+'.aligned')
		data = f.read()
		f.close()
		consensus = None
		for cluster_aln in data.split('===========================================\n')[1:]:
			tmpseqs = {}
			tmpcons = ''
			for part in cluster_aln.split('\n\n'):
				for line in part.split('\n'):
					line=line.rstrip()
					if not line: continue
					if line[0] == 'A':
						consensusid=line.split(' ')[-1].split(':')[0];
						continue
					read_id = line.split(':  +  ')[0]
					seq     = line.split(':  +  ')[1]
					if read_id == 'Consensus':
						tmpcons += seq.split(' ')[0];
					else:
						try: 		 tmpseqs[read_id] += seq
						except KeyError: tmpseqs[read_id]  = seq
			self.consensuses[consensusid].alignmentStr = tmpcons
			for readid, seq in tmpseqs.iteritems():
				try: self.consensuses[consensusid].readpairs[int(readid)].alignmentStr = seq
				except KeyError:
					config.logfile = open(config.logfile.name,'a',1)
					config.logfile.write('WARNING: cluster '+str(self.id)+', read '+str(readid).rstrip()+' is not supposed to be in consensus '+str(consensusid)+'.\n')
					config.logfile.close()
					print 'WARNING: cluster '+str(self.id)+', read '+str(readid)+' is not supposed to be in consensus '+str(consensusid)+'.'
	
	def loadconsensussequences(self, config, indata):
		##load singelton consensus sequences
		if indata.tempFileFolder: f = open(indata.tempFileFolder+'/SEAseqtemp/cluster.'+str(self.id)+'.consensus.fasta')
		else: f = open(config.path+'/sortedReads/cluster.'+str(self.id)+'.consensus.fasta')
		data = f.read()
		f.close()
		for consensus_identities in data.split('>')[1:]:
			if   consensus_identities[0] == 'c':
				consensusid = consensus_identities.split(' ')[0].split('_')[-1]
				readcount  = int(consensus_identities.split(' ')[2])
				assert self.consensuses[consensusid].readcount == readcount,	'ERRORSHMERROR 1'
			elif consensus_identities[0] == 's':
				consensusid = consensus_identities.split(' ')[0].split('_')[-1]
				read_id     = int(consensus_identities.split(' ')[1].split('\n')[0])
				assert self.consensuses[consensusid].readcount  == 1,		'ERRORSHMERROR 2'
				assert self.consensuses[consensusid].seedpairid == read_id,	'ERRORSHMERROR 3, '+str(self.consensuses[consensusid].seedpairid)+' != '+str(read_id)
			consensus = self.consensuses[consensusid]
			tmpseq = ''
			for line in consensus_identities.split('\n')[1:]:
			    line=line.rstrip()
			    tmpseq += line
			consensus.sequence = sequence(consensus.id,tmpseq,tmpseq)
			if self.consensuses[consensusid].readcount  == 1:
				self.consensuses[consensusid].alignmentStr = tmpseq
				self.consensuses[consensusid].readpairs[self.consensuses[consensusid].seedpairid].alignmentStr = tmpseq

	def getDefinedAmplicons(self, ):
		for amplicon in self.amplicons.values():
			if amplicon.allelecount >= 1: self.definedamplicons[amplicon.type] = amplicon	

	def createConsensusFasta(self,config,indata):
		output = ''
		if self.definedampliconcount == 0:
		    output += '\t0 defined amplicons, cluster '+str(self.id)+'\n'
		    #return [output, cluster]
		    return (output ,1)

		# creating the temporary fasta file
		if indata.tempFileFolder:
		    import os
		    try: os.mkdir(indata.tempFileFolder+'/SEAseqtemp')
		    except: pass
		    self.blastfile = open(indata.tempFileFolder+'/SEAseqtemp/blastinput.'+str(self.id)+'.fa','w')
		    self.rdpfile   = open(indata.tempFileFolder+'/SEAseqtemp/rdpinput.'+str(self.id)+'.fa','w')
		else:
		    self.blastfile= open(         config.path+'/sortedReads/blastinput.'+str(self.id)+'.fa','w')
		    self.rdpfile  = open(         config.path+'/sortedReads/rdpinput.'+str(self.id)+'.fa','w')
		output+= '\tmaking fasta for cluster '+str(self.id)+'\n'
		fastaentries = 0
		self.blastamplicons = {}
		self.rdpAmplicons = {}
		self.rdpConfidence = {}
		for amplicon in self.definedamplicons.values():
		    self.blastamplicons[amplicon.type] = {}
		    self.rdpAmplicons[amplicon.type] = {}
		    self.rdpConfidence[amplicon.type] = {}
		    for consensus in amplicon.goodalleles:
			try:
			    r1 = consensus.sequence.seq.split('NNNNNNNNNN')[0]
			    r2 = consensus.sequence.seq.split('NNNNNNNNNN')[1]
			except IndexError:
			    config.logfile = open(config.logfile.name,'a')
			    config.logfile.write('WARNING: Skipping amplicon '+amplicon.type+', consensus '+str(consensus.id)+' for cluster '+str(self.id)+' beacause splitting consensus sequence on N10 failed.\n')
			    output += 'WARNING: Skipping amplicon '+amplicon.type+', consensus '+str(consensus.id)+' for cluster '+str(self.id)+' beacause splitting consensus sequence on N10 failed.\n'
			    print 'WARNING: Skipping amplicon '+amplicon.type+', consensus '+str(consensus.id)+' for cluster '+str(self.id)+' beacause splitting consensus sequence on N10 failed.'
			    config.logfile.close()
			    continue
			self.blastfile.write('>'+amplicon.type+'|tempSep|'+str(consensus.id)+'|tempSep|r1\n'+r1+'\n'+
					'>'+amplicon.type+'|tempSep|'+str(consensus.id)+'|tempSep|r2\n'+r2+'\n')
			self.rdpfile.write('>'+amplicon.type+'|tempSep|'+str(consensus.id)+'|tempSep|\n'+r1+r2+'\n')
			self.blastamplicons[amplicon.type][str(consensus.id)] = {'r1':None,'r2':None}
			self.rdpAmplicons[amplicon.type][str(consensus.id)] = None
			self.rdpConfidence[amplicon.type][str(consensus.id)] = None
			fastaentries +=1
		self.blastfile.close()
		self.rdpfile.close()
		output+=  '\tfasta file created (cluster '+str(self.id)+')\n'
		
		# check that there were reads to align
		if fastaentries == 0:
		    self.blastHits = 'No fasta produced.'
		    import os
		    os.remove(self.blastfile.name)
		    os.remove(self.rdpfile.name)
		    output += '\tNo reads in fasta file, cluster '+str(self.id)+'\n'
		    #return [output, cluster]
		    return (output, 2)
		else:
			return (output, 0)

	def blastAllAmplicons(self,config,indata):
		output = ''
		from Bio.Blast.Applications import NcbiblastnCommandline
		from Bio.Blast import NCBIXML
		from cStringIO import StringIO
		import time
		#setting up blast
		BLAST_identity = config.minBlastIdentity
		database = config.blastDb
		max_target_seqs = 5000
		cline = NcbiblastnCommandline(query=self.blastfile.name, db=database ,evalue=0.001, outfmt=5, num_threads=8,perc_identity=BLAST_identity,max_target_seqs=max_target_seqs)#, out=infile+'.blastout')
		#cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn', out=infile+'.'+config.blastid+'.blastout')
		import time
		output+='\tStarting BLAST, cluster '+str(self.id)+'\n'
		starttime = time.time()
		blast_handle = cline.__call__()
		#output+=str(self.id)+' '+ str(blast_handle)[0:100].replace('\n','--NEWLINE--')+'\n'
		output+='\tBLAST search finished, cluster '+str(self.id)+', run time was '+str(round(time.time()-starttime))+'s.\n\n'
		
		import subprocess, sys
		from cStringIO import StringIO
		rdpProcess = subprocess.Popen(['java','-Xmx1g','-jar','/proj/b2014005/SEAseq/bin/RDP/RDPTools/classifier.jar','classify',self.rdpfile.name,'--outputFile',''+self.rdpfile.name+'.out','--format','allrank'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		rdpProcessOut, errdata = rdpProcess.communicate()
		if rdpProcess.returncode != 0: print 'RDP view Error code', rdpProcess.returncode, errdata; sys.exit()
		rdpProcessOut = StringIO(rdpProcessOut)
		rdpDataTmp = open(self.rdpfile.name+'.out','r')
		rdpData = rdpDataTmp.read()
		rdpDataTmp.close()
		
		#convert blast output to parsable handle
		blast_handle = StringIO(blast_handle[0])
		blast_handle.seek(0)
		records = NCBIXML.parse(blast_handle)
		del blast_handle
			    		
		# parse through the blast records and sort them by amplicon, allele and readnumber/part
		for blast_record in records:
	    
		    # get the information from the query header
		    amptype = blast_record.query.split('|tempSep|')[0]
		    allele =blast_record.query.split('|tempSep|')[1]
		    readnumber =blast_record.query.split('|tempSep|')[2]
		    
		    self.blastamplicons[amptype][allele][readnumber] = blast_record
		    if blast_record == None: output += amptype+' '+allele+' '+readnumber+'ajajajaj\n'; print 'ERROR: '+str(amptype)+' '+str(allele)+' '+str(readnumber)+'ajajajaj\n'
		
		#print 'CLUSTER',self.id
		for line in rdpData.rstrip().split('\n'):
		    #print line
		    line    = line.rstrip().split('\t')
		    header  = line[0]
		    empty   = line[1]
		    tmpData = line[2:]
		    #print header, tmpData
		    assert empty == '', 'ERROR, the first column for '+header+' is not empty (empty="'+empty+'"):\n'+rdpData
		    amptype   = header.split('|tempSep|')[0]
		    allele    = header.split('|tempSep|')[1]
		    rdpClassification = {}
		    rdpConfidence = {}
		    for i in range(0,len(tmpData),3):
			name = tmpData[i]
			rank = tmpData[i+1]
			confidence = float(tmpData[i+2])
			#print str(i)+'\t'+str(rank)+'\t'+str(name)+'\t'+str(confidence)
			if float(confidence) >= config.rdpConfidence/100.0: rdpClassification[rank] = name
			else: rdpClassification[rank] = 'LowConfidence'
			rdpConfidence[rank] = confidence
		    self.rdpAmplicons[amptype][allele] = rdpClassification
		    self.rdpConfidence[amptype][allele] = rdpClassification

		# remove temporary files
		import os
		os.remove(self.blastfile.name)
		os.remove(self.rdpfile.name)
		os.remove(self.rdpfile.name+'.out')
		
		return output
	
	def parseRdpAmplicons(self,config,indata):
		output = '\nRDP classification:\n'
		#store rdp hits for the random match estimation
		self.rdpHits = {}
		orgToClass = {}

		# parse through the blast result one amplicon at the time
		for amplicon in self.rdpAmplicons:
	    
		    broken = False
		    self.rdpHits[amplicon] = {}
		    consensuses = self.rdpAmplicons[amplicon]
		    output += '\tAmplicon: '+amplicon     +'\n'
	    
		    # for each consensus get the blast records for the two sequences (reads) and check for completeness of file
		    for consensus, rdpClassification in consensuses.iteritems():
			output += '\t\tConsensus/Allele/Variant '+str(consensus)     +'\n'
			for rank, value in rdpClassification.iteritems(): output +='\t\t\t'+rank+' is '+str(value)+'\n'# with '+str(round(100*values[1],2))+'%conf\n'
	    	    
		#if cluster is monoclonal for all defined and there are more than one amplicon check overlap and make a classification
		monoclonalAmpliconsArray = [self.amplicons[amplicon].monoclonal for amplicon in self.definedamplicons]
		if (self.definedampliconcount == monoclonalAmpliconsArray.count(True)):
		    
			classificationsInAllAmps = {}
			for amplicon, consensuses in self.rdpAmplicons.iteritems():
			    if len(consensuses) == 0: output += amplicon+' has no RDP-classifications.\n';continue
			if self.rdpAmplicons:
			    #print self.rdpAmplicons.keys(),self.rdpAmplicons[self.rdpAmplicons.keys()[0]].keys()
			    for firstAmpRank in self.rdpAmplicons[self.rdpAmplicons.keys()[0]][ self.rdpAmplicons[self.rdpAmplicons.keys()[0]].keys()[0] ]:
				firstAmpRankInAll = True
				for amplicon in self.rdpAmplicons:
					consensus = self.rdpAmplicons[amplicon].keys()[0]
					if firstAmpRank not in self.rdpAmplicons[amplicon][ consensus ]:
					    firstAmpRankInAll = False
					    print firstAmpRank, 'not in',self.rdpAmplicons[amplicon].keys()
				if firstAmpRankInAll:
					classificationsInAllAmps[firstAmpRank] = []
					rankValue = self.rdpAmplicons[amplicon][consensus][firstAmpRank]
					rankValueInAll = True
					lowConfidence = False
					for amplicon in self.rdpAmplicons:
						consensus = self.rdpAmplicons[amplicon].keys()[0]
						if rankValue != self.rdpAmplicons[amplicon][consensus][firstAmpRank]: rankValueInAll = False
						if rankValue == 'LowConfidence' or self.rdpAmplicons[amplicon][consensus][firstAmpRank] == 'LowConfidence' and 'LowConfidence' not in classificationsInAllAmps[firstAmpRank]: lowConfidence = True
					if rankValueInAll: classificationsInAllAmps[firstAmpRank].append(rankValue)
					elif lowConfidence: classificationsInAllAmps[firstAmpRank].append('LowConfidence')
					else: classificationsInAllAmps[firstAmpRank].append('MissMatch')
			toNext = {}
			if classificationsInAllAmps and self.definedampliconcount >= 2:
			    #print 'classificationsInAllAmps:',classificationsInAllAmps
			    output += '\tRDP-Classification overlaps:\n'
			    tmp1 = []
			    tmp2 = []
			    for rank in ['domain','phylum','class','order','family','genus']:
			    #for rank in ['domain','phylum','class','subclass','order','suborder','family','genus']:
				if rank in classificationsInAllAmps and classificationsInAllAmps[rank]:
					toNext[rank] = classificationsInAllAmps[rank]
					if len(classificationsInAllAmps[rank]) > 1:
						output += '\t\t'+rank +' is '+ ', '.join(classificationsInAllAmps[rank][:-1])+' or '+classificationsInAllAmps[rank][-1]+'.\n'
					else:  output += '\t\t'+rank +' is '+classificationsInAllAmps[rank][0]+'.\n'
				else: break
			    output += 'RDPtaxdata\t'+str(toNext)+'\n'
			elif self.definedampliconcount == 1:
			    output += 'Only one defined amplicon.\n'
			else:
			    output += 'No classification overlap found.\n'
			    #output += str([])+'\t'+str([])+'\n'
			    output += 'RDPtaxdata\t'+str(toNext)+'\n'

		return output

	def parseBlastAmplicons(self,config,indata):
		output = ''
		#store blast hits for the random match estimation
		self.blastHits = {}
		orgToClass = {}

		#load the mapping to organism name dictionary
		#from SEAseqLib.mainLibrary import gi2orgname
		#local_gi2org = {}
		#if config.gidatabase:
		#    tmp_file = open(config.gidatabase)
		#    tmp_string = tmp_file.read()
		#    tmp_file.close()
		#    local_gi2org = eval(tmp_string)

		# parse through the blast result one amplicon at the time
		for amplicon in self.blastamplicons:
	    
		    broken = False
		    self.blastHits[amplicon] = {}
		    consensuses = self.blastamplicons[amplicon]
		    output += '\tAmplicon: '+amplicon     +'\n'
	    
		    # for each consensus get the blast records for the two sequences (reads) and check for completeness of file
		    for consensus in consensuses:
			output += '\t\tConsensus/Allele/Variant '+str(consensus)     +'\n'
			try:
			    r1 = consensuses[consensus]['r1']
			    if r1 == None:
				output += '\t\tno read1 matches\n'
				continue
			except KeyError:
			    output += '\t\tConsensus cannot be read as read 1 sequence are not available from infile\n'
			    output += 'WARNING: Running on a uncomplete input file! Maybe the SEAseq meta program is still runnning?\n'
			    continue
			try:
			    r2 = consensuses[consensus]['r2']
			    if r2 == None:
				output += '\t\tno read2 matches\n'
				continue
			except KeyError:
			    output += '\t\tConsensus cannot be read as read 2 sequence are not available from infile\n'
			    output += 'WARNING: Running on a uncomplete input file! Maybe the SEAseq meta program is still runnning?\n'
			    continue
	    
			# get all organisms that part one maps to
			hitInfo = {}
			in_r1 = {}
			in_r2 = {}
			org2gi = {}
			for in_read, read, readStr in [[in_r1,r1,'r1'],[in_r2,r2,'r2']]:
				for alignment in read.alignments:
				    for hsp in alignment.hsps:
					perc_identity = float(hsp.identities) 	/	float(hsp.align_length)	*100
					perc_coverage = float(hsp.align_length)	/	float(read.query_letters)*100
					try: gi_number = alignment.title.split(' ')[1].split('|')[1]
					except IndexError:
					    import sys
					    gi_number = alignment.title.split(' ')[0].split('|')[1]
					    sys.stderr.write( 'This alignment title is not good: '+alignment.title +'\tsplitting on zero gi => '+gi_number+'.\n')
					#try:
					#    organism = local_gi2org[gi_number]
					#except KeyError:
					organism = gi2orgname(gi_number,lock=config.dbLock)
					    #try: organism = gi2orgname(gi_number)
					    #except urllib2.URLError: organism = alignment.title.split(' ')[1:]
					#    local_gi2org[gi_number] = organism
					if not config.subSpecies: organism = ' '.join(organism.split(' ')[:2])
					org2gi[organism] = gi_number
					import re
					if re.match('Prevotella',organism) and config.skipPrevotella: continue
					if perc_identity >= config.minBlastIdentity and perc_coverage >= config.minBlastCoverage and organism not in in_read:
					    in_read[organism] = alignment
					    try:
						hitInfo[organism][readStr]['pi'].append(perc_identity)
						hitInfo[organism][readStr]['pc'].append(perc_coverage)
						hitInfo[organism][readStr]['ss'].append(str(hsp.sbjct_start)+'-'+str(hsp.sbjct_start+hsp.align_length))
					    except KeyError:
						if readStr == 'r1': tmp = 'r2'
						else: tmp = 'r1'
						hitInfo[organism] = {readStr:{'pi':[perc_identity],'pc':[perc_coverage],'ss':[str(hsp.sbjct_start)+'-'+str(hsp.sbjct_start+hsp.align_length)]},tmp:{'pi':[],'pc':[],'ss':[]}}

			# get all organisms that both parts map to and print them also print if no overlap is found
			in_both_reads = []
			for organism in in_r1:
			    if organism in in_r2:
				in_both_reads.append(organism)
			for organism in in_both_reads:
			    output +=       '\t\t\t'+organism    +'\n'
			    output += '\t\t\t\tpart1:\tHIT#'+'\n\t\t\t\t\tHIT#'.join([str(i)+': identity='+str(hitInfo[organism]['r1']['pi'][i])+'%, coverage='+str(hitInfo[organism]['r1']['pc'][i])+'%, pos='+str(hitInfo[organism]['r1']['ss'][i]) for i in range(len(hitInfo[organism]['r1']['ss']))]) + '\n'
			    output += '\t\t\t\tpart2:\tHIT#'+'\n\t\t\t\t\tHIT#'.join([str(i)+': identity='+str(hitInfo[organism]['r2']['pi'][i])+'%, coverage='+str(hitInfo[organism]['r2']['pc'][i])+'%, pos='+str(hitInfo[organism]['r2']['ss'][i]) for i in range(len(hitInfo[organism]['r2']['ss']))]) + '\n'
			    self.blastHits[amplicon][organism] = hitInfo[organism]
			    orgToClass[organism] = getClassification(gi=org2gi[organism],lock=config.dbLock)
			if not in_both_reads:
			    output +=       '\t\t\tNo alignment supported by both reads with >='+str(config.minBlastIdentity)+'% identity and '+str(config.minBlastCoverage)+'% alignment length coverage'     +'\n'
	    
		#if cluster is monoclonal for all defined and there are more than one amplicon check overlap and make a classification
		monoclonalAmpliconsArray = [self.amplicons[amplicon].monoclonal for amplicon in self.definedamplicons]
		if len(self.blastHits) > 1 and (self.definedampliconcount == monoclonalAmpliconsArray.count(True)):
		    
		    self.hitReduction = { amplicon:len(organisms) for amplicon, organisms in self.blastHits.iteritems() }
		    
		    # find organisms present in all amplicons hitlists
		    self.organismsInAllAmplicons = []
		    for organism in self.blastHits[self.blastHits.keys()[0]]:
			inAll = True
			for amplicon, organisms in self.blastHits.iteritems():
				if organism not in organisms: inAll = False
			if inAll: self.organismsInAllAmplicons.append(organism)

		    # find classifications present within all organisms within hitlist for each amplicon (=consensus as they are monoclonal)
		    if self.hitReduction and self.hitReduction[self.hitReduction.keys()[0]] > 0 and self.hitReduction[self.hitReduction.keys()[0]] < sum([value for value in self.hitReduction.values()]):
			classifications = {}
			classificationsInAllAmps = {}
			for amplicon, organisms in self.blastHits.iteritems():
			    if len(organisms) == 0: output += amplicon+' has no BLAST-hits.\n';continue
			    classifications[amplicon] = {}
			    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
				    #try: 		value = orgToClass[organisms.keys()[0]][rank]
				    #except KeyError:value = 'NotSet'
				    #except IndexError:value = 'NotSet'
				    classifications[amplicon][rank] = []
				    for organism in organisms:
					try: value = orgToClass[organism][rank]
					except KeyError:value = 'NoInformationAvailable'
					if value not in classifications[amplicon][rank]: classifications[amplicon][rank].append(value)
			    #for organism in organisms:
			    #	for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
			    #		try: value = orgToClass[organism][rank]
			    #		except KeyError: value = 'NotSet'
			    #		if value != classifications[amplicon][rank]: classifications[amplicon][rank] = False
			    #output +=       '\tAmplicon '+amplicon+' have the following classifications for all hits:'     +'\n'
			    #for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
			    #	if classifications[amplicon][rank]: output += '\t\t'+rank + '\t' + classifications[amplicon][rank]+'\n'
			    #	else: output += '\t\t'+rank + '\tnot same for all\n'
			if classifications:
			    for firstAmpRank in classifications[ classifications.keys()[0] ]:
					    firstAmpRankInAll = True
					    for amplicon in classifications:
						    if firstAmpRank not in classifications[amplicon]:
							firstAmpRankInAll = False
							print firstAmpRank, 'not in',classifications[amplicon].keys()
					    if firstAmpRankInAll:
						    classificationsInAllAmps[firstAmpRank] = []
						    for rankValue in classifications[classifications.keys()[0]][firstAmpRank]:
							    rankValueInAll = True
							    for amplicon in classifications:
								    if rankValue not in classifications[amplicon][firstAmpRank]: rankValueInAll = False
							    if rankValueInAll: classificationsInAllAmps[firstAmpRank].append(rankValue)
			if classificationsInAllAmps:
			    output += '\tClassification overlaps:\n'
			    tmp1 = []
			    tmp2 = []
			    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
				if classificationsInAllAmps[rank]:
					if len(classificationsInAllAmps[rank]) > 1:
						output += '\t\t'+rank +' is '+ ', '.join(classificationsInAllAmps[rank][:-1])+' or '+classificationsInAllAmps[rank][-1]+'.\n'
					else:  output += '\t\t'+rank +' is '+classificationsInAllAmps[rank][0]+'.\n'
					#tmp1.append(classificationsInAllAmps[rank]); tmp2.append(rank);
				else: classificationsInAllAmps[rank] = ['NoMatchFound']
			    #output += '\t\t'+':\n'.join([', '.join(rankValues[:-1])+' or '+rankValues[-1] for rankValues in tmp1])+'\n'
			    #output += '('+':'.join(tmp2)+')\n'
			    #output += str(tmp1)+'\t'+str(tmp2)+'\n'
			    output += 'taxdata\t'+str(classificationsInAllAmps)+'\n'
			else:
			    output += 'No classification overlap found.\n'
			    #output += str([])+'\t'+str([])+'\n'
			    output += 'taxdata\t'+str(classificationsInAllAmps)+'\n'
		    
		#    # find level classification that are the same for all amplicons
		#    classInAllAmps = {}
		#    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
		#	try: value = classifications[self.blastHits.keys()[0]][rank]
		#	except KeyError: value = 'NotSet'
		#	classInAllAmps[rank] = value
		#    for amplicon in self.blastHits.keys():
		#	for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
		#		try:
		#			if classifications[amplicon][rank] != classInAllAmps[rank]:classInAllAmps[rank] = False
		#		except KeyError:classInAllAmps[rank] = False
		#    output +='\tAll amplicons have the following classifications:'     +'\n'
		#    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
		#	if classInAllAmps[rank]: output += '\t\t'+rank + '\t' + classInAllAmps[rank]+'\n'
		#	else: output += '\t\t'+rank + '\tnot same for all\n'	
		    
		    
		    self.hitReduction = { amplicon:len(organisms) for amplicon, organisms in self.blastHits.iteritems() }
		    self.hitReduction['inAll'] = len(self.organismsInAllAmplicons)
		    
		    output +=       '\tOrganims in all amplicons:'     +'\n'
		    for organism in self.organismsInAllAmplicons:
			output += '\t\t'+str(organism)     +'\n'
			
		    output += str(self.hitReduction)+'\n'
		    if   len(self.organismsInAllAmplicons) >1:   output += 'cluster is '+', '.join(self.organismsInAllAmplicons[:-1])+' or ' +self.organismsInAllAmplicons[-1] +'\n'
		    elif len(self.organismsInAllAmplicons) ==1:  output += 'cluster is '+self.organismsInAllAmplicons[-1] +'\n'
		#    if len(classInAllAmps) > 0:
		#	output += 'matching clasification: '
		#	tmp1 = []
		#	tmp2 = []
		#	for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
		#	    if classInAllAmps[rank]: tmp1.append(classInAllAmps[rank]); tmp2.append(rank);
		#	    else: break
		#	output += ':'.join(tmp1)+' '
		#	output += '('+':'.join(tmp2)+')\n'
		#	output += str(tmp1)+'\t'+str(tmp2)
		
		return output

class Amplicon(object):

	def __init__(self, amplicon_type, primerpair):
		self.allels = []#hold the consensus sequences produced
		self.type = amplicon_type
		self.primer = primerpair
		self.allelecount = None
		self.goodalleles = []

	@property
	def readcount(self):
		return sum([consensus.readcount for consensus in self.allels])
		
	def addallele(self,consensus ):
		if consensus.type == self.type: self.allels.append(consensus)
		else:
			import sys
			sys.stderr.write('ERROR: Amplicon type does not match the Consensus type.\n')
			raise ValueError

	def sortalleles(self):
		#sort alleles by number of reads/consesussequence

		import operator

		temporaryDict = {}
		for consensus in self.allels:
			try:		temporaryDict[consensus.readcount].append(consensus)
			except KeyError:temporaryDict[consensus.readcount] = [consensus]
		
		sortedAlleles = []
		for count, consensuses in sorted(temporaryDict.iteritems(), key=operator.itemgetter(0))[::-1]:
			for consensus in consensuses: sortedAlleles.append(consensus)
		
		self.allels = sortedAlleles
			
	def checkmono(self, config):
		output = ''
		output += '\t'+self.type+' '+str(self.readcount)+' reads in total.\n'
		self.allelecount = 0
		self.monoclonal = None
		
		self.sortalleles()
		
		mostRepresentedConsensus = self.allels[0]
		allowedAllelLevelVariation = config.allowedAllelLevelVariation #0.6 = 60% variation ie second allele most have readcount >= 40% of "mostRepresentedConsensus" allele readcount
		
		for consensus in self.allels:
			consensus.percentagesupport = 100*float(consensus.readcount)/float(self.readcount)
			if consensus.readcount > 1 and consensus.percentagesupport > 1:
				output += '\t\tConsensus '+consensus.id+' supported by '+str(round(consensus.percentagesupport,2))+'% of readpop ('+str(consensus.readcount)+' reads)\t'+str(round(100*float(consensus.readcount)/float(mostRepresentedConsensus.readcount),0))+'%\t'+consensus.sequence.seq+'\n'
			if consensus.readcount >= config.minReadCountPerConsensus:
				if  consensus.percentagesupport >= config.minReadPopSupportConsensus:
					if allowedAllelLevelVariation:
						if (float(consensus.readcount)/float(mostRepresentedConsensus.readcount)) >= (1.0-allowedAllelLevelVariation):
							self.allelecount += 1
							if consensus not in self.goodalleles: self.goodalleles.append(consensus)
					else:
						self.allelecount += 1
						if consensus not in self.goodalleles: self.goodalleles.append(consensus)
		
		if self.allelecount == 1:
			self.monoclonal = True
			output += '\tMonoclonal for '+self.type + '\n'
		elif self.allelecount >1:
			self.monoclonal = False
			output += '\tPolyclonal for '+self.type + '\n'
		else:   output += '\tNo "good allels" found for '+self.type + '\n'
		
		return output

class Consensus(object):

	def __init__(self, idnumber, amplicon_type=None, primerpair=None, sequence=None):
		self.readpairs = {}
		self.type = amplicon_type
		self.primer = primerpair
		self.id = idnumber
		self.sequence = sequence
		self.alignmentStr = None

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

	def alignmentoutput(self, config, verb = 1):
		#make output for alnignments
		output = ''
		output += 'Consensus number '+self.id+' from '+str(self.readcount)+' read pairs\t'+self.alignmentStr+'\n'
		withalnstr = ''
		withoutstr = ''
		pair = self.readpairs[self.seedpairid]
		try:			withalnstr += 'Read pair id = '+str(self.seedpairid)+'    \t'+pair.consensusidentity+'\t'+pair.p1+'\t'+pair.alignmentStr +'\n'
		except AttributeError:	withoutstr += 'Read pair id = '+str(self.seedpairid)+'    \t'+pair.consensusidentity+'\t'+pair.p1+'\tAlignmentStr not available.\n'
		for readid, pair in self.readpairs.iteritems():
			if readid == self.seedpairid: continue
			try:			withalnstr += 'Read pair id = '+str(pair.id)+'    \t'+pair.consensusidentity+'\t'+pair.p1+'\t'+pair.alignmentStr +'\n'
			except AttributeError:	withoutstr += 'Read pair id = '+str(pair.id)+'    \t'+pair.consensusidentity+'\t'+pair.p1+'\tAlignmentStr not available.\n'
		if verb < 2: withoutstr = ''
		if verb: return output + withalnstr + withoutstr + '\n'
	

if __name__ == "__main__": lib_main()
