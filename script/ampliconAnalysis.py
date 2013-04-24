#! /usr/bin/env python2.7
#! /bin/env python

from tnt.lib import *
import os
import time

MASTER = os.getppid()
C_HANDLE = sequence('c handle',"CTAAGTCCATCCGCACTCCT","CTAAGTCCATCCGCACTCCT")

def main():
#--------------------------MAIN PROGRAM-----------------------------
    indata = getindata()
    
    indata.logfile.write('Starting.\n')
    indata.logfile.write('Running script '+time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+'.\n')
    indata.logfile.write('Master process id='+str(MASTER)+'\n')
    indata.logfile.write('cmd: '+' '.join(sys.argv)+'\n')

    #deciding if to run multiprocessing or single process for debugging
    import multiprocessing
    if indata.debug: #single process // serial

	results=[] # create holder for processed reads
	indata.logfile.write('Part1: Running in debug mode identifying handles ...\n')
	
	progress = Progress(indata.reads2process, logfile=indata.logfile) # creates a progress "bar" thingy
	#itarate through reads and do the "magicFunction"
	with progress:
	    for tmp in getPairs(indata):
		progress.update()
		results.append(foreachread(tmp))
	
	indata.logfile.write('Part1: search finished\n')
	
    else: # multiple processes in parallel
	#create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=1000)
	results = WorkerPool.imap(foreachread,getPairs(indata),chunksize=100)

    # output log message about whats happening
    if not indata.debug: indata.logfile.write('Part1: running in multiproccessing mode using '+str(indata.cpus)+' processes  ...\n')
    else: indata.logfile.write('Part2: running the multiprocessing results handeling in serial\n')

    progress = Progress(indata.reads2process, logfile=indata.logfile)
    summary = SEAseqSummary()
    with progress:
	for pair in results:
	    progress.update()
	    summary.add(pair)
    if not indata.debug:WorkerPool.close()
    if not indata.debug:WorkerPool.join()
    
    print summary.part1()

    # Part1 - Parallel - each read pair:
    # check all reads and tryt to find "C handle"
    # and thereby also the barcode also try to identify the amplicon parts and any adapter sequences if present
    # output some initial statistics

    # Part2 - Serial:
    # some kind of clustering of barcodes with predetermined number of missmatches
    
    if True:
	summary.reducebarcodes(indata)
	f= open('clusters.tempfile','w')
	f.write(str(summary.clusters))
	f.close()
    else:
	summary.loadclusters('clusters.tempfile')    

    indata.logfile.write('Preparing for sorting of reads to clusters\n')
    cid_by_bc = {}
    outfiles = {}
    import os
    try: os.makedirs('temporary.cluster.files')
    except OSError:pass
    progress = Progress(len(summary.clusters),verb='minimal')
    min_reads = 100
    with progress:
	for cluster_id, infodist in summary.clusters.iteritems():
	    progress.update()
	    if int(infodist['total']) > min_reads:
		#outfiles[cluster_id] = open('temporary.cluster.files/'+str(cluster_id)+'.reads','w')
		#outfiles[cluster_id].close()
		for barcode in infodist['barcodes']:
		    cid_by_bc[barcode] = cluster_id
	    else: pass#print 'low read cluster'

    # Part3 - Parallel? - each barocde group: (assumes that we have monoclonal beads)
    # calculate statistics, level of clonality of amplicons,
    
    del summary
    
    if indata.debug: #single process // serial

	results=[] # create holder for processed reads
	indata.logfile.write('Part3: Running in debug mode identifying handles ...\n')
	
	progress = Progress(indata.reads2process, logfile=indata.logfile) # creates a progress "bar" thingy
	#itarate through reads and do the "magicFunction"
	with progress:
	    for tmp in getPairs(indata):
		progress.update()
		results.append(foreachread2(tmp))
	
	indata.logfile.write('Part3: search finished\n')
	
    else: # multiple processes in parallel
	#create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=1000)
	results = WorkerPool.imap(foreachread,getPairs(indata),chunksize=100)

    # output log message about whats happening
    if not indata.debug: indata.logfile.write('Part2: running in multiproccessing mode using '+str(indata.cpus)+' processes  ...\n')
    else: indata.logfile.write('Part4: running the multiprocessing results handeling in serial\n')

    tempfiles = {}
    progress = Progress(indata.reads2process, logfile=indata.logfile)
    with progress:
	for pair in results:
	    progress.update()
	    if pair.n15 and pair.n15.len == 15:# and not (pair.r1.illuminaadapter or pair.r2.illuminaadapter):
		try:
		    f = open('temporary.cluster.files/'+str(cid_by_bc[pair.n15.seq])+'.reads','a')
		    f.write('>'+pair.r1.header+'_r1\n'+pair.r1.seq+'\n>'+pair.r2.header+'_r2\n'+pair.r2.seq+'\n')
		    f.close()
		    del f
		    tempfiles[cid_by_bc[pair.n15.seq]] = True
		except KeyError: pass#print 'low read cluster read ie not printed'
    if not indata.debug:WorkerPool.close()
    if not indata.debug:WorkerPool.join()
    indata.logfile.write('reads sorted into cluster tempfiles\n')

    #Go through CLUSTERS

    if indata.debug: #single process // serial

	results=[] # create holder for processed reads
	for infile in tempfiles.keys():
	    results.append(foreach3(infile))
	
    else: # multiple processes in parallel
	#create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=100)
	results = WorkerPool.imap(foreach3,tempfiles.keys(),chunksize=5)

    monoclonal_clusts = 0
    tot_clust  = 0

    beadtypes = {}
    for [tot_reads, output, monoclonal,genome] in results:
	if tot_reads > min_reads:
	    print output
	    tot_clust+=1
	    if monoclonal: monoclonal_clusts+=1
	    try: beadtypes[genome] += 1
	    except KeyError: beadtypes[genome] = 1

    print 'We found',tot_clust,'clusters, with atleast',min_reads,'reads per cluster, out of theese were',str(round(100*float(monoclonal_clusts)/tot_clust,2)),'% monoclonal'
    total = sum([ count for beadtype,count in beadtypes.iteritems()])
    for beadtype,count in beadtypes.iteritems(): print beadtype, str(round(100*float(count)/total,2))+'% ('+str(count)+' beads)'

    import shutil
    shutil.rmtree('temporary.cluster.files/')
    import os
    os.remove('clusters.tempfile')
    os.remove('seed_bcs.tempfile')
    os.remove('raw_bcs.tempfile')
    # Part4 - Serial:
    # statistics summary
    
#--------------------------MAIN PROGRAM END-------------------------

#--------------------- Functions, Subroutines and classes --------------------	
def foreach3(infile):
    infile = 'temporary.cluster.files/'+str(infile)+'.reads'
    return classify_cluser(infile=infile)

def foreachread(tmp):

    # unpack info
    [pair, indata] = tmp
    
    # convert to SEAseq readpair
    pair = SEAseqpair(pair.header, pair.r1, pair.r2)
    
    pair.identify(C_HANDLE, indata)
    pair.getN15()
    
    return pair

def foreachread2(tmp):

    # unpack info
    [pair, indata] = tmp
    
    # convert to SEAseq readpair
    pair = SEAseqpair(pair.header, pair.r1, pair.r2)
    
    pair.identify(C_HANDLE, indata)
    pair.getN15()
    pair.identifyIllumina(indata)
    #if not (pair.r1.illuminaadapter or pair.r2.illuminaadapter): pass#pair.getSpecific(indata)
    
    return pair
    
def getindata():
    import argparse
    argparser = argparse.ArgumentParser(description='Analysis of SEAseq amplicon data.', formatter_class=argparse.RawTextHelpFormatter,)
    argparser.add_argument(	'--debug',		dest='debug', 			action='store_true', 			required=False,	default=False,	help='debug (run as regular single process python script).')
    argparser.add_argument(	'-skip',		dest='skip',	metavar='N',				type=int,	required=False,	default=0,	help='skip the first N read pairs in files (default 0).')
    argparser.add_argument(	'-stop',		dest='stop',	metavar='N',				type=int,	required=False,	default=0,	help='stop after N read pairs, set to 0 to disable (default 0).')
    argparser.add_argument(	'-r1',			dest='reads1',	metavar='FILE',				type=file,	required=True, 			help='indata "fastq"-file read1.')
    argparser.add_argument(	'-r2',			dest='reads2',	metavar='FILE',				type=file,	required=True,	default=None,	help='indata "fastq"-file read2.')
    argparser.add_argument(	'-p',			dest='cpus',	metavar='N',				type=int,	required=False,	default=1,	help='The number of processes to start (default 1).')
    argparser.add_argument(	'-o',			dest='outfile',	metavar='outfile',			type=str,	required=False,	default=False,	help='Print output to outfile (default stdout).')
    argparser.add_argument(	'-l',			dest='logfile',	metavar='logfile',			type=str,	required=False,	default=False,	help='Print log messages to logfile (default stderr).')
    argparser.add_argument(	'-random',		dest='n',	metavar='N',				type=int,	required=False,	default=0,	help='Use a random subset of N read pairs, this option is slower (default 0 = off). Can not be used in combination with "-skip" or "-stop"')
    argparser.add_argument(	'-hm',			dest='handlemm',metavar='N',				type=int,	required=False,	default=0,	help='Number off missmatches allowed in handle sequence (default 0)')
    indata = argparser.parse_args(sys.argv[1:])
    indata.selftest = False
    
    if indata.outfile: indata.outfile = open(indata.outfile, 'w',1)
    else: indata.outfile = sys.stdout
    
    if indata.logfile: indata.logfile = open(indata.logfile, 'w',1)
    else: indata.logfile = sys.stderr
    
    # get the readcount
    indata.logfile.write('Getting readcount ... ')
    indata.numreads=bufcount(indata.reads1.name)/4
    #indata.numreads=10000
    indata.logfile.write(str(indata.numreads)+' read pairs in fastq files.\n');

    # calculate the number of reads to process
    indata.reads2process = indata.numreads
    if indata.skip: indata.reads2process -= indata.skip
    if indata.stop: indata.reads2process = indata.stop
    if indata.n:    indata.reads2process = indata.n

    return indata



#####
#check if run or imported // call main() or not
#####
if __name__ == "__main__":
    main()
#END of script
