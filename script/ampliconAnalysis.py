#! /usr/bin/env python2.7
#! /bin/env python

from tnt.lib import *
import os
import time

MASTER = os.getppid()
C_HANDLE = sequence('c handle',"CTAAGTCCATCCGCACTCCT","CTAAGTCCATCCGCACTCCT")

def main():
#--------------------------MAIN PROGRAM-----------------------------
    import gc
    indata = getindata()
    
    indata.logfile.write('Starting.\n')
    indata.logfile.write('Running script '+time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+'.\n')
    indata.logfile.write('Master process id='+str(MASTER)+'\n')
    indata.logfile.write('cmd: '+' '.join(sys.argv)+'\n')

    if not indata.analyzeclust:

	
	if not indata.onlysort:
	    summary = cluster_barcodes(indata)

	## SORT READS TO BARCODE CLUSTERS
	if indata.onlysort:
	    summary = SEAseqSummary()
	    summary.loadclusters(indata.outfolder+'/clusters.tempfile')
	sortreads(summary.clusters,indata)
        if indata.onlysort: sys.exit()

    ## GO THROUGH EACH CLUSTER AND IDENTIFY THE READS
    if not indata.analyzeclust: indata.analyzeclust = indata.outfolder+'/temporary.cluster.files'
    
    blast_it(indata)

    ## CLEANUP
    if not (indata.keeptemp):
	import shutil
	shutil.rmtree(indata.outfolder+'/'+'temporary.cluster.files/')
	import os
	os.remove(indata.outfolder+'/'+'clusters.tempfile')
	os.remove(indata.outfolder+'/'+'seed_bcs.tempfile')
	os.remove(indata.outfolder+'/'+'raw_bcs.tempfile')

    indata.outfile.close()
    indata.logfile.write('All Done.\n')
    indata.logfile.close()
#--------------------------MAIN PROGRAM END-------------------------

#--------------------- Functions, Subroutines and classes --------------------	
def cluster_barcodes(indata):
    	    ##
	    ## FIND CHANDLE AND BARCODE SEQUENCE IN EACH READ PAIR
	    ##
	    indata.logfile.write('Part1: Clustering barcodes START\n')
	    #deciding if to run multiprocessing or single process for debugging
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
		import multiprocessing
		#create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
		WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=1000000)
		results = WorkerPool.imap_unordered(foreachread,getPairs(indata),chunksize=100)
	
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
	    
	    indata.outfile.write(str( summary.part1() )+'\n')
    
	    ##
	    ## CLUSTER BARCODES
	    ##
	
	    indata.logfile.write('Clustering the barcodes.\n')
	    summary.reducebarcodes(indata)
	    if indata.onlycluster:
		f= open(indata.outfolder+'/clusters.tempfile','w')
		f.write(str(summary.clusters))
		f.close()
		indata.logfile.write('Part1: Clustering barcodes END\n')
		sys.exit()
	    indata.logfile.write('Part1: Clustering barcodes END\n')
	    return summary	

def sortreads(clusters,indata):
	indata.logfile.write('Part2: Sorting reads to clusters START\n')
	indata.logfile.write('\nPreparing for sorting of reads to clusters\n')

	import os
	try: os.makedirs(indata.outfolder+'/temporary.cluster.files')
	except OSError:pass
	try: os.makedirs(indata.outfolder+'/temporary.cluster.files/barcodes')
	except OSError:pass

	cid_by_bc = {}
	progress = Progress(len(clusters),verb='minimal', logfile=indata.logfile)
	with progress:
	    for cluster_id, infodist in clusters.iteritems():
		progress.update()
		if int(infodist['total']) >= indata.mrc:
		    for barcode in infodist['barcodes']:
			cid_by_bc[barcode] = cluster_id
		else: pass#print 'low read cluster'

	clusters.clear()
	del clusters;del cluster_id; del infodist; # delete summary object to save memory, this does not reduce memory usage why?!?
	import gc
	gc.collect()
	indata.cid_by_bc = cid_by_bc

	if indata.debug: #single process // serial
	    results=[] # create holder for processed reads
	    indata.logfile.write('Part3: Running in debug mode sorting reads to clusters ...\n')
	    progress = Progress(indata.reads2process, logfile=indata.logfile) # creates a progress "bar" thingy
	    #itarate through reads and do the "magicFunction"
	    with progress:
		for tmp in getPairs(indata):
		    progress.update()
		    results.append(foreachread(tmp))
	    indata.logfile.write('Part3: sort finished start writing to files.\n')
	else: # multiple processes in parallel
	    import multiprocessing
	    #create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	    WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10000)
	    results = WorkerPool.imap_unordered(foreachread,getPairs(indata),chunksize=1000)

	# output log message about whats happening
	if not indata.debug: indata.logfile.write('Part2: Sorting reads to cluster '+str(indata.cpus)+' processes  ...\n')
	else: indata.logfile.write('Part4: writing to cluster tempfiles\n')

	progress = Progress(indata.reads2process, logfile=indata.logfile)
	outstrs = {};
	barcodes2files={}
	with progress:
	    if indata.sortfmt == 'fq':
		f1 = open(indata.outfolder+'/temporary.cluster.files/sorted.reads.1.fq','w')
		f2 = open(indata.outfolder+'/temporary.cluster.files/sorted.reads.2.fq','w')
	    elif indata.sortfmt == 'fa':
		f1 = open(indata.outfolder+'/temporary.cluster.files/sorted.reads.fa','w')
	    for pair in results:
		progress.update()
		if pair.cid:# and not (pair.r1.illuminaadapter or pair.r2.illuminaadapter):
		    if indata.sortfmt == 'fq':
			f1.write('_'.join(pair.r1.header.split(' ')) + '_' + str(pair.cid) + '_' + pair.n15.seq + '\n'
				+pair.r1.seq+'\n+\n'
				+pair.r1.qual+'\n')
			f2.write('_'.join(pair.r2.header.split(' ')) + '_' + str(pair.cid) + '_' + pair.n15.seq + '\n'
				+pair.r2.seq+'\n+\n'
				+pair.r2.qual+'\n')
		    elif indata.sortfmt == 'fa':
			f1.write('>' + '_'.join(pair.r1.header.split(' ')) + '_r1_' + str(pair.cid) + '_' + pair.n15.seq + '\n'
				+pair.r1.seq+'\n'+
				'>' + '_'.join(pair.r2.header.split(' ')) + '_r2_' + str(pair.cid) + '_' + pair.n15.seq + '\n'
				+pair.r2.seq+'\n')
	    f1.close()
	    if indata.sortfmt == 'fq': f2.close()

	if not indata.debug:WorkerPool.close()
	if not indata.debug:WorkerPool.join()
	indata.logfile.write('Reads sorted into cluster tempfiles\n')
	indata.logfile.write('Part2: Sorting reads to clusters END\n')
	gc.collect()
        return 0

def blast_it(indata):
    indata.logfile.write('Part3: Identifying clusters by Blast START\n')
    indata.blastid = time.strftime("%A_%d_%b_%Y_%H_%M_%S",time.localtime())

    infile = indata.analyzeclust+'/sorted.reads.fa'
    lines = bufcount(infile)
    reads = lines/2
    pairs = reads/2

    indata.logfile.write('Splitting infiles ...\n')
    tempfiles = {}
    f = open(infile)
    sub=1;rc=0.00;f_sub = open(infile+'.'+str(sub)+'_sub.fa','w');tempfiles[f_sub.name]=True
    progress = Progress(lines, logfile=indata.logfile ,unit='line')
    with progress: # Split fasta to smaller parts for parallel blasting and results handelling
	for line in f:
	    rc+=0.25
	    progress.update()
	    f_sub.write(line)
	    if rc >= indata.readsinsub:
		sub+=1;
		f_sub.close();
		#if sub > 100: break
		f_sub = open(infile+'.'+str(sub)+'_sub.fa','w');
		rc=0.00;tempfiles[f_sub.name]=True
    f_sub.close()
    indata.logfile.write('Done. Resulted in '+str(sub)+' parts in total.\n')

    if indata.debug: #single process // serial
	results=[] # create holder for processed reads
	for tmp in yielder(tempfiles.keys(), indata):
	    results.append(foreach3(tmp))
    else: # multiple processes in parallel
	import multiprocessing
	#create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10)
	results = WorkerPool.imap_unordered(foreach3,yielder(tempfiles.keys(), indata),chunksize=1)

    indata.logfile.write('Now doing blasting and xml parsing ... \n')
    progress = Progress(len(tempfiles.keys()), logfile=indata.logfile ,unit='file')
    info_dict = {'total':0}
    with progress: # blast and prse results in parallel store info in dict
	for tmp in results:
	    progress.update()
	
	    info_dict['total'] += tmp['total']
	    for cluster_id, data in tmp.iteritems():

		if cluster_id != 'total':

		    try: info_dict[cluster_id]['total'] += tmp[cluster_id]['total']
		    except KeyError: info_dict[cluster_id] = {'total':tmp[cluster_id]['total']}

		    for genome, count in data.iteritems():
			if genome != 'total':
			    try:info_dict[cluster_id][genome] += count
			    except KeyError:info_dict[cluster_id][genome] = count

    indata.logfile.write('Parsing info_dict ... \n')
    
    if indata.printblast:
	import os
	try: os.makedirs(indata.outfolder+'/blastReports')
	except OSError:pass
    orgname = indata.outfile.name
    graph_info = {}
    import numpy as np
    for indata.mrc in [2,10,20,30,40,50,60,70,80,90,100,125,150,175,200,225,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000]:
	monoclonal_clusts = 0
	non_class = 0
	tot_clust  = 0
	removed_primer_dimer=0
	pdpercent=90
	beadtypes = {}
	indata.logfile.write('Runnning mrc '+str(indata.mrc)+' ...\n')
	indata.outfile.close()
	indata.outfile = open(orgname+'.mrc_'+str(indata.mrc),'w')
	if indata.mrc == 2 and indata.printblast:
	    temp = 0
#	    for i in info_dict.keys():
#		if info_dict[i]['total'] >= 2:
#		    temp+=1
#	    progress = Progress(temp, logfile=indata.logfile ,unit='cluster')
	    progress = Progress(len(info_dict), logfile=indata.logfile ,unit='cluster')
	    progress.__enter__()
	for cluster_id, data in info_dict.iteritems(): # go through info dict and summarize each cluster then summarize all clusters #could be done in parallel
	    if cluster_id == 'total':continue
	    if info_dict[cluster_id]['total'] >= indata.mrc:
		if indata.mrc == 2 and indata.printblast:
		    progress.update()
		    fx = open(indata.outfolder+'/blastReports/c'+str(cluster_id)+'.txt','w')
		    fx.write( '######### CLUSTER '+str(cluster_id)+' ##########'+data['output'])
		    fx.close()
		# print some info
		output = 'Cluster number '+str(cluster_id)+':\n'
		output += 'Total number of read pairs '+str(data['total'])+'\n'
		
		# check what amplicons were found in the results and determine cluster genome
		amplicons = 0
		genome = 'Unknown'
		max_rc = 0
		monoclonal = None
		
		try: nohitperc = round(100*float(data['No Hits'])/data['total'],2)
		except KeyError: nohitperc = 0
		try: disagreeperc = round(100*float(data['Pair Disagree'])/data['total'],2)
		except KeyError: disagreeperc = 0
		beforeremovalreads = data['total']
		postremtotal = data['total']
		if indata.skipnohits:
		    try:postremtotal=postremtotal-data['No Hits']
		    except KeyError: pass
		    try:postremtotal=postremtotal-data['Pair Disagree']
		    except KeyError: pass
		    output += 'Total number of (SE) reads after removing "NoHits" and "Pair Disagree"-read pairs '+str(postremtotal)+'\n'# (='+str(results['total']/2)+'pairs)\n'
		
		for hit,count in data.iteritems():
		    if indata.skipnohits and postremtotal == 0:
			try:output += 'No Hits'+'\t'+str('NA ')+'% ('+str(data['No Hits'])+' read pairs) (originally '+str(nohitperc)+'% before no hits removal)\n';
			except KeyError:pass
			try:output += 'Pair Disagree'+'\t'+str('NA ')+'% ('+str(data['Pair Disagree'])+' read pairs) (originally '+str(disagreeperc)+'% before no hits removal)\n';
			except KeyError:pass
			break
		    
		    if hit != 'output': percentage = round(100*float(count)/postremtotal,2)
		    if hit == 'total' or hit == 'output': continue
		    elif hit == 'No Hits' or hit == 'Pair Disagree':
			    if hit == 'No Hits':	beforeperc = nohitperc
			    elif hit == 'Pair Disagree':beforeperc = disagreeperc
			    if not indata.skipnohits:	output += hit+'\t'+str(percentage)+'% ('+str(count)+' read pairs)\n'
			    else:			output += hit+'\t'+str('NA ')+'% ('+str(count)+' read pairs) (originally '+str(beforeperc)+'% before no hits removal)\n'
			    continue
		    else: output += hit+'\t'+str(percentage)+'% ('+str(count)+' read pairs)\n'

		    if percentage >= indata.rqpct and (hit != 'No Hits' and hit != 'Pair Disagree'):
			monoclonal = True
			genome = hit
		    if count > max_rc:
			    max_rc = count;
			    genome = hit
		    if percentage >= indata.gpct: amplicons += 1 # if more than 2% of read pop else disregard
		
		#output += str(amplicons)+'amplicons found genome thought to be '+genome +'\n'
		if amplicons == 0 and not monoclonal:
		    monoclonal = None
		    genome = 'Unknown'
		elif amplicons == 1:
		    monoclonal = True
		    if round(100*float(data[genome])/postremtotal,2) < indata.rqpct: monoclonal = False; genome = 'Mixed'
		elif amplicons > 1 and not monoclonal:
		    monoclonal = False
		    genome = 'Mixed'
		elif (amplicons > 1) and round(100*float(data[genome])/postremtotal,2) >= indata.rqpct and monoclonal:pass
		else: indata.logfile.write('WARNING: something is really odd!!!\n')
		output += 'Cluster classified as monoclonal='+str(monoclonal)+' and genome is '+str(genome)+'.\n'
		output += '\n'
	
		#for [tot_reads, output, monoclonal,genome,nohitperc,beforeremovalreads] in results:
		#	progress.update()
		if beforeremovalreads >= indata.mrc:
		    if nohitperc > pdpercent or postremtotal == 0:
			removed_primer_dimer+=1;
			output += 'This cluster has to many pairs classified as primerdimer or similar and will therefore be discarded.\n'
			continue
		    indata.outfile.write(str(output)+'\n')
		    tot_clust+=1
		    if monoclonal: monoclonal_clusts+=1
		    elif monoclonal== None: non_class+=1
		    try: beadtypes[genome] += 1
		    except KeyError: beadtypes[genome] = 1
	if indata.mrc == 2 and indata.printblast:progress.__exit__()
	    
	if tot_clust > 0:
	    moncperc = str(round(100*float(monoclonal_clusts)/tot_clust,2))
	    noclassperc=str(round(100*float(non_class)/tot_clust,2))
	else: moncperc=noclassperc='0'
	indata.outfile.write(str( 'We found '+str(tot_clust)+' clusters, with atleast '+str(indata.mrc)+' reads per cluster, out of theese were '+moncperc+'% monoclonal ('+str(monoclonal_clusts)+')')+' and '+noclassperc+'% not classifiable ('+str(non_class)+'), '+str(removed_primer_dimer)+' clusters were classified as primer dimer (>='+str(pdpercent)+'% nohits) and removed\n')
	graph_info[indata.mrc] = {'monoclonal':moncperc,'clustercount':tot_clust,'nonclass':noclassperc,'primerdimer':removed_primer_dimer}
    
	total = sum([ count for beadtype,count in beadtypes.iteritems()])
	graph_info[indata.mrc]['beadtype_total']=total;
	for beadtype,count in beadtypes.iteritems():
	    indata.outfile.write(str( beadtype +' '+ str(round(100*float(count)/total,2))+'% ('+str(count)+' clusters)')+'\n');
	    graph_info[indata.mrc][beadtype]=count;
#    if indata.printblast:fx.close()

    indata.logfile.write('Making graphs ... \n')
    temp_x=graph_info.keys()
    temp_x.sort()
    y1=[];x =[];y2=[]
    for i in temp_x:
	    x.append(i)
	    y1.append(graph_info[i]['monoclonal'])
	    y2.append(graph_info[i]['clustercount'])
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rc
    for scale in [[0,5000,0,20],[0,500,0,20],[0,1000,0,20]]:
	    fig = plt.figure(figsize=(20, 15), dpi=100)
	    ax = fig.add_subplot(111)
	    ax.plot(x, y1, '-', label = 'Percentage Monoclonal')
	    #ax.plot(x, Rn, '-', label = 'Rn')
	    ax2 = ax.twinx()
	    ax2.plot(x, y2, '-r', label = 'Number of clusters')
	    lines, labels = ax.get_legend_handles_labels()
	    lines2, labels2 = ax2.get_legend_handles_labels()
	    ax2.legend(lines + lines2, labels + labels2, loc=7)

	    ax.grid(b=True, which='both')
	    ax.set_xlabel('Read pairs per Barcode Cluster')
	    ax.set_ylabel('Percentage Monoclonal')
	    ax2.set_ylabel('Number of Clusters')
	    
	    ax.set_ylim(0,100)
	    ax2.set_ylim(0, 1000)
	    
	    ax.set_xlim(scale[0],scale[1])
	    ax2.set_xlim(scale[0],scale[1])

	    ax.set_xticks(np.arange(scale[0],scale[1]+1,scale[1]/10))
	    ax.set_yticks(np.arange(0,101,10))
	    ax2.set_yticks(np.arange(0,1001,100))

	    plt.savefig(indata.outfolder+'/finalresult.x_'+str(scale[0])+'-'+str(scale[1])+'.y_'+str(scale[2])+'-'+str(scale[3])+'.pdf')
	    plt.close()
	    
    indata.logfile.write( 'done\n')

    indata.logfile.write('Part3: Identifying clusters by Blast END\n')
    
def yielder(lista, indata):
    for entry in lista: yield [entry, indata]
    
def foreach3(tmp):
    [infile,indata] = tmp
    tmp = classify_cluser(indata,infile=infile)
    return tmp

def foreachread(tmp):

    # unpack info
    pair, indata = tmp
    del tmp
    
    # convert to SEAseq readpair
    pair = SEAseqpair(pair.header, pair.r1, pair.r2)
    
    pair.identify(C_HANDLE, indata)
    pair.getN15()
    
    if indata.cid_by_bc: pair.get_cid(indata)
    
    return pair

def getindata():
    import argparse
    argparser = argparse.ArgumentParser(description='Analysis of SEAseq amplicon data.', formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument(	'--debug',		dest='debug', 			action='store_true', 			required=False,	default=False,	help='Debug (run as regular single process python script).')
    argparser.add_argument(	'-skip',		dest='skip',	metavar='N',				type=int,	required=False,	default=0,	help='Skip the first N read pairs in files (default 0).')
    argparser.add_argument(	'-stop',		dest='stop',	metavar='N',				type=int,	required=False,	default=0,	help='Stop after N read pairs, set to 0 to disable (default 0).')
    argparser.add_argument(	'-r1',			dest='reads1',	metavar='FILE',				type=file,	required=True, 			help='Indata "fastq"-file read1.')
    argparser.add_argument(	'-r2',			dest='reads2',	metavar='FILE',				type=file,	required=True,	default=None,	help='Indata "fastq"-file read2.')
    argparser.add_argument(	'-p',			dest='cpus',	metavar='N',				type=int,	required=False,	default=1,	help='The number of processes to start (default 1).')
    argparser.add_argument(	'-o',			dest='outfile',	metavar='outfile',			type=str,	required=True,	default=False,	help='Print output to outfile.')
    argparser.add_argument(	'-l',			dest='logfile',	metavar='logfile',			type=str,	required=False,	default=False,	help='Print log messages to logfile (default stderr).')
    argparser.add_argument(	'-random',		dest='n',	metavar='N',				type=int,	required=False,	default=0,	help='Use a random subset of N read pairs, this option is slower (default 0 = off). Can not be used in combination with "-skip" or "-stop"')
    argparser.add_argument(	'-hm',			dest='handlemm',metavar='N',				type=int,	required=False,	default=0,	help='Number off missmatches allowed in handle sequence (default 0)')
    argparser.add_argument(	'-bm',			dest='bcmm',	metavar='N',				type=int,	required=False,	default=0,	help='Number off missmatches allowed in barcode sequence during clustering (default 0)')
    argparser.add_argument(	'-outfolder',		dest='outfolder',metavar='FOLDER',			type=str,	required=True,	default='.',	help='Folder where temporary outputfiles are stored (default current working dir).')
    argparser.add_argument(	'--keeptemp',		dest='keeptemp', 		action='store_true', 			required=False,	default=False,	help='Keep temporary files (default is to delete tempfiles).')
    argparser.add_argument(	'-analyzeclusters',	dest='analyzeclust',metavar='FOLDER',			type=str,	required=False,	default=False,	help='Only analyze cluster files (default False)')
    argparser.add_argument(	'-analysis',		dest='analyzeclust',metavar='[blast/bowtie]',		type=str,	required=False,	default='blast',help='Type of analysis mapping by bowtie or use blast (default blast) NOTE: currently only does blast!')
    argparser.add_argument(	'-blastsettting',	dest='blastsetting',metavar='\["strict"\|"sloppy"\]',	type=str,	required=False,	default='strict',help='Setting for the blast either "strict" or "sloppy" (default False)')
    argparser.add_argument(	'-mrc',			dest='mrc',	metavar='N',				type=int,	required=False,	default=1,	help='Minimum number of reads per cluster to consider it (default 1) DISABLED: tests from 10 to 1000')
    argparser.add_argument(	'-seed',		dest='seed',	metavar='N',				type=int,	required=False,	default=100,	help='Number of top barcodes (with most reads) to use as seeds in clustering(default 100)')
    argparser.add_argument(	'-gpct',		dest='gpct',	metavar='N',				type=float,	required=False,	default=2.0,	help='Disreagard genomes with less than N percent of read population (default 2)')
    argparser.add_argument(	'-rqpct',		dest='rqpct',	metavar='N',				type=float,	required=False,	default=95.0,	help='Require at least N percent of read population for the major amplicon type to count the cluster as monoclonal (default 95, overrides the "-gpct" option)')
    argparser.add_argument(	'-ris',			dest='readsinsub',metavar='N',				type=float,	required=False,	default=1000,	help='Read pairs to be placed in each sub part for blasting (default 1000)')
    argparser.add_argument(	'--skipnohits',		dest='skipnohits', 		action='store_true', 			required=False,	default=False,	help='Skip the No Hits reads in all calculations (default False).')
    argparser.add_argument(	'--printblast',		dest='printblast', 		action='store_true', 			required=False,	default=False,	help='print details of the BLAST searcch for each cluster (default False).NOT WORKING!')
    argparser.add_argument(	'--onlycluster',	dest='onlycluster',		action='store_true',			required=False,	default=False,	help='Only do clustering of n15 sequences (default False)')
    argparser.add_argument(	'--onlysort',		dest='onlysort',		action='store_true',			required=False,	default=False,	help='Only sort reads according to n15 sequences (default False)')
    argparser.add_argument(	'-sortfmt',		dest='sortfmt',	metavar='[fa/fq]',			type=str,	required=False,	default='fa',	help='Format to output reads to fa=fasta or fq=fastq (default fastq)')
    argparser.add_argument(	'-blastdb',		dest='blastdb',	metavar='STR',			type=str,	required=False,	default="/bubo/proj/b2011011/SEAseq/reference/4amplicons/4amplicons.fasta",	help='blast database (default "/bubo/proj/b2011011/SEAseq/reference/4amplicons/4amplicons.fasta")')
    indata = argparser.parse_args(sys.argv[1:])
    indata.selftest = False

    assert indata.reads1 != indata.reads2, 'Error: read 1 and read 2 cannot be same file!\n'
    
    import os
    try: os.makedirs(indata.outfolder)
    except OSError:pass
    
    if indata.outfile: indata.outfile = open(indata.outfile, 'w',1)
    else: indata.outfile = sys.stdout
    
    if indata.logfile: indata.logfile = open(indata.logfile, 'w',1)
    else: indata.logfile = sys.stderr
    
    # get the readcount
    if not indata.analyzeclust:
	indata.logfile.write('Getting readcount ... ')
	indata.numreads=bufcount(indata.reads1.name)/4
	#indata.numreads=10000
	indata.logfile.write(str(indata.numreads)+' read pairs in fastq files.\n');

	# calculate the number of reads to process
	indata.reads2process = indata.numreads
	if indata.skip: indata.reads2process -= indata.skip
	if indata.stop: indata.reads2process = indata.stop
	if indata.n:    indata.reads2process = indata.n
   
    indata.cid_by_bc = False

    return indata



#####
#check if run or imported // call main() or not
#####
if __name__ == "__main__":
    main()
#END of script
