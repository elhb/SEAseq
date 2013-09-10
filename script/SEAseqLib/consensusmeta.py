import sys

class RunStatCounter(object):
    def __init__(self, config):
        self.config = config
        self.ampliconcombinations = {}
        self.clustercount = 0
        self.definedclusters = 0
        self.definedclustersMono = 0
        self.undefinedclusters = 0
        self.junkclusters = 0
        self.lowreadclusters = 0
        
    def addcluster(self, cluster):
        
        self.clustercount += 1
        
        if cluster.lowread:
            self.lowreadclusters += 1
            return
        if cluster.ampliconpairs == 0:
            self.junkclusters +=1
            return
        if cluster.definedampliconcount < 1:
            self.undefinedclusters += 1
            return
        
        self.definedclusters += 1
        
        # count combo of defined amplicons
        ampliconnames = cluster.definedamplicons.keys()
        ampliconnames.sort()
        ampliconcombo = '/'.join(ampliconnames)        
        try:
            self.ampliconcombinations[ampliconcombo]['count'] += 1
        except KeyError:
            self.ampliconcombinations[ampliconcombo] = {'count':1}
            #self.ampliconcombinations[ampliconcombo]['monos'] = {'All':0}
            self.ampliconcombinations[ampliconcombo]['monos'] = {}
            self.ampliconcombinations[ampliconcombo]['poly'] = 0
        
        # count combo of monoclonal amplicons
        monoAmps = {}
        for ampname, amplicon in cluster.definedamplicons.iteritems():
            if amplicon.monoclonal: monoAmps[amplicon.type] = amplicon.monoclonal
        
        if not monoAmps.keys():
            self.ampliconcombinations[ampliconcombo]['poly'] += 1
        else:
            mononames = monoAmps.keys()
            mononames.sort()
            monoCombo = '/'.join(mononames)
            try:
                self.ampliconcombinations[ampliconcombo]['monos'][monoCombo] += 1
            except KeyError:
                self.ampliconcombinations[ampliconcombo]['monos'][monoCombo] = 1
            if monoCombo == ampliconcombo:
                #self.ampliconcombinations[ampliconcombo]['monos']['All'] += 1
                self.definedclustersMono += 1
    
    def createsummary(self, indata):
        output = ''
        output += ('##### SUMMARY #####'+'\n')
        output += (  str(self.clustercount)   +' clusters processed, out of these were:'+'\n')
        output += (  str(self.junkclusters)   +' ('+str(round(100*float(self.junkclusters   )/float(self.clustercount),2))+'%) clusters of only adapter sequences or faulty primers'+'\n')
        output += (  str(self.lowreadclusters)+' ('+str(round(100*float(self.lowreadclusters)/float(self.clustercount),2))+'%) clusters with to few reads to be analysed'+'\n')
        output += (  str(self.definedclusters)+' ('+str(round(100*float(self.definedclusters)/float(self.clustercount),2))+'%) clusters had defined amplicons'+'\n')

        for ampliconcombo, data in self.ampliconcombinations.iteritems():
            count = data['count']
            output += (  str(count)+' ('+str(round(100*float(count)/float(self.clustercount),2))+'%) '+ ampliconcombo+' '+ 'whereof:'+'\n')
            for monoCombo, count2 in data['monos'].iteritems():
                output += (  '\t'+' '+str(count2)+' ('+str(round(100*float(count2)/float(count),2))+'%) '+'were monoclonal for'+' '+monoCombo+'\n')
            output += (  '\t'+' '+str(data['poly'])+' ('+str(round(100*float(data['poly'])/float(count),2))+'%) '+'were polyclonal\n')
        
        tmppercentage = 0
        if self.definedclusters: tmppercentage = round(100*float(self.definedclustersMono)/float(self.definedclusters),2)
        output += (
            str(self.definedclusters)+' ('+str(round(100*float(self.definedclusters)/float(self.clustercount),2))+'%) has at least one defined amplicon out of these are '+str(self.definedclustersMono)+' ('+str(tmppercentage)+'%) monoclonal for the defined amplicon(s)\n'+
            '(ie there is only one consensus sequence of that type with more than '+str(indata.minimum_reads)+' reads and '+str(indata.minimum_support)+'% support, clustering done with '+str(indata.clustering_identity)+'% identity cutoff)\n'
            )
        return output

def meta(indata):

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd, skip=indata.skip ,stop=indata.stop ,random=indata.n)
    import os
    import sys
    if os.path.exists(config.outfile): os.remove(config.outfile)
    config.openconnections()
    
    writelogheader(config.logfile)

    # settings
    config.logfile.write('Get infiles from config-file.\n')
    config.load()
    config.getreads2process()
    
    # set primerpairs
    from SEAseqLib.mainLibrary import PrimerPair
    config.primerpairs = {}
    config.primerpairs['16s']      = PrimerPair('GTGBCAGCMGCCGCGGTAA',         'ACAHBTCACRRCACGAGCTGACGAC',    '16s')
    config.primerpairs['its']      = PrimerPair('G?GBCTTBTACWCACYGCCCGTC',     'CTCYDRNWGCCVRGGCATCCACC',      'its')
    config.primerpairs['ecoli']    = PrimerPair('TGCGAACGCGCGAATCAACTGG',      'AAGCGCGCGGCTGAATTACTGG',       'ecoli')
    config.primerpairs['m13']      = PrimerPair('GCCTCGTTCCGGCTAAGTAACATGGAG', 'AGTTGCGCCGACAATGACAACAACC',    'm13')
    config.primerpairs['myco']     = PrimerPair('ATGCCGCAGCCAAGAACGCATC',      'TTCGTGGCACTTGCCGAACTGG',       'myco')
    config.primerpairs['lambda']   = PrimerPair('TCAGCTATGCGCCGACCAGAACAC',    'TTCCATGACCGCACCAACAGGCTC',     'lambda')

    
    import multiprocessing as mp
    man = mp.Manager()
    clusterq = man.Queue(1000)
    # Run subprocess that get read pars and adds clusters to queue last add a END
    reader = mp.Process(target=getClustersAndPairs,args=(config,clusterq))
    reader.start()
    
    # make a worker pool that works with the clusters returned from the generator
    from SEAseqLib.mainLibrary import Progress
    if indata.debug: #single process // serial
	config.logfile.write('debugging: ');
	sys.stdout.write('debugging: ')
	config.logfile.write('Running in debug mode ')
	results=[] # create holder for processed reads
	progress = Progress(config.clustercount, logfile=config.logfile) # creates a progress "bar" thingy
	with progress:
	    for cluster_pairs in clusteriterator(clusterq, indata):
		progress.update()
		results.append(foreachcluster_meta(cluster_pairs))
	config.logfile.write('finished, making summary ... \n')
    else: # multiple processes in parallel
	import multiprocessing
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10000)
	results = WorkerPool.imap_unordered(foreachcluster_meta,clusteriterator(clusterq,indata),chunksize=1)
	#results = WorkerPool.imap(foreachcluster_meta,clusteriterator(clusterq, indata),chunksize=10)

    if not indata.debug: config.logfile.write('Part1: Per cluster action '+str(indata.cpus)+' processes  ...\n')
    progress = Progress(config.clustercount, logfile=config.logfile, unit='cluster',mem=True)
    counter = RunStatCounter(config)
    
    statstable = open(config.path+'/meta.statstable','w',1)
    statsheader = ['clusterid','number of reads in total','number of adaper reads','number of strange primers','its reads','16s reads','its','16s','its monoclonal','16s monoclonal','number of consensus types','number of consensus types with good support']
    for amptype in ['ecoli','myco','lambda','m13']: statsheader.append(amptype+' reads');statsheader.append(amptype+' monoclonal');statsheader.append(amptype)
    statstable.write('\t'.join(statsheader))

    seqdump = open(config.path+'/meta.sequences.fa','w',1)

    with progress:
	for tmp in results:
	    progress.update()
            [output, cluster] = tmp
            
            counter.addcluster(cluster)

            config.outfile.write(output)

            # SEQDICT printing to file ONLY for both monoclonal
            #if _its and _16s:
            #    for cons_type in seqdict:
            #	for consensus in seqdict[cons_type]:
            #	    seqdump.write(seqdict[cons_type][consensus])

	    #if return_info:
		#statstable.write( '\n'+'\t'.join([str(return_info[stat]) for stat in statsheader]) )
		#assert len(return_info) == len(statsheader)
	seqdump.close()
	statstable.close()

    reader.join()
    	
    # create a nice run summary
    config.outfile.write(counter.createsummary(indata))
    
    config.logfile.write('Done.\n')
    return 0

def foreachcluster_meta(cluster_pairs):
    indata = cluster_pairs[1]
    cluster = cluster_pairs[0][0]
    config = cluster_pairs[0][1]
    verb = 3
    
    lowreadcutoff = 1
    output =  '\n###--- Cluster number '+str(cluster.id)+' -> '+str(cluster.readcount)+' pairs. ---###\n'

    if cluster.readcount < lowreadcutoff:
	output += 'Less than '+str(lowreadcutoff)+' read pairs.\n'
        cluster.lowread = True
	return [output,cluster]
    
    else:
        cluster.lowread = False
        
        pairsOut = cluster.createtempfile(config)
        filesOut = cluster.clusterreadpairs(config, indata)
        
        aligmnmentsOut = ''
        perAmpOut = ''
        
        if cluster.ampliconpairs < 0:
            cluster.loadconsensuses(config)
            cluster.loadconsensusalignemnts(config)
            cluster.loadconsensussequences(config)
            cluster.consensusesToAmplicons(config)

            for amplicon in cluster.amplicons.values():
                for consensus in amplicon.allels: aligmnmentsOut += consensus.alignmentoutput(config)
    
            for amplicon in cluster.amplicons.values(): perAmpOut += amplicon.checkmono(indata)
            cluster.getDefinedAmplicons()
            
            cluster.removetempfiles(config)
        
        output += 'There are '+str(cluster.adaptercount)+' illumina adapter reads.\n'
	output += 'There are '+str(cluster.primererrors)+' primer missmatch reads.\n'
        if cluster.ampliconpairs == 0:
            output += '0 amplicon(s) have enough data (>=1 cons with >= '+str(indata.minimum_support)+'% support and >= '+str(indata.minimum_reads)+' reads)\n'
        if cluster.ampliconpairs > 0:
            output += str(cluster.definedampliconcount)+' amplicon(s) have enough data (>=1 cons with >= '+str(indata.minimum_support)+'% support and >= '+str(indata.minimum_reads)+' reads):\n'
            output += perAmpOut + '\n'
            output += '# Details:\n'
            output += aligmnmentsOut + '\n'
            output += pairsOut + '\n'
            output += filesOut + '\n'

	return [output, cluster]

def clusteriterator(clusterq, indata):  # a generator that yields clusters from queue until it finds a cluster = END
    while True:
	pairs = clusterq.get()
	if pairs == 'END': break
	yield [pairs, indata]
    
def getClustersAndPairs(config,clusterq):
    
    import os
    import time
    config.logfile.write('Reader initiated pid='+str(os.getpid())+'.\n')

    currentclusterid = 1
    from SEAseqLib.mainLibrary import BarcodeCluster
    cluster = BarcodeCluster(currentclusterid)
    #pairs = []
    config.infiles['r1'] = [config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq']
    config.infiles['r2'] = [config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq']
    
    from SEAseqLib.mainLibrary import getPairs
    for tmp in getPairs(config):

	pair, config = tmp
	del tmp
	
	pair.cid = int(pair.header.split(':')[-1].split('_')[1])
	
	if pair.cid == currentclusterid:
            cluster.addreadpair(pair)
	elif pair.cid == currentclusterid+1:
	    while clusterq.qsize() > 160 or clusterq.full(): time.sleep(1); #config.logfile.write('waiting for queue ...\n')
            clusterq.put([cluster,config])
	    currentclusterid = pair.cid
            cluster = BarcodeCluster(currentclusterid)
            cluster.addreadpair(pair)
	else: sys.stdout.write('ERROR 1979 in consensus creation get clusters and pairs.\n')

    clusterq.put('END')
    config.logfile.write('Reader exiting.\n')