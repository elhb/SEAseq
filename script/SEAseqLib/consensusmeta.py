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
        
        import os
	if os.path.islink(config.path+'/meta.statstable'): os.unlink(config.path+'/meta.statstable')
	self.statstable = open(config.path+'/meta.statstable','w',1)
        self.statsheader = ['clusterid','number of reads in total','number of adaper reads','number of strange primers','its reads','16s reads','its','16s','its monoclonal','16s monoclonal','number of consensus types','number of consensus types with good support','monoclonal for all defined amplicons']
        #for amptype in ['ecoli','myco','lambda','m13']: statsheader.append(amptype+' reads');statsheader.append(amptype+' monoclonal');statsheader.append(amptype)
        self.statstable.write('\t'.join(self.statsheader)+'\n')

        
    def addcluster(self, cluster):
        
        self.clustercount += 1
        monoForAllDefined = None
        if cluster.lowread:
            self.lowreadclusters += 1
            return
        if cluster.ampliconpairs == 0:
            self.junkclusters +=1
            return
        if cluster.definedampliconcount < 1:
            self.undefinedclusters += 1
            #return
        else:
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
            
            monoForAllDefined = False
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
                    monoForAllDefined = True
    
        #print to stats info file
        #NOTE: ONLY for 16s its stuff!!!
        #statstable.write( '\n')
        self.statstable.write(str(cluster.id                                     )+'\t')#'clusterid'
        self.statstable.write(str(cluster.readcount                              )+'\t')#'number of reads in total'
        self.statstable.write(str(cluster.adaptercount                           )+'\t')#'number of adaper reads'
        self.statstable.write(str(cluster.primererrors                           )+'\t')#'number of strange primers'
        try:            self.statstable.write(str(cluster.amplicons['its'].readcount             )+'\t')#'its reads',
        except KeyError:self.statstable.write(str(0                                              )+'\t')#'its reads',
        try:            self.statstable.write(str(cluster.amplicons['16s'].readcount             )+'\t')#'16s reads',
        except KeyError:self.statstable.write(str(0                                              )+'\t')#'16s reads',
        try:            self.statstable.write(str(bool(cluster.amplicons['its'].allelecount)     )+'\t')#'its'
        except KeyError:self.statstable.write(str(False                                          )+'\t')#'its'
        try:            self.statstable.write(str(bool(cluster.amplicons['16s'].allelecount)     )+'\t')#'16s'
        except KeyError:self.statstable.write(str(False                                          )+'\t')#'16s'
        try:            self.statstable.write(str(cluster.amplicons['its'].monoclonal            )+'\t')#'its monoclonal'
        except KeyError:self.statstable.write(str(False                                          )+'\t')#'its monoclonal'
        try:            self.statstable.write(str(cluster.amplicons['16s'].monoclonal            )+'\t')#'16s monoclonal'
        except KeyError:self.statstable.write(str(False                                          )+'\t')#'16s monoclonal'
        for name, primerpair in self.config.primerpairs.iteritems(): pass
        self.statstable.write(str(cluster.ampliconcount                          )+'\t')#'number of consensus types'
        self.statstable.write(str(cluster.definedampliconcount                   )+'\t')#'number of consensus types with good support'
        self.statstable.write(str(monoForAllDefined                              )+'\n')#'cluster is mono for all defined amplicons'

    
    def createsummary(self, config):
        output = ''
        output += ('##### SUMMARY #####'+'\n')
        output += (  str(self.clustercount)   +' clusters processed, out of these were:'+'\n')
        output += (  str(self.junkclusters)   +' ('+str(round(100*float(self.junkclusters   )/float(self.clustercount),2))+'%) clusters of only adapter sequences or faulty primers'+'\n')
        output += (  str(self.lowreadclusters)+' ('+str(round(100*float(self.lowreadclusters)/float(self.clustercount),2))+'%) clusters with to few reads to be analysed'+'\n')
        output += (  str(self.definedclusters)+' ('+str(round(100*float(self.definedclusters)/float(self.clustercount),2))+'%) clusters had defined amplicons'+'\n')

        for ampliconcombo, data in self.ampliconcombinations.iteritems():
            count = data['count']
            output += (  str(count)+' ('+str(round(100*float(count)/float(self.clustercount),2))+'% of total, '+str(round(100*float(count)/float(self.definedclusters),2))+'% of defined) '+ ampliconcombo+' '+ 'whereof:'+'\n')
            for monoCombo, count2 in data['monos'].iteritems():
                output += (  '\t'+' '+str(count2)+' ('+str(round(100*float(count2)/float(count),2))+'%) '+'were monoclonal for'+' '+monoCombo+'\n')
            output += (  '\t'+' '+str(data['poly'])+' ('+str(round(100*float(data['poly'])/float(count),2))+'%) '+'were polyclonal for the defined amplicon(s)\n')
        output += (  str(self.undefinedclusters)+' ('+str(round(100*float(self.undefinedclusters)/float(self.clustercount),2))+'%) were undefined.\n')
        
        tmppercentage = 0
        if self.definedclusters: tmppercentage = round(100*float(self.definedclustersMono)/float(self.definedclusters),2)
        output += (
            str(self.definedclusters)+' ('+str(round(100*float(self.definedclusters)/float(self.clustercount),2))+'%) has at least one defined amplicon out of these are '+str(self.definedclustersMono)+' ('+str(tmppercentage)+'%) monoclonal for the defined amplicon(s)\n'+
            '(ie there is only one consensus sequence of that type with more than '+str(config.minReadCountPerConsensus)+' reads and '+str(config.minReadPopSupportConsensus)+'% support, clustering done with '+str(config.minConsensusClusteringIdentity)+'% identity cutoff and '+str(100*config.allowedAllelLevelVariation)+'% allowed allele level variation)\n'
            )
        self.statstable.close()
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
    
    config.loadPrimers()
    ## set primerpairs
    #from SEAseqLib.mainLibrary import PrimerPair
    #config.primerpairs = {}
    #for line in config.primerset:
    #    if line[0] != "#":
    #        line = line.rstrip().split('\t')
    #        try:config.primerpairs[line[0]] = PrimerPair(line[1],line[2],line[3])
    #        except IndexError: print line;print 'ERROR: Wrong number of columns in primerset file.'; sys.exit()
    
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
	progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile) # creates a progress "bar" thingy
	with progress:
	    for cluster_pairs in clusteriterator(clusterq, indata):
		progress.update()
		results.append(foreachcluster_meta(cluster_pairs))
	config.logfile.write('finished, making summary ... \n')
    else: # multiple processes in parallel
	import multiprocessing
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=1000)
	results = WorkerPool.imap_unordered(foreachcluster_meta,clusteriterator(clusterq,indata),chunksize=1)
	#results = WorkerPool.imap(          foreachcluster_meta,clusteriterator(clusterq, indata),chunksize=1)

    if not indata.debug: config.logfile.write('Part1: Per cluster action '+str(indata.cpus)+' processes  ...\n')
    progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile, unit='cluster',mem=True)
    counter = RunStatCounter(config)
    
    #compressing takes forever skip this and do later if needed
    #import gzip
    #clusterdump = gzip.open(config.path+'/clusters.pickle.gz','wb',9)
    if indata.tempFileFolder:
	clusterdump = open(indata.tempFileFolder+'/SEAseqtemp/clusters.pickle','wb')
    else:
	import os
	if os.path.exists(config.path+'/clusters.pickle') and os.path.islink(config.path+'/clusters.pickle'): os.unlink(config.path+'/clusters.pickle')
	clusterdump = open(config.path+'/clusters.pickle','wb')
    
    #import cPickle

    with progress:
	for tmp in results:
	    progress.update()
            [output, cluster,picklestring] = tmp
            
            # Dump clusters to file
            #clusterdump.write(repr(cPickle.dumps(cluster))+'\n')
            #cPickle.dump(cluster,clusterdump)
            clusterdump.write(picklestring)

            # do the pivkling on subprocess process only print here ...
            counter.addcluster(cluster)
            config.outfile.write(output)

	clusterdump.close()
        if indata.tempFileFolder:
            import shutil
	    import os
	    if os.path.exists(config.path+'/clusters.pickle') and os.path.islink(config.path+'/clusters.pickle'): os.unlink(config.path+'/clusters.pickle')
            shutil.move(indata.tempFileFolder+'/SEAseqtemp/clusters.pickle',config.path+'/clusters.pickle')

    reader.join()
    	
    # create a nice run summary
    config.outfile.write(counter.createsummary(config))
    
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
        import cPickle
	return [output, cluster, cPickle.dumps(cluster,-1)]
    
    else:
        cluster.lowread = False
        
        pairsOut = cluster.createtempfile(config, indata)
        #import os;
        #os.remove(config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa')
        #return [output,cluster]
        filesOut = cluster.clusterreadpairs(config, indata)
        
        aligmnmentsOut = ''
        perAmpOut = ''
        
        if cluster.ampliconpairs > 0:
            cluster.loadconsensuses(config, indata)
            cluster.loadconsensusalignemnts(config, indata)
            cluster.loadconsensussequences(config, indata)
            cluster.consensusesToAmplicons(config)

            for amplicon in cluster.amplicons.values():
                for consensus in amplicon.allels: aligmnmentsOut += consensus.alignmentoutput(config)
    
            for amplicon in cluster.amplicons.values(): perAmpOut += amplicon.checkmono(config,cluster.readcount)
            cluster.getDefinedAmplicons()
            
            cluster.removetempfiles(config, indata)
        
        output += 'There are '+str(cluster.adaptercount)+' illumina adapter reads.\n'
	output += 'There are '+str(cluster.primererrors)+' primer missmatch reads.\n'
        if cluster.ampliconpairs == 0:
            output += '0 amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% of readpop support, '+str(config.minPercentagePerConsensus)+'% of cluster total readcount and >= '+str(config.minReadCountPerConsensus)+' reads)\n'
        if cluster.ampliconpairs > 0:
            output += str(cluster.definedampliconcount)+' amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% of readpop support, '+str(config.minPercentagePerConsensus)+'% of cluster total readcount and >= '+str(config.minReadCountPerConsensus)+' reads)\n'
        output += perAmpOut + '\n'
        output += '# Details:\n'
        output += '# Alignemnts:\n'
        output += aligmnmentsOut + '\n'
        output += '# Pairs:\n'
        output += pairsOut + '\n'
        output += '# Original cdhit data:\n'
        output += filesOut + '\n'

        import cPickle
	return [output, cluster, cPickle.dumps(cluster,-1)]

def clusteriterator(clusterq, indata):  # a generator that yields clusters from queue until it finds a cluster = END
    while True:
	pairs = clusterq.get()
	if pairs == 'END': break
	yield [pairs, indata]
    
def getClustersAndPairs(config,clusterq):
    
    import os
    import time
    import sys
    config.logfile.write('Reader initiated pid='+str(os.getpid())+'.\n')
    tmpcounter = 0

    from SEAseqLib.mainLibrary import BarcodeCluster
    cluster = BarcodeCluster(1)
    #pairs = []
    if   os.path.exists(config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq'):    config.infilesDictionary['r1'] = [config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq']
    elif os.path.exists(config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq.gz'): config.infilesDictionary['r1'] = [config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq.gz']
    else: print "Error: infiles are not around." ; sys.exit()
    if   os.path.exists(config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq'):    config.infilesDictionary['r2'] = [config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq']
    elif os.path.exists(config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq.gz'): config.infilesDictionary['r2'] = [config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq.gz']
    else:  print "Error: infiles are not around." ; sys.exit()
    missingClustersFlag = False
    missingClusters = []
    
    from SEAseqLib.mainLibrary import getPairs
    for tmp in getPairs(config):

	pair, config = tmp
	del tmp
	
	pair.cid = int(pair.header.split(':')[-1].split('_')[1])
	
	if pair.cid == cluster.id:
            cluster.addreadpair(pair)
	elif pair.cid == cluster.id+1:
	    while clusterq.qsize() > 160 or clusterq.full(): time.sleep(1); #config.logfile.write('waiting for queue ...\n')
            clusterq.put([cluster,config])
            tmpcounter += 1
            cluster = BarcodeCluster(pair.cid)
            cluster.addreadpair(pair)
	elif pair.cid > cluster.id+1:
	    while clusterq.qsize() > 160 or clusterq.full(): time.sleep(1);
            clusterq.put([cluster,config])
            tmpcounter += 1
	    if not missingClustersFlag:
		config.logfile.write('WARNING: missing cluster(s), are you running on a subset of data?\n')
		missingClustersFlag = True
	    for i in xrange(cluster.id+1,pair.cid): missingClusters.append(str(i))
	    tmp = pair.cid
            cluster = BarcodeCluster(tmp)
            cluster.addreadpair(pair)
	else: sys.stdout.write('ERROR 1979 in consensus creation get clusters and pairs.\nCluster id in pair is '+str(pair.cid)+' and the current cluster id is '+str(currentclusterid)+'\n')
    if missingClustersFlag: config.logfile.write('WARNING: Cluster(s) '+', '.join(missingClusters[:-1])+' & '+missingClusters[-1]+' were missing ('+str(len(missingClusters))+' clusters), are you running on a subset of data?\n')

    clusterq.put([cluster,config]) # add the last cluster to queue
    clusterq.put('END')
    config.logfile.write('Reader exiting after adding '+str(tmpcounter+1)+' clusters to queue.\n')
