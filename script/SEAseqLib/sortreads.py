def foreachread_sort(tmp):

    # unpack info
    pair, config = tmp
    del tmp
    
    # convert to SEAseq readpair
    from SEAseqLib.mainLibrary import sequence
    #pair = SEAseqpair(pair.header, pair.r1, pair.r2)
    
    C_HANDLE = sequence('c handle',"CTAAGTCCATCCGCACTCCT","CTAAGTCCATCCGCACTCCT")
    pair.identify(C_HANDLE, config)
    pair.getN15()
    
    pair.get_cid(config)
    
    return pair

def load_clusters(config,queue):
    import multiprocessing
    import os
    config.logfile.write('\nLoading clustering data ...\n')
    config.logfile.write( 'Loader pid='+str(os.getpid())+'\n')
    from SEAseqLib.mainLibrary import SEAseqSummary
    summary = SEAseqSummary()
    summary.loadclusters(config.clusters_file)
    clustercount = len(summary.clusters)
    config.logfile.write('Data Loaded, reformatting ...\n')
    from SEAseqLib.mainLibrary import Progress, getPairs
    progress = Progress(len(summary.clusters), logfile=config.logfile, unit = 'cluster', mem = True)
    chunksize = 5000
    with progress:
	for cluster_id, infodist in summary.clusters.iteritems():
	    progress.update()
	    if int(infodist['total']) >= config.read_count_per_barcode_cluster_cutoff:
		for barcode in infodist['barcodes']:
		    config.cid_by_bc[barcode] = cluster_id
	    else: pass#print 'low read cluster'
    summary.clusters.clear()
    queue.put(clustercount)
    config.logfile.write('Done returning to master process.\n')

def sortreads(indata):

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd, skip=indata.skip ,stop=indata.stop ,random=indata.n)
    config.openconnections()
    writelogheader(config.logfile)

    #settings
    config.logfile.write('Get infiles from config-file.\n')
    config.load()
    config.getreads2process()

    import multiprocessing
    man = multiprocessing.Manager()
    config.cid_by_bc = man.dict()
#    chunk_list = man.list([[[cluster+chunk*chunksize,[],[], 0, 0] for cluster in xrange(chunksize)] for chunk in xrange(config.clustercount/chunksize+1)])
    queue = multiprocessing.Queue()
    p = multiprocessing.Process(target=load_clusters,args=(config,queue))
    p.start()
    clustercount = queue.get()
    p.join()
    
    import os
    try: os.makedirs(config.path+'/sortedReads')
    except OSError:pass

    from SEAseqLib.mainLibrary import Progress
    if indata.debug: #single process // serial
	results=[] # create holder for processed reads
	config.logfile.write('Part1: Sorting reads to clusters ...\n')
	config.logfile.write('Running in debug mode processing reads ...\n')
	progress = Progress(config.reads2process, logfile=config.logfile, mem=True) # creates a progress "bar" thingy
	#itarate through reads and do the "magicFunction"
	with progress:
	    for tmp in getPairs(config):
		progress.update()
		results.append(foreachread_sort(tmp))
	config.logfile.write('finished, now sorting reads to clusters in memory ... \n')
    else: # multiple processes in parallel
	import multiprocessing
        from SEAseqLib.mainLibrary import getPairs
	#create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10000)
	results = WorkerPool.imap_unordered(foreachread_sort,getPairs(config),chunksize=1000)

    config.logfile.write('Allocating sorting memory  ...\n')
    chunksize = 5000
    chunk_list = [[[cluster+chunk*chunksize,[],[]] for cluster in xrange(chunksize)] for chunk in xrange(config.clustercount/chunksize+1)]
    config.logfile.write(' done.\n')

    if not indata.debug: config.logfile.write('Sorting reads to cluster '+str(indata.cpus)+' processes  ...\n')
    progress = Progress(config.reads2process, logfile=config.logfile, mem = True)
    with progress:
	if config.sortformat == 'fq':
	    f1 = open(config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq','w')
	    f2 = open(config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq','w')
	elif config.sortformat == 'fa':
	    f1 = open(config.path+'/sortedReads/sorted_by_barcode_cluster.fa','w')
	for pair in results:
	    progress.update()
	    if pair.cid:
		if config.sortformat == 'fq':
                    if indata.trimr1:   chunk_list[pair.cid/chunksize][pair.cid%chunksize][1].append('_'.join(pair.r1.header.split(' ')) + '_' + str(pair.cid) + '_' + pair.n15.seq + '\n'+pair.r1.seq[:-indata.trimr1]+'\n+\n'+pair.r1.qual[:-indata.trimr1]+'\n')
                    else:               chunk_list[pair.cid/chunksize][pair.cid%chunksize][1].append('_'.join(pair.r1.header.split(' ')) + '_' + str(pair.cid) + '_' + pair.n15.seq + '\n'+pair.r1.seq+'\n+\n'+pair.r1.qual+'\n')
                    if indata.trimr2:   chunk_list[pair.cid/chunksize][pair.cid%chunksize][2].append('_'.join(pair.r2.header.split(' ')) + '_' + str(pair.cid) + '_' + pair.n15.seq + '\n'+pair.r2.seq[:-indata.trimr2]+'\n+\n'+pair.r2.qual[:-indata.trimr2]+'\n')
                    else:	        chunk_list[pair.cid/chunksize][pair.cid%chunksize][2].append('_'.join(pair.r2.header.split(' ')) + '_' + str(pair.cid) + '_' + pair.n15.seq + '\n'+pair.r2.seq+'\n+\n'+pair.r2.qual+'\n')
		elif config.sortformat == 'fa':
		    f1.write('>' + '_'.join(pair.r1.header.split(' ')) + '_r1_' + str(pair.cid) + '_' + pair.n15.seq + '\n'+pair.r1.seq+'\n'+
			     '>' + '_'.join(pair.r2.header.split(' ')) + '_r2_' + str(pair.cid) + '_' + pair.n15.seq + '\n'+pair.r2.seq+'\n')

    if config.sortformat == 'fq':
	config.logfile.write('\nWriting to fastq files sorted by cluster id ...\n')
	progress = Progress(clustercount, logfile=config.logfile, unit = 'cluster')
	with progress:
	    for chunk in chunk_list:
		for cluster in chunk:
			if cluster[0] <= config.clustercount: progress.update()
			for i in xrange(len(cluster[1])):
			    f1.write(cluster[1][i])
			    f2.write(cluster[2][i])	

    #close connections
    if config.sortformat in ['fa','fq']:	f1.close()
    if config.sortformat == 'fq':		f2.close()

    if not indata.debug: WorkerPool.close()
    if not indata.debug: WorkerPool.join()
    config.logfile.write('Reads sorted into sortedReads/sorted_by_barcode_cluster.1.fq\n')
    config.logfile.write('Part1: Sorting reads to clusters END\n----------\n')
    config.logfile.close()
    config.outfile.write('Finished as expected\n')
    config.outfile.close()
    return 0

