def foreachread_cluster(tmp):

    from SEAseqLib.mainLibrary import sequence

    # unpack info
    pair, config = tmp
    del tmp
    
    # convert to SEAseq readpair
    #pair = SEAseqpair(pair.header, pair.r1, pair.r2)

    C_HANDLE = sequence('c handle',"CTAAGTCCATCCGCACTCCT","CTAAGTCCATCCGCACTCCT")
    pair.identify(C_HANDLE, config)
    pair.getN15()
    
    return pair

def clusterbarcodes(indata):
    
    from SEAseqLib.mainLibrary import Configuration, writelogheader

    config = Configuration(indata.path, indata.cmd, skip=indata.skip ,stop=indata.stop ,random=indata.n)
    config.openconnections()
    
    writelogheader(config.logfile)

    # settings
    config.logfile.write('Get infiles and settings from config-file.\n')
    config.load()
    config.getreads2process()
    config.save()

    config.logfile.write('Part1: identifying barcode sequences in reads.\n')
    
    #deciding if to run multiprocessing or single process for debugging
    from SEAseqLib.mainLibrary import getPairs, Progress

    if indata.debug: #single process // serial
	results=[] # create holder for processed reads
	config.logfile.write('Running in debug mode identifying handles ...\n')
	progress = Progress(config.reads2process, logfile=config.logfile) # creates a progress "bar" thingy
	with progress:
	    for tmp in getPairs(config): #itarate through reads and do the "magicFunction"
		progress.update()
		results.append(foreachread_cluster(tmp))
	config.logfile.write('handle identification finished.\n')

    else: # multiple processes in parallel, create worker pool that iterates through the reads and does the "magicFunction" and sends results to "results"
	import multiprocessing
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=1000000) 
	results = WorkerPool.imap_unordered(foreachread_cluster,getPairs(config),chunksize=100)
    if not indata.debug: config.logfile.write('Running in multiproccessing mode using '+str(indata.cpus)+' processes  ...\n')
    else: config.logfile.write('Running the multiprocessing results handeling in serial ... \n')
    
    from SEAseqLib.mainLibrary import SEAseqSummary
    progress = Progress(config.reads2process, logfile=config.logfile)
    summary = SEAseqSummary()
    f1 = open(config.absolutePath+'/nonCreads.1.fq','w')
    f2 = open(config.absolutePath+'/nonCreads.2.fq','w')
    totalMMs = 0
    totalHandles = 0
    with progress:
	for pair in results:
	    progress.update()
	    summary.add(pair)
	    if pair.handle_start==None or pair.handle_end==None:
		f1.write(pair.r1.header + '\n'+pair.r1.seq+'\n+\n'+pair.r1.qual+'\n')
		f2.write(pair.r2.header + '\n'+pair.r2.seq+'\n+\n'+pair.r2.qual+'\n')
            else:
                totalHandles+=1
                totalMMs += pair.missMatchesInTheHandle
    if not indata.debug:WorkerPool.close()
    if not indata.debug:WorkerPool.join()
    config.outfile.write(str( summary.part1() )+'\n')
    config.outfile.write('\nFound '+str(round(100*float(totalMMs)/float(totalHandles*20),2))+'% missmatches in the C handle sequences.\n')
    config.outfile.write('Total MMs: '+str(totalMMs)+', total handles: '+str(totalHandles)+'.\n')
    config.logfile.write('Part1: finished barcode sequences identified.\n')
    f1.close()
    f2.close()
    
    config.logfile.write('Part2: Clustering the identified barcodes.\n')
    summary.reducebarcodes(config)
    config.logfile.write('Writing cluster info to file.\n')
    f = open(config.clusters_file,'w')
    f.write(str(summary.clusters))
    f.close()
    config.load()
    config.set('clustercount', len(summary.clusters))
    config.save()
    config.logfile.write('Part2: Clustering barcodes END\n----------\n')
    return summary

