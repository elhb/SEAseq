import sys


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
		#config.logfile.write('.');
		#sys.stdout.write('.')
	config.logfile.write('finished, making summary ... \n')
    else: # multiple processes in parallel
	import multiprocessing
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10000)
	results = WorkerPool.imap_unordered(foreachcluster_meta,clusteriterator(clusterq,indata),chunksize=1)
	#results = WorkerPool.imap(foreachcluster_meta,clusteriterator(clusterq, indata),chunksize=10)

    nonecount = 0
    processed = 0
    lowread = 0
    onlyjunk = 0
    typecounter = {'ITS':0,'16S':0,'undefined':0,'both ITS and 16S':0}
    monoclonal = {'ITS':0,'16S':0,'undefined':'NA','both':{'only ITS':0,'only 16S':0,'None':0,'both':0}}
    for amptype in ['ecoli','myco','lambda','m13']:
        typecounter[amptype] = 0
        monoclonal[amptype] = 0
        
    if not indata.debug: config.logfile.write('Part1: Per cluster action '+str(indata.cpus)+' processes  ...\n')
    progress = Progress(config.clustercount, logfile=config.logfile, unit='cluster',mem=True)
    statstable = open(config.path+'/meta.statstable','w',1)
    statsheader = ['clusterid','number of reads in total','number of adaper reads','number of strange primers','its reads','16s reads','its','16s','its monoclonal','16s monoclonal','number of consensus types','number of consensus types with good support']
    for amptype in ['ecoli','myco','lambda','m13']: statsheader.append(amptype+' reads');statsheader.append(amptype+' monoclonal');statsheader.append(amptype)
    statstable.write('\t'.join(statsheader))
    seqdump = open(config.path+'/meta.sequences.fa','w',1)
    with progress:
	for tmp in results:
	    seqdict = {}
	    progress.update()
	    if tmp == 'INITIAL':
		nonecount+=1
		if nonecount > 1: print 'Error!!! to many none clusters!!'; raise ValueError
		return_info = False
	    elif tmp[0] == 'LOW READ':
		processed += 1
		lowread += 1
		return_info = tmp[1]
                output = tmp[2]
                config.outfile.write( output )
	    elif tmp[0] == 'ONLY JUNK':
		processed += 1
		onlyjunk += 1
		return_info = tmp[1]
	    else:
		[output, _its, _16s, return_info, seqdict] = tmp
		config.outfile.write( output )
		processed += 1

		#print tmp[1:]		
		#if _its != None and _16s != None:
                if return_info['its'] and return_info['16s']:
		    typecounter['both ITS and 16S']+=1
		#    if _its and _16s: monoclonal['both']['both'] += 1
		#    if _its and not _16s: monoclonal['both']['only ITS'] += 1
		#    if _16s and not _its: monoclonal['both']['only 16S'] += 1
		#    if not _16s and not _its: monoclonal['both']['undefined'] += 1
		#elif _its != None:
                    if return_info['its monoclonal'] and return_info['16s monoclonal']: monoclonal['both']['both'] += 1
		    if return_info['its monoclonal'] and not return_info['16s monoclonal']: monoclonal['both']['only ITS'] += 1
		    if return_info['16s monoclonal'] and not return_info['its monoclonal']: monoclonal['both']['only 16S'] += 1
		    if not return_info['16s monoclonal'] and not return_info['its monoclonal']: monoclonal['both']['None'] += 1
		elif return_info['its']:
		    typecounter['ITS'] += 1
		    if return_info['its monoclonal']: monoclonal['ITS'] += 1
		elif return_info['16s']:
		    typecounter['16S'] += 1
		    if return_info['16s monoclonal']: monoclonal['16S'] += 1
		else:
                    nonefound = True
                    for amptype in ['ecoli','myco','lambda','m13']:
                        if return_info[amptype]:
                            typecounter[amptype] +=1
                            if return_info[amptype+' monoclonal']: monoclonal[amptype]+=1
                            nonefound = False;break
		    if nonefound: typecounter['undefined'] += 1

		# SEQDICT printing to file ONLY for both monoclonal
		if _its and _16s:
		    for cons_type in seqdict:
			for consensus in seqdict[cons_type]:
			    seqdump.write(seqdict[cons_type][consensus])

	    if return_info:
		statstable.write( '\n'+'\t'.join([str(return_info[stat]) for stat in statsheader]) )
		#assert len(return_info) == len(statsheader)
	seqdump.close()
	statstable.close()

    reader.join()
    
    defined_clust = 0
    defined_clust_mono = 0
    config.outfile.write('##### SUMMARY #####'+'\n')
    config.outfile.write(  str(processed)+ ' clusters processed, out of these were:'+'\n')
    config.outfile.write(  str(onlyjunk) + ' ('+str(round(100*float(onlyjunk)/float(processed),2))+'%) clusters of only adapter sequences or faulty primers'+'\n')
    config.outfile.write(  str(lowread)  + ' ('+str(round(100*float(lowread)/ float(processed),2))+'%) clusters with to few reads to be analysed'+'\n')
    for name,count in typecounter.iteritems():
	if name != None and name != 'undefined':defined_clust += count
        if name == 'both ITS and 16S':
	    config.outfile.write(  str(count)+' ('+str(round(100*float(count)/float(processed),2))+'%) '+ name+' '+ 'whereof:'+'\n')
	    for name2,count2 in monoclonal['both'].iteritems():
                percentage2 = 0
                if count != 0: percentage2 = round(100*float(count2)/float(count),2)
		config.outfile.write(  '\t'+' '+str(count2)+' ('+str(percentage2)+'%) '+'were monoclonal for'+' '+name2+'\n')
                if name2 == 'both': defined_clust_mono += count2
	    continue
        if name != None and name != 'undefined':
            defined_clust_mono += monoclonal[name]
            percentage = 0
            if count: percentage = round(100*float(monoclonal[name])/float(count),2)
            config.outfile.write(  str(count)+' ('+str(round(100*float(count)/float(processed),2))+'%) '+ name+' '+ 'whereof'+' '+ str(monoclonal[name])+' ('+str(percentage)+'%) were monoclonal'+'\n')
        else:config.outfile.write(  str(count)+' ('+str(round(100*float(count)/float(processed),2))+'%) '+ name+' '+ 'whereof'+' '+ str(monoclonal[name])+' ( NA %) were monoclonal'+'\n')

    
    tmppercentage = 0
    if defined_clust: tmppercentage = round(100*float(defined_clust_mono)/float(defined_clust),2)
    config.outfile.write(
        str(defined_clust)+' ('+str(round(100*float(defined_clust)/float(processed),2))+'%) has at least one defined amplicon out of these are '+str(defined_clust_mono)+' ('+str(tmppercentage)+'%) monoclonal for the defined amplicon(s)\n'+
        '(ie there is only one consensus sequence of that type with more than '+str(indata.minimum_reads)+' reads and '+str(indata.minimum_support)+'% support, clustering done with '+str(indata.clustering_identity)+'% identity cutoff)\n'
        )
	
    # create a nice run summary
    config.logfile.write('Done.\n')
    return 0

def foreachcluster_meta(cluster_pairs):
    indata = cluster_pairs[1]
    cluster = cluster_pairs[0][0]
    config = cluster_pairs[0][1]
    verb = 3
    
    return_info = {
	'clusterid':cluster.id,
	'number of reads in total':cluster.readcount,
	'number of adaper reads':None,
	'number of strange primers':None,
	'number of consensus types':None,
	'number of consensus types with good support':None
    }
    for amptype in ['ecoli','myco','lambda','m13','its','16s']:
        return_info[amptype+' reads'] = None
        return_info[amptype] = None
        return_info[amptype+' monoclonal'] = None    
    
    from SEAseqLib.mainLibrary import readpair, sequence, UIPAC2REGEXP

    adaptercount = 0
    primererror = 0
    lowreadcutoff = 1
    output = ''
    output =  '\n###--- Cluster number '+str(cluster.id)+' -> '+str(cluster.readcount)+' pairs. ---###\n'
    if cluster.readcount < lowreadcutoff:
	output += 'Less than '+str(lowreadcutoff)+' read pairs.\n'
	return ['LOW READ',return_info,output]
    else:
        output += cluster.createtempfile(config)
        cluster.clusterreadpairs(config, indata)
        cluster.loadconsensuses(config)
        cluster.loadconsensusalignemnts(config)

	perccutoff  = indata.minimum_support#5.0
	countcutoff = indata.minimum_reads#5
	#remove less than 5% (should I remove singletons aswell?)
	tmp = {}
	typecounter={}
	for cons_type, consensuses in types.iteritems():
	    tmp[cons_type] = {}
	    tmp[cons_type]['mono'] = False
	    types[cons_type]['mono'] = False
	    for consensus, data in consensuses.iteritems():
		if consensus == 'total' or consensus == 'mono': tmp[cons_type][consensus] = data;continue
		#print cons_type,consensus,data, consensuses['total']
		percentage = round(100*float(data['support'])/consensuses['total'],2)
		if percentage >= perccutoff and data['support'] >= countcutoff:
		    tmp[cons_type][consensus] = data
		    typecounter[cons_type] = True
	types = tmp

	seqdict = {}
	# make cluster summary and seqdict for dumping to fasta
	output += 'There are '+str(adaptercount)+' illumina adapter reads.\n'
	output += 'There are '+str(primererror)+' primer missmatch reads.\n'
	#output += 'There are '+str(len(types))+' type(s) of amplicon(s) in cluster:\n'
	output += str(len(typecounter))+' amplicon(s) have enough data (>=1 cons with >= '+str(perccutoff)+'% support and >= '+str(countcutoff)+' reads):\n'
	for cons_type, consensuses in types.iteritems():
	    types[cons_type]['mono'] = False
	    output += '\t'+cons_type+' '+str(consensuses['total'])+' reads in total.\n'
	    seqdict[cons_type] = {}
	    for consensus, data in consensuses.iteritems():
		if consensus == 'total' or consensus == 'mono':continue
		percentage = round(100*float(data['support'])/consensuses['total'],2)
		output += '\t\tConsensus '+consensus+' supported by '+str(percentage)+'% of readpop ('+str(data['support'])+' reads)\t'+data['sequence'].replace('-','')+'\n'
		#if percentage >= 95.0: types[cons_type]['mono'] = True
		try:
		    seqdict[cons_type][consensus] = '\n>cluster='+str(cluster.id)+'.amplicon='+str(cons_type)+'.consensus='+str(consensus)+'_r1.'+str(percentage)+'%_of_'+str(consensuses['total'])+'reads\n'+data['sequence'].replace('-','').split('NNNNNNNNNN')[0]
		    seqdict[cons_type][consensus] += '\n>cluster='+str(cluster.id)+'.amplicon='+str(cons_type)+'.consensus='+str(consensus)+'_r2.'+str(percentage)+'%_of_'+str(consensuses['total'])+'reads\n'+data['sequence'].replace('-','').split('NNNNNNNNNN')[1]
		except IndexError:
		    output += 'WARNING: consensus sequence not properly splittable into r1 and r2, not dumping this cluster ('+str(cluster.id)+').\n'
		if len(consensuses)-2 == 1: types[cons_type]['mono'] = True
	output += '\n'

        return_info['number of adaper reads'] = cluster.adaptercount
	return_info['number of strange primers'] = primererror
	return_info['number of consensus types'] = len(types)
	return_info['number of consensus types with good support'] = len(typecounter)
	return_info['its reads'] = types['ITS']['total']
	return_info['16s reads'] = types['16S']['total']
	return_info['its'] = bool('ITS' in typecounter)
	return_info['16s'] = bool('16S' in typecounter)
	return_info['its monoclonal'] = types['ITS']['mono']
	return_info['16s monoclonal'] = types['16S']['mono']
        for amptype in ['ecoli','myco','lambda','m13']:
            return_info[amptype+' reads'] = types[amptype]['total']
            return_info[amptype] = bool(amptype in typecounter)
            return_info[amptype+' monoclonal'] = types[amptype]['mono']

        cluster.removetempfiles(config)

	_its = None
	_16s = None
        _amptype = {'its':None,'16s':None,'ecoli':None,'myco':None,'m13':None,'lambda':None}
	if 'ITS' in typecounter:
	    _its = types['ITS']['mono']
	    if _its: output+='Cluster is monoclonal for ITS.\n'
	if '16S' in typecounter:
	    _16s = types['16S']['mono']
	    if _16s:output+='Cluster is monoclonal for 16S.\n'
        for amptype in ['ecoli','myco','lambda','m13']:
            if amptype in typecounter:
                _amptype[amptype] = types[amptype]['mono']
                if _amptype[amptype]:output+='Cluster is monoclonal for '+amptype+'.\n'
        

	return [output, _its,_16s,return_info,seqdict]#str(cluster_pairs[1][0].cid)+' has '+str(len(cluster_pairs[1]))+' read pairs'

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