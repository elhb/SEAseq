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
    if cluster.readcount<lowreadcutoff:
	output += 'Less than '+str(lowreadcutoff)+' read pairs.\n'
	return ['LOW READ',return_info,output]
    else:
	
        import re
        readsbyheader = {}
	tmpcounter = 0
	int2header = {}
        
        
        ############################################################################################
        #          GO TROUGH ALL READS PAIRS IN CLUSTER NAD PREPARE FOR CONSENSUS CLUSTERING       #
        ############################################################################################
	if verb >=3: output += 'PAIRS:\n'
	f = open(config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa','w')
	for pair in cluster.readpairs:
	    
	    tmpcounter += 1
	    if verb >=3: output += str(tmpcounter)+'\t'
            pair.primererror = None
	    
	    #identify subparts of read
	    C_HANDLE = sequence('c handle',"CTAAGTCCATCCGCACTCCT","CTAAGTCCATCCGCACTCCT")
	    pair.identify(C_HANDLE, config)
	    pair.getN15()
	    pair.identifyIllumina(config)
	    
	    #match adaptersequence
	    if pair.isillumina:
		adaptercount +=1;
		if verb >=3: output+='ADAPTER\n';
		continue
	    
            # Match the primer pairs
            for name, primerpair in config.primerpairs.iteritems():
                pair.matchprimerpair(primerpair)

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
            if verb >=3: output += pair.p1+'\t'+pair.p2+'\t';
            
	    # check that primers match
	    if pair.p1 == '???' or pair.p1 != pair.p2 or len(pair.matchingprimerpairs) != 1:
		primererror+=1;
                if pair.p1 == '???':                        pair.primererror = 'primers-not-identifiable'
                elif pair.p1 != pair.p2:                    pair.primererror = 'fwd-rev-pair-missmatch'
                elif len(pair.matchingprimerpairs) != 1:    pair.primererror = 'more-than-one-pair-match'
                if verb >=3:
		    output+=pair.primererror;
		    output += pair.r1.seq +' '+ pair.r2.seq+'\n'
		continue

	    #look for unexpected primer sequences
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
		    if verb >=3:
			output+='PRIMER ODD COMBO\t';
			output += pair.r1.seq +' '+ pair.r2.seq+'\n'
		    continue

	    # Add sequences to output
	    if verb >=3: output += pair.r1.seq +' '+ pair.r2.seq+'\n'
	    
	    # prepare for clustering by writing read pair to tempfile and saving temporary id number and mappinf header to pair-object
            int2header[tmpcounter] = pair.header
	    readsbyheader[pair.header] = pair
            tem_seq = pair.r1.seq[pair.handle_end:][len( config.primerpairs[pair.p1].fwd )+1:]+'NNNNNNNNNN'+pair.r2.revcomp().seq[:-(len(    config.primerpairs[pair.p1].rev   )+1)]
	    f.write('>'+str(tmpcounter)+'\n'+ tem_seq +'\n')
	f.close()
	return_info['number of adaper reads'] = cluster.adaptercount
	return_info['number of strange primers'] = primererror

        assert adaptercount == cluster.adaptercount, '\n\n\t\t\t######### ERROR-SCHMERROR_1!!!! ##########\n\n'
        assert primererror  == cluster.primererrors, '\n\n\t\t\t######### ERROR-SCHMERROR_2!!!! ##########\n\n'

        ############################################################################################
        #                                DO CONSENSUS CLUSTERING                                   #
        ############################################################################################

	#check that there is data to work with
	if adaptercount+primererror == cluster.readcount:
	    output += 'All adapter and/or primer error.\n'
	    import os
	    os.remove( f.name )
	    return ['ONLY JUNK',return_info]

	# Cluster Read pairs
	import subprocess
	from cStringIO import StringIO
	import time
	import multiprocessing
	tempo = time.time()
	cdhit = subprocess.Popen( ['cd-hit-454','-i',config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa','-o',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa','-g','1','-c',str(indata.clustering_identity/100.0)], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	cdhit_out, errdata = cdhit.communicate()
	if cdhit.returncode != 0:
                print 'cmd: '+' '.join( ['cd-hit-454','-i',config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa','-o',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa','-g','1','-c',str(indata.clustering_identity/100.0)])
		print 'cd-hit cluster='+str(cluster.id)+' view Error code', cdhit.returncode, errdata
		sys.exit()
	seconds = round(time.time()-tempo,2)

	# Build consensus sequences for read pair clusters
	ccc = subprocess.Popen( ['cdhit-cluster-consensus',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr',config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned'], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	ccc_out, errdata = ccc.communicate()
	if ccc.returncode != 0:
                print 'cmd: '+' '.join( ['cdhit-cluster-consensus',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr',config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus',config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned'])
		print 'cluster='+str(cluster.id)+' cdhit-cluster-consensus view Error code', ccc.returncode, errdata
		print ccc_out
		sys.exit()

	# output info from temporary files
	if verb >=5: 
	    for info in [
		[config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus.fasta','\nCD-HIT consensus sequences:\n'],
		[config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr','\nClustering details:\n'],
		[config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned','\nAlignment details:\n']
		]:
		filename , message =info
		f = open(filename)
		tmp = f.read()
		output += message
		output += tmp
		f.close()
	
	# get alignments from file
	f = open(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned')
	data = f.read()
	f.close()
	consensuses = {}
	for cluster_aln in data.split('===========================================\n'):
	    for part in cluster_aln.split('\n\n'):
		for line in part.split('\n'):
		    line=line.rstrip()
		    if not line: continue
		    if line[0] == 'A':
			#print '### HERE: '+line
			consensusid=line.split(' ')[-1].split(':')[0];
			consensuses[consensusid] = {}
			continue
		    read_id = line.split(':  +  ')[0].rstrip()
		    try: seq = line.split(':  +  ')[1]
		    except: print line+'<- Problematic';raise ValueError
		    if read_id == 'Consensus': seq = seq.split(' ')[0]
		    try:		consensuses[consensusid][read_id][0] += seq
		    except KeyError:	consensuses[consensusid][read_id] = [seq,'NotSet']

#	get identity from file
	f = open(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr')
	data = f.read()
	f.close()
	for consensus_identities in data.split('>Cluster '):
		consensusid = consensus_identities.split('\n')[0]
		if consensusid not in consensuses: consensuses[consensusid] = {}
		for line in consensus_identities.split('\n')[1:]:
		    line=line.rstrip()
		    if not line: continue
		    read_id     = line.split('>')[1].split('.')[0]
		    identity = line.split('/')[-1]
		    try:		consensuses[consensusid][read_id][1] = identity
		    except KeyError:	consensuses[consensusid][read_id] = ['SEQUENCE NOT LOADED',identity]

	##load singelton consensus sequences
	f = open(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus.fasta')
	data = f.read()
	f.close()
	for consensus_identities in data.split('>')[1:]:
		if consensus_identities[0] == 'c': continue
		consensusid = consensus_identities.split(' ')[0].split('_')[-1]
		read_id     = consensus_identities.split(' ')[1].split('\n')[0]
		if consensusid not in consensuses: consensuses[consensusid] = {}
		if 'Consensus' not in consensuses[consensusid]: consensuses[consensusid]['Consensus'] = ['', 'NotSet']
		if consensuses[consensusid][read_id][0] == 'SEQUENCE NOT LOADED': consensuses[consensusid][read_id][0] = ''
		for line in consensus_identities.split('\n')[1:]:
		    line=line.rstrip()
		    consensuses[consensusid][read_id][0] += line
		    consensuses[consensusid]['Consensus'][0] += line

	#make output for alnignments
	if verb >=2: output += '\n'
	types = {'ITS':{'total':0},'16S':{'total':0},'ecoli':{'total':0},'myco':{'total':0},'m13':{'total':0},'lambda':{'total':0}}
	for consensusid in consensuses:
	    if  not consensusid: continue
	    primerpairs = []
	    if 'Consensus' in consensuses[consensusid]:
		if verb >=2: output += 'Consensus number '+consensusid+' from '+str(len(consensuses[consensusid])-1)+' read pairs'
		if verb >=2: output += ':\t'+consensuses[consensusid]['Consensus'][0]+'\n'
	    for read_id, tmp in consensuses[consensusid].iteritems():
		[seq, identity] = tmp
		if identity[-1] == '*': identity = 'SEED'
	        if read_id == 'Consensus': continue#output +=read_id+'\t'+seq+'\n'
	        else:
		    #if verb >=2: output += str(read_id)+'\t'+int2header[int(read_id)].split('_')[0]+'   \t'+identity+'\t'+reads[int2header[int(read_id)]].p1+'\t'+seq
		    if verb >=2: output += 'Read pair id = '+str(read_id)+'    \t'+identity+'\t'+readsbyheader[int2header[int(read_id)]].p1+'\t'+seq
		    if verb >=2: output += '\n'
		    primerpairs.append(readsbyheader[int2header[int(read_id)]].p1)
	    if   primerpairs and primerpairs.count(primerpairs[0]) == len(primerpairs):
		if verb >=2: output += 'consensus is '+primerpairs[0]
	    elif primerpairs and primerpairs.count(primerpairs[0]) != len(primerpairs):
		if verb >=2: output += 'WARNING: mixed consensus clustering!'
	    try:
		types[primerpairs[0]][consensusid]={'sequence':consensuses[consensusid]['Consensus'][0], 'support':len(consensuses[consensusid])-1}
		types[primerpairs[0]]['total']+=len(consensuses[consensusid])-1
	    except KeyError: raise ERIKERROR
#		types[primerpairs[0]] = {consensusid:{'sequence':consensuses[consensusid]['Consensus'][0], 'support':len(consensuses[consensusid])-1}}
#		types[primerpairs[0]]['total']=len(consensuses[consensusid])-1
	    if verb >=2: output += '\n\n'

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
		#    try:
		#	seqdict[cons_type][consensus] = '\n>cluster='+str(cluster.id)+'.amplicon='+str(cons_type)+'.consensus='+str(consensus)+'_r1.'+str(percentage)+'%_of_'+str(consensuses['total'])+'reads\n'+data['sequence'].replace('-','').split('NNNNNNNNN')[0]
		#	seqdict[cons_type][consensus] += '\n>cluster='+str(cluster.id)+'.amplicon='+str(cons_type)+'.consensus='+str(consensus)+'_r2.'+str(percentage)+'%_of_'+str(consensuses['total'])+'reads\n'+data['sequence'].replace('-','').split('NNNNNNNNN')[1]
		#    except IndexError:
		#	print '\n\n###Cluster ==',cid;
		#	print 'Erro when splitting sequence:'
		#	print 'cluster='+str(cluster.id)+'.amplicon='+str(cons_type)+'.consensus='+str(consensus)+'.'+str(percentage)+'%_of_'+str(consensuses['total'])+'reads'
		#	print data['sequence']
		#	raise ValueError
		if len(consensuses)-2 == 1: types[cons_type]['mono'] = True
	output += '\n'

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


	import os
	os.remove(config.path+'/sortedReads/temporary.'+str(cluster.id)+'.fa')
	os.remove(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa')
	os.remove(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.fa.clstr')
	os.remove(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.consensus.fasta')
	os.remove(config.path+'/sortedReads/cluster.'+str(cluster.id)+'.aligned')

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