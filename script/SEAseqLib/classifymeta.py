def clusterGenerator(config,indata):
    
    import cPickle
    import gzip
    
    #filename=config.path+'/meta.clusters.pickle.gz'
    filename=config.path+'/meta.clusters.pickle'
    #clusterdump = gzip.open(filename,'rb')
    clusterdump = open(filename)
    
    while True:
        try:
            cluster = cPickle.load(clusterdump)
            yield [cluster,config,indata]
        except EOFError:
            config.logfile.write('All clusters read from file.\n')
            break

def foreachCluster(tmp):
    cluster,config,indata = tmp
    perAmpOut = ''
    output = '######## Cluster '+str(cluster.id)+' ###################\n'
    indata.minimum_reads = 5
    indata.minimum_support = 5
    for amplicon in cluster.amplicons.values(): perAmpOut += amplicon.checkmono(indata)
    output += 'There are '+str(cluster.adaptercount)+' illumina adapter reads.\n'
    output += 'There are '+str(cluster.primererrors)+' primer missmatch reads.\n'
    if cluster.ampliconpairs == 0:
        output += '0 amplicon(s) have enough data (>=1 cons with >= '+str(indata.minimum_support)+'% support and >= '+str(indata.minimum_reads)+' reads)\n'
    if cluster.ampliconpairs > 0:
        output += str(cluster.definedampliconcount)+' amplicon(s) have enough data (>=1 cons with >= '+str(indata.minimum_support)+'% support and >= '+str(indata.minimum_reads)+' reads):\n'
    output += perAmpOut + '\n'

    if indata.tempfilefolder:
        import os
        try: os.mkdir(indata.tempfilefolder+'/SEAseqtemp')
        except: pass
        blastfile = open(indata.tempfilefolder+'/SEAseqtemp/blastinput.'+str(cluster.id)+'.fa','w')
    else:blastfile= open(         config.path+'/sortedReads/blastinput.'+str(cluster.id)+'.fa','w')
    
    output+= 'making fasta '+str(cluster.id)+'\n'
    fastaentries = 0
    amplicons = {}
    for amplicon in cluster.definedamplicons.values():
        amplicons[amplicon.type] = {}
        for consensus in amplicon.goodalleles:
            amplicons[amplicon.type][str(consensus.id)] = {'1':None,'2':None}
            r1 = consensus.sequence.seq.split('NNNNNNNNNN')[0]
            r2 = consensus.sequence.seq.split('NNNNNNNNNN')[1]
            blastfile.write('>'+amplicon.type+'.'+str(consensus.id)+'.r1\n'+r1+'\n'+
                            '>'+amplicon.type+'.'+str(consensus.id)+'.r2\n'+r2+'\n')
            fastaentries +=1
    blastfile.close()
    output+=  'done '+str(cluster.id)+'\n'
    
    if fastaentries == 0:
        import os
        os.remove(blastfile.name)
        output += 'No reads in fasta file '+str(cluster.id)+'\n'
        return output


    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
    from cStringIO import StringIO
    import time

    #setting up blast
    BLAST_identity = indata.identity

    database = indata.database
    cline = NcbiblastnCommandline(query=blastfile.name, db=database ,evalue=0.001, outfmt=5, num_threads=8,perc_identity=BLAST_identity)#, out=infile+'.blastout')
    #cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn', out=infile+'.'+config.blastid+'.blastout')
    
    output+='Starting BLAST '+str(cluster.id)+'\n'
    blast_handle = cline.__call__()
    output+=str(cluster.id)+' '+ str(blast_handle)[0:100].replace('\n','--NEWLINE--')+'\n'
    output+='BLAST search done '+str(cluster.id)+'\n'

    import os
    os.remove(blastfile.name)

    blast_handle = StringIO(blast_handle[0])
    blast_handle.seek(0)
    records = NCBIXML.parse(blast_handle)
    
    from SEAseqLib.mainLibrary import gi2orgname
    local_gi2org = {}

    #for the random match estimation
    tmp_matches = {}

    if indata.gidatabase:
	tmp_file = open(indata.gidatabase)
	tmp_string = tmp_file.read()
	tmp_file.close()
	local_gi2org = eval(tmp_string)
    else:local_gi2org = {}
    
    identity_cutoff = 99
    alignment_length_cutoff = 95
    for blast_record in records:

        amptype = blast_record.query.split('.')[0]
        allele =blast_record.query.split('.')[1]
        readnumber =blast_record.query.split('.')[2]
        
        amplicons[amptype][allele][readnumber] = blast_record
        if blast_record == None: output += amptype+' '+allele+' '+readnumber+'ajajajaj\n'
        
        #for alignment in blast_record.alignments[:10]:
        #    for hsp in alignment.hsps[:1]:
        #        perc_identity = float(hsp.identities) 	/	float(hsp.align_length)	*100
        #        perc_coverage = float(hsp.align_length)	/	float(blast_record.query_letters)	*100
        #        gi_number = alignment.title.split(' ')[1].split('|')[1]
        #        try:
        #            organism = local_gi2org[gi_number]
        #        except KeyError:
        #            organism = gi2orgname(gi_number)
        #            local_gi2org[gi_number] = organism
        #        if perc_identity >= identity_cutoff and perc_coverage >= alignment_length_cutoff:
        #            print perc_identity,perc_coverage,organism
    
    for amplicon in amplicons:
        broken = False
        consensuses = amplicons[amplicon]

        output += '\t'+amplicon     +'\n'

        for consensus in consensuses:
            try:
                r1 = consensuses[consensus]['1']
                if r1 == None:
                    output += '\t\tno read1 matches\n'
                    continue
            except KeyError:
                output += '\t\tConsensus cannot be read as read 1 sequence are not available from infile\n'
                output += 'WARNING: Running on a uncomplete input file! Maybe the SEAseq meta program is still runnning?\n'
                continue
            try:
                r2 = consensuses[consensus]['2']
                if r2 == None:
                    output += '\t\tno read2 matches\n'
                    continue
            except KeyError:
                output += '\t\tConsensus cannot be read as read 2 sequence are not available from infile\n'
                output += 'WARNING: Running on a uncomplete input file! Maybe the SEAseq meta program is still runnning?\n'
                continue

            #amplicon	= 	r1.query.split('.')[1].split('=')[1]
            #support	= float(r1.query.split('r1.')[-1].split('_')[0].replace('%','')	)
            #readpop = int(r1.query.split('r1.')[-1].split('_')[-1].replace('reads','')	)

            output +=       '\t\tConsensus '+str(consensus)+': '+amplicon+' supported by '+str(support) +'% ('+str(readpop)+')'     +'\n'

            in_r1 = {}
            for alignment in r1.alignments:
                for hsp in alignment.hsps:
                    perc_identity = float(hsp.identities) 	/	float(hsp.align_length)	*100
                    perc_coverage = float(hsp.align_length)	/	float(r1.query_letters)	*100
                    gi_number = alignment.title.split(' ')[1].split('|')[1]
                    try:
                        organism = local_gi2org[gi_number]
                    except KeyError:
                        organism = gi2orgname(gi_number)
                        local_gi2org[gi_number] = organism
                    #print organism
                    output += organism
                    if perc_identity >= identity_cutoff and perc_coverage >= alignment_length_cutoff: in_r1[organism] = alignment

            in_r2 = {}
            for alignment in r2.alignments:
                for hsp in alignment.hsps:
                    perc_identity = float(hsp.identities)	/	float(hsp.align_length)	*100
                    perc_coverage = float(hsp.align_length)	/	float(r2.query_letters)	*100
                    gi_number = alignment.title.split(' ')[1].split('|')[1]
                    try:
                        organism = local_gi2org[gi_number]
                    except KeyError:
                        organism = gi2orgname(gi_number)
                        local_gi2org[gi_number] = organism
                    #print organism
                    output += organism
                    if perc_identity >= identity_cutoff and perc_coverage >= alignment_length_cutoff: in_r2[organism] = alignment


            for organism in in_r1:
                if organism in in_r2:
                    in_both_reads[amplicon].append(organism)

            for organism in in_both_reads[amplicon]:
                output +=       '\t\t\t'+organism     +'\n'
                tmp_matches[amplicon][organism] = True
            if not in_both_reads[amplicon]:
                output +=       '\t\t\tNo alignment supported by both reads with >='+str(identity_cutoff)+'% identity and '+str(alignment_length_cutoff)+'% alignment length coverage'     +'\n'
                if amplicon == '16S':
                    noblasthit_16s +=1
                    broken = True;
                elif amplicon == 'ITS':
                    noblasthit_its += 1
                    broken = True;
                    break;
        if broken: break
    
    
    
    return output


def classifymeta(indata):

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    import os
    if os.path.exists(config.outfile): os.remove(config.outfile)
    config.openconnections()
    writelogheader(config.logfile)

    # settings
    config.load()

    if indata.debug: #single process // serial
	config.logfile.write('debugging: ');
        import sys
	sys.stdout.write('debugging: ')
	config.logfile.write('Running in debug mode ')
	results=[] # create holder for processed reads
	#progress = Progress(config.clustercount, logfile=config.logfile) # creates a progress "bar" thingy
	#with progress:
        tmcounter = 0
	for cluster in clusterGenerator(config, indata):
		#progress.update()
		results.append(foreachCluster(cluster))
                tmcounter +=1
                if tmcounter == 20: break
	config.logfile.write('finished, making summary ... \n')
    else: # multiple processes in parallel
	import multiprocessing
	WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10000)
	results = WorkerPool.imap_unordered(foreachCluster,clusterGenerator(config,indata),chunksize=1)
	#results = WorkerPool.imap(          foreachcluster,clusterGenerator(config),chunksize=1)

    # will paralellise this when stuff works
    tmcounter = 0
    for cluster in results:
        print cluster
        tmcounter +=1
        if tmcounter == 20: break

    #for the random match estimation
    matches = {}

    #from SEAseqLib.mainLibrary import Progress
    #config.logfile.write('Parsing BLASTreport ....\n')
    #progress = Progress(len(data),unit='cluster',logfile=config.logfile)
    #with progress:

    import sys
    sys.exit(0)

    #CUTTNIG AND PASTING FROMHERE:::

    for pelleplut in list:
	for cluster_id, amplicons in data.iteritems():
	    
	    
	    progress.update()

	    config.outfile.write(      '### Cluster '+str(cluster_id)+' ###'     +'\n')
	    identity_cutoff = BLAST_identity#%
	    alignment_length_cutoff = indata.length#95#%

	    #for amplicon, consensuses in amplicons.iteritems():
	    in_both_reads = {'ITS':[],'16S':[]}


	    if tmp_matches['ITS'] and tmp_matches['16S']:
		matches[cluster_id] = tmp_matches

	    hitagree = False
	    config.outfile.write(      '\nSupported by both:'     +'\n')
	    for organism in in_both_reads['ITS']:
		if organism in in_both_reads['16S']:
		    config.outfile.write(      '\t'+organism+' is present in both'     +'\n')
		    hitagree = True

	    if hitagree:	blasthitsagree += 1
	    elif not broken:
		hitsdonotagree += 1
		config.outfile.write(      '\tno hit is present in both'     +'\n')
	    config.outfile.write(      ''     +'\n')

	assert len(data) - noblasthit_its - noblasthit_16s - blasthitsagree - hitsdonotagree == 0, '\n\nError '+str(len(data))+' - '+str(noblasthit_its)+' - '+str(noblasthit_16s)+' - '+str(blasthitsagree)+' - '+str(hitsdonotagree)+' != 0\n\n'
	config.outfile.write(      'out of '+str(len(data))+' analyzed clusters (with both amplicons present and monoclonal for both) were:'     +'\n')
	config.outfile.write(      str(noblasthit_its) +' ('+str(round(100*float(noblasthit_its)/float(len(data)),2))+'%) removed because the ITS sequence gave no hits supported by both reads with >='+str(identity_cutoff)+'% identity and '+str(alignment_length_cutoff)+'% alignment length coverage.'     +'\n')
	config.outfile.write(      str(noblasthit_16s) +' ('+str(round(100*float(noblasthit_16s)/float(len(data)),2))+'%) removed because the 16S sequence gave no hits supported by both reads with >='+str(identity_cutoff)+'% identity and '+str(alignment_length_cutoff)+'% alignment length coverage.'     +'\n')
	config.outfile.write(      str(hitsdonotagree+blasthitsagree) +' ('+str(round(100*float(hitsdonotagree+blasthitsagree)/float(len(data)),2))+'%) have hits for both ITS and 16S.'     +'\n')
        if hitsdonotagree+blasthitsagree > 0:
            config.outfile.write(  str(hitsdonotagree) +' ('+str(round(100*float(hitsdonotagree)/float(len(data)),2))+'%) removed because none of the ITS and 16S hits did agree.'+	 ' ('+str(round(100*float(hitsdonotagree)/float(hitsdonotagree+blasthitsagree),2))+'% out of clusters with hits for both)'     +'\n')
            config.outfile.write(  'For '+str(blasthitsagree) +' ('+str(round(100*float(blasthitsagree)/float(len(data)),2))+'%) clusters did the ITS and 16S hits agree at least once.'+' ('+str(round(100*float(blasthitsagree)/float(hitsdonotagree+blasthitsagree),2))+'% out of clusters with hits for both)'     +'\n')
	
            groups_its = [] 
            groups_16s = []
            for cluster in matches:
                groups_its.append(matches[cluster]['ITS'])
                groups_16s.append(matches[cluster]['16S'])
            
            import random
            total_tries = 10000
            random_hits = 0
            for i in xrange(total_tries):
                hits_its = random.choice(groups_its)
                hits_16s = random.choice(groups_16s)
                
                at_least_one_match = False
                for hit in hits_its:
                    if hit in hits_16s: at_least_one_match = True
                
                if at_least_one_match: random_hits += 1
            
            config.outfile.write('Random matches between hits in 16S/ITS hit groups, made '+str(total_tries)+' tries and got '+str(random_hits)+' hits ('+str(round(100*float(random_hits)/float(total_tries),2))+'%).\n')
	else:
            config.outfile.write('No clusters found with hits for both were found'     +'\n')
        
    config.logfile.write('Classification done.\n')
    return 0