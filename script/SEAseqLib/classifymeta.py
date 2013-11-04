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
    
    # unpack input information
    cluster,config,indata = tmp

    # initiate the output string and set some variables
    output = '######## Cluster '+str(cluster.id)+' ###################\n'
    output += 'Cluster contains '+str(cluster.readcount)+' raw reads.\n'
    #indata.minReadCountPerConsensus = 5
    #indata.minReadPopSupportConsensus = 5
    #config.minBlastIdentity = 99
    #config.minBlastCoverage = 95

    # get the cluster and per amplicon information and append it to output
    perAmpOut = ''
    for amplicon in cluster.amplicons.values(): perAmpOut += amplicon.checkmono(config)
    output += 'There are '+str(cluster.adaptercount)+' illumina adapter reads.\n'
    output += 'There are '+str(cluster.primererrors)+' primer missmatch reads.\n'
    output += 'There are '+str(cluster.readcount-cluster.primererrors-cluster.adaptercount)+' "good reads".\n'
    if cluster.ampliconpairs == 0:
        output += '0 amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% support and >= '+str(config.minReadCountPerConsensus)+' reads)\n'
    if cluster.ampliconpairs > 0:
        output += str(cluster.definedampliconcount)+' amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% support and >= '+str(config.minReadCountPerConsensus)+' reads):\n'
    output += perAmpOut + '\n'
    
    if cluster.definedampliconcount == 0:
        output += '\t0 defined amplicons, cluster '+str(cluster.id)+'\n'
        return [output, cluster]

    # creating the temporary fasta file
    if indata.tempFileFolder:
        import os
        try: os.mkdir(indata.tempFileFolder+'/SEAseqtemp')
        except: pass
        blastfile = open(indata.tempFileFolder+'/SEAseqtemp/blastinput.'+str(cluster.id)+'.fa','w')
    else:blastfile= open(         config.path+'/sortedReads/blastinput.'+str(cluster.id)+'.fa','w')    
    output+= '\tmaking fasta for cluster '+str(cluster.id)+'\n'
    fastaentries = 0
    amplicons = {}
    for amplicon in cluster.definedamplicons.values():
        amplicons[amplicon.type] = {}
        for consensus in amplicon.goodalleles:
            amplicons[amplicon.type][str(consensus.id)] = {'r1':None,'r2':None}
            r1 = consensus.sequence.seq.split('NNNNNNNNNN')[0]
            r2 = consensus.sequence.seq.split('NNNNNNNNNN')[1]
            blastfile.write('>'+amplicon.type+'|tempSep|'+str(consensus.id)+'|tempSep|r1\n'+r1+'\n'+
                            '>'+amplicon.type+'|tempSep|'+str(consensus.id)+'|tempSep|r2\n'+r2+'\n')
            fastaentries +=1
    blastfile.close()
    output+=  '\tfasta file created (cluster '+str(cluster.id)+')\n'
    
    # check that there were reads to align
    if fastaentries == 0:
        import os
        os.remove(blastfile.name)
        output += '\tNo reads in fasta file, cluster '+str(cluster.id)+'\n'
        return [output, cluster]


    # Align consensus sequences by blast
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
    from cStringIO import StringIO
    import time
    #setting up blast
    BLAST_identity = config.minBlastIdentity
    database = config.blastDb
    cline = NcbiblastnCommandline(query=blastfile.name, db=database ,evalue=0.001, outfmt=5, num_threads=8,perc_identity=BLAST_identity)#, out=infile+'.blastout')
    #cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn', out=infile+'.'+config.blastid+'.blastout')
    import time
    output+='\tStarting BLAST, cluster '+str(cluster.id)+'\n'
    starttime = time.time()
    blast_handle = cline.__call__()
    #output+=str(cluster.id)+' '+ str(blast_handle)[0:100].replace('\n','--NEWLINE--')+'\n'
    output+='\tBLAST search finished, cluster '+str(cluster.id)+', run time was '+str(round(time.time()-starttime))+'s.\n\n'
    # remove temporary file
    import os
    os.remove(blastfile.name) 
    #convert blast output to parsable handle
    blast_handle = StringIO(blast_handle[0])
    blast_handle.seek(0)
    records = NCBIXML.parse(blast_handle)
    
    #load the mapping to organism name dictionary
    from SEAseqLib.mainLibrary import gi2orgname
    local_gi2org = {}
    if config.gidatabase:
        tmp_file = open(config.gidatabase)
        tmp_string = tmp_file.read()
        tmp_file.close()
        local_gi2org = eval(tmp_string)
    else:local_gi2org = {}

    #for the random match estimation
    cluster.blastHits = {}
    
    # parse through the blast records and sort them by amplicon, allele and readnumber/part
    for blast_record in records:

        # get the information from the query header
        amptype = blast_record.query.split('|tempSep|')[0]
        allele =blast_record.query.split('|tempSep|')[1]
        readnumber =blast_record.query.split('|tempSep|')[2]
        
        amplicons[amptype][allele][readnumber] = blast_record
        if blast_record == None: output += amptype+' '+allele+' '+readnumber+'ajajajaj\n'
    
    # parse through the blast result one amplicon at the time
    for amplicon in amplicons:

        broken = False
        cluster.blastHits[amplicon] = {}
        consensuses = amplicons[amplicon]
        output += '\tAmplicon: '+amplicon     +'\n'

        # for each consensus get the blasst records for the two sequences and check for completeness of file
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
                    if perc_identity >= config.minBlastIdentity and perc_coverage >= config.minBlastCoverage:
                        in_r1[organism] = alignment
                        try:
                            hitInfo[organism]['r1']['pi'].append(perc_identity)
                            hitInfo[organism]['r1']['pc'].append(perc_coverage)
                            hitInfo[organism]['r1']['ss'].append(str(hsp.sbjct_start)+'-'+str(hsp.sbjct_start+hsp.align_length))
                        except KeyError:
                            hitInfo[organism] = {'r1':{'pi':[perc_identity],'pc':[perc_coverage],'ss':[str(hsp.sbjct_start)+'-'+str(hsp.sbjct_start+hsp.align_length)]},'r2':{'pi':[],'pc':[],'ss':[]}}

            # get all organisms that part two maps to
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
                    if perc_identity >= config.minBlastIdentity and perc_coverage >= config.minBlastCoverage:
                        in_r2[organism] = alignment
                        try:
                            hitInfo[organism]['r2']['pi'].append(perc_identity)
                            hitInfo[organism]['r2']['pc'].append(perc_coverage)
                            hitInfo[organism]['r2']['ss'].append(str(hsp.sbjct_start)+'-'+str(hsp.sbjct_start+hsp.align_length))
                        except KeyError:
                            hitInfo[organism] = {'r2':{'pi':[perc_identity],'pc':[perc_coverage],'ss':[str(hsp.sbjct_start)+'-'+str(hsp.sbjct_start+hsp.align_length)]},'r1':{'pi':[],'pc':[],'ss':[]}}

            # get all organisms that both parts map to and print them also ptint if no overlap is found
            in_both_reads = []
            for organism in in_r1:
                if organism in in_r2:
                    in_both_reads.append(organism)
            for organism in in_both_reads:
                output +=       '\t\t\t'+organism    +'\n'
                output += '\t\t\t\tpart1:\tHIT#'+'\n\t\t\t\t\tHIT#'.join([str(i)+': identity='+str(hitInfo[organism]['r1']['pi'][i])+'%, coverage='+str(hitInfo[organism]['r1']['pc'][i])+'%, pos='+str(hitInfo[organism]['r1']['ss'][i]) for i in range(len(hitInfo[organism]['r1']['ss']))]) + '\n'
                output += '\t\t\t\tpart2:\tHIT#'+'\n\t\t\t\t\tHIT#'.join([str(i)+': identity='+str(hitInfo[organism]['r2']['pi'][i])+'%, coverage='+str(hitInfo[organism]['r2']['pc'][i])+'%, pos='+str(hitInfo[organism]['r2']['ss'][i]) for i in range(len(hitInfo[organism]['r2']['ss']))]) + '\n'
                cluster.blastHits[amplicon][organism] = hitInfo[organism]
            if not in_both_reads:
                output +=       '\t\t\tNo alignment supported by both reads with >='+str(config.minBlastIdentity)+'% identity and '+str(config.minBlastCoverage)+'% alignment length coverage'     +'\n'
    
    return [output, cluster]

def classifymeta(indata):

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    import os
    if os.path.exists(config.outfile): os.remove(config.outfile)
    config.openconnections()
    writelogheader(config.logfile)

    # settings
    config.load()
    
    if indata.tempFileFolder:
        config.logfile.write('Copying database to temporary location for fast access:\n ');
        import os, shutil, glob
        pid = ''#+os.getpid()
        dbBaseName = config.blastDb.split('/')[-1]
        for filename in glob.glob(config.blastDb+'*'):
                noPathFileName = filename.split('/')[-1]
                src = filename
                dst = indata.tempFileFolder+'/SEAseqtemp/'+str(pid)+''+noPathFileName
                if not os.path.exists(dst):
                    config.logfile.write('\tcopying: '+src+' to '+dst+ '\n');
                    shutil.copyfile(src, dst)
                    config.logfile.write('\tDone.\n');
        config.blastDb = indata.tempFileFolder+'/SEAseqtemp/'+str(pid)+''+dbBaseName
        if config.gidatabase:
            filename = config.gidatabase
            noPathFileName = filename.split('/')[-1]
            src = filename
            dst = indata.tempFileFolder+'/SEAseqtemp/'+str(pid)+''+noPathFileName
            if not os.path.exists(dst):
                config.logfile.write('\tcopying: '+src+' to '+dst+ '\n');
                shutil.copyfile(src, dst)
                config.logfile.write('\tDone.\n');
            config.gidatabase = indata.tempFileFolder+'/SEAseqtemp/'+str(pid)+''+noPathFileName

    config.logfile.write('Starting to align clusters:\n ');
    if indata.debug: #single process // serial
	config.logfile.write('debugging:\n ');
        import sys
	sys.stdout.write('debugging:\n ')
	config.logfile.write('Running in debug mode ')
	results=[] # create holder for processed reads
	#progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile) # creates a progress "bar" thingy
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
    moreThanOneAmpAndMono4All = 0
    moreThanOneAmpMono4AllAndAllHaveHits = 0
    orgInAllAmpsCounter = 0
    noMatchAmp = {}
    matches = {}
    from SEAseqLib.mainLibrary import Progress
    progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile, unit='cluster',mem=True, printint = 1)
    if indata.stop: progress = Progress(indata.stop, logfile=config.logfile, unit='cluster',mem=True, printint = 1)
    with progress:
        for tmp in results:
            tmcounter +=1
            if tmcounter <= indata.skip:
                progress.update(); continue
            [output, cluster] = tmp
            config.outfile.write( output+'\n')
            
            if cluster.definedampliconcount > 1:
                
                # get combo
                ampliconnames = cluster.definedamplicons.keys()
                ampliconnames.sort()
                ampliconcombo = '/'.join(ampliconnames)
                # get monocombo
                monoAmps = {}
                for ampname, amplicon in cluster.definedamplicons.iteritems():
                    if amplicon.monoclonal: monoAmps[amplicon.type] = amplicon.monoclonal
                mononames = monoAmps.keys()
                mononames.sort()
                monoCombo = '/'.join(mononames)
                
                atLeastOneOrgInAllAmps = False
                if ampliconcombo == monoCombo: # monoclonal for all defined amplicons

                    moreThanOneAmpAndMono4All += 1

                    allHaveHits = True
                    noneHaveHits = True
                    for amplicon in ampliconnames:
                        if not cluster.blastHits[amplicon]:
                            allHaveHits = False
                            try:            noMatchAmp[amplicon].append(cluster.id)
                            except KeyError:noMatchAmp[amplicon] = [cluster.id]
                        else: noneHaveHits = False
                    if noneHaveHits:
                        try:            noMatchAmp['any defined amplicon'].append(cluster.id)
                        except KeyError:noMatchAmp['any defined amplicon'] = [cluster.id]


                    for organism in cluster.blastHits[ampliconnames[0]]:
                        notInAll = None
                        for amplicon in ampliconnames[1:]:
                            try:
                                if organism not in cluster.blastHits[amplicon]:
                                    notInAll = True
                            except KeyError: pass
                        if not notInAll: # ie organism inAll amplicons
                            atLeastOneOrgInAllAmps = True
                    
                    if allHaveHits:
                        moreThanOneAmpMono4AllAndAllHaveHits += 1
                        matches[cluster.id] = {}
                        for amplicon in ampliconnames:
                            matches[cluster.id][amplicon] = [organism for organism in cluster.blastHits[amplicon]]
                
                if atLeastOneOrgInAllAmps: orgInAllAmpsCounter += 1
    
            progress.update()
            if indata.stop and tmcounter >= indata.stop: break

    config.outfile.write('##### SUMMARY #####'+'\n')
    #assert len(data) - noblasthit_its - noblasthit_16s - blasthitsagree - hitsdonotagree == 0, '\n\nError '+str(len(data))+' - '+str(noblasthit_its)+' - '+str(noblasthit_16s)+' - '+str(blasthitsagree)+' - '+str(hitsdonotagree)+' != 0\n\n'
    config.outfile.write(      'out of '+str(moreThanOneAmpAndMono4All)+' analyzed clusters with more than one amplicon defined and monoclonal for all the defined amplicon, were:'     +'\n')
    config.outfile.write(      '\t'+str(moreThanOneAmpMono4AllAndAllHaveHits) +' ('+str(round(100*float(moreThanOneAmpMono4AllAndAllHaveHits)/float(moreThanOneAmpAndMono4All),2))+'%) clusters that had at least one blast hit for all defined amplicons, out of these were:\n')
    config.outfile.write(      '\t\t'+str(orgInAllAmpsCounter) +' ('+str(round(100*float(orgInAllAmpsCounter)/float(moreThanOneAmpMono4AllAndAllHaveHits),2))+'%) clusters where atleast one organism were identified within all the hitLists of all defined amplicons.\n')
    for amplicon in noMatchAmp:  
        config.outfile.write(  '\t'+str(len(noMatchAmp[amplicon])) +' ('+str(round(100*float(len(noMatchAmp[amplicon]))/float(moreThanOneAmpAndMono4All),2))+'%) clusters had no BLAST hits for '+amplicon+' supported by both reads with >='+str(config.minBlastIdentity)+'% identity and '+str(config.minBlastCoverage)+'% alignment length coverage.'     +'\n')
    #config.outfile.write(      str(hitsdonotagree+blasthitsagree) +' ('+str(round(100*float(hitsdonotagree+blasthitsagree)/float(len(data)),2))+'%) have hits for both ITS and 16S.'     +'\n')
    #if hitsdonotagree+blasthitsagree > 0:
    #    config.outfile.write(  str(hitsdonotagree) +' ('+str(round(100*float(hitsdonotagree)/float(len(data)),2))+'%) removed because none of the ITS and 16S hits did agree.'+         ' ('+str(round(100*float(hitsdonotagree)/float(hitsdonotagree+blasthitsagree),2))+'% out of clusters with hits for both)'     +'\n')
    #    config.outfile.write(  'For '+str(blasthitsagree) +' ('+str(round(100*float(blasthitsagree)/float(len(data)),2))+'%) clusters did the ITS and 16S hits agree at least once.'+' ('+str(round(100*float(blasthitsagree)/float(hitsdonotagree+blasthitsagree),2))+'% out of clusters with hits for both)'     +'\n')
    #
    config.logfile.write('All clusters aligned, starting random match estimation ...\n')
    if matches:
        groups = {}
        for cluster, amplicons in matches.iteritems():
            for amplicon, hitList in amplicons.iteritems():
                try:
                    groups[amplicon].append(hitList)
                except KeyError:
                    groups[amplicon]=[hitList]
        
        import random
        total_tries = 10000
        random_hits = 0
        for i in xrange(total_tries):
            currentHitList = {}
            for amplicon, listOfHitLists in groups.iteritems():
                currentHitList[amplicon] = random.choice(listOfHitLists)
            
            atLeastOneOrgInAllHitLists = False
            for organismHit in currentHitList[groups.keys()[0]]:
                notInAllHitLists = None
                for amplicon in groups.keys()[1:]:
                    if organismHit not in currentHitList[amplicon]:
                        notInAllHitLists = True
                if not notInAllHitLists: # ie organism in All hitLists
                    atLeastOneOrgInAllHitLists = True
            if atLeastOneOrgInAllHitLists: random_hits += 1

        

        config.outfile.write(
            '\nThere was '+str(random_hits)+' ('+str(round(100*float(random_hits)/float(total_tries),2))+'%) occasion when at least one organism was present in all BLAST-hitList, when trying '+
            str(total_tries)+' times to find overlaps between amplicon '+', '.join([amplicon for amplicon in groups.keys()[:-1]])+' and '+groups.keys()[-1]+' hitLists each one randomly chosen from a population of BarcodeClusters.\n'
        )
        config.outfile.write('Rand match was based on:\n')
        for amplicon, listOfHitLists in groups.iteritems():
            config.outfile.write( str(len(listOfHitLists)) +' Hitlists for amplicon '+amplicon+'.\n')
        
        mostCommonOrganism = {}
        for amplicon, listOfHitLists in groups.iteritems():
            mostCommonOrganism[amplicon] = {'numberOfHitlists':0,'organisms':{}}
            for hitList in listOfHitLists:
                mostCommonOrganism[amplicon]['numberOfHitlists']+=1
                for organism in hitList:
                    try: mostCommonOrganism[amplicon]['organisms'][organism]+=1
                    except KeyError:mostCommonOrganism[amplicon]['organisms'][organism]=1
        
        mostCommonCountToShow = config.mostCommonToShow
        import operator
        for amplicon, info in mostCommonOrganism.iteritems():
            config.outfile.write( 'Amplicon '+amplicon+' had '+str(mostCommonOrganism[amplicon]['numberOfHitlists'])+' hitLists, the '+str(mostCommonCountToShow)+' most common organisms were:\n')
            if mostCommonCountToShow > len(mostCommonOrganism[amplicon]['organisms']):
                tmpList = sorted(mostCommonOrganism[amplicon]['organisms'].iteritems(), key=operator.itemgetter(1))[::-1]
            else:
                tmpList = sorted(mostCommonOrganism[amplicon]['organisms'].iteritems(), key=operator.itemgetter(1))[::-1][:mostCommonCountToShow]
            for organism, count in tmpList:
                config.outfile.write('\t'+str(count) +'\t('+str(round(100*float(count)/float(mostCommonOrganism[amplicon]['numberOfHitlists']),2))+'%) of the hitlists had hits towards\t'+organism+ '.\n')

        
    else:
        config.outfile.write('No clusters found with hits for both were found.\n')
        
    config.logfile.write('Classification done.\n')

    return 0