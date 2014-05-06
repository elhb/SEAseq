def classify_main():
    pass

def doBLASTrandMatch(config,blastclassmatches):
    randmatchrounds = open(config.path+'/blastRandMatchRounds.txt','w')
    if not blastclassmatches: randmatchrounds.write('# No matches found exiting.\n');return
    from SEAseqLib.mainLibrary import Progress, popRandom, matchBLASTclass
    output = ''
    blastclassmatchesBackup = blastclassmatches.copy()
#    print blastclassmatchesBackup
    blastRandomCounter = {}
    total_fakeClusters = 0
    repeats = 100
    randmatchrounds.write('# HERE WE GO #\n')
    progress = Progress(repeats, logfile=config.logfile, unit='round',mem=True, printint = 5)
    with progress:
        for i in range(repeats):
            progress.update()
            blastRandomCounterRound = {}
            randmatchrounds.write('\n##############################\nRandom matching repeat '+str(i+1)+'\n')
            blastclassmatches = blastclassmatchesBackup.copy()
            emptyList = False
            tmcounter2 = 0
            while not emptyList:
                tmcounter2 += 1
                total_fakeClusters += 1 
                blastAmplicons = {}
    #            print blastclassmatches.keys()
                for amplicon in blastclassmatches:
                    poppedItem, blastclassmatches[amplicon] = popRandom(blastclassmatches[amplicon])
                    if blastclassmatches[amplicon] == [] :
                        randmatchrounds.write('Amplicon '+amplicon+' list is empty after '+str(tmcounter2)+' pops.\n')
                        emptyList = True
                    blastAmplicons[amplicon] = poppedItem
    
    #            print 'blastAmplicons =',blastAmplicons
                classificationsInAllAmps = matchBLASTclass(blastAmplicons)
    
                if classificationsInAllAmps:
                    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
                        if rank in classificationsInAllAmps and classificationsInAllAmps[rank]:
    
                                if len(classificationsInAllAmps[rank]) > 1:
                                        overlaps = ', '.join(classificationsInAllAmps[rank][:-1])+' or '+classificationsInAllAmps[rank][-1]
                                else:   overlaps = classificationsInAllAmps[rank][0]
    
                        else:
                            classificationsInAllAmps[rank] = ['NoMatchFound']
                            overlaps = 'NoMatchFound'
    
                        try: blastRandomCounter[rank][overlaps] +=1
                        except KeyError:
                            if rank not in blastRandomCounter: blastRandomCounter[rank] = {overlaps:1}
                            else:                            blastRandomCounter[rank][overlaps] = 1
    
                        try: blastRandomCounterRound[rank][overlaps] +=1
                        except KeyError:
                            if rank not in blastRandomCounterRound: blastRandomCounterRound[rank] = {overlaps:1}
                            else:                            blastRandomCounterRound[rank][overlaps] = 1
    
    #                output += 'taxdata\t'+str(classificationsInAllAmps)+'\n'
                else:
                    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
                        overlaps = 'No classification overlap found for any rank.'
                        try: blastRandomCounter[rank][overlaps] +=1
                        except KeyError:
                            if rank not in blastRandomCounter: blastRandomCounter[rank] = {overlaps:1}
                            else:                            blastRandomCounter[rank][overlaps] = 1
    
                        try: blastRandomCounterRound[rank][overlaps] +=1
                        except KeyError:
                            if rank not in blastRandomCounterRound: blastRandomCounterRound[rank] = {overlaps:1}
                            else:                            blastRandomCounterRound[rank][overlaps] = 1                    
    #                output += 'No classification overlap found.\n'
    #                output += 'taxdata\t'+str(classificationsInAllAmps)+'\n'
    
            randmatchrounds.write('\nRandom BLAST matching repeat '+str(i+1)+' gave '+str(tmcounter2)+' "fakeclusters" where the classifications were matched:\n')
            blastrandomprint(randmatchrounds,blastRandomCounterRound)
    randmatchrounds.close()
    config.outfile.write('\nRandom BLAST matching ('+str(repeats)+' repeats) gave '+str(total_fakeClusters)+' "fakeclusters" where the classifications were matched:\n')
    blastrandomprint(config.outfile,blastRandomCounter)

def blastrandomprint(fileHandle,dictionary,outType='New'):
    if outType == 'Old':
        for rank, names in dictionary.iteritems():
            fileHandle.write('\tfor rank '+rank+':\n')
            for name, count in names.iteritems():
                fileHandle.write('\t\t'+str(count)+' ('+str(round(100*float(count)/float(total_fakeClusters),2))+'%) were '+name+'\n')
    elif outType == 'New':
        for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
            names = dictionary[rank]
            fileHandle.write('\n\t---------- '+rank)
            
            total = sum(names.values())
            try:count = total-names['NoInformationAvailable']
            except KeyError: count = total
            
            try: count2 = names['NoInformationAvailable']
            except KeyError:count2 = 0
            try:percentage = round(100*float(count2)/float(total),2)
            except ZeroDivisionError:percentage=0.00
            #fileHandle.write('\n\t'+str(count2)+' are excluded due to '+'NoInformationAvailable'+' ('+str(total-count2)+' clusters have info at this rank).')
            fileHandle.write('\n\t'+str(count2)+' clusters that are excluded as NoInformationAvailable for any BlastHit at this rank ('+str(total-count2)+' clusters have info at this rank).')
            
            try: count2 = names['NoMatchFound']
            except KeyError: count2 = 0
            try:percentage = round(100*float(count2)/float(count),2)
            except ZeroDivisionError:percentage=0.00
            #fileHandle.write('\n\t'+str(count2)+' ('+str(percentage)+'%) are MissMatch.')
            fileHandle.write('\n\t'+str(count2)+' ('+str(percentage)+'%) clusters where no overlap is found between the amplicon hitlists for this rank.')
            
            try: matchCount = count-names['NoMatchFound']
            except KeyError: matchCount = count
            if count: percentage = round(100*float(matchCount)/float(count),2)
            else: percentage = 0.00
            #fileHandle.write('\n\t'+str(matchCount)+' clusters ('+str(percentage)+'%) where the amplicons RDPclassinfo has >=1 match at the '+rank+' level, out of theese are:')
            fileHandle.write('\n\t'+str(matchCount)+' clusters ('+str(percentage)+'%) where the amplicons hitlist has >=1 match at the '+rank+' level, out of theese are:')
            
            if count:
                for name, count2 in names.iteritems():
                    if name == 'NoMatchFound' or name == 'NoInformationAvailable': continue
                    percentage = round(100*float(count2)/float(matchCount),2)                        
                    fileHandle.write('\n\t\t'+str(count2)+' ('+str(percentage)+'%) were '+name+'')

def rdprandomprint(output2,dictionary):
    assert type(output2) == file or type(output2) == str, 'Error in function "rdprandomprint()": the output2 variable for the function must be a string or file handle.\n'
    output = ''
    for rank in ['domain','phylum','class','order','family','genus']:

        try: names = dictionary[rank]
        except KeyError: output += 'Warning: '+rank+' is not present in dictionary.\n'; continue
        output += '\n\t---------- '+rank
        
        total = sum(names.values())
        try:count = total-names['LowConfidence']
        except KeyError: count = total
        
        try: count2 = names['LowConfidence']
        except KeyError:count2 = 0
        try:percentage = round(100*float(count2)/float(total),2)
        except ZeroDivisionError:percentage=0.00
        output += '\n\t'+str(count2)+' are excluded due to LowConfidence ('+str(total-count2)+' clusters have info at this rank).'
        
        try: count2 = names['MissMatch']
        except KeyError: count2 = 0
        try:percentage = round(100*float(count2)/float(count),2)
        except ZeroDivisionError:percentage=0.00
        output += '\n\t'+str(count2)+' ('+str(percentage)+'%) are MissMatch.'
        
        try: matchCount = count-names['MissMatch']
        except KeyError: matchCount = count
        if count: percentage = round(100*float(matchCount)/float(count),2)
        else: percentage = 0.00
        output += '\n\t'+str(matchCount)+' clusters ('+str(percentage)+'%) where the amplicons RDPclassinfo has >=1 match at the '+rank+' level, out of theese are:'
        
        if count:
            for name, count2 in names.iteritems():
                if name == 'MissMatch' or name == 'LowConfidence': continue
                percentage = round(100*float(count2)/float(matchCount),2)                        
                output += '\n\t\t'+str(count2)+' ('+str(percentage)+'%) were '+name+''
    
    if   type(output2) == file: output2.write(output)
    elif type(output2) == str:  return output

def rdpMatchRound(rdpclassmatches):
    from SEAseqLib.mainLibrary import popRandom, matchRDPclass
    rdpRandomCounterRound = {}
    output = ''
    emptyList = False
    tmcounter2 = 0
    while not emptyList:
        tmcounter2 += 1
        rdpAmplicons = {}
        #print 'classmaatchkeys->',rdpclassmatches.keys()
        for amplicon in rdpclassmatches:
            #print amplicon
            prepoplen = len(rdpclassmatches[amplicon])
            if len(rdpclassmatches[amplicon]) == 0:
                output += 'Amplicon '+amplicon+' list is empty after '+str(tmcounter2)+' pops.\n'
                emptyList = True
                continue
            poppedItem, rdpclassmatches[amplicon] = popRandom(rdpclassmatches[amplicon])
            assert prepoplen == len(rdpclassmatches[amplicon])+1, 'poping Error: '+str(prepoplen)+' != '+str(len(rdpclassmatches[amplicon])+1)
            #print '## HERE '+amplicon+'-> ',len(rdpclassmatches[amplicon])
            if rdpclassmatches[amplicon] == [] or len(rdpclassmatches[amplicon]) == 0:
                output += 'Amplicon '+amplicon+' list is empty after '+str(tmcounter2)+' pops.\n'
                emptyList = True
            rdpAmplicons[amplicon] = {0:poppedItem}
        #print emptyList
        classificationsInAllAmps = matchRDPclass(rdpAmplicons)
        #print 'fakecluster',tmcounter2,
        if classificationsInAllAmps:
            for rank in ['domain','phylum','class','order','family','genus']:
                if rank in classificationsInAllAmps and classificationsInAllAmps[rank]:
                    #print rank, classificationsInAllAmps[rank]
                    try: rdpRandomCounterRound[rank][classificationsInAllAmps[rank][0]] +=1
                    except KeyError:
                        if rank not in rdpRandomCounterRound: rdpRandomCounterRound[rank] = {classificationsInAllAmps[rank][0]:1}
                        else:                            rdpRandomCounterRound[rank][classificationsInAllAmps[rank][0]] = 1
                else: print rank, 'not in comparison'
        else:
            pass#print 'No classification overlap found.'
        if tmcounter2 > 10000000: output += '\nWARNING: Now at 10M pops and list is still not empty, BREAKING loop.\n';break
    output += rdprandomprint('A string:\n',rdpRandomCounterRound)
    return (rdpRandomCounterRound, output, tmcounter2)

def yieldThis(timesToYield,whatToYield):
    if type(whatToYield) == dict: whatToYieldBackup = whatToYield.copy()
    for i in  xrange(timesToYield):
        if type(whatToYield) == dict: whatToYield = whatToYieldBackup.copy()
        yield whatToYield

def doRDPrandMatch(config,rdpclassmatches):
    
    #
    # Initiialize vars and import functionality
    #
    from SEAseqLib.mainLibrary import Progress, popRandom, matchRDPclass
    randmatchrounds = open(config.path+'/rdpRandMatchRounds.txt','w')
    if not rdpclassmatches: randmatchrounds.write('# No matches found exiting.\n');return
    rdpRandomCounter = {}
    total_fakeClusters = 0
    repeats = 100
    
    #
    # do rounds of random matching
    #
    randmatchrounds.write('# HERE WE GO #\n')
    debug = False
    
    if debug: #single process // serial
        config.logfile.write('debugging:\n ');
        import sys
        config.logfile.write('Running in debug mode ...\n')
        results=[] # create holder for processed reads
        progress = Progress(repeats, logfile=config.logfile, unit='round',mem=True, printint = 10)
        i = 0
        with progress:
            for rdpclassmatchesWorkingCopy in yieldThis(repeats,rdpclassmatches):
                i += 1
                # update progress and do round
                progress.update()
                (rdpRandomCounterRound, output, tmcounter2) = rdpMatchRound(rdpclassmatchesWorkingCopy)
                randmatchrounds.write('\n\n##############################\nRandom matching repeat '+str(i)+'\n')
                randmatchrounds.write('\nRandom RDP matching repeat '+str(i)+' gave '+str(tmcounter2)+' "fakeclusters" where the classifications were matched:\n')
                randmatchrounds.write(output)
    
                # add round to total counters
                total_fakeClusters += tmcounter2
                for rank, names in rdpRandomCounterRound.iteritems():
                    for name, count in names.iteritems():
                        try: rdpRandomCounter[rank][name] += count
                        except KeyError:
                            if rank not in rdpRandomCounter: rdpRandomCounter[rank] = {name:count}
                            else:                            rdpRandomCounter[rank][name] = count
        config.logfile.write('finished, making summary ... \n')

    else: # multiple processes in parallel
        import multiprocessing
        WorkerPool = multiprocessing.Pool(multiprocessing.cpu_count(),maxtasksperchild=10000)
        results = WorkerPool.imap_unordered(rdpMatchRound,yieldThis(repeats,rdpclassmatches),chunksize=1)
        #results = WorkerPool.imap(         rdpMatchRound,yieldThis(repeats,rdpclassmatches),chunksize=1)

        progress = Progress(repeats, logfile=config.logfile, unit='round',mem=True, printint = 10)
        i = 0
        with progress:
            for (rdpRandomCounterRound, output, tmcounter2) in results:
                i += 1
        
                # update progress and do round
                progress.update()
                #rdpRandomCounterRound, output, tmcounter2 = rdpMatchRound(rdpclassmatches)
                randmatchrounds.write('\n\n##############################\nRandom matching repeat '+str(i)+'\n')
                randmatchrounds.write('\nRandom RDP matching repeat '+str(i)+' gave '+str(tmcounter2)+' "fakeclusters" where the classifications were matched:\n')
                randmatchrounds.write(output)
    
                # add round to total counters
                total_fakeClusters += tmcounter2
                for rank, names in rdpRandomCounterRound.iteritems():
                    for name, count in names.iteritems():
                        try: rdpRandomCounter[rank][name] += count
                        except KeyError:
                            if rank not in rdpRandomCounter: rdpRandomCounter[rank] = {name:count}
                            else:                            rdpRandomCounter[rank][name] = count
        WorkerPool.close()
        WorkerPool.join()
    randmatchrounds.close()

    # print total output
    config.outfile.write('\nRandom RDP matching ('+str(repeats)+' repeats) gave '+str(total_fakeClusters)+' "fakeclusters" where the classifications were matched:\n')
    rdprandomprint(config.outfile, rdpRandomCounter)

def resultsHandling(tmp,config,clusterdump,counter,moreThanOneAmpAndMono4All,moreThanOneAmpMono4AllAndAllHaveHits,orgInAllAmpsCounter,noMatchAmp,matches,singleAmpliconHitlists,blastclassmatches,rdpclassmatches):
    [output, cluster,picklestring] = tmp
    config.outfile.write( output+'\n')
    clusterdump.write(picklestring)
    counter.addcluster(cluster)
    
    if cluster.definedampliconcount == 1:
        amplicon = cluster.definedamplicons.keys()[0]
        if cluster.blastHits != 'No fasta produced.' and cluster.blastHits[amplicon]:
            singleAmpliconHitlists[cluster.id] = {}
            singleAmpliconHitlists[cluster.id][amplicon] = [organism for organism in cluster.blastHits[amplicon]]
        
    elif cluster.definedampliconcount > 1 and cluster.blastHits != 'No fasta produced.':
        
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
                    #for classification in cluster.classifications[amplicon]:
                    try: blastclassmatches[amplicon].append(cluster.classifications[amplicon])
                    except KeyError: blastclassmatches[amplicon] = [cluster.classifications[amplicon]]
            for amplicon in ampliconnames:
                for consensus, rdpclass in cluster.rdpAmplicons[amplicon].iteritems():
                    try: rdpclassmatches[amplicon].append(rdpclass)
                    except KeyError: rdpclassmatches[amplicon] = [rdpclass]
        
        if atLeastOneOrgInAllAmps: orgInAllAmpsCounter += 1
    return [cluster,config,counter,moreThanOneAmpAndMono4All,moreThanOneAmpMono4AllAndAllHaveHits,orgInAllAmpsCounter,noMatchAmp,matches,singleAmpliconHitlists,blastclassmatches,rdpclassmatches]

def foreachCluster(tmp):
    
    import cPickle
    
    # unpack input information
    cluster,config,indata = tmp

    # initiate the output string and set some variables
    output = '######## Cluster '+str(cluster.id)+' ###################\n'
    output += 'Cluster contains '+str(cluster.readcount)+' raw reads.\n'
    if config == 'skipped':
        output += '# SKIPPING CLUSTER!\n\n'
        cluster.blastHits = 'No fasta produced.'
        return [output, cluster, cPickle.dumps(cluster,-1)]

    # get the cluster and per amplicon information and append it to output
    perAmpOut = ''
    for amplicon in cluster.amplicons.values(): perAmpOut += amplicon.checkmono(config,cluster.readcount)
    output += 'There are '+str(cluster.adaptercount)+' illumina adapter reads.\n'
    output += 'There are '+str(cluster.primererrors)+' primer missmatch reads.\n'
    output += 'There are '+str(cluster.readcount-cluster.primererrors-cluster.adaptercount)+' "good reads".\n'
    if cluster.ampliconpairs == 0:
        output += '0 amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% support and >= '+str(config.minReadCountPerConsensus)+' reads)\n'
    if cluster.ampliconpairs > 0:
        output += str(cluster.definedampliconcount)+' amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% support and >= '+str(config.minReadCountPerConsensus)+' reads):\n'
    output += perAmpOut + '\n'
    
    if indata.reFilter:
        output += 'Skipping all classsification as reFilter=True.\n'
        cluster.blastHits = 'No fasta produced.'
        return [output, cluster, cPickle.dumps(cluster,-1)]

######### DO NOT BLAST AMPLICONS THAT ARE NOT MONOCLONAL FOR ALL DEFINED OR HAVE <= 1 DEFINED AMPLICON
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
        
        if ampliconcombo != monoCombo: # not monoclonal for all defined amplicons
            output += 'Not monoclonal for all defined amplicons, skipping BLAST and RDP classification.\n'
            cluster.blastHits = 'No fasta produced.'
            return [output, cluster, cPickle.dumps(cluster,-1)]
    else:
        output += 'Less than two defined amplicons, skipping BLAST and RDP classification.\n'
        cluster.blastHits = 'No fasta produced.'
        return [output, cluster, cPickle.dumps(cluster,-1)]
#########

    # make fasta file with consensus sequences to classify
    tmpOut, errorCode = cluster.createConsensusFasta(config,indata)
    output += tmpOut
    if errorCode:
        return [output, cluster, cPickle.dumps(cluster,-1)]

    # Align consensus sequences by blast
    output += cluster.blastAllAmplicons(config,indata)

    output += cluster.parseBlastAmplicons(config,indata)
    output += cluster.parseRdpAmplicons(config,indata)

    return [output, cluster, cPickle.dumps(cluster,-1)]

def classifymeta(indata):

    #import multiprocessing, logging
    #logger = multiprocessing.log_to_stderr()
    #logger.setLevel(multiprocessing.SUBDEBUG)

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    import os
    if os.path.exists(config.outfile): os.remove(config.outfile)
    config.openconnections()
    writelogheader(config.logfile)

    # settings
    config.load()
    
    if indata.tempFileFolder and not indata.debug and indata.tempFileFolder != config.path and config.blastDb != '/sw/data/uppnex/blast_databases/nt':
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
        if indata.tempFileFolder:
            database='/proj/b2011011/SEAseq/reference/NCBItaxonomy.db'
            filename = database
            noPathFileName = filename.split('/')[-1]
            src = filename
            dst = indata.tempFileFolder+'/SEAseqtemp/'+noPathFileName
            # database=indata.tempFileFolder+'/SEAseqtemp/NCBItaxonomy.db'
            if not os.path.exists(dst):
                config.logfile.write('\tcopying: '+src+' to '+dst+ '\n');
                shutil.copyfile(src, dst)
                config.logfile.write('\tDone.\n');

    config.loadPrimers()

    from multiprocessing import Manager
    man = Manager()
    config.dbLock = man.Lock()

    tmcounter = 0
    counter = RunStatCounter(config)
    moreThanOneAmpAndMono4All = 0
    moreThanOneAmpMono4AllAndAllHaveHits = 0
    orgInAllAmpsCounter = 0
    noMatchAmp = {}
    matches = {}
    blastclassmatches = {}
    rdpclassmatches = {}
    singleAmpliconHitlists = {}

    if indata.tempFileFolder:
        clusterdump = open(indata.tempFileFolder+'/SEAseqtemp/clusters.pickle','wb')
    else:
        import os
        if os.path.exists(config.path+'/clusters.pickle') and os.path.islink(config.path+'/clusters.pickle'): os.unlink(config.path+'/clusters.pickle')
        clusterdump = open(config.path+'/clusters.pickle','wb')

    from SEAseqLib.mainLibrary import Progress, clusterGenerator
    config.logfile.write('Starting to align clusters:\n ');
    if indata.debug: #single process // serial
        config.logfile.write('debugging:\n ');
        import sys
        sys.stdout.write('debugging:\n ')
        config.logfile.write('Running in debug mode ...\n')
        results=[] # create holder for processed reads
        progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile, unit='cluster',mem=True, printint = 1)
        if indata.stop: progress = Progress(indata.stop, logfile=config.logfile, unit='cluster',mem=True, printint = 1)
        with progress:
            tmcounter = 0
            for cluster in clusterGenerator(config, indata):
                tmcounter +=1
                progress.update()
                if tmcounter <= indata.skip: continue
                #results.append(foreachCluster(cluster))
                tmp = foreachCluster(cluster)
                [cluster,config,counter,moreThanOneAmpAndMono4All,moreThanOneAmpMono4AllAndAllHaveHits,orgInAllAmpsCounter,noMatchAmp,matches,singleAmpliconHitlists,blastclassmatches,rdpclassmatches] = resultsHandling(tmp,config,clusterdump,counter,moreThanOneAmpAndMono4All,moreThanOneAmpMono4AllAndAllHaveHits,orgInAllAmpsCounter,noMatchAmp,matches,singleAmpliconHitlists,blastclassmatches,rdpclassmatches)
                if indata.stop and tmcounter >= indata.stop: break
        config.logfile.write('finished, making summary ... \n')
    else: # multiple processes in parallel
        import multiprocessing
        WorkerPool = multiprocessing.Pool(indata.cpus,maxtasksperchild=10000)
        results = WorkerPool.imap_unordered(foreachCluster,clusterGenerator(config,indata),chunksize=10)
        #results = WorkerPool.imap(           foreachCluster,clusterGenerator(config,indata),chunksize=10)

    # will paralellise this when stuff works
    progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile, unit='cluster',mem=True, printint = 1)
    if indata.stop: progress = Progress(indata.stop, logfile=config.logfile, unit='cluster',mem=True, printint = 1)
    with progress:
        for tmp in results:
            if indata.debug: break
            tmcounter +=1
            if tmcounter <= indata.skip: progress.update(); continue
            [cluster,config,counter,moreThanOneAmpAndMono4All,moreThanOneAmpMono4AllAndAllHaveHits,orgInAllAmpsCounter,noMatchAmp,matches,singleAmpliconHitlists,blastclassmatches,rdpclassmatches] = resultsHandling(tmp,config,clusterdump,counter,moreThanOneAmpAndMono4All,moreThanOneAmpMono4AllAndAllHaveHits,orgInAllAmpsCounter,noMatchAmp,matches,singleAmpliconHitlists,blastclassmatches,rdpclassmatches)
    
            progress.update()
            if indata.stop and tmcounter >= indata.stop: break
    WorkerPool.close()
    WorkerPool.join()

    config.outfile.write(counter.createsummary(config))

    config.outfile.write('##### SUMMARY #####'+'\n')
    config.outfile.write(      'out of '+str(moreThanOneAmpAndMono4All)+' analyzed clusters with more than one amplicon defined and monoclonal for all the defined amplicon, were:'     +'\n')
    
    if moreThanOneAmpAndMono4All:
        config.outfile.write(      '\t'+str(moreThanOneAmpMono4AllAndAllHaveHits) +' ('+str(round(100*float(moreThanOneAmpMono4AllAndAllHaveHits)/float(moreThanOneAmpAndMono4All),2))+'%) clusters that had at least one blast hit for all defined amplicons, out of these were:\n')
    else:
        config.outfile.write(      '\t'+str(moreThanOneAmpMono4AllAndAllHaveHits) +' ('+str(round(100*0,2))+'%) clusters that had at least one blast hit for all defined amplicons, out of these were:\n')
    
    if moreThanOneAmpMono4AllAndAllHaveHits:
        config.outfile.write(      '\t\t'+str(orgInAllAmpsCounter) +' ('+str(round(100*float(orgInAllAmpsCounter)/float(moreThanOneAmpMono4AllAndAllHaveHits),2))+'%) clusters where atleast one organism were identified within all the hitLists of all defined amplicons.\n')
    else:
        config.outfile.write(      '\t\t'+str(orgInAllAmpsCounter) +' ('+str(round(100*0,2))+'%) clusters where atleast one organism were identified within all the hitLists of all defined amplicons.\n')
    
    for amplicon in noMatchAmp:  
        config.outfile.write(  '\t'+str(len(noMatchAmp[amplicon])) +' ('+str(round(100*float(len(noMatchAmp[amplicon]))/float(moreThanOneAmpAndMono4All),2))+'%) clusters had no BLAST hits for '+amplicon+' supported by both reads with >='+str(config.minBlastIdentity)+'% identity and '+str(config.minBlastCoverage)+'% alignment length coverage.'     +'\n')

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
    config.outfile.write('\n')
    if singleAmpliconHitlists:
        config.outfile.write('Checking organism distribution for clusters with single amplicons:\n')
        groups = {}
        for cluster, amplicons in singleAmpliconHitlists.iteritems():
            for amplicon, hitList in amplicons.iteritems():
                try:
                    groups[amplicon].append(hitList)
                except KeyError:
                    groups[amplicon]=[hitList]
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
        config.outfile.write('No clusters found with single amplicons and hits for that amp were found.\n')

    config.logfile.write('Here goes the new random matching ...\n')
    
    config.logfile.write('RDP random matching ...\n')
    doRDPrandMatch(config,rdpclassmatches)
    config.logfile.write('RDP random matching done.\n')
   
    config.logfile.write('BLAST random matching ...\n')
    doBLASTrandMatch(config,blastclassmatches)
    config.logfile.write('BLAST random matching done.\n')
    config.logfile.write('Classification done.\n')

    clusterdump.close()
    #if indata.tempFileFolder:
    #    import shutil
    #    import os
    #    if os.path.exists(config.path+'/clusters.pickle') and os.path.islink(config.path+'/clusters.pickle'): os.unlink(config.path+'/clusters.pickle')
    #    shutil.move(indata.tempFileFolder+'/SEAseqtemp/clusters.pickle',config.path+'/clusters.pickle')

    return 0

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
        
        self.statstable = open(config.path+'/classify.statstable','w',1)
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
    
if __name__ == "__main__": classify_main()