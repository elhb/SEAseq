def classify_main():
    pass

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
                    blastclassmatches[amplicon] = [classification for classification in cluster.classifications[amplicon]]
            for amplicon in ampliconnames:
                rdpclassmatches[amplicon] = [rdpclass for rdpclass in cluster.rdpAmplicons[amplicon]]
        
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
    for amplicon in cluster.amplicons.values(): perAmpOut += amplicon.checkmono(config)
    output += 'There are '+str(cluster.adaptercount)+' illumina adapter reads.\n'
    output += 'There are '+str(cluster.primererrors)+' primer missmatch reads.\n'
    output += 'There are '+str(cluster.readcount-cluster.primererrors-cluster.adaptercount)+' "good reads".\n'
    if cluster.ampliconpairs == 0:
        output += '0 amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% support and >= '+str(config.minReadCountPerConsensus)+' reads)\n'
    if cluster.ampliconpairs > 0:
        output += str(cluster.definedampliconcount)+' amplicon(s) have enough data (>=1 cons with >= '+str(config.minReadPopSupportConsensus)+'% support and >= '+str(config.minReadCountPerConsensus)+' reads):\n'
    output += perAmpOut + '\n'
    
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