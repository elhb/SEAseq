def compare(indata):
    
    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    import os
    import sys
    if os.path.exists(config.outfile): os.remove(config.outfile)
    config.openconnections()
    writelogheader(config.logfile)

    # settings
    config.load()

    infile = None
    
    if os.path.exists(config.path+'/classify.out.txt'):
        infile = open(config.path+'/classify.out.txt','r',1024*100)
    elif os.path.exists(config.path+'/classify.out.txt.gz'):
        import gzip
        infile = gzip.open(config.path+'/classify.out.txt.gz')
    else:
        config.logfile.write(config.path+'/classify.out.txt(.gz) does not exist please complete the classify step before trying to run SEAseq compare.\n')
        config.outfile.write(config.path+'/classify.out.txt(.gz) does not exist please complete the classify step before trying to run SEAseq compare.\n')
        sys.exit()
    
    config.logfile.write('Reading indata ...\n')
    temp0 = 0
    temp1 = 0
    temp2 = 0
    temp3 = 0
    clusterOrganisms = {}
    overlaps = {}
    classifications = {}
    ranks = {}
    rdpClassifications = {}
    rdpRanks = {}
    rdpRanksWithInfo = {}
    summaryHeader = False
    import re
    for line in infile:
        clusterisLine   = re.search('cluster is',line)
        overlapLine     = re.match("{'",line)
        classificationline = re.match("taxdata",line)#[u'Bacteria', 'NotSet'] ['superkingdom', 'kingdom']
        rdpLine = re.match("RDPtaxdata",line)
        if not summaryHeader: summaryHeader = re.match('##### SUMMARY #####',line)
        else: break
        
        if clusterisLine:
            temp0 +=1 
            try: clusterOrganisms[ line.rstrip() ] += 1
            except KeyError: clusterOrganisms[ line.rstrip() ] = 1
        elif overlapLine:
            temp1 += 1
            try: overlaps[line.rstrip()]+=1
            except KeyError:overlaps[line.rstrip()]=1
        elif classificationline:
            temp2 += 1
            #tmp1 = eval(line.split('\t')[0])
            classificationsInAllAmps = eval(line.split('\t')[1])
            #if not tmp1 and not tmp2: temp2 -= 1
            #for i in range(len(tmp2)):
            for rank, rankValues in classificationsInAllAmps.iteritems():
                if not rankValues: print 'WARNING!!';break
                try:
                    ranks[rank] += 1
                except KeyError:
                    ranks[rank]  = 1
                    classifications[rank] = {"['NoInformationAvailable']":0,"['NoMatchFound']":0}
#                for rankValue in rankValues:
                try:             classifications[rank][str(rankValues)] += 1
                except KeyError: classifications[rank][str(rankValues)]  = 1
        elif rdpLine:
            temp3 += 1
            rdpClassificationsInAllAmps = eval(line.split('\t')[1])
            for rank, rankValues in rdpClassificationsInAllAmps.iteritems():
                if not rankValues: break
                try:
                    rdpRanks[rank] += 1
                except KeyError:
                    rdpRanks[rank]  = 1
                    rdpClassifications[rank] = {"['MissMatch']":0,"['LowConfidence']":0}
                try:             rdpClassifications[rank][str(rankValues)] += 1
                except KeyError: rdpClassifications[rank][str(rankValues)]  = 1

    if not summaryHeader:
        config.outfile.write('# WARNING: The classification outfile is not complete!\n')
        config.logfile.write('# WARNING: The classification outfile is not complete!\n')
    config.logfile.write('Analyzing organism distribution ...\n')    
    config.outfile.write('Organism distribution:\nTotal '+str(temp0)+' clusters, '+str(len(clusterOrganisms))+' organism combinations/classifications.\n')
    config.outfile.write('#\t%\tclassification\n')
    #total = sum([count for defenition, count in clusterOrganisms.iteritems()])
    import operator
    for defenition, count in sorted(clusterOrganisms.iteritems(), key=operator.itemgetter(1))[::-1]:
    #for defenition, count in clusterOrganisms.iteritems():
        config.outfile.write( str(count) + '\t'+str(round(100*float(count/float(temp0)),2)) + '%\t' + defenition.replace('cluster is ','') +'\n')
    
    def getAmplicons(overlaps):
        amplicons = {}
        for data in overlaps:
            data = eval(data)
            for key in data:
                amplicons[key] = True
        del amplicons['inAll']
        return amplicons
    
    def printTable(overlaps,outfile,total):
        amplicons = getAmplicons(overlaps)
        wReductions = []
        outfile.write('\n#clusters\t%\t'+'\t'.join([amplicon for amplicon in amplicons])+'\toverlap\treduction\tVIKTADReduction\tpercentage reduction\n')
        import operator
        percentageReductions = []
        for data,count in sorted(overlaps.iteritems(), key=operator.itemgetter(1))[::-1]:
        #for data,count in overlaps.iteritems():
            data = eval(data)
            shortestList = min( [data[amplicon] for amplicon in amplicons] )
            reduction = shortestList - data['inAll']
            wReduction = reduction*count
            wReductions.append(wReduction)
            if shortestList > 0: percentageReduction = 100*float(reduction)/float(shortestList)
            else:                percentageReduction = 0
            for i in xrange(count): percentageReductions.append(percentageReduction)
            outfile.write( str(count) + '\t'+str(round(100*float(count)/float(total),2)) + '%\t' + '\t'.join([str(data[amplicon]) for amplicon in amplicons]) +'\t'+str(data['inAll'])+'\t'+str(reduction)+'\t'+str(wReduction)+'\t'+str(round(percentageReduction,2))+'%\n')
        outfile.write('Average reduction: '+str(float(sum(wReductions))/float(total))+'\n')
        outfile.write('Average % reduction: '+str(round(float(sum(percentageReductions))/float(len(percentageReductions)),2))+'%\n')

    config.logfile.write('Analyzing hitlist overlaps ...\n')
    config.outfile.write('\nOverlap analysis:\n')
    config.logfile.write('All clusters ...\n')
    config.outfile.write('\nTotal '+str(temp1)+' clusters:')
    #printTable(overlaps,config.outfile,temp1)
    if not overlaps: tmp = '\n\n# WARNING: No overlaps found! maybe the classification is still running.\n'; config.logfile.write(tmp); config.outfile.write(tmp); sys.exit()
    config.logfile.write('Removing all cluster without hits for all amplicons ...\n')
    remove = []
    amplicons = getAmplicons(overlaps)
    for data,count in overlaps.iteritems():
        data_str = data
        data = eval(data)
        if 0 in [data[amplicon] for amplicon in amplicons]: remove.append(str(data_str))
    for data in remove:
        temp1 -= overlaps[data]
        del overlaps[data]

    config.outfile.write('\n'+str(temp1)+' clusters with hits:')
    #printTable(overlaps,config.outfile,temp1)

    config.logfile.write('Removing clusters without overlap ...\n')
    remove = []
    amplicons = getAmplicons(overlaps)
    for data,count in overlaps.iteritems():
        data_str = data
        data = eval(data)
        bothOk = True
        for amplicon in amplicons:
            if data['inAll'] == 0: bothOk = False
        if not bothOk:
            remove.append(str(data_str))
    for data in remove:
        temp1 -= overlaps[data]
        del overlaps[data]
    withOvelap = temp1

    if overlaps:
        config.outfile.write('\n'+str(temp1)+' clusters with overlap:')
        printTable(overlaps,config.outfile,temp1)
    
    config.logfile.write('Removing clusters without reduction ...\n')
    remove = []
    amplicons = getAmplicons(overlaps)
    for data,count in overlaps.iteritems():
        data_str = data
        data = eval(data)
        bothOk = True
        for amplicon in amplicons:
            if data[amplicon] == data['inAll'] or data['inAll'] == 0: bothOk = False
        if not bothOk:
            remove.append(str(data_str))
    for data in remove:
        temp1 -= overlaps[data]
        del overlaps[data]
    withReduction = temp1

    if overlaps:
        config.outfile.write('\n'+str(temp1)+' clusters with hitList reduction:')
        printTable(overlaps,config.outfile,temp1)
    
    config.outfile.write('\nFor '+str(round(100*float(withReduction)/float(withOvelap),2)) + '% of the clusters with hitlist overlap, is the overlap list shorter than the "best" of the defined amplicon lists.\n')

    config.outfile.write('\nOut of a total of '+str(temp2)+' cluster with BLAST hits for both amplicons there are:')
    for rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species','subspecies']:
        config.outfile.write('\n\t---------- '+rank)
        try:count = ranks[rank]-classifications[rank]["['NoInformationAvailable']"]
        except KeyError: count = 0
        count2 = classifications[rank]["['NoInformationAvailable']"]
        try:percentage = round(100*float(count2)/float(count),2)
        except ZeroDivisionError:percentage=0.00
        config.outfile.write('\n\t'+str(count2)+' clusters that are excluded as NoInformationAvailable for any BlastHit at this rank ('+str(ranks[rank]-count2)+' clusters have info at this rank).')
        count2 = classifications[rank]["['NoMatchFound']"]
        try:percentage = round(100*float(count2)/float(count),2)
        except ZeroDivisionError:percentage=0.00
        config.outfile.write('\n\t'+str(count2)+' ('+str(percentage)+'%) clusters where no overlap is found between the amplicon hitlists for this rank.')
        if count: percentage = round(100*float(count-classifications[rank]["['NoMatchFound']"])/float(temp2-classifications[rank]["['NoInformationAvailable']"]),2)
        else: percentage = 0.00
        config.outfile.write('\n\t'+str(count-classifications[rank]["['NoMatchFound']"])+' clusters ('+str(percentage)+'%) where the amplicons hitlist has >=1 match at the '+rank+' level, out of theese are:')
        #print '\n\t'+str(count)+' clusters ('+str(percentage)+'%) where the amplicons hitlist has >=1 match at the '+rank+' level, out of theese are:'
        if count:
            #values = classifications[rank]
            for value, count2 in classifications[rank].iteritems():
                if value == "['NoMatchFound']" or value == "['NoInformationAvailable']": continue
                value = eval(value)
                #print value
                percentage = round(100*float(count2)/float(count-classifications[rank]["['NoMatchFound']"]),2)
                if len(value) > 1: config.outfile.write('\n\t\t'+str(count2)+' ('+str(percentage)+'%) are '+', '.join(value[:-1])+' or '+value[-1]+'.')
                else: config.outfile.write('\n\t\t'+str(count2)+' ('+str(percentage)+'%) are '+value[0]+'.')
    config.outfile.write('\n\n')

    config.outfile.write('\nOut of a total of '+str(temp3)+' cluster with RDP info for both amplicons there are:')
    for rank in ['domain','phylum','class','order','family','genus']:
        config.outfile.write('\n\t---------- '+rank)
        try:count = rdpRanks[rank]-rdpClassifications[rank]["['LowConfidence']"]
        except KeyError: count = 0
        count2 = rdpClassifications[rank]["['LowConfidence']"]
        try:percentage = round(100*float(count2)/float(count),2)
        except ZeroDivisionError:percentage=0.00
        config.outfile.write('\n\t'+str(count2)+' are excluded due to LowConfidence ('+str(rdpRanks[rank]-count2)+' clusters have info at this rank).')
        count2 = rdpClassifications[rank]["['MissMatch']"]
        try:percentage = round(100*float(count2)/float(count),2)
        except ZeroDivisionError:percentage=0.00
        config.outfile.write('\n\t'+str(count2)+' ('+str(percentage)+'%) are MissMatch.')
        if count: percentage = round(100*float(count-rdpClassifications[rank]["['MissMatch']"])/float(temp3-rdpClassifications[rank]["['LowConfidence']"]),2)
        else: percentage = 0.00
        config.outfile.write('\n\t'+str(count-rdpClassifications[rank]["['MissMatch']"])+' clusters ('+str(percentage)+'%) where the amplicons RDPclassinfo has >=1 match at the '+rank+' level, out of theese are:')
        if count:
            for value, count2 in rdpClassifications[rank].iteritems():
                if value == "['MissMatch']" or value == "['LowConfidence']": continue
                value = eval(value)
                percentage = round(100*float(count2)/float(count-rdpClassifications[rank]["['MissMatch']"]),2)
                if len(value) > 1: config.outfile.write('\n\t\t'+str(count2)+' ('+str(percentage)+'%) are '+', '.join(value[:-1])+' or '+value[-1]+'.')
                else: config.outfile.write('\n\t\t'+str(count2)+' ('+str(percentage)+'%) are '+value[0]+'.')
    config.outfile.write('\n\n')

    config.logfile.write('Done.\n')
    infile.close()
    
    return 0