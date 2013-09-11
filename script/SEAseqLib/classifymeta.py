def clusterGenerator(config):
    
    import cPickle
    import gzip
    
    filename=config.path+'/meta.clusters.pickle.gz'
    clusterdump = gzip.open(filename,'rb')
    
    while True:
        try:
            cluster = cPickle.load(clusterdump)
            yield cluster
        except EOFError:
            config.logfile.write('All clusters read from file.\n')
            break

def foreachCluster(cluster,indata,config):
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
    
    print 'making fasta'
    for amplicon in cluster.definedamplicons.values():
        for consensus in amplicon.goodalleles:
            r1 = consensus.sequence.seq.split('NNNNNNNNNN')[0]
            r2 = consensus.sequence.seq.split('NNNNNNNNNN')[1]
            blastfile.write('>'+amplicon.type+'.'+str(consensus.id)+'.r1\n'+r1+'\n'+
                            '>'+amplicon.type+'.'+str(consensus.id)+'.r2\n'+r2+'\n')
    blastfile.close()
    print 'done'

    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
    from cStringIO import StringIO
    import time

    #setting up blast
    BLAST_identity = indata.identity

    database = indata.database
    cline = NcbiblastnCommandline(query=blastfile.name, db=database ,evalue=0.001, outfmt=5, num_threads=8,perc_identity=BLAST_identity)#, out=infile+'.blastout')
    #cline = NcbiblastnCommandline(query=infile, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn', out=infile+'.'+config.blastid+'.blastout')

    print ( 'Starting BLAST')
    blast_handle = cline.__call__()
    print cluster.id, str(blast_handle)[0:100].replace('\n','--NEWLINE--')
    print ('BLAST search done')    

    blast_handle = StringIO(blast_handle[0])
    blast_handle.seek(0)
    records = NCBIXML.parse(blast_handle)
    
    from SEAseqLib.mainLibrary import gi2orgname
    local_gi2org = {}
    
    identity_cutoff = 99
    alignment_length_cutoff = 95
    for blast_record in records:
        print blast_record.query
        for alignment in blast_record.alignments[:10]:
            for hsp in alignment.hsps[:1]:
                perc_identity = float(hsp.identities) 	/	float(hsp.align_length)	*100
                perc_coverage = float(hsp.align_length)	/	float(blast_record.query_letters)	*100
                gi_number = alignment.title.split(' ')[1].split('|')[1]
                try:
                    organism = local_gi2org[gi_number]
                except KeyError:
                    organism = gi2orgname(gi_number)
                    local_gi2org[gi_number] = organism
                if perc_identity >= identity_cutoff and perc_coverage >= alignment_length_cutoff:
                    print perc_identity,perc_coverage,organism
        
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
    
    for cluster in clusterGenerator(config):
            print foreachCluster(cluster,indata,config)

    import sys
    sys.exit(0)

    #CUTTNIG AND PASTING FROMHERE:::


    #f = open(infile+'.'+config.blastid+'.blastout')
    #records = NCBIXML.parse(f)
    records = NCBIXML.parse(blast_handle)

    config.printblast = True
    if config.printblast: pass
    
    config.logfile.write('Reformatting BLAST output to enable clusterwise parsing\n')
    data = {}
    for blast_record in records:

	#>cluster=1.amplicon=ITS.consensus=83_r1.27.08%_of_2626reads
	splitted = blast_record.query.split('.')
	cluster	= int(	splitted[0].split('=')[1]	)
	amplicon	= splitted[1].split('=')[1]
	consensus	= int(	splitted[2].split('=')[1].split('_')[0]	)
	readnumber	= int(	splitted[2].split('=')[1].split('_')[1].split('r')[1]	)
	
	if cluster not in data: data[cluster] = {'ITS':{},'16S':{}}

	try:			data[cluster][amplicon][consensus][readnumber] = blast_record
	except KeyError:	data[cluster][amplicon][consensus] = { readnumber : blast_record }
    config.logfile.write('Done' +'\n')

    #for the random match estimation
    matches = {}
    
    # initiating vars
    noblasthit_its = 0
    noblasthit_16s = 0
    blasthitsagree = 0
    hitsdonotagree = 0
    
    if indata.gidatabase:
	tmp_file = open(indata.gidatabase)
	tmp_string = tmp_file.read()
	tmp_file.close()
	local_gi2org = eval(tmp_string)
    else:local_gi2org = {}
    
    from SEAseqLib.mainLibrary import Progress
    config.logfile.write('Parsing BLASTreport ....\n')
    progress = Progress(len(data),unit='cluster',logfile=config.logfile)
    with progress:
	for cluster_id, amplicons in data.iteritems():
	    
	    #for the random match estimation
	    tmp_matches = {'ITS':{},'16S':{}}
	    
	    progress.update()

	    config.outfile.write(      '### Cluster '+str(cluster_id)+' ###'     +'\n')
	    identity_cutoff = BLAST_identity#%
	    alignment_length_cutoff = indata.length#95#%

	    #for amplicon, consensuses in amplicons.iteritems():
	    in_both_reads = {'ITS':[],'16S':[]}
	    for amplicon in ['ITS','16S']:
		broken = False
		consensuses = amplicons[amplicon]
		config.outfile.write(      '\t'+amplicon     +'\n')
		if not consensuses:
		    config.outfile.write(      '\t\tNo support'     +'\n')
		    if amplicon == 'ITS': break
		if len(consensuses) != 1:
		    config.outfile.write('\t\tNot monoclonal'+'\n');
		    break

		for consensus in consensuses:

		    try: r1 = consensuses[consensus][1]
		    except KeyError:
			config.outfile.write('\t\tConsensus cannot be read as read 1 sequence are not available from infile\n')
			config.logfile.write('WARNING: Running on a uncomplete input file! Maybe the SEAseq meta program is still runnning?\n')
			continue
		    try: r2 = consensuses[consensus][2]
		    except KeyError:
			config.outfile.write('\t\tConsensus cannot be read as read 2 sequence are not available from infile\n')
			config.logfile.write('WARNING: Running on a uncomplete input file! Maybe the SEAseq meta program is still runnning?\n')
			continue

		    amplicon	= 	r1.query.split('.')[1].split('=')[1]
		    support	= float(r1.query.split('r1.')[-1].split('_')[0].replace('%','')	)
		    readpop = int(r1.query.split('r1.')[-1].split('_')[-1].replace('reads','')	)

		    config.outfile.write(      '\t\tConsensus '+str(consensus)+': '+amplicon+' supported by '+str(support) +'% ('+str(readpop)+')'     +'\n')

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
			    if perc_identity >= identity_cutoff and perc_coverage >= alignment_length_cutoff: in_r2[organism] = alignment


		    for organism in in_r1:
			if organism in in_r2:
			    in_both_reads[amplicon].append(organism)

		    for organism in in_both_reads[amplicon]:
			config.outfile.write(      '\t\t\t'+organism     +'\n')
			tmp_matches[amplicon][organism] = True
		    if not in_both_reads[amplicon]:
			config.outfile.write(      '\t\t\tNo alignment supported by both reads with >='+str(identity_cutoff)+'% identity and '+str(alignment_length_cutoff)+'% alignment length coverage'     +'\n')
			if amplicon == '16S':
			    noblasthit_16s +=1
			    broken = True;
			elif amplicon == 'ITS':
			    noblasthit_its += 1
			    broken = True;
			    break;
		if broken: break

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