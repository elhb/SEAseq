def makegraphs(indata):

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    config.openconnections()
    writelogheader(config.logfile)

    import os
    try:
	os.mkdir(config.path+'/graphs')
    except: pass


    # settings
    config.load()

#    if indata.graphs == 'all': indata.graphs = 'abcdefghijklmnopqrstuvxyz';
#    indata.graphs = list(indata.graphs)
    
    #read indata file
    if os.path.exists(config.path+'/meta.statstable'):
	config.logfile.write('Loading data ... \n')
	f = open(config.path+'/meta.statstable','r')
	data = f.read()
	f.close()
	header = data.split('\n')[0].split('\t')
	data = data.split('\n')[1:]
	tmp = {}
	incomplete = False
	for cluster in data:
	    if not cluster: continue
	    cluster = cluster.split('\t')
	    cid = cluster[0]
	    #assert len(cluster) == len(header), '\nhead='+ str(header) +'\nclust='+ str(cluster)+'\n'
	    if len(cluster) == len(header): tmp[cid] = {}
	    else:
		incomplete = True
		print 'Warning errounus row in stat table.'
		config.logfile.write('Warning errounus row in stat table, you are creating graphs for a incomplete dataset, SEAseq meta is probably still running.\n')
		continue
	    for i in xrange(len(header)):
		name = header[i]
		value = cluster[i]
		tmp[cid][name] = eval(value)
		#try: tmp[cid][name] = eval(value)
		#except TypeError: print cid,name,value
	data = tmp
	config.logfile.write('Loaded.\n')

    config.logfile.write('Setting scale and update intervall ... \n')
    xscale = [int(i) for i in indata.xscale.split('-')]
    yscale = [int(i) for i in indata.yscale.split('-')]
    
    #['clusterid','number of reads in total','number of adaper reads','number of strange primers','its reads','16s reads','its','16s','its monoclonal','16s monoclonal','number of consensus types','number of consensus types with good support','monoclonal for all defined amplicons']    
   
    if not indata.highres:
	step = xscale[1]/200
	if indata.step:
	    step = indata.step
	step = max([step,1])
	x_range = range(xscale[0],xscale[1]+1, step)
	config.logfile.write('Updates y every '+str(step)+'th x value (x = '+', '.join([	str(i) for i in x_range[ : min( [ 5,len(x_range) ] ) ]	])+' ... '+str(xscale[1])+').\n')
    else:
	x_range = range(xscale[0],xscale[1]+1)
	config.logfile.write('Highres: updates y every x value (step = 1).\n')


    if os.path.exists(config.path+'/meta.statstable'): config.logfile.write('Using a total of '+str(len(data))+' clusters.\n')
    
    graph_info = {'good':{},'total':{}}
    graph_info_non_cumulative = {'good':{},'total':{}}
    for x_current in x_range:
	for rc_type in ['total','good']:
	    graph_info[rc_type][ x_current ] = {'totalClusterCount':0,'16s':0,'its':0,'both':0,'undefinedClusterCount':0,'definedClusterCount':0,'16s_mono':0,'its_mono':0,'both_mono':0,'both_16s_mono':0,'both_its_mono':0,'monoForAllDefinedAmps':0}
	    graph_info_non_cumulative[rc_type][ x_current ] = {'totalClusterCount':0,'undefinedClusterCount':0,'definedClusterCount':0,'monoForAllDefinedAmps':0}

    config.logfile.write('Counting clusters ... \n')
    if os.path.exists(config.path+'/meta.statstable'):
	#for cumulative count
	for cid in data:
	    compare_pairs = [	['total',data[cid]['number of reads in total']],    ['good',data[cid]['number of reads in total']-data[cid]['number of adaper reads']-data[cid]['number of strange primers']]    ]
	    for x_current in x_range:
		breaker = []
		for tmp in compare_pairs:
		    [rc_type, comp_value] = tmp
		
		    if comp_value > x_current:
			graph_info[rc_type][ x_current ]['totalClusterCount'] += 1
    
			if data[cid]['number of consensus types with good support'] >= 1:
			    graph_info[rc_type][ x_current ]['definedClusterCount'] += 1
			elif data[cid]['number of consensus types with good support'] == 0:
			    graph_info[rc_type][ x_current ]['undefinedClusterCount'] += 1
			else: print 'ERROR: this should not be possible'
			
			if data[cid]['monoclonal for all defined amplicons']:
			    graph_info[rc_type][ x_current ]['monoForAllDefinedAmps'] += 1
    
			
		    else: breaker.append(True)
		if len(breaker) == len(compare_pairs): break

	#for noncumulative count
	for cid in data:
	    compare_pairs = [	['total',data[cid]['number of reads in total']],    ['good',data[cid]['number of reads in total']-data[cid]['number of adaper reads']-data[cid]['number of strange primers']]    ]
	    for tmp in compare_pairs:
		[rc_type, comp_value] = tmp
		
		x_current = comp_value
		if x_current > xscale[1]: continue

		try:
		    graph_info_non_cumulative[rc_type][ x_current ]['totalClusterCount'] += 1
		
		    if data[cid]['number of consensus types with good support'] >= 1:
			graph_info_non_cumulative[rc_type][ x_current ]['definedClusterCount'] += 1
		    elif data[cid]['number of consensus types with good support'] == 0:
			graph_info_non_cumulative[rc_type][ x_current ]['undefinedClusterCount'] += 1
		    else: print 'ERROR: this should not be possible'
		    
		    if data[cid]['monoclonal for all defined amplicons']:
			graph_info_non_cumulative[rc_type][ x_current ]['monoForAllDefinedAmps'] += 1
		
		except KeyError: pass
		

	config.logfile.write('Calculating percentages ... \n')
	for x_current in x_range:
	    for rc_type in ['total','good']:
    
		for [total_id,count_id] in [['definedClusterCount','monoForAllDefinedAmps']]:
		    #for cumulative count
		    total  = graph_info[rc_type][ x_current ][total_id]
		    count = graph_info[rc_type][ x_current ][count_id]
		    if total:
			tmp_percentage = round(100*float(count)/float(total),2)
			graph_info[rc_type][ x_current ][count_id] = tmp_percentage
		    else:
			graph_info[rc_type][ x_current ][count_id] = 0.0
		    
		    #for noncumulative count
		    total  = graph_info_non_cumulative[rc_type][ x_current ][total_id]
		    count = graph_info_non_cumulative[rc_type][ x_current ][count_id]		    
		    if total:
			tmp_percentage = round(100*float(count)/float(total),2)
			graph_info_non_cumulative[rc_type][ x_current ][count_id] = tmp_percentage
		    else:
			graph_info_non_cumulative[rc_type][ x_current ][count_id] = 0.0

    #for cumulative count
    if os.path.exists(config.path+'/meta.statstable'):
	config.logfile.write('Preparing variables for plotting ... \n')
	for rc_type in ['total','good']:

	    temp_x=graph_info[rc_type].keys()
	    temp_x.sort()
	    x =[];y1=[];y2=[];y3=[];y4=[];
	    for i in temp_x:
		    x.append(i)
		    y1.append(graph_info[rc_type][i]['totalClusterCount'])
		    y2.append(graph_info[rc_type][i]['definedClusterCount'])
		    y3.append(graph_info[rc_type][i]['undefinedClusterCount'])
		    y4.append(graph_info[rc_type][i]['monoForAllDefinedAmps'])

	    config.logfile.write('Creating graphics ... \n')
	    import numpy as np
	    import matplotlib.pyplot as plt
	    from matplotlib import rc
		    
	    fig = plt.figure(figsize=(20, 15), dpi=100)
	    ax = fig.add_subplot(111)
	    if incomplete: ax.set_title('WARNING: incomplete dataset! '+ config.jobName+' ' +config.path)
	    else : ax.set_title(config.jobName+' ' +config.path)
	    ax.plot(x, y4, '--b', label = 'Perc. mono. for all defined amplicons')
	    ax2 = ax.twinx()
	    ax2.plot(x, y1, '-b', label = 'Total number of clusters')
	    ax2.plot(x, y2, '-g', label = 'Number of defined clusters')
	    ax2.plot(x, y3, '-r', label = 'Number of undefined clusters')
	    
	    lines, labels   = ax.get_legend_handles_labels()
	    lines2, labels2 = ax2.get_legend_handles_labels()
	    ax2.legend(lines + lines2, labels + labels2, loc=7)
	    #ax.legend(lines, labels, loc=7)

	    ax.grid(b=True, which='both')
	    ax.set_xlabel('Read pairs per Barcode Cluster ( '+rc_type+' reads )')
	    ax.set_ylabel('Percentage Monoclonal')
	    ax2.set_ylabel('Number of Clusters')
	    
	    ax.set_ylim(0,100)
	    ax2.set_ylim(yscale[0],yscale[1])
	    
	    ax.set_xlim(xscale[0],xscale[1])
	    ax2.set_xlim(xscale[0],xscale[1])

	    ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/20))
	    ax2.set_yticks(np.arange(yscale[0],yscale[1]+1,yscale[1]/20))
	    ax.set_yticks(np.arange(0,101,5))

	    plt.savefig(                     config.path+'/graphs/'+rc_type+'_pairs_per_barcode_with_amp.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+config.jobName+'.pdf')
	    config.logfile.write('Created: '+config.path+'/graphs/'+rc_type+'_pairs_per_barcode_with_amp.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+config.jobName+' (.pdf and .values)'+'.\n')
	    plt.close()

    #for noncumulative count
    if os.path.exists(config.path+'/meta.statstable'):
	config.logfile.write('Preparing variables for plotting ... \n')
	for rc_type in ['total','good']:

	    temp_x=graph_info_non_cumulative[rc_type].keys()
	    temp_x.sort()
	    x =[];y1=[];y2=[];y3=[];y4=[];
	    for i in temp_x:
		    x.append(i)
		    y1.append(graph_info_non_cumulative[rc_type][i]['totalClusterCount'])
		    y2.append(graph_info_non_cumulative[rc_type][i]['definedClusterCount'])
		    y3.append(graph_info_non_cumulative[rc_type][i]['undefinedClusterCount'])
		    y4.append(graph_info_non_cumulative[rc_type][i]['monoForAllDefinedAmps'])

	    config.logfile.write('Creating graphics ... \n')
	    import numpy as np
	    import matplotlib.pyplot as plt
	    from matplotlib import rc
		    
	    fig = plt.figure(figsize=(20, 15), dpi=100)
	    ax = fig.add_subplot(111)
	    if incomplete: ax.set_title('WARNING: incomplete dataset! '+ config.jobName+' ' +config.path)
	    else : ax.set_title(config.jobName+' ' +config.path)
	    ax.plot(x, y4, '--b', label = 'Perc. mono. for all defined amplicons')
	    ax2 = ax.twinx()
	    ax2.plot(x, y1, '-b', label = 'Total number of clusters')
	    ax2.plot(x, y2, '-g', label = 'Number of defined clusters')
	    ax2.plot(x, y3, '-r', label = 'Number of undefined clusters')
	    
	    lines, labels   = ax.get_legend_handles_labels()
	    lines2, labels2 = ax2.get_legend_handles_labels()
	    ax2.legend(lines + lines2, labels + labels2, loc=7)
	    #ax.legend(lines, labels, loc=7)

	    ax.grid(b=True, which='both')
	    ax.set_xlabel('Read pairs per Barcode Cluster ( '+rc_type+' reads )')
	    ax.set_ylabel('Percentage Monoclonal')
	    ax2.set_ylabel('Number of Clusters')
	    
	    ax.set_ylim(0,100)
	    ax2.set_ylim(yscale[0],yscale[1])
	    
	    ax.set_xlim(xscale[0],xscale[1])
	    ax2.set_xlim(xscale[0],xscale[1])

	    ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/20))
	    ax2.set_yticks(np.arange(yscale[0],yscale[1]+1,yscale[1]/20))
	    ax.set_yticks(np.arange(0,101,5))

	    plt.savefig(                     config.path+'/graphs/'+rc_type+'_pairs_per_barcode_with_amp.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+config.jobName+'.noncumulative.pdf')
	    config.logfile.write('Created: '+config.path+'/graphs/'+rc_type+'_pairs_per_barcode_with_amp.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+config.jobName+' (.pdf and .values)'+'.\n')
	    plt.close()

    if os.path.exists(config.path+'/cluster.graphStats'):
	
	f = open(config.path+'/cluster.graphStats','r')
	reads_in_clusters = eval(f.read())
	f.close()
	
	temp_x=reads_in_clusters.keys()
	temp_x.sort()
	y=[];x =[];y2=[];
	for i in x_range:
		x.append(i)
		cumulative = []
		for i2 in range(i,max(reads_in_clusters.keys()),1):
		    try: cumulative.append( reads_in_clusters[i2] )
		    except KeyError: pass
		y.append( sum(cumulative) )
		try: y2.append(reads_in_clusters[i])
		except KeyError: y2.append(0)

	x=x
	y=y

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib import rc
		
	fig = plt.figure(figsize=(20, 15), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_title(config.jobName+' Reads Pairs per Barcode Cluster (Raw reads directly after clustering). ' +config.path)
	ax.plot(x, y, '-b', label = 'Cumulative Cluster Count')
	ax.plot(x, y2, '-r', label = 'Non-cumulative Cluster Count')
	
	lines, labels   = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc=7)

	ax.grid(b=True, which='both')
	ax.set_xlabel('Read pairs per Barcode Cluster (raw reads)')
	ax.set_ylabel('Number of Clusters')
	
	ax.set_ylim(yscale[0],yscale[1])
	ax.set_xlim(xscale[0],xscale[1])

	ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/20))
	ax.set_yticks(np.arange(yscale[0],yscale[1]+1,yscale[1]/20))
	plt.savefig(config.path+'/graphs/rawReadPairsPerBarcodeCluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+'.pdf')
	config.logfile.write('Created: '+config.path+'/graphs/rawReadPairsPerBarcodeCluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+'.pdf'+'.\n')
	plt.close()

	config.logfile.write( 'done\n')

    if os.path.exists(config.path+'/clusters.pickle') or os.path.exists(config.path+'/clusters.pickle.gz'):
	#from SEAseqLib.classifymeta import clusterGenerator
	
	incomplete2 = False
	
	config.logfile.write('Loading all amplicon alleles:\n ');
        import sys
        
	from SEAseqLib.mainLibrary import Progress, clusterGenerator
        progress = Progress(config.numberOfBarcodeClustersIdentified, logfile=config.logfile, unit='cluster',mem=True, printint = 5)
        if indata.stop: progress = Progress(indata.stop, logfile=config.logfile, unit='cluster',mem=True, printint = 5)
        
	readsPerAmpliconAlleleDict = {'total':{'total':0}}
	
	with progress:
            tmcounter = 0
            for cluster in clusterGenerator(config, indata):
                [cluster,config,indata] = cluster
		tmcounter +=1
                progress.update()
                if tmcounter <= indata.skip: continue
                
		# do something wit the cluster
		for amplicon in cluster.amplicons.values():
		    readsPerAmpliconAlleleDict['total']['total'] +=1
		    try:             readsPerAmpliconAlleleDict[amplicon.type]['total'] +=1
		    except KeyError: readsPerAmpliconAlleleDict[amplicon.type]= {'total':1}
		    
		    for consensus in amplicon.allels:
			try:             readsPerAmpliconAlleleDict['total'][consensus.readcount] += 1
			except KeyError: readsPerAmpliconAlleleDict['total'][consensus.readcount]  = 1
			try:             readsPerAmpliconAlleleDict[amplicon.type][consensus.readcount] += 1
			except KeyError: readsPerAmpliconAlleleDict[amplicon.type][consensus.readcount]  = 1
		
                if indata.stop and tmcounter >= indata.stop: break
        config.logfile.write('all alleles loaded.\n')

	config.logfile.write('Preparing allele readcounts for plotting ... \n')

	temp_x = x_range
	plots = readsPerAmpliconAlleleDict.keys()
	
	xs = {};ys = {};
	for i2 in range(len(plots)): xs[i2] = [];ys[i2] = [];
	for i in temp_x:
		
		for i2 in range(len(plots)):
		    xs[i2].append(i)
		    try:   	     ys[i2].append(readsPerAmpliconAlleleDict[ plots[i2] ][i])
		    except KeyError: ys[i2].append(0)

	config.logfile.write('Creating graphics ... \n')
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib import rc
		
	fig = plt.figure(figsize=(20, 15), dpi=100)
	ax = fig.add_subplot(111)
	if incomplete2: ax.set_title('WARNING: incomplete dataset! '+ config.jobName+' ' +config.path)
	else : ax.set_title(config.jobName+' ' +config.path)
	for i2 in range(len(plots)): ax.plot(xs[i2], ys[i2], label = 'Reads per Amplicon variant -- '+plots[i2])
	#ax.plot(x, y4, '--b', label = 'Reads per Amplicon')
	#ax2 = ax.twinx()
	#ax2.plot(x, y1, '-b', label = 'Total number of clusters')
	#ax2.plot(x, y2, '-g', label = 'Number of defined clusters')
	#ax2.plot(x, y3, '-r', label = 'Number of undefined clusters')
	
	lines, labels   = ax.get_legend_handles_labels()
	#lines2, labels2 = ax2.get_legend_handles_labels()
	#ax2.legend(lines + lines2, labels + labels2, loc=7)
	ax.legend(lines, labels, loc=7)

	ax.grid(b=True, which='both')
	ax.set_xlabel('Read pairs per amplicon allele/variant/consensus')
	ax.set_ylabel('Number of Amplicons alleles (or variants or consensuses)')
	#ax2.set_ylabel('Number of Clusters')
	
	#ax.set_ylim(0,100)
	ax.set_ylim(yscale[0],yscale[1])
	
	ax.set_xlim(xscale[0],xscale[1])
	#ax2.set_xlim(xscale[0],xscale[1])

	#x.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/20))
	ax.set_yticks(np.arange(yscale[0],yscale[1]+1,yscale[1]/20))
	#ax.set_yticks(np.arange(0,101,5))

	plt.savefig(                     config.path+'/graphs/pairs_per_ampAllele.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+config.jobName+'.noncumulative.pdf')
	config.logfile.write('Created: '+config.path+'/graphs/pairs_per_ampAllele.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+config.jobName+'.noncumulative.pdf'+'.\n')
	plt.close()


    config.logfile.write( 'done\n')
    return 0


