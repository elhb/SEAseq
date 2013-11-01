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

    if indata.graphs == 'all': indata.graphs = 'abcdefghijklmnopqrstuvxyz';
    indata.graphs = list(indata.graphs)
    
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
    
    #['clusterid','number of reads in total','number of adaper reads','number of strange primers','its reads','16s reads','its','16s','its monoclonal','16s monoclonal','number of consensus types','number of consensus types with good support']
    
   
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
    for x_current in x_range:
	for rc_type in ['total','good']:
	    graph_info[rc_type][ x_current ] = {'all':0,'16s':0,'its':0,'both':0,'None':0,'any':0,'16s_mono':0,'its_mono':0,'both_mono':0,'both_16s_mono':0,'both_its_mono':0,'any_mono':0}

    config.logfile.write('Counting clusters ... \n')
    if os.path.exists(config.path+'/meta.statstable'):
	for cid in data:
	    for x_current in x_range:
		breaker = []
		compare_pairs = [	['total',data[cid]['number of reads in total']],    ['good',data[cid]['number of reads in total']-data[cid]['number of adaper reads']-data[cid]['number of strange primers']]    ]
		for tmp in compare_pairs:
		    [rc_type, comp_value] = tmp
		
		    if comp_value > x_current:
			graph_info[rc_type][ x_current ]['all'] += 1
    
			if data[cid]['16s'] and data[cid]['its']:
			    graph_info[rc_type][ x_current ]['both'] += 1
			    graph_info[rc_type][ x_current ]['any'] += 1
			    if data[cid]['16s monoclonal'] and data[cid]['its monoclonal']:
				graph_info[rc_type][ x_current ]['both_mono'] += 1
				graph_info[rc_type][ x_current ]['any_mono'] += 1
			    if not data[cid]['16s monoclonal'] and data[cid]['its monoclonal']:
				graph_info[rc_type][ x_current ]['both_its_mono'] += 1
			    if data[cid]['16s monoclonal'] and not data[cid]['its monoclonal']:
				graph_info[rc_type][ x_current ]['both_16s_mono'] += 1
    
			elif data[cid]['16s'] and not data[cid]['its']:
			    graph_info[rc_type][ x_current ]['16s'] += 1
			    graph_info[rc_type][ x_current ]['any'] += 1
			    if data[cid]['16s monoclonal']:
				graph_info[rc_type][ x_current ]['16s_mono'] += 1
				graph_info[rc_type][ x_current ]['any_mono'] += 1
			
			elif data[cid]['its'] and not data[cid]['16s']:
			    graph_info[rc_type][ x_current ]['its'] += 1
			    graph_info[rc_type][ x_current ]['any'] += 1
			    if data[cid]['its monoclonal']:
				graph_info[rc_type][ x_current ]['its_mono'] += 1
				graph_info[rc_type][ x_current ]['any_mono'] += 1
			    
			elif not data[cid]['16s'] and not data[cid]['its']:
			    graph_info[rc_type][ x_current ]['None'] += 1
			
			else: print 'ERROR: this should not be possible'
			
		    else: breaker.append(True)
		if len(breaker) == len(compare_pairs): break

	config.logfile.write('Calculating percentages ... \n')
	for rc_type in ['total','good']:
	    f = open( config.path+'/graphs/'+rc_type+'_read_pairs_per_barcode_cluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+'.values' ,'w' )
	    f.write(
		'x'			+'\t'+ 
		'all'		+'\t'+ 
		'None'		+'\t'+ 
		'16s'		+'\t'+ 
		'16s_mono'		+'\t'+ 
		'16s_mono %'	+'\t'+ 
		'its'		+'\t'+ 
		'its_mono'		+'\t'+ 
		'its_mono %'	+'\t'+ 
		'both'		+'\t'+ 
		'both_mono'		+'\t'+ 
		'both_mono %'	+'\t'+ 
		'both'		+'\t'+ 
		'both_16s_mono'	+'\t'+ 
		'both_16s_mono %'	+'\t'+ 
		'both'		+'\t'+ 
		'both_its_mono'	+'\t'+ 
		'both_its_mono %'	+'\t'+
		'any'		+'\t'+ 
		'any_mono'		+'\t'+ 
		'any_mono %'	+'\n'
	    )
	    f.close()
	for x_current in x_range:
	    for rc_type in ['total','good']:
		f = open( config.path+'/graphs/'+rc_type+'_read_pairs_per_barcode_cluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+'.values' ,'a' )
		f.write(str(x_current) +'\t'+     str(graph_info[rc_type][ x_current ]['all'])	+'\t'+	str(graph_info[rc_type][ x_current ]['None'])	+'\t'	)
    
		for [total_id,count_id] in [['16s','16s_mono'],['its','its_mono'],['both','both_mono'],['both','both_16s_mono'],['both','both_its_mono'],['any','any_mono']]:
		    total  = graph_info[rc_type][ x_current ][total_id]
		    count = graph_info[rc_type][ x_current ][count_id]
		    
		    f.write(	str(total)+'\t'+	str(count)+'\t'    )
		    if total:
			tmp_percentage = round(100*float(count)/float(total),2)
			graph_info[rc_type][ x_current ][count_id] = tmp_percentage
			f.write(str(tmp_percentage)+'\t')
		    else:
			graph_info[rc_type][ x_current ][count_id] = 0.0
			f.write('0.0\t')
		f.write('\t'+str(x_current)+'\n')
		f.close()

    if 'a' in indata.graphs and os.path.exists(config.path+'/meta.statstable'):
	config.logfile.write('Preparing variables for plotting ... \n')
	for rc_type in ['total','good']:

	    temp_x=graph_info[rc_type].keys()
	    temp_x.sort()
	    x =[];y1=[];y2=[];y3=[];y4=[];y5=[];y6=[];y7=[];y8=[];y9=[];y10=[]
	    for i in temp_x:
		    x.append(i)
		    y1.append(graph_info[rc_type][i]['all'])
		    y2.append(graph_info[rc_type][i]['both'])
		    y3.append(graph_info[rc_type][i]['16s'])
		    y4.append(graph_info[rc_type][i]['its'])
		    y5.append(graph_info[rc_type][i]['None'])
		    y6.append(graph_info[rc_type][i]['16s_mono'])
		    y7.append(graph_info[rc_type][i]['its_mono'])
		    y8.append(graph_info[rc_type][i]['both_mono'])
		    y9.append(graph_info[rc_type][i]['both_its_mono'])
		    y10.append(graph_info[rc_type][i]['both_16s_mono'])

	    config.logfile.write('Creating graphics ... \n')
	    import numpy as np
	    import matplotlib.pyplot as plt
	    from matplotlib import rc
		    
	    fig = plt.figure(figsize=(20, 15), dpi=100)
	    ax = fig.add_subplot(111)
	    if incomplete: ax.set_title('WARNING: incomplete dataset! '+ indata.sample+' ' +config.path)
	    else : ax.set_title(indata.sample+' ' +config.path)
	    ax.plot(x, y6, '--g', label = 'Percentage monoclonal 16S cluster')
	    ax.plot(x, y7, '--y', label = 'Percentage monoclonal ITS cluster')
	    ax.plot(x, y8, '--b', label = 'Perc. both mono. both in cluster')
	    ax.plot(x, y9, '--r', label = 'Perc. ITS mono. both in cluster')
	    ax.plot(x, y10, '--k', label = 'Perc. 16S mono. both in cluster')
	    ax2 = ax.twinx()
	    ax2.plot(x, y1, '-r', label = 'Total number of clusters')
	    ax2.plot(x, y2, '-b', label = 'Number of 16S/ITS clusters')
	    ax2.plot(x, y3, '-g', label = 'Number of 16S only clusters')
	    ax2.plot(x, y4, '-y', label = 'Number of ITS only clusters')
	    ax2.plot(x, y5, '-k', label = 'Number of undefined clusters')
	    
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

	    ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/10))
	    ax2.set_yticks(np.arange(yscale[0],yscale[1]+1,min(100,yscale[1]/10)))
	    ax.set_yticks(np.arange(0,101,5))

	    plt.savefig(config.path+'/graphs/'+rc_type+'_read_pairs_per_barcode_cluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+'.pdf')
	    config.logfile.write('Created: '+config.path+'/graphs/'+rc_type+'_read_pairs_per_barcode_cluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+' (.pdf and .values)'+'.\n')
	    plt.close()

    if 'b' in indata.graphs and os.path.exists(config.path+'/meta.statstable'):
	config.logfile.write('Preparing variables for plotting ... \n')
	for cluster_type in ['its','16s','both']:

	    temp_x=graph_info['total'].keys()
	    temp_x.sort()
	    
	    x =[];y1=[];y2=[];y3=[];y4=[];y5=[];y6=[];y7=[];y8=[];

	    for i in temp_x:
		    x.append(i)
		    y1.append(graph_info['good' ][i][cluster_type])
		    y2.append(graph_info['total'][i][cluster_type])
		    y3.append(graph_info['good' ][i][cluster_type+'_mono'])
		    y4.append(graph_info['total'][i][cluster_type+'_mono'])
		    if cluster_type == 'both':
			y5.append(graph_info['good' ][i][cluster_type+'_its_mono'])
			y6.append(graph_info['good' ][i][cluster_type+'_16s_mono'])
			y7.append(graph_info['total'][i][cluster_type+'_its_mono'])
			y8.append(graph_info['total'][i][cluster_type+'_16s_mono'])

	    config.logfile.write('Creating graphics ... \n')
	    import numpy as np
	    import matplotlib.pyplot as plt
	    from matplotlib import rc
		    
	    fig = plt.figure(figsize=(20, 15), dpi=100)
	    ax = fig.add_subplot(111)
	    
	    if incomplete:	ax.set_title('WARNING: incomplete dataset! '+'Clusters with '+cluster_type+' consensus sequence(s)'+'. '+indata.sample+' '+config.path)
	    else:		ax.set_title(				     'Clusters with '+cluster_type+' consensus sequence(s)'+'. '+indata.sample+' '+config.path)

	    ax.plot(x, y3, '--b', label =     'Perc. mono. good pairs')
	    ax.plot(x, y4, '--r', label =     'Perc. mono. totalpairs')
	    if cluster_type == 'both':
		ax.plot(x, y5, '--g', label = '% mono.its good pairs')
		ax.plot(x, y7, '--m', label = '% mono.its totalpairs')
		ax.plot(x, y6, '--c', label = '% mono.16s good pairs')
		ax.plot(x, y8, '--y', label = '% mono.16s totalpairs')

	    ax2 = ax.twinx()
	    ax2.plot(x, y1, '-b', label =     'Number of good pairs')
	    ax2.plot(x, y2, '-r', label =     'Number of totalpairs')
	    
	    lines, labels   = ax.get_legend_handles_labels()
	    lines2, labels2 = ax2.get_legend_handles_labels()
	    ax2.legend(lines + lines2, labels + labels2, loc=7)
	    #ax.legend(lines, labels, loc=7)

	    ax.grid(b=True, which='both')
	    ax.set_xlabel('Read pairs per Barcode Cluster')
	    ax.set_ylabel('Percentage Monoclonal')
	    ax2.set_ylabel('Number of Clusters')
	    
	    ax.set_ylim(0,100)
	    ax2.set_ylim(yscale[0],yscale[1])
	    
	    ax.set_xlim(xscale[0],xscale[1])
	    ax2.set_xlim(xscale[0],xscale[1])

	    ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/10))
	    ax2.set_yticks(np.arange(yscale[0],yscale[1]+1,min(100,yscale[1]/10)))
	    ax.set_yticks(np.arange(0,101,5))

	    plt.savefig(                     config.path+'/graphs/'+cluster_type+'_read_pairs_per_barcode_cluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+'.pdf')
	    config.logfile.write('Created: '+config.path+'/graphs/'+cluster_type+'_read_pairs_per_barcode_cluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+' (.pdf and .values)'+'.\n')
	    plt.close()

    if 'c' in indata.graphs and os.path.exists(config.path+'/meta.statstable'):
	config.logfile.write('Preparing variables for plotting ... \n')
	for rc_type in ['total','good']:

	    temp_x=graph_info[rc_type].keys()
	    temp_x.sort()
	    x =[];y1=[];y2=[];y3=[];y4=[];y5=[];y6=[];y7=[];y8=[];y9=[];y10=[]
	    for i in temp_x:
		    x.append(i)
		    y1.append(graph_info[rc_type][i]['all'])
		    y2.append(graph_info[rc_type][i]['any'])
		    y3.append(graph_info[rc_type][i]['None'])
		    y4.append(graph_info[rc_type][i]['any_mono'])

	    config.logfile.write('Creating graphics ... \n')
	    import numpy as np
	    import matplotlib.pyplot as plt
	    from matplotlib import rc
		    
	    fig = plt.figure(figsize=(20, 15), dpi=100)
	    ax = fig.add_subplot(111)
	    if incomplete: ax.set_title('WARNING: incomplete dataset! '+ indata.sample+' ' +config.path)
	    else : ax.set_title(indata.sample+' ' +config.path)
	    ax.plot(x, y4, '--b', label = 'Perc. mono. for defined amplicons')
	    ax2 = ax.twinx()
	    ax2.plot(x, y1, '-r', label = 'Total number of clusters')
	    ax2.plot(x, y2, '-b', label = 'Number of defined clusters')
	    ax2.plot(x, y3, '-k', label = 'Number of undefined clusters')
	    
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

	    ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/10))
	    ax2.set_yticks(np.arange(yscale[0],yscale[1]+1,min(100,yscale[1]/10)))
	    ax.set_yticks(np.arange(0,101,5))

	    plt.savefig(                     config.path+'/graphs/'+rc_type+'_pairs_per_barcode_with_amp.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+'.pdf')
	    config.logfile.write('Created: '+config.path+'/graphs/'+rc_type+'_pairs_per_barcode_with_amp.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+indata.sample+' (.pdf and .values)'+'.\n')
	    plt.close()

    if os.path.exists(config.path+'/cluster.graphStats'):
	
	f = open(config.path+'/cluster.graphStats','r')
	reads_in_clusters = eval(f.read())
	f.close()
	
	temp_x=reads_in_clusters.keys()
	temp_x.sort()
	y=[];x =[]
	for i in x_range:
		x.append(i)
		try: y.append(reads_in_clusters[i])
		except KeyError: y.append(0)
	x=x
	y=y

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib import rc
		
	fig = plt.figure(figsize=(20, 15), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_title(indata.sample+' Reads Pairs per Barcode Cluster (Raw reads directly after clustering). ' +config.path)
	ax.plot(x, y, '-b', label = 'Total number of clusters')
	
	lines, labels   = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc=7)

	ax.grid(b=True, which='both')
	ax.set_xlabel('Read pairs per Barcode Cluster (raw reads)')
	ax.set_ylabel('Number of Clusters')
	
	ax.set_ylim(yscale[0],yscale[1])
	ax.set_xlim(xscale[0],xscale[1])

	ax.set_xticks(np.arange(xscale[0],xscale[1]+1,xscale[1]/10))
	ax.set_yticks(np.arange(yscale[0],yscale[1]+1,min(100,yscale[1]/10)))
	plt.savefig(config.path+'/rawReadPairsPerBarcodeCluster.x_scale_'+str(xscale[0])+'-'+str(xscale[1])+'.y_scale_'+str(yscale[0])+'-'+str(yscale[1])+'.pdf')
	plt.close()

	config.logfile.write( 'done\n')

    config.logfile.write( 'done\n')
    return 0


