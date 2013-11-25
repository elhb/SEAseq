def clone(indata):

    # indata has:
    # '-clonePath'	    dest='clonePath'
    # '--clusterAndSort'    dest='cloneClusterAndSort'
    # '--meta'	            dest='cloneMeta'

    import os
    import shutil
    from SEAseqLib.init import init
    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = init(indata)
    if config == 1:
        import sys
        sys.exit()

    indata.clonePath = os.path.abspath(indata.clonePath)

    config.logfile.write('Cloning from '+indata.clonePath+'.\n')
    for src, dst in [
        [indata.clonePath+'/config',indata.path+'/config'],
        [indata.clonePath+'/addfqs.log.txt',indata.path+'/addfqs.log.txt'],
        [indata.clonePath+'/addfqs.out.txt',indata.path+'/addfqs.out.txt'],
        [indata.clonePath+'/setVariables.log.txt',indata.path+'/setVariables.log.txt'],
        [indata.clonePath+'/setVariables.out.txt',indata.path+'/setVariables.out.txt']
        ]:
        config.logfile.write('Copying from source: '+src+' -> '+dst+'.\n')
        shutil.copy(src,dst)
    toLink = []
    if indata.cloneClusterAndSort or indata.cloneMeta:
        toLink += [[indata.clonePath+'/barcode_clusters_dictionary',indata.path+'/barcode_clusters_dictionary']]
        toLink += [[indata.clonePath+'/cluster.graphStats',indata.path+'/cluster.graphStats']]
        toLink += [[indata.clonePath+'/cluster.log.txt',indata.path+'/cluster.log.txt']]
        toLink += [[indata.clonePath+'/cluster.out.txt',indata.path+'/cluster.out.txt']]
        toLink += [[indata.clonePath+'/predetermined_cluster_centers.fa',indata.path+'/predetermined_cluster_centers.fa']]
        toLink += [[indata.clonePath+'/raw_barcode_sequences.fa',indata.path+'/raw_barcode_sequences.fa']]
        config.logfile.write('creating sortedReads folder.\n')
        os.mkdir(indata.path+'/sortedReads')
        if os.path.exists(indata.clonePath+'/sortedReads/sorted_by_barcode_cluster.1.fq'):
            toLink += [[indata.clonePath+'/sortedReads/sorted_by_barcode_cluster.1.fq',indata.path+'/sortedReads/sorted_by_barcode_cluster.1.fq']]
            toLink += [[indata.clonePath+'/sortedReads/sorted_by_barcode_cluster.2.fq',indata.path+'/sortedReads/sorted_by_barcode_cluster.2.fq']]
        elif os.path.exists(indata.clonePath+'/sortedReads/sorted_by_barcode_cluster.1.fq.gz'):
            toLink += [[indata.clonePath+'/sortedReads/sorted_by_barcode_cluster.1.fq.gz',indata.path+'/sortedReads/sorted_by_barcode_cluster.1.fq.gz']]
            toLink += [[indata.clonePath+'/sortedReads/sorted_by_barcode_cluster.2.fq.gz',indata.path+'/sortedReads/sorted_by_barcode_cluster.2.fq.gz']]
        toLink += [[indata.clonePath+'/sort.log.txt',indata.path+'/sort.log.txt']]
        toLink += [[indata.clonePath+'/sort.out.txt',indata.path+'/sort.out.txt']]
    if indata.cloneMeta:
        if os.path.exists(indata.clonePath+'/meta.clusters.pickle'):
            toLink += [[indata.clonePath+'/meta.clusters.pickle',indata.path+'/meta.clusters.pickle']]
            toLink += [[indata.clonePath+'/meta.smaller.out.txt',indata.path+'/meta.smaller.out.txt']]
            toLink += [[indata.clonePath+'/meta.statstable',indata.path+'/meta.statstable']]
        elif os.path.exists(indata.clonePath+'/meta.clusters.pickle.gz'):
            toLink += [[indata.clonePath+'/meta.clusters.pickle',indata.path+'/meta.clusters.pickle.gz']]
            toLink += [[indata.clonePath+'/meta.smaller.out.txt',indata.path+'/meta.smaller.out.txt.gz']]
            toLink += [[indata.clonePath+'/meta.statstable',indata.path+'/meta.statstable.gz']]
        toLink += [[indata.clonePath+'/meta.log.txt',indata.path+'/meta.log.txt']]
        toLink += [[indata.clonePath+'/meta.out.txt',indata.path+'/meta.out.txt']]
    for src, dst in toLink:
        config.logfile.write('Linking from source: '+src+' -> '+dst+'.\n')
        os.symlink(src,dst)

    #settings
    config.load()

    oldValues = {}
    for varName,varValue in config.__dict__.iteritems():
        oldValues[varName] = varValue
    
    for varName,varValue in [['path',indata.path],['absolutePath',os.path.abspath(indata.path)]]:
        try:
            if oldValues[varName] != varValue:
                config.__dict__[varName] = varValue
                config.logfile.write('Changeing varible '+varName+' from '+str(oldValues[varName])+' to '+str(varValue)+'.\n')
        except KeyError: config.logfile.write('Varible '+varName+' is runspecific and not stored in config.\n')
    
    config.save()
    
    config.logfile.write('Finished cloning.\n')
    config.outfile.write('Finished cloning.\n')