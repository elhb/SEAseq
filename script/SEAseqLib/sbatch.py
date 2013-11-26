def sbatch(indata):

    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    config.openconnections()
    writelogheader(config.logfile)

    #settings
    config.logfile.write('Get infiles from config-file.\n')
    config.load()
    config.getreads2process()
    
    
    import os
    import sys
    uppmaxprojectid = indata.project
    from socket import gethostname
    if gethostname().split('.')[0][:5] == 'milou':corenumber = 16
    else: corenumber = 8

    if not indata.sendonly:
        config.logfile.write('Creating sbatch scripts.\n')
        
        if not indata.skipCluster:
            f = open(config.path +'/sbatch.cluster.sh','w')
            f.write(
                '#! /bin/bash -l'+'\n'+
                '#SBATCH -A '+uppmaxprojectid+''+'\n'+
                '#SBATCH -n '+str(corenumber)+' -p node'+'\n')
            if not indata.small:
                    f.write('#SBATCH -t 24:00:00'+'\n')
            else:
                    f.write('#SBATCH -t 1:00:00'+'\n')
            f.write('#SBATCH -J clust_'+config.jobName+'_'+config.path+'\n')
            f.write('#SBATCH -e '+config.absolutePath+'/sbatch.cluster.stderr.txt'+'\n')
            f.write('#SBATCH -o '+config.absolutePath+'/sbatch.cluster.stdout.txt'+'\n')
            f.write('#SBATCH --mail-type=All'+'\n')
            f.write('#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n')
            f.write('echo "$(date) Running on: $(hostname)"'+'\n')
            f.write('cd '+os.getcwd()+'\n')
            f.write('module load python/2.7'+'\n')
            f.write(sys.argv[0]+' clusterbarcodes -path '+config.path+' -p '+str(indata.cpus)+' -stop '+str(indata.stop)+' -skip '+str(indata.skip))
            f.write('\n')
            f.close()
        
        if not indata.skipSort:
            f = open( config.path +'/sbatch.sortreads.sh','w')
            f.write(
                '#! /bin/bash -l'+'\n'+
                '#SBATCH -A '+uppmaxprojectid+''+'\n'+
                '#SBATCH -n '+str(corenumber)+' -p node'+'\n')
            if not indata.small:
                    f.write('#SBATCH -t 24:00:00'+'\n')
                    if gethostname().split('.')[0][:5] != 'milou':f.write('#SBATCH -C fat'+'\n')
            else:
                    f.write('#SBATCH -t 1:00:00'+'\n')
            f.write(
                '#SBATCH -J sort_'+config.jobName+'_'+config.path+'\n'+
                '#SBATCH -e '+config.absolutePath+'/sbatch.sortreads.stderr.txt'+'\n'+
                '#SBATCH -o '+config.absolutePath+'/sbatch.sortreads.stdout.txt'+'\n'+
                '#SBATCH --mail-type=All'+'\n'+
                '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
                'echo "$(date) Running on: $(hostname)"'+'\n'+
                'cd '+os.getcwd()+'\n'+
                'module load python/2.7'+'\n'+
                sys.argv[0]+' sortreads -path '+config.path+' -p '+str(indata.cpus)+' -stop '+str(indata.stop)+' -skip '+str(indata.skip))
            f.write('\n')
            f.close()
        
        if not indata.skipMeta:
            f = open( config.path +'/sbatch.meta.sh','w')
            f.write(
                '#! /bin/bash -l'+'\n'+
                '#SBATCH -A '+uppmaxprojectid+''+'\n'+
                '#SBATCH -n '+str(corenumber)+' -p node'+'\n')
            if not indata.small:
                    f.write('#SBATCH -t 24:00:00'+'\n')
            else:
                    f.write('#SBATCH -t 1:00:00'+'\n')
            f.write(
                '#SBATCH -J meta_'+config.jobName+'_'+config.path+'\n'+
                '#SBATCH -e '+config.absolutePath+'/sbatch.meta.stderr.txt'+'\n'+
                '#SBATCH -o '+config.absolutePath+'/sbatch.meta.stdout.txt'+'\n'+
                '#SBATCH --mail-type=All'+'\n'+
                '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
                'echo "$(date) Running on: $(hostname)"'+'\n'+
                'cd '+os.getcwd()+'\n'+
                'module load python/2.7'+'\n'+
                sys.argv[0]+' meta -path '+config.path+' -p '+str(indata.cpus))
            f.write(
                '\n'+
                'grep -vP "^((Read)|([0-9]+\t)|P|c|(Co))" '+config.path+'/meta.out.txt | cat -s > '+config.path+'/meta.smaller.out.txt\n'
            )
            f.close()

        if not indata.skipClassify:    
            f = open( config.path +'/sbatch.classify.sh','w')
            f.write(
                '#! /bin/bash -l'+'\n'+
                '#SBATCH -A '+uppmaxprojectid+''+'\n'+
                '#SBATCH -n '+str(corenumber)+' -p node'+'\n')
            if not indata.small:
                    f.write('#SBATCH -t 15:00:00'+'\n')
            else:
                    f.write('#SBATCH -t 1:00:00'+'\n')
            f.write(
                '#SBATCH -J classify_'+config.jobName+'_'+config.path+'\n'+
                '#SBATCH -e '+config.absolutePath+'/sbatch.classify.stderr.txt'+'\n'+
                '#SBATCH -o '+config.absolutePath+'/sbatch.classify.stdout.txt'+'\n'+
                '#SBATCH --mail-type=All'+'\n'+
                '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
                'echo "$(date) Running on: $(hostname)"'+'\n'+
                'cd '+os.getcwd()+'\n'+
                'module load bioinfo-tools blast/2.2.28+ biopython'+'\n'+
                'module unload python/2.6.6'+'\n'+
                'module load python/2.7'+'\n'+
                sys.argv[0]+' classifymeta -path '+config.path+' -p '+str(indata.cpus)+'\n'+
                'grep "##### SUMMARY #####" '+config.path+'/classify.out.txt -A 10000000 > '+config.path+'/classify.summary.txt'+'\n'
            )
            f.close()
    
        f = open( config.path +'/sbatch.gzip.sh','w')
        f.write(
            '#! /bin/bash -l'+'\n'+
            '#SBATCH -A '+uppmaxprojectid+''+'\n'+
            '#SBATCH -n '+str(corenumber)+' -p node'+'\n'+
            '#SBATCH -t 72:00:00'+'\n'+
            '#SBATCH -J gzip_'+config.jobName+'_'+config.path+'\n'+
            '#SBATCH -e '+config.absolutePath+'/sbatch.gzip.stderr.txt'+'\n'+
            '#SBATCH -o '+config.absolutePath+'/sbatch.gzip.stdout.txt'+'\n'+
            '#SBATCH --mail-type=All'+'\n'+
            '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
            'echo "$(date) Running on: $(hostname)"'+'\n'+
            'cd '+os.getcwd()+'\n'+
            'gzip -v9 '+config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq &'+'\n'
            'gzip -v9 '+config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq &'+'\n'
            'gzip -v9 '+config.path+'/nonCreads.1.fq &'+'\n'
            'gzip -v9 '+config.path+'/nonCreads.2.fq &'+'\n'
            'gzip -v9 '+config.path+'/meta.out.txt &'+'\n'
            'gzip -v9 '+config.path+'/meta.smaller.out.txt &'+'\n'
            'gzip -v9 '+config.path+'/meta.statstable &'+'\n'
            'gzip -v9 '+config.path+'/meta.clusters.pickle &'+'\n')
        if corenumber == 8: f.write('wait'+'\n')
        f.write('gzip -v9 '+config.path+'/classify.out.txt &'+'\n'
            'gzip -v9 '+config.path+'/classify.statstable &'+'\n'
            'gzip -v9 '+config.path+'/classify.clusters.pickle &'+'\n'            
            'wait'+'\n'
        )
        f.close()
    
    if indata.send or indata.sendonly:
	config.logfile.write('Placing scripts in jobqueue.\n')
	import subprocess

	if not indata.skipCluster:
            command = ['sbatch',config.path +'/sbatch.cluster.sh']
            sbatch = subprocess.Popen( command, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            sbatch_out, errdata = sbatch.communicate()
            if sbatch.returncode != 0:
                    print 'sbatch view Error code', sbatch.returncode, errdata
                    print sbatch_out
                    sys.exit()
            cluster_jobid = sbatch_out.split('\n')[0].split(' ')[3]
            config.logfile.write('Queued barcode clustering with JOBID '+cluster_jobid+'.\n')
        else: cluster_jobid = 0

	if not indata.skipSort:
            if cluster_jobid:   command = ['sbatch','--dependency=afterok:'+str(cluster_jobid),config.path +'/sbatch.sortreads.sh']
            else:               command = ['sbatch',                                           config.path +'/sbatch.sortreads.sh']
            sbatch = subprocess.Popen( command, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            sbatch_out, errdata = sbatch.communicate()
            if sbatch.returncode != 0:
                    print 'sbatch view Error code', sbatch.returncode, errdata
                    print sbatch_out
                    sys.exit()
            sort_jobid = sbatch_out.split('\n')[0].split(' ')[3]
            config.logfile.write('Queued sorting of reads with JOBID '+sort_jobid+'.\n')
        else: sort_jobid = 0

	if not indata.skipMeta:
            if sort_jobid: command = ['sbatch','--dependency=afterok:'+str(sort_jobid),config.path +'/sbatch.meta.sh']
            else:          command = ['sbatch',                                        config.path +'/sbatch.meta.sh']
            sbatch = subprocess.Popen( command, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            sbatch_out, errdata = sbatch.communicate()
            if sbatch.returncode != 0:
                    print 'sbatch view Error code', sbatch.returncode, errdata
                    print sbatch_out
                    sys.exit()
            meta_jobid = sbatch_out.split('\n')[0].split(' ')[3]
            config.logfile.write('Queued metagenomics analysis with JOBID '+meta_jobid+'.\n')
        else: meta_jobid = 0

	if not indata.skipClassify:
            if meta_jobid:  command = ['sbatch','--dependency=afterok:'+str(meta_jobid),config.path +'/sbatch.classify.sh']
            else:           command = ['sbatch',                                        config.path +'/sbatch.classify.sh']
            sbatch = subprocess.Popen( command, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            sbatch_out, errdata = sbatch.communicate()
            if sbatch.returncode != 0:
                    print 'sbatch view Error code', sbatch.returncode, errdata
                    print sbatch_out
                    sys.exit()
            classify_jobid = sbatch_out.split('\n')[0].split(' ')[3]
            config.logfile.write('Queued classifygenomics analysis with JOBID '+classify_jobid+'.\n')
        else: classify_jobid = 0

    config.logfile.write('Done.\n')
    return 0