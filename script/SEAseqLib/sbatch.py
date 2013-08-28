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

    config.logfile.write('Creating sbatch scripts.\n')
    
    f = open(config.path +'/sbatch.cluster.sh','w')
    f.write(
	'#! /bin/bash -l'+'\n'+
	'#SBATCH -A b2011011'+'\n'+
	'#SBATCH -n 8 -p node'+'\n'+
	'#SBATCH -t 24:00:00'+'\n'+
	'#SBATCH -J clust_'+config.path+'\n'+
	'#SBATCH -e '+config.abspath+'/sbatch.cluster.stderr.txt'+'\n'+
	'#SBATCH -o '+config.abspath+'/sbatch.cluster.stdout.txt'+'\n'+
	'#SBATCH --mail-type=All'+'\n'+
	'#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
	'echo "$(date) Running on: $(hostname)"'+'\n'+
	'cd '+os.getcwd()+'\n'+
	'module load python/2.7'+'\n'+
	sys.argv[0]+' clusterbarcodes -path '+config.path+' -bm 2 -hm 4 -seed 2000 -p 8'+'\n'
    )
    f.close()

    f = open( config.path +'/sbatch.sortreads.sh','w')
    f.write(
	'#! /bin/bash -l'+'\n'+
	'#SBATCH -A b2011011'+'\n'+
	'#SBATCH -n 8 -p node'+'\n'+
	'#SBATCH -C fat'+'\n'+
	'#SBATCH -t 24:00:00'+'\n'+
	'#SBATCH -J sort_'+config.path+'\n'+
	'#SBATCH -e '+config.abspath+'/sbatch.sortreads.stderr.txt'+'\n'+
	'#SBATCH -o '+config.abspath+'/sbatch.sortreads.stdout.txt'+'\n'+
	'#SBATCH --mail-type=All'+'\n'+
	'#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
	'echo "$(date) Running on: $(hostname)"'+'\n'+
	'cd '+os.getcwd()+'\n'+
	'module load python/2.7'+'\n'+
	sys.argv[0]+' sortreads -path '+config.path+' -p8'+'\n'
    )
    f.close()
    
    f = open( config.path +'/sbatch.meta.sh','w')
    f.write(
	'#! /bin/bash -l'+'\n'+
	'#SBATCH -A b2011011'+'\n'+
	'#SBATCH -n 8 -p node'+'\n'+
	'#SBATCH -t 24:00:00'+'\n'+
	'#SBATCH -J meta_'+config.path+'\n'+
	'#SBATCH -e '+config.abspath+'/sbatch.meta.stderr.txt'+'\n'+
	'#SBATCH -o '+config.abspath+'/sbatch.meta.stdout.txt'+'\n'+
	'#SBATCH --mail-type=All'+'\n'+
	'#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
	'echo "$(date) Running on: $(hostname)"'+'\n'+
	'cd '+os.getcwd()+'\n'+
	'module load python/2.7'+'\n'+
	sys.argv[0]+' meta -path '+config.path+' -p8'+'\n'+
        'grep -vP "^((Read)|([0-9]+\t)|P|c|(Co))" '+config.path+'/meta.out.txt | cat -s > '+config.path+'/meta.smaller.out.txt\n'
    )
    f.close()

    f = open( config.path +'/sbatch.classify.sh','w')
    f.write(
	'#! /bin/bash -l'+'\n'+
	'#SBATCH -A b2011011'+'\n'+
	'#SBATCH -n 8 -p node'+'\n'+
	'#SBATCH -t 05:00:00'+'\n'+
	'#SBATCH -J classify_'+config.path+'\n'+
	'#SBATCH -e '+config.abspath+'/sbatch.classify.stderr.txt'+'\n'+
	'#SBATCH -o '+config.abspath+'/sbatch.classify.stdout.txt'+'\n'+
	'#SBATCH --mail-type=All'+'\n'+
	'#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
	'echo "$(date) Running on: $(hostname)"'+'\n'+
	'cd '+os.getcwd()+'\n'+
	'module load bioinfo-tools blast/2.2.28+ biopython'+'\n'+
	'module unload python/2.6.6'+'\n'+
	'module load python/2.7'+'\n'+
	sys.argv[0]+' classifymeta -path '+config.path+''+'\n'
    )
    f.close()

    f = open( config.path +'/sbatch.gzip.sh','w')
    f.write(
	'#! /bin/bash -l'+'\n'+
	'#SBATCH -A b2011011'+'\n'+
	'#SBATCH -n 8 -p node'+'\n'+
	'#SBATCH -t 72:00:00'+'\n'+
	'#SBATCH -J gzip_'+config.path+'\n'+
	'#SBATCH -e '+config.abspath+'/sbatch.gzip.stderr.txt'+'\n'+
	'#SBATCH -o '+config.abspath+'/sbatch.gzip.stdout.txt'+'\n'+
	'#SBATCH --mail-type=All'+'\n'+
	'#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'+
	'echo "$(date) Running on: $(hostname)"'+'\n'+
	'cd '+os.getcwd()+'\n'+
	'gzip -v9 '+config.path+'/sortedReads/sorted_by_barcode_cluster.1.fq &'+'\n'
	'gzip -v9 '+config.path+'/sortedReads/sorted_by_barcode_cluster.2.fq &'+'\n'
	'gzip -v9 '+config.path+'/nonCreads.1.fq &'+'\n'
	'gzip -v9 '+config.path+'/nonCreads.2.fq &'+'\n'
	'gzip -v9 '+config.path+'/meta.out.txt &'+'\n'
	'gzip -v9 '+config.path+'/classify.out.txt &'+'\n'
	'wait'+'\n'
    )
    f.close()
    
    if indata.send:
	config.logfile.write('Placing scripts in jobqueue.\n')
	import subprocess

	sbatch = subprocess.Popen( ['sbatch',config.path +'/sbatch.cluster.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	sbatch_out, errdata = sbatch.communicate()
	if sbatch.returncode != 0:
		print 'sbatch view Error code', sbatch.returncode, errdata
		print sbatch_out
		sys.exit()
	cluster_jobid = sbatch_out.split('\n')[0].split(' ')[3]
	config.logfile.write('Queued barcode clustering with JOBID '+cluster_jobid+'.\n')

	sbatch = subprocess.Popen( ['sbatch','--dependency=afterok:'+cluster_jobid,config.path +'/sbatch.sortreads.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	sbatch_out, errdata = sbatch.communicate()
	if sbatch.returncode != 0:
		print 'sbatch view Error code', sbatch.returncode, errdata
		print sbatch_out
		sys.exit()
	sort_jobid = sbatch_out.split('\n')[0].split(' ')[3]
	config.logfile.write('Queued sorting of reads with JOBID '+sort_jobid+'.\n')

	sbatch = subprocess.Popen( ['sbatch','--dependency=afterok:'+sort_jobid,config.path +'/sbatch.meta.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	sbatch_out, errdata = sbatch.communicate()
	if sbatch.returncode != 0:
		print 'sbatch view Error code', sbatch.returncode, errdata
		print sbatch_out
		sys.exit()
	meta_jobid = sbatch_out.split('\n')[0].split(' ')[3]
	config.logfile.write('Queued metagenomics analysis with JOBID '+meta_jobid+'.\n')

	sbatch = subprocess.Popen( ['sbatch','--dependency=afterok:'+meta_jobid,config.path +'/sbatch.classify.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	sbatch_out, errdata = sbatch.communicate()
	if sbatch.returncode != 0:
		print 'sbatch view Error code', sbatch.returncode, errdata
		print sbatch_out
		sys.exit()
	classify_jobid = sbatch_out.split('\n')[0].split(' ')[3]
	config.logfile.write('Queued classifygenomics analysis with JOBID '+classify_jobid+'.\n')

    
    config.logfile.write('Done.\n')
    return 0

