#! /bin/env python
#! /usr/bin/env python2.7

MASTER = os.getppid()

def main():
#--------------------------MAIN PROGRAM-----------------------------
    indata = getindata()
    
    indata.logfile.write('Starting.\n')
    indata.logfile.write('Running script '+time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+'.\n')
    indata.logfile.write('Master process id='+str(MASTER)+'\n')
    indata.logfile.write('cmd: '+' '.join(sys.argv)+'\n')
    
    # Part1 - Parallel - each read pair:
    # check all reads and tryt to find "C handle"
    # and thereby also the barcode also try to identify the amplicon parts and any adapter sequences if present
    # output some initial statistics
    
    # Part2 - Serial:
    # some kind of clustering of barcodes with predetermined number of missmatches
    
    # Part3 - Parallel? - each barocde group: (assumes that we have monoclonal beads)
    # calculate statistics, level of clonality of amplicons, 

    # Part4 - Serial:
    # statistics summary
    
#--------------------------MAIN PROGRAM END-------------------------

#--------------------- Functions, Subroutines and classes --------------------	
def getindata():
    import argparse
    argparser = argparse.ArgumentParser(description='Analysis of SEAseq amplicon data.', formatter_class=argparse.RawTextHelpFormatter,)
    argparser.add_argument(	'--debug',		dest='debug', 			action='store_true', 			required=False,	default=False,	help='debug (run as regular single process python script).')
    argparser.add_argument(	'-skip',		dest='skip',	metavar='N',				type=int,	required=False,	default=0,	help='skip the first N read pairs in files (default 0).')
    argparser.add_argument(	'-stop',		dest='stop',	metavar='N',				type=int,	required=False,	default=0,	help='stop after N read pairs, set to 0 to disable (default 0).')
    argparser.add_argument(	'-r1',			dest='reads1',	metavar='FILE',				type=file,	required=True, 			help='indata "fastq"-file read1.')
    argparser.add_argument(	'-r2',			dest='reads2',	metavar='FILE',				type=file,	required=True,	default=None,	help='indata "fastq"-file read2.')
    argparser.add_argument(	'-p',			dest='cpus',	metavar='N',				type=int,	required=False,	default=1,	help='The number of processes to start (default 1).')
    argparser.add_argument(	'-o',			dest='outfile',	metavar='outfile',			type=str,	required=False,	default=False,	help='Print output to outfile (default stdout).')
    argparser.add_argument(	'-l',			dest='logfile',	metavar='logfile',			type=str,	required=False,	default=False,	help='Print log messages to logfile (default stderr).')
    argparser.add_argument(	'-random',		dest='n',	metavar='N',				type=int,	required=False,	default=0,	help='Use a random subset of N read pairs, this option is slower (default 0 = off). Can not be used in combination with "-skip" or "-stop"')
    indata = argparser.parse_args(sys.argv[1:])
    
    if indata.outfile: indata.outfile = open(indata.outfile, 'w',1)
    else: indata.outfile = sys.stdout
    
    if indata.logfile: indata.logfile = open(indata.logfile, 'w',1)
    else: indata.logfile = sys.stderr
    
    # get the readcount
    indata.logfile.write('Getting readcount ... ')
    indata.numreads=bufcount(indata.reads1.name)/4
    indata.logfile.write(str(numreads)+' read pairs in fastq files.\n');

    # calculate the number of reads to process
    indata.reads2process = indata.numreads
    if indata.skip: indata.reads2process -= indata.skip
    if indata.stop: indata.reads2process = indata.stop
    if indata.n:    indata.reads2process = indata.n

    return indata



#####
#check if run or imported // call main() or not
#####
if __name__ == "__main__":
    main()
#END of script
