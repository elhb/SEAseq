def init(indata):
    
    from SEAseqLib.mainLibrary import Configuration, writelogheader
    import os
    import sys
    tmp_log = ''

    if os.path.exists(indata.path+'/'+'config'):
	sys.stderr.write('This analysis has already been initiated try another command.\n')
	return 1
    
    config = Configuration(indata.path, indata.cmd)
    
    try:
	os.mkdir(indata.path)
	tmp_log += 'Folder '+indata.path+' created sucessfully.\n'
    except OSError as inst:
	if inst[0] != 17: print inst; return
	_continue = 'yes'
	while not _continue or _continue[0] not in ['Y','y', 'N','n']:
	    _continue = raw_input('WARNING: the folder '+indata.path+' already excists. Continue anyway? (yes/no) ')
	tmp_log += 'WARNING: the folder '+indata.path+' already excists. Continue anyway? (yes/no) '+_continue+'\nUsing already exsisting folder ...\n'
	if _continue[0] in ['Y','y']: pass
	elif _continue[0] in ['N','n']: return
	else:
	    sys.stderr.write('Error 1\n.'); return 1

    config.openconnections()

    import time
    writelogheader(config.logfile)
    config.logfile.write(tmp_log)

    config.logfile.write('Creating config file:\n')
    config.save()

    config.logfile.write('Analysis '+config.absolutePath+' sucesfully initiated.\n')
    return config