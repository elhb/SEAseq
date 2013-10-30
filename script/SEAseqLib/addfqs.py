def addfqs(indata):
    
    from SEAseqLib.mainLibrary import Configuration, writelogheader, bufcount

    config = Configuration(indata.path, indata.cmd)
    config.openconnections()
    
    writelogheader(config.logfile)

    config.logfile.write('Reading current config ...\n')
    config.load()

    config.logfile.write('Adding filenames to config file.\n')
    import os
    import sys
    r1 = os.path.abspath(indata.reads1.name)
    r2 = os.path.abspath(indata.reads2.name)
    if r1 in config.infilesDictionary['r1'] + config.infilesDictionary['r2'] or r2 in config.infilesDictionary['r1'] + config.infilesDictionary['r2']:
	sys.stderr.write('ERROR:\nat least one of the files:\n'+r1+'\n'+r2+'\nare already in the config file.\n');
	config.logfile.write('ERROR:\nat least one of the files:\n'+r1+'\n'+r2+'\nare already in the config file.\nExiting withiout changes.\n');
	return 1
    config.infilesDictionary['r1'].append(r1)
    config.infilesDictionary['r2'].append(r2)
    
    config.logfile.write('Getting readcount ...\n')
    config.readCountsList.append(bufcount(r1)/4)
    
    config.save()
    
    config.logfile.write('Files '+r1+' and '+r2+' sucessfully added to infiles dist in config.\n\n')
    config.logfile.close()
    return 0