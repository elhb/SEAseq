def setVariables(indata):
    from SEAseqLib.mainLibrary import Configuration, writelogheader
    config = Configuration(indata.path, indata.cmd)
    import os
    if os.path.exists(config.outfile): os.remove(config.outfile)
    config.openconnections()
    writelogheader(config.logfile)

    # settings
    config.load()
    oldValues = {}
    for varName,varValue in config.__dict__.iteritems():
        oldValues[varName] = varValue

    indata.absolutePath = os.path.abspath(indata.path)
    for varName,varValue in indata.__dict__.iteritems():
        try:
            if oldValues[varName] != varValue:
                config.__dict__[varName] = varValue
                config.logfile.write('Changeing varible '+varName+' from '+str(oldValues[varName])+' to '+str(varValue)+'.\n')
        except KeyError: config.logfile.write('Varible '+varName+' is runspecific and not stored in config.\n')
    
    config.save()
    
    config.logfile.write('Finished setting variables.\n')