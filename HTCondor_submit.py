#!/usr/bin/env python

'''Add important stuff'''
import os.path
import optparse
import subprocess

os.system("mkdir -p cond")

''' Inputs for the skim code '''
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--infiles', dest='infiles', help='input list of files', default='', type='string')
parser.add_option('-v', '--infilesv0', dest='infilesv0', help='V0 input list of files', default='', type='string')
parser.add_option('-o', '--outfiles', dest='outfiles', help='output files', default='', type='string')
parser.add_option('-s', '--subfiles', dest='subfiles', help='HTCondor submission file', default='', type='string')
(opt, args) = parser.parse_args()
inFiles = opt.infiles
inFilesV0 = opt.infilesv0
outFiles = opt.outfiles
subFiles = opt.subfiles

''' Read list of files '''
listOfFiles = open(inFiles+'.txt', 'r')
Lines = listOfFiles.readlines()
print ("Number of files/jobs: "+str(len(Lines)))


''' Start the write submission file '''
fsubfile = open(subFiles+".sub", "w")
command_lines = '''universe   = vanilla
getenv     = True
executable = sub_skim.sh
+JobFlavour           = "tomorrow"
requirements = (OpSysAndVer =?= "CentOS7")
RequestCpus = 2
'''

''' Loop over files '''
i=0
for line in Lines:
    outtempfiles = open(inFiles+"_part"+str(i)+".txt", "w")
    outtempfiles.write(line)
    outtempfiles.close()
    temp = '''
log        = cond/'''+subFiles+'''_part_'''+str(i)+'''.log
output     = cond/'''+subFiles+'''_part_'''+str(i)+'''.out
error      = cond/'''+subFiles+'''_part_'''+str(i)+'''.err
arguments = '''+inFiles+'''_part'''+str(i)+'''.txt '''+inFilesV0+'''.txt '''+outFiles+'''_'''+str(i)+'''.root ''' ''' 
queue
'''
    command_lines += temp
    i=i+1

fsubfile.write(command_lines)
fsubfile.close()
subprocess.call(["condor_submit", subFiles+".sub"])
