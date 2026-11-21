#!/usr/bin/python

import sys
from gifutils import *

args=sys.argv

inputdirectory=args[1]
filelist=glob.glob(inputdirectory+"/*.nii.gz")

masksdirectory = args[2]
cpplinksdirectory = args[3]
xml = args[4]
outputdirectory = args[5]
taskid = int(os.environ['SGE_TASK_ID'])

input = filelist[taskid-1]

inputfilename=os.path.basename(input)
inputname=re.sub(r'\.nii.gz$', '', inputfilename)

mask=os.path.join(masksdirectory, inputfilename)

cppsdirectory=os.path.join(cpplinksdirectory, inputname)

print("************************************************")
print("gifclusterjob.py ("+ str(taskid) +"):")
print("************************************************")
print("input: \t" + input)
print("mask: \t" + mask)
print("xml: \t" + xml)
print("cpps: \t" + cppsdirectory)
print("outputdirectory: \t" + outputdirectory)

GIF(input, mask, xml, cppsdirectory, outputdirectory)
