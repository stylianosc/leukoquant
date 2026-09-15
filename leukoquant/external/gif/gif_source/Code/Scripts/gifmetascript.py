#! /usr/bin/env python

import os, glob, sys, getopt, re, time

from gifclusterutils import *
from gifutils import *

def main(argv):
  inputdirectory = ''
  allcppsdirectory = ''
  workingdirectory = ''
  initiallabelsdirectory = ''
  clustersegmentation=0
  helpmessage = "\n\nUsage: gifmetascript.py -i <inputs> -c <cpps> -w <working_dir> -l <labels_dir> -s <clustersegmentation>\n\n"
  
  try:
    opts, args = getopt.getopt(argv,"hi:c:w:l:r:s:",["inputdirectory=","allcppsdirectory=","workingdirectory=","initiallabelsdirectory=","clustersegmentation="])
  except getopt.GetoptError:
    print(helpmessage)
    sys.exit(2)

  if len(opts) < 5:
    print(helpmessage)
    print("ERROR: not enough arguments: ", opts)
    sys.exit()

  for opt, arg in opts:
    if opt == '-h':
      print(helpmessage)
      sys.exit()
    elif opt in ("-i", "--inputdirectory"):
      inputdirectory = os.path.abspath(arg)
    elif opt in ("-c", "--allcppsdirectory"):
      allcppsdirectory = os.path.abspath(arg)
    elif opt in ("-w", "--workingdirectory"):
      workingdirectory = os.path.abspath(arg)
    elif opt in ("-l", "--initiallabelsdirectory"):
      initiallabelsdirectory = os.path.abspath(arg)
    elif opt in ("-s", "--clustersegmentation"):
      clustersegmentation = int(arg)

  gif_numberofiterations=13
  
  print('Script path: \t\t\t', os.path.abspath(__file__))
  print('Input directory: \t\t', inputdirectory)
  print('Cpps directory: \t\t', allcppsdirectory)
  print('Working directory: \t\t', workingdirectory)
  print('Initial Labels directory: \t', initiallabelsdirectory)
  print('Use cluster-segmentation: \t', clustersegmentation)
  print('Number of Iterations: \t\t', gif_numberofiterations)
  time.sleep(2)

  scriptdirectory = os.path.dirname(os.path.abspath(__file__))

  os.chdir(workingdirectory)
  
  filelist=glob.glob(inputdirectory+"/*.nii.gz")
  
  numberofimages = len(filelist)
  print("Number of images found: \t%d" % numberofimages)
  time.sleep(2)

  print("************************************************")
  print("masking inputs...")
  print("************************************************")
  print("the script will create (if inexistent) a directory in " + os.path.join(workingdirectory,'masks'))
  print("and will roughly mask each input (if not already done)")
  print("************************************************")
  time.sleep(2)
  
  masksdirectory=os.path.join(workingdirectory,'masks')
  ensuredir(masksdirectory)
  
  for file in filelist:
    
    filename=os.path.basename(file)
    maskfile=os.path.join(masksdirectory,filename)
    createmask(file,maskfile)
  
  print("************************************************")
  print("make cpp file links...")
  print("************************************************")
  print("the script will create displacementfields links in " + os.path.join(workingdirectory,'cpp-links'))
  print("************************************************")
  time.sleep(2)

  cpplinksdirectory=os.path.join(workingdirectory,'cpp-links')
  ensuredir(cpplinksdirectory)

  for file1 in filelist:
    name1=os.path.basename(file1)
    name1=re.sub(r'\.nii.gz$', '', name1)
    
    print("cpp links for reference " + name1 + "...")
    
    directory1=os.path.join(cpplinksdirectory,name1)
    ensuredir(directory1)
    
    for file2 in filelist:
      name2=os.path.basename(file2)
      name2=re.sub(r'\.nii.gz$', '', name2)
      targetfile=name1 + "_VS_" + name2 + ".nii.gz"
      targetfile=os.path.join(allcppsdirectory, targetfile)
      linkfile=name2 + ".nii.gz"
      linkfile=os.path.join(directory1, linkfile)
      createlink(targetfile, linkfile)

      
  print("************************************************")
  print("creating databases...")
  print("************************************************")
  print("the script will create a folder per GIF iteration in " + workingdirectory)
  print("************************************************")
  time.sleep(2)
  
  for gif_iteration in range(1, gif_numberofiterations+2):
    name="db"+str(gif_iteration)
    database=os.path.join(workingdirectory, name);
    createdatabase(inputdirectory, initiallabelsdirectory, database)
    

  print("************************************************")
  print("Geodesic Image Flow Algorithm...")
  print("************************************************")
  time.sleep(2)

  for gif_iteration in range(10, gif_numberofiterations+1):
    
    print("************************************************")
    print("GIF ITERATION " + str(gif_iteration) + "...")
    print("************************************************")
    time.sleep(2)
    
    databasename="db"+str(gif_iteration)
    nextdatabasename="db"+str(gif_iteration+1)
    databasepath=os.path.join(workingdirectory, databasename);
    xml=os.path.join(databasepath, "db.xml")
    nextdatabasepath=os.path.join(workingdirectory, nextdatabasename);
    outputdirectory=nextdatabasepath
    
    if clustersegmentation >= 1:
      scriptpath=os.path.join(scriptdirectory,"gifclusterjob.py")
      submit_GIF_vector(scriptpath, inputdirectory, masksdirectory, cpplinksdirectory, xml, outputdirectory)
      continue
    
    for input in filelist:
      inputfilename=os.path.basename(input)
      inputname=re.sub(r'\.nii.gz$', '', inputfilename)

      mask=os.path.join(masksdirectory, inputfilename)
      cppsdirectory=os.path.join(cpplinksdirectory, inputname)
      
      GIF(input, mask, xml, cppsdirectory, outputdirectory)


if __name__ == "__main__":
   main(sys.argv[1:])
