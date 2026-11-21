#! /usr/bin/env python

import os, glob, re, subprocess, neuromorphometricslabels
from imageprocessingutils import *


def createdatabase(inputs, labels, database):

  print("creating database at " + database)

  labelxmlfile=os.path.join(database, "labels.xml")
  databasexmlfile=os.path.join(database, "db.xml")
  labelsdir=os.path.join(database, "labels")
  priorsdir=os.path.join(database, "priors")
  tissuesdir=os.path.join(database, "tissues")
  geosdir=os.path.join(database, "labels_geo")
  databasexmlcontent='<?xml version="1.0"?> \n\
<document>\n\
<data>\n\
<fname>T1s</fname>\n\
<path>' + os.path.basename(inputs) + '</path>\n\
<descr>T1 MRI Data</descr>\n\
<sform>1</sform>\n\
</data>\n\
\n\
<info>\n\
<fname>labels</fname>\n\
<path>labels</path>\n\
<gpath>labels_geo</gpath>\n\
<descr>Segmentation</descr>\n\
<extra>labels.xml</extra>\n\
<type>LabelSeg</type>\n\
</info>\n\
\n\
</document>'
  
  labelxmlcontent=neuromorphometricslabels.getlabels()
  
  ensuredir(database)

  if not os.path.exists(databasexmlfile):
    f=open(databasexmlfile, 'w+')
    f.write(databasexmlcontent)
  if not os.path.exists(labelxmlfile):
    f=open(labelxmlfile, 'w+')
    f.write(labelxmlcontent)
  
  ensuredir(labelsdir)
  ensuredir(priorsdir)
  ensuredir(tissuesdir)
  ensuredir(geosdir)
  createlink(inputs, os.path.join(database,os.path.basename(inputs)))

  labelslist=glob.glob(labels+"/*.nii.gz")
  
  for file in labelslist:
    filename=os.path.basename(file)
    link=os.path.join(labelsdir,filename)
    createlink(file, link)




def GIF(input, mask, xml, cpps, outputdirectory):

  currentdir=os.getcwd()
  
  dbdirectory=os.path.dirname(xml)
  inputfilename=os.path.basename(input)
  inputname=re.sub(r'\.nii.gz$', '', inputfilename)

  print("Geodesic Information Flow Propagation with arguments:")
  print("target \t\t" + input)
  print("mask \t\t" + mask)
  print("xml \t\t" + xml)
  print("cpps \t\t" + cpps)
  print("outputdir \t" + outputdirectory)

  outputfilename1=inputname + "_labels_label.nii.gz"
  outputfilename2=inputname + ".nii.gz";
  outputfilename1=os.path.join(outputdirectory, outputfilename1)
  outputfilename2=os.path.join(outputdirectory, "labels", outputfilename2)
  outputfilename3=os.path.join(outputdirectory, inputname + "-labels.nii.gz")

  if os.path.exists(outputfilename1) or os.path.exists(outputfilename2) or os.path.exists(outputfilename3):
    print("output file exists, skipping")
    return None
    
  gif="seg_GIF"
  gif_args="-in " + input + " -mask " + mask + " -db " + xml + " -cpp " + cpps + " -out " + outputdirectory + " -v 2"

  gif_command=gif + " " + gif_args
  print("gif command: " + gif_command)

  os.system(gif_command)

  # renaming in a consistent way the output files...
  labelfile   = os.path.join(outputdirectory, inputname + "_labels_label.nii.gz")
  priorfile   = os.path.join(outputdirectory, inputname + "_labels_prior.nii.gz")
  tissuesfile = os.path.join(outputdirectory, inputname + "_labels_tissue.nii.gz")
  geofile     = os.path.join(outputdirectory, inputname + "_labels_geo.nii.gz")

  if not os.path.exists(labelfile):
    labelfile = os.path.join(outputdirectory, inputname + "_labels.nii.gz")

  if os.path.exists(labelfile):

    resetrange_cmd = "seg_maths " + labelfile + " -range " + labelfile
    os.system(resetrange_cmd)
    
    if os.path.exists(os.path.join(outputdirectory,"labels")):
      os.rename(labelfile, os.path.join(outputdirectory,"labels",inputfilename))
    else:
      os.rename(labelfile, os.path.join(outputdirectory,inputname + "-labels.nii.gz"))

  if os.path.exists(priorfile):
    if os.path.exists(os.path.join(outputdirectory,"priors")):
      os.rename(priorfile, os.path.join(outputdirectory,"priors",inputfilename))
    else:
      os.rename(priorfile, os.path.join(outputdirectory,inputname + "-priors.nii.gz"))

  if os.path.exists(tissuesfile):
    if os.path.exists(os.path.join(outputdirectory,"tissues")):
      os.rename(tissuesfile, os.path.join(outputdirectory,"tissues",inputfilename))
    else:
      os.rename(tissuesfile, os.path.join(outputdirectory,inputname + "-tissues.nii.gz"))

  if os.path.exists(geofile):
    if os.path.exists(os.path.join(outputdirectory,"labels_geo")):
      os.rename(geofile, os.path.join(outputdirectory,"labels_geo",inputfilename))
    else:
      os.rename(geofile, os.path.join(outputdirectory,inputname + "-geo.nii.gz"))
  

def mergelabelsandmask(label, mask, listoflabelstomerge, output):
  print("masking file " + label)

  if not os.path.exists(label) or not os.path.exists(mask):
    print("At least one of the input files do not exist, skipping")
    return None

  inputfilename=os.path.basename(label)
  inputname=re.sub(r'\.nii.gz$', '', inputfilename)
  temporarydirectory="/tmp/"

  firstlabelfile = temporarydirectory + inputname + "_" + str(listoflabelstomerge[0]) + ".nii.gz"
  addanddilatecommand="fslmaths " + firstlabelfile + " "

  outputlabelfile = output 

  if os.path.exists(outputlabelfile):
    print("file " + outputlabelfile + " exists, skipping")
    return None

  for i in range(0,len(listoflabelstomerge)):
    l = listoflabelstomerge[i]
    print("extracting label l=" + str(l))
    labelfile = temporarydirectory + inputname + "_" + str(l) + ".nii.gz"
    thresholdcommand="fslmaths " + label + " -thr " + str(l-1) + ".5 -uthr " + str(l) + ".5 " + labelfile
    print("extraction command: " + thresholdcommand)
    if not os.path.exists(labelfile):
      os.system(thresholdcommand)

    if os.path.exists(labelfile):
      addanddilatecommand = addanddilatecommand + " -add " + labelfile
  
  alllabelsfile=temporarydirectory + inputname + "_alllabels.nii.gz"
  addanddilatecommand = addanddilatecommand + " " + alllabelsfile
  print("addition command: " + addanddilatecommand)
  if not os.path.exists(alllabelsfile):
    os.system(addanddilatecommand)

  dilatedlabelfile=temporarydirectory + inputname + "_alllabels_dil.nii.gz"
  dilatecommand="seg_maths " + alllabelsfile + " -dil 1 " + dilatedlabelfile
  print("dilate command: " + dilatecommand)
  if not os.path.exists(dilatedlabelfile):
    os.system(dilatecommand)

  resampledmaskfile=temporarydirectory + inputname + "_resampledmask.nii.gz"
  resampledmaskcommand="reg_resample -ref " + label + " -flo " + mask + " -inter 0 -res " + resampledmaskfile + " -platf 0"
  print("resample mask command: " + resampledmaskcommand)
  if not os.path.exists(resampledmaskfile):
    os.system(resampledmaskcommand)

  labelsmaskedfile=temporarydirectory + inputname + "_initialmask.nii.gz"
  initialmaskcommand="fslmaths " + label + " -mas " + resampledmaskfile + " " + labelsmaskedfile
#  initialmaskcommand="fslmaths " + label + " -mas " + mask + " " + labelsmaskedfile
  print("initial mask command: " + initialmaskcommand)
  if not os.path.exists(labelsmaskedfile):
    os.system(initialmaskcommand)

  maskfile=temporarydirectory + inputname + "_mask.nii.gz"
  createmaskcommand="fslmaths " + dilatedlabelfile + " -add " + labelsmaskedfile + " " + maskfile
  print("mask creation command: " + createmaskcommand)
  if not os.path.exists(maskfile):
    os.system(createmaskcommand)

  almostlabelfile1=temporarydirectory + inputname + "_almost1.nii.gz"
  maskcommand="fslmaths " + label + " -mas " + maskfile + " " + almostlabelfile1
  print("mask command: " + maskcommand)
  if not os.path.exists(almostlabelfile1):
    os.system(maskcommand)

  almostlabelfile2=temporarydirectory + inputname + "_almost2.nii.gz"
  alsmost1command="seg_maths " + almostlabelfile1 + " -sub 1 -thr 0.5 " + almostlabelfile2
  print("alsmost1 command: " + alsmost1command)
  if not os.path.exists(almostlabelfile2):
    os.system(alsmost1command)

  alsmost2command="seg_maths " + maskfile + " -dil 6 -fill -ero 2 -add " + almostlabelfile2 + " " + outputlabelfile
  print("alsmost2 command: " + alsmost2command)
  if not os.path.exists(outputlabelfile):
    os.system(alsmost2command)
  

  print("cleaning up")
  for l in listoflabelstomerge:
    labelfile = temporarydirectory + inputname + "_" + str(l) + ".nii.gz"
    os.remove(labelfile)
  
  os.remove(alllabelsfile)
  os.remove(dilatedlabelfile)
  os.remove(maskfile)
  os.remove(almostlabelfile1)
  os.remove(almostlabelfile2)
  os.remove(labelsmaskedfile)
  os.remove(resampledmaskfile)
  
  
