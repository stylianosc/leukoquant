#! /usr/bin/env python

import os, sys, glob, re

def submit_pairwiseregistration_vector(regf3dclusterjob, file, inputdirectory, allcppsdirectory):
  filelist=glob.glob(inputdirectory+"/*.nii.gz")
  numberofjobs=len(filelist)
  inputfilename=os.path.basename(file)
  inputname=re.sub(r'\.nii.gz$', '', inputfilename)
  qsub_array_cmd="qsub \
	-l h_rt=36:00:00 \
	-l tmem=7.8G \
	-l h_vmem=7.8G \
	-l vf=7.8G \
	-l s_stack=10240  \
	-j y \
	-S /bin/csh \
	-b y \
	-V \
	-o " + allcppsdirectory + " \
	-N regf3d \
	-t 1:" + str(numberofjobs) + " \
	python " + regf3dclusterjob + " " + file + " " + inputdirectory + " " + allcppsdirectory + ";"


  print("submitting jobs with command \n" + qsub_array_cmd)
  os.system(qsub_array_cmd)

def submit_GIF_vector(gifclusterjob, inputdirectory, masksdirectory, cpplinksdirectory, xml, outputdirectory):
  filelist=glob.glob(inputdirectory+"/*.nii.gz")
  numberofjobs=len(filelist)
  databasename=os.path.basename(outputdirectory)
  formerdatabasename=os.path.basename(os.path.dirname(xml))
  gifprocedure = "qstat | awk '{print $3}' | grep GIF_" + formerdatabasename + " | tr '\\n' ','"
  gif_jobs_to_wait_for=os.popen(gifprocedure).read()
  regprocedure = "qstat | awk '{print $3}' | sort | uniq | grep regf3d | tr '\\n' ','"
  reg_jobs_to_wait_for=""
  if formerdatabasename == 'db1':
    reg_jobs_to_wait_for=os.popen(regprocedure).read()
  maskprocedure = "qstat | awk '{print $3}' | sort | uniq | grep MASK | tr '\\n' ','"
  mask_jobs_to_wait_for=""
  if formerdatabasename == 'db6':
    mask_jobs_to_wait_for=os.popen(maskprocedure).read()

  qsub_array_cmd="qsub \
	-l h_rt=36:00:00 \
	-l tmem=7.8G \
	-l h_vmem=7.8G \
	-l vf=7.8G \
	-l s_stack=10240  \
	-j y \
	-S /bin/csh \
	-b y \
	-V \
	-o " + outputdirectory + " \
	-N GIF_" + databasename + " \
	-t 1:" + str(numberofjobs) + " \
	-hold_jid " + str(gif_jobs_to_wait_for) + str(reg_jobs_to_wait_for) + str(mask_jobs_to_wait_for) + " \
	python " + gifclusterjob + " " + inputdirectory + " " + masksdirectory + " " + cpplinksdirectory + " " + xml + " " + outputdirectory + ";"

  print("submitting jobs with command \n" + qsub_array_cmd)
  os.system(qsub_array_cmd)


def submit_mergemasks_vector(maskclusterjob, labelsdirectory, brainmasksdirectory, outputdirectory):
  filelist=glob.glob(brainmasksdirectory+"/*.nii.gz")
  numberofjobs=len(filelist)
  formerdatabasename="db14"
  newdatabasename="db15"
  gifprocedure = "qstat | awk '{print $3}' | grep GIF_" + formerdatabasename + " | tr '\\n' ','"
  gif_jobs_to_wait_for=os.popen(gifprocedure).read()
  
  qsub_array_cmd="qsub \
	-l h_rt=12:00:00 \
	-l tmem=3.7G \
	-l h_vmem=3.7G \
	-l vf=3.7G \
	-l s_stack=10240  \
	-j y \
	-S /bin/csh \
	-b y \
	-V \
	-o " + outputdirectory + " \
	-N MASK_" + newdatabasename + " \
	-t 1:" + str(numberofjobs) + " \
	-hold_jid " + str(gif_jobs_to_wait_for) + " \
	python " + maskclusterjob + " " + labelsdirectory + " " + brainmasksdirectory + " " + outputdirectory + ";"

  print("submitting jobs with command \n" + qsub_array_cmd)
  os.system(qsub_array_cmd)
