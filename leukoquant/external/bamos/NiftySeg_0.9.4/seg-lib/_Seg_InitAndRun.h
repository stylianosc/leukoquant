//
//  Seg_InitAndRun.h
//  TreeSeg
//
//  Created by Carole Sudre on 15/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 07/05/2026.
//
#ifndef _SEG_IR_H
#define _SEG_IR_H

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include "_TreeEM_new.h"
#include "_seg_tools.h"
#include "_DirichletPriors.h"


using namespace std;

vector<nifti_image *> ReadFromFilenamesVector(vector<string> FilenamesVector);
nifti_image* ReadFromFilename(string Filename);
nifti_image * CreateDataImage(vector<nifti_image *> ImagesToSegment);
TreeEM * CreateTreeSegmentationInit(SEG_PARAMETERS * segment_param);


int FirstNotNULL(vector<nifti_image *> ImagesVector);

int LastNotNULL(vector<nifti_image *> ImagesVector);
typedef std::map<int, int>MAP;
bool IsNormalised(nifti_image * ImageToCheck);
int * ReadCountfromFile(string CountFile);
vector<int *>ReadCountfromFiles(vector<string> CountFiles);
vector<int *>ReadDistributionFromTextFile(string DistributionTextFile);
vector<int *> ReadCountFromDistributionTextFile(string DistributionTextFile);
TreeEM * ReadTreeFromTextFileWithAdapt(string TreeTextFile, SEG_PARAMETERS * segment_param, vector<nifti_image *> PriorsGCVector, vector<nifti_image*>PriorsIOVector, vector<string> AdaptTextFile);
TreeEM * ReadTreeFromTextFile(string TreeTextFile, SEG_PARAMETERS * segment_param,char * ChangePath);
TreeEM * ReadChildLevel(istream & text, TreeEM * TreeResult, SEG_PARAMETERS * segment_param);
map<int,int> ReadReclassifDecisionFromFile(string DecisionTextFile);
nifti_image * NoiseMaskOtsu(nifti_image * CreatedImage );
float * NoiseCovPriorsOtsu(nifti_image * CreatedImage);
vector<nifti_image *> SmoothAndAverageImagesForPriors(nifti_image * Summarised1, nifti_image * Summarised2, float Weight, SEG_PARAMETERS * segment_param,int Level);
nifti_image * LinearCombWeighted(vector<nifti_image *> VectorImagesComb, vector<float> VectorWeightComb);
vector<nifti_image *> Demerge4DImage(nifti_image * Image4D);
nifti_image * CreateNewPrior(vector<nifti_image*> DoubledVectorPriors,vector<int> ListToCombine);
nifti_image * CopyFloatNiiImage(nifti_image * FloatNiiToCopy);
vector<string> ReadListFilenamesFromTextFile(string TextFile);
vector<int> ReadReclassifLesionFromFile(string ClassifTextFile);
template <typename DTYPE >nifti_image * SavePartialResult(DTYPE * PartialResult, nifti_image * DataImage, string filename){
    if (DataImage==NULL) {
        cout<<"No original image given"<<endl;
        return NULL;
    }
    nifti_image * Result=nifti_copy_nim_info(DataImage);

    // Saving Partial Result means by definition that it is only on one dimension thus need to change DataImage if multispectral
    Result->dim[0]=3;
    Result->dim[4]=1;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    //cout<<"the number of voxels in Result is "<<Result->nvox;
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    // Have to check beforehand that PartialResult has really the proper number of elements
    float * FloatPartialResult=new float[Result->nvox];
    float * FloatPartialResult_PTR=FloatPartialResult;
    DTYPE * PartialResult_PTR=PartialResult;
    for (int i=0; i<Result->nvox; i++,FloatPartialResult_PTR++,PartialResult_PTR++) {
        *FloatPartialResult_PTR=(float)(*PartialResult_PTR);
    }
    float * ResultData_PTR=static_cast<float *>(Result->data);
    int nonZeroCount=0;
    FloatPartialResult_PTR=FloatPartialResult;
    //cout<<"the number of voxels is "<<Result->nvox;
    for (int i=0; i<Result->nvox; i++,FloatPartialResult_PTR++,ResultData_PTR++) {
        (*ResultData_PTR)=(float)(*FloatPartialResult_PTR);
        if (*ResultData_PTR>0) {
            nonZeroCount++;
        }
    }
    delete [] FloatPartialResult;
    //cout<<" the number of non zero values in Result is "<<nonZeroCount;
    //cout<<nifti_validfilename(filename)<<endl;
    nifti_set_filenames(Result, filename.c_str(), 0, 0);
    string fnametest=Result->fname;
    //cout<<fnametest<<endl;
    return Result;
}

#endif


