//
//  TreeEM.cpp
//  TreeSeg
//
//  Created by Carole Sudre on 08/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 07/05/2026.
//

#include "_TreeEM.h"
#include <fstream>
#include <vector>
#include <iostream>
#include "_seg_matrix.h"
#include "_Seg_InitAndRun.h"
#include <math.h>
#include "_SVDCov.h"
#include "_EigenvaluesCov.h"
//#include "Seg_EMOperations.h"
//#include "Seg_SizeAndMask.h"
//using namespace std;

TreeEM::TreeEM(){
    
    this->Parent = NULL;
    this->DataImage=NULL;
    this->Mask = NULL;
    this->Priors = NULL;
    //this->Short_2_Long_Indices = NULL;
    //this->Long_2_Short_Indices;
    this->L2S=NULL;
    this->S2L=NULL;
    
    this->Distribution=NULL;
    this->NonNormResp=NULL;
    this->NormResp=NULL;
    this->NormWeight=1;
    
    this->ParametersDistribution=NULL;
    this->SplitCheck=NULL;
    this->MergeCheck=NULL;
   // this->BFBasisFunctions=NULL;
    this->BFCorrection=NULL;
    this->BFCoeffs=NULL;

}

static int numbbins=100;
static int BForder=3;
static bool BFFlag=1;

TreeEM::TreeEM(nifti_image * DataInput){
    
    TreeEM * Parent=NULL;
    nifti_image * DataImage=DataInput;
    this->Mask=NULL;
    this->L2S=NULL;
    this->S2L=NULL;
    this->Priors=NULL;
    
}

/************* COPY METHODS ***********************/

TreeEM * TreeEM::CopyTree(TreeEM * ParentNew){
    TreeEM * CopiedTree=new TreeEM();
    
    // Copy all elements in the node according to preestablished copying methods
    CopiedTree->Parent=ParentNew;
    if (this->GetDataDirect()==NULL && ParentNew==NULL) {
        cout<<"Improper beginning of copy"<<endl;
        CopiedTree->DataImage=this->FindRoot()->CopyImage();
        CopiedTree->Mask=this->FindRoot()->CopyMask();
    }
    else{
        //cout<<"Proper beginning of copy"<<endl;
        CopiedTree->DataImage=this->CopyImage();
        CopiedTree->Mask=this->CopyMask();
    }
    CopiedTree->Priors=this->CopyPriors();
    CopiedTree->L2S=this->CopyL2S();
    CopiedTree->S2L=this->CopyS2L();
    CopiedTree->Distribution=this->CopyDistribution();
    CopiedTree->NonNormResp=this->CopyNonNormResp();
    CopiedTree->NormResp=this->CopyNormResp();
    CopiedTree->NonNormWeight=this->GetNonNormWeight();
    CopiedTree->NormWeight=this->GetNormWeight();
    CopiedTree->ParametersDistribution=this->CopyParameters();
    CopiedTree->SplitCheck=this->CopySplitCheck();
    CopiedTree->MergeCheck=this->CopyMergeCheck();
    CopiedTree->BFCorrection=this->CopyBFCorrection();
    CopiedTree->BFCoeffs=this->CopyBFCoeffs();
    
    
    // Copying for all children of this but with the proper parent pointer
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        CopiedTree->Children.push_back((this->GetChild(c))->CopyTree(CopiedTree));
    }
    //cout<<"Node copied"<<endl;
    return CopiedTree;
}

// Returns a pointer to a nifti_image that is a strict copy from the data image directly obtained in this
nifti_image * TreeEM::CopyImage(){
    if (this->GetDataDirect()==NULL) {// Meaning that it is not the root
        //cout<<"Direct Data is NULL"<<endl;
        return NULL;
    }
    
    else{ // DataImage is directly available so we have to copy it;
    if(this->GetDataImage()->datatype!=DT_FLOAT){
        this->ConvertDataImageToFloat();
    }
    //cout<<this->GetDataImage()<<" and "<<this->GetDataDirect()<<endl;
    nifti_image * CopiedImage=nifti_copy_nim_info(this->GetDataImage());
        // CopiedImage is at the moment with a NULL pointer to the data. According memory must then be allocated and filled with the values
    int numbvox=this->GetDataImage()->nvox;
    CopiedImage->data=(void *)calloc(this->GetDataImage()->nvox,sizeof(float));
    float * CopiedImageData=static_cast<float *>(CopiedImage->data);
    float * CopiedImageData_PTR=CopiedImageData;
    float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
    for (int i=0; i<numbvox; i++,Data_PTR++,CopiedImageData_PTR++) {
        *CopiedImageData_PTR=*Data_PTR;
    }
    //cout<<"Image copied"<<endl;
    return CopiedImage;
    }
}

// Returns as a nifti_image binary converted the mask used if a mask is used. Similar to CopyImage since the Mask can only be at the root also
nifti_image * TreeEM::CopyMask(){
    if (this->GetMaskDirect()==NULL) {
        return NULL;
    }
    else{
    if(this->GetMask()->datatype!=DT_BINARY){
        this->MakeMaskBinary();
    }
    nifti_image * CopiedMask=nifti_copy_nim_info(this->GetMask());
    int numbvox=this->GetMask()->nvox;
    CopiedMask->data=(void *)calloc(this->GetMask()->nvox,sizeof(bool));
    bool * CopiedMaskData=static_cast<bool *>(CopiedMask->data);
    bool * CopiedMaskData_PTR=CopiedMaskData;
    bool * Mask_PTR=static_cast<bool *>(this->GetMask()->data);
    for (int i=0; i<numbvox; i++,Mask_PTR++,CopiedMaskData_PTR++) {
        *CopiedMaskData_PTR=*Mask_PTR;
    }
    //cout<<" Mask copied "<<endl;
    return CopiedMask;
    }
}

// Returns a pointer to a nifti_image float converted that is the copy of the priors used.
nifti_image * TreeEM::CopyPriors(){
    if (this->GetPriorsDirect()==NULL) {
        return NULL;
    }
    else{
    if (this->GetPriorsDirect()->datatype!=DT_FLOAT) {
        this->MakePriorsFloat();
    }
    nifti_image * CopiedPriors=nifti_copy_nim_info(this->GetPriorsDirect());
    int numbvox=this->GetPriors()->nvox;
    CopiedPriors->data=(void *)calloc(this->GetPriors()->nvox,sizeof(float));
    float * CopiedPriorsData=static_cast<float *>(CopiedPriors->data);
    float * CopiedPriorsData_PTR=CopiedPriorsData;
    float * Priors_PTR=static_cast<float *>(this->GetPriors()->data);
    for (int i=0; i<numbvox; i++,Priors_PTR++,CopiedPriorsData_PTR++) {
        *CopiedPriorsData_PTR=*Priors_PTR;
    }
    //cout<<" Priors copied"<< endl;
    return CopiedPriors;
    }
}

// Similar to previous ones since it is only stored at the Root. If exists makes the copy of L2S and returns a int pointer
int * TreeEM::CopyL2S(){
    if (this->GetL2SDirect()==NULL) {
        return NULL;
    }
    else {
        int numel=this->GetNumberElements();
        int * CopiedL2S=new int[numel];
        int * CopiedL2S_PTR=CopiedL2S;
        int * L2S_PTR=this->GetL2S();
        for (int i=0; i<numel; i++,CopiedL2S_PTR ++,L2S_PTR++) {
            *CopiedL2S_PTR=*L2S_PTR;
        }
        return CopiedL2S;
    }
}

int * TreeEM::CopyS2L(){
    if (this->GetS2LDirect()==NULL) {
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        int * CopiedS2L=new int[numelmasked];
        int * CopiedS2L_PTR=CopiedS2L;
        int * S2L_PTR=this->GetS2L();
        for (int i=0; i<numelmasked; i++,CopiedS2L_PTR++,S2L_PTR++) {
            *CopiedS2L_PTR=*S2L_PTR;
        }
        return CopiedS2L;
    }
}

float * TreeEM::CopyDistribution(){
    if (this->GetS2L()==NULL) {
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        float * CopiedDistribution=new float[numelmasked];
        float * CopiedDistribution_PTR=CopiedDistribution;
        float * Distribution_PTR=this->GetDistribution();
        for (int i=0; i<numelmasked; i++,CopiedDistribution_PTR++,Distribution_PTR++) {
            *CopiedDistribution_PTR=*Distribution_PTR;
        }
        return CopiedDistribution;
    }
}

float * TreeEM::CopyNonNormResp(){
    if (this->GetS2L()==NULL) {
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        float * CopiedNonNormResp=new float[numelmasked];
        float * CopiedNonNormResp_PTR=CopiedNonNormResp;
        float * NonNormResp_PTR=this->GetNonNormResp();
        for (int i=0; i<numelmasked; i++,CopiedNonNormResp_PTR++,NonNormResp_PTR++) {
            *CopiedNonNormResp_PTR=*NonNormResp_PTR;
        }
        return CopiedNonNormResp;
    }
}

float * TreeEM::CopyNormResp(){
    if (this->GetS2L()==NULL) {
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        float * CopiedNormResp=new float[numelmasked];
        float * CopiedNormResp_PTR=CopiedNormResp;
        float * NormResp_PTR=this->GetNormResp();
        for (int i=0; i<numelmasked; i++,CopiedNormResp_PTR++,NormResp_PTR++) {
            *CopiedNormResp_PTR=*NormResp_PTR;
        }
        return CopiedNormResp;
    }
}

Parameters * TreeEM::CopyParameters(){
    if (this->GetParameters()==NULL) {
        return NULL;
    }
    else{
        Parameters * CopiedParameters=new Parameters;
        CopiedParameters->DistributionType=this->GetDistributionType();
        CopiedParameters->SizeParameters=this->GetSizeParameters();
        int sizeParams=CopiedParameters->SizeParameters;
        if (sizeParams==0) { // Distinction done if the size of the Parameters is more than 0 or not to allocate proper amount of memory
            CopiedParameters->ValueParameters=NULL;
        }
        else{
        CopiedParameters->ValueParameters=new float[CopiedParameters->SizeParameters];
        float * CopiedValueParameters_PTR=CopiedParameters->ValueParameters;
        
        float * ValueParameters_PTR=this->GetParametersValue();
        for (int i=0; i<sizeParams; i++,CopiedValueParameters_PTR++,ValueParameters_PTR++) {
            *CopiedValueParameters_PTR=*ValueParameters_PTR;
        }
        }
        return CopiedParameters;
    }
}

bool * TreeEM::CopyMergeCheck(){
    if (this->GetMergeCheck()==NULL) {
        return NULL;
    }
    else{
    int numbDirectLeaves=this->GetNumberDirectLeaves();
    if (numbDirectLeaves==0) {
        return NULL;
    }
    bool * MergeCheck_PTR=this->GetMergeCheck();
    bool * CopiedMergeCheck=new bool[numbDirectLeaves*numbDirectLeaves];
    for(int i=0;i<numbDirectLeaves*numbDirectLeaves;i++){
        CopiedMergeCheck[i]=false;
    }
    bool * CopiedMergeCheck_PTR=CopiedMergeCheck;
    for (int i=0; i<numbDirectLeaves*numbDirectLeaves; i++,CopiedMergeCheck_PTR++,MergeCheck_PTR++) {
        *CopiedMergeCheck_PTR=*MergeCheck_PTR;
    }
    return CopiedMergeCheck;
    }
}

bool * TreeEM::CopySplitCheck(){
    if (this->GetSplitCheck()==NULL) {
        return NULL;
    }
    else{
        int numbDirectLeaves=this->GetNumberDirectLeaves();
        if (numbDirectLeaves==0) {
            return NULL;
        }
        bool * SplitCheck_PTR=this->GetSplitCheck();
        bool * CopiedSplitCheck=new bool[numbDirectLeaves];
        for(int i=0;i<numbDirectLeaves;i++){
            CopiedSplitCheck[i]=false;
        }
        bool * CopiedSplitCheck_PTR=CopiedSplitCheck;
        for (int i=0; i<numbDirectLeaves; i++,CopiedSplitCheck_PTR++,SplitCheck_PTR++) {
            *CopiedSplitCheck_PTR=*SplitCheck_PTR;
        }
        return CopiedSplitCheck;
    }
}

// Returns a float pointer to an array where the BFBasisFunctions have been copied
float * TreeEM::CopyBFCorrection(){
    // Normally the BF are at the root and only stored there, elsewhere it is a NULL pointer
    if (this->GetBFCorrectionDirect()==NULL) {
        return NULL;
    }
    else{
        int numelmasked=this->GetNumberMaskedElements();
        int numbmodal=this->GetNumberModalities();
        //int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
        float * CopiedBFCorrection=new float[numbmodal*numelmasked];//{0};
        for (int i=0; i<numbmodal*numelmasked; i++) {
            CopiedBFCorrection[i]=0;
        }
        float * CopiedBFCorrection_PTR=CopiedBFCorrection;
        float * BFCorrection_PTR=this->GetBFCorrectionDirect();
        for (int i=0; i<numbmodal*numelmasked; i++,CopiedBFCorrection_PTR++,BFCorrection_PTR++) {
            *CopiedBFCorrection_PTR=*BFCorrection_PTR;
        }
        return CopiedBFCorrection;
    }
}

// Returns a vector of pointers to the arrays containing the coefficients corresponding to the different modalities
float * TreeEM::CopyBFCoeffs(){
    //Declaration of the element to return;
    int numbmodal=this->GetNumberModalities();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * BFCoeffsToCopy=this->GetBFCoeffsDirect();


    // 1st case where there is no available coeffs, the size of the vector to copy is 0;
    if (this->GetBFCoeffsDirect()==NULL) {
        return NULL;
    }

    else{
        float * CopiedBFCoeffs=new float[numbmodal*numbBF];
        for(int i=0;i<numbmodal*numbBF;i++){
            CopiedBFCoeffs[i]=0;
        }
        for(int l=0;l<numbBF*numbmodal;l++){
            CopiedBFCoeffs[l]=BFCoeffsToCopy[l];
        }
        return CopiedBFCoeffs;
    }
}

TreeEM::~TreeEM (){
    MakeEmpty ();
}

void TreeEM::MakeEmpty (){
    int numbchild=this->GetNumberChildren();
    for (int c = 0; c < numbchild; c++)
    {
        if(this->GetChild(c)!=NULL){
            delete this->Children [c];
            this->Children[c] = NULL;
        }
    }
    this->Children.clear();
    if (this->GetBFCoeffsDirect()!=NULL) {
        delete[] this->BFCoeffs;
        this->BFCoeffs=NULL;
    }
    if (this->GetBFCorrectionDirect()!=NULL) {
        delete[] this->BFCorrection;
        this->BFCorrection=NULL;
    }
    if (this->GetDataDirect()!=NULL) {
        nifti_image_free(this->GetDataDirect());
        this->DataImage=NULL;
    }
    if (this->GetMaskDirect()!=NULL) {
        nifti_image_free(this->GetMaskDirect());
        this->Mask=NULL;
    }
    if (this->GetL2SDirect()!=NULL) {
        delete [] this->GetL2SDirect();
        this->L2S=NULL;
    }
    if (this->GetS2LDirect()!=NULL) {
        delete [] this->GetS2LDirect();
        this->S2L=NULL;
    }
    if (this->GetPriorsDirect()!=NULL) {
        nifti_image_free(this->GetPriorsDirect());
        this->Priors=NULL;
    }
    if (this->NonNormResp!=NULL){
        delete this->NonNormResp;
    }
    if (this->NormResp!=NULL){
        delete [] this->NormResp;
        this->NormResp=NULL;
    }
    if (this->Distribution!=NULL){
        delete []this->Distribution;
        this->Distribution=NULL;
    }
    if (this->ParametersDistribution!=NULL){
        delete this->ParametersDistribution;
        this->ParametersDistribution=NULL;
    }
    if (this->GetSplitCheck()!=NULL) {
        delete [] this->GetSplitCheck();
        this->SplitCheck=NULL;
    }
    if (this->GetMergeCheck()!=NULL) {
        delete [] this->GetMergeCheck();
        this->MergeCheck=NULL;
    }

    this->Parent=NULL;
    return;
}

// Get Functions
TreeEM * TreeEM::GetParent(){
    return this->Parent;
}

vector<TreeEM *> TreeEM::GetChildren(){
    return this->Children;
}

// Returns the number of children
int TreeEM::GetNumberChildren(){
    return this->GetChildren().size();
}

// Gets the child of index c after checking that we have at least c+1 children
TreeEM * TreeEM::GetChild(int c){
    if (c>=this->GetNumberChildren()) {
        cout<<"Getting over bounds of number of children"<<endl;
        return NULL;
    }
    else {
        return this->GetChildren()[c];
    }
}

nifti_image * TreeEM::GetDataDirect(){
    return this->DataImage;
}

nifti_image * TreeEM::GetMaskDirect(){
    return this->Mask;
}

nifti_image * TreeEM::GetMask(){
    if (this->GetMaskDirect()!=NULL) {
        //cout<<"Has direct Mask reference"<<endl;
        return this->GetMaskDirect();
    }
    else if (this->GetParent()==NULL){ // Means it is a root or an initial node
        //cout<<"Root or IN"<<endl;
        return NULL;
    }
    else {
        return this->GetParent()->GetMask();
    }
}

nifti_image * TreeEM::GetPriors(){
    if (this->Priors !=NULL){
        return this->Priors;
    }
    else if(this->GetParent()==NULL){// Means it is a root or an initial node
        return NULL;
    }
    else {
     return this->GetParent()->GetPriors();
    }
}

nifti_image * TreeEM::GetPriorsDirect(){
    return this->Priors;
}

float * TreeEM::GetNonNormResp(){
    return this->NonNormResp;
}

float * TreeEM::GetNormResp(){
    return this->NormResp;
}

float TreeEM::GetNonNormWeight(){
    return this->NonNormWeight;
}

float TreeEM::GetNormWeight(){
    return this->NormWeight;
}

float TreeEM::GetPartNormWeight(){
    if (this->Parent==NULL || this->GetParent()->GetPriorsDirect()!=NULL) {
        return this->GetNormWeight();
    }
    else return this->GetNormWeight()*(this->GetParent())->GetPartNormWeight();
}

Parameters * TreeEM::GetParameters(){
    return this->ParametersDistribution;
}

int TreeEM::GetDistributionType(){
    return this->GetParameters()->DistributionType;
}

float * TreeEM::GetDistribution(){
    return this->Distribution;
}

int * TreeEM::GetL2SDirect(){
    return this->L2S;
}

int * TreeEM::GetS2LDirect(){
    return this->S2L;
}

bool * TreeEM::GetSplitCheck(){
    return this->SplitCheck;
}

bool * TreeEM::GetMergeCheck(){
    return this->MergeCheck;
}

/*Returns the number of bins used for the histogram. Has to be the same for each member. Therefore only stored at the root and put to 0 everywhere else*/
int TreeEM::GetNumbbins(){
    return numbbins;
}

// Returns the pointer to the beginning of the values for the means in the parameters values
/* WARNING : do not create a copy of the mean values in another array*/
float * TreeEM::GetMean(){
    if (this->GetNumberChildren()!=0) {
        return NULL;
    }
    else {
        return this->GetParameters()->ValueParameters;
    }
}

// Creates the float array and returns the pointer to it which will hold the mean for the considered NormResp given the values (no need for the considered tree to be a leaf)
float * TreeEM::GetMeanDirect(){
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    float * MeanResult=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        MeanResult[i]=0;
    }
    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * Data_PTR=Data;
    int * L2S_PTR=this->GetL2S();
    
    for (int m=0; m<numbmodal; m++) {
        Data_PTR=&Data[m*numel];
        L2S_PTR=this->GetL2S();
        float * NormResp_PTR=this->GetNormResp();
        float SumNormResp=0;
        for (int i=0; i<numel; i++,Data_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                MeanResult[m]+=*Data_PTR*(*NormResp_PTR);
                SumNormResp+=*NormResp_PTR;
                NormResp_PTR++;
            }
        }
        MeanResult[m]/=SumNormResp;
    }
    return MeanResult;
}

float * TreeEM::GetDiagVarianceDirect(){
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    float * MeanUsed=this->GetMeanDirect();
    float * Data=static_cast<float *>(this->GetDataImage()->data);
    float * Data_PTR=Data;
    int * L2S_PTR=this->GetL2S();
    float * NormResp_PTR=this->GetNormResp();
    float * DiagVarianceResult=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        DiagVarianceResult[i]=0;
    }
    float SumNormResp=0;
    
    for (int m=0; m<numbmodal; m++) {
        Data_PTR=&Data[m*numel];
        L2S_PTR=this->GetL2S();
        NormResp_PTR=this->GetNormResp();
        SumNormResp=0;
        for (int i=0; i<numel; i++,Data_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                DiagVarianceResult[m]+=*NormResp_PTR*powf((*Data_PTR-MeanUsed[m]), 2);
                SumNormResp+=*NormResp_PTR;
                NormResp_PTR++;
            }
        }
        DiagVarianceResult[m]/=SumNormResp;
    }
    delete [] MeanUsed;
    return DiagVarianceResult;
}

// Returns the pointer to the beginning of the parameters values concerned by the variance
/*WARNING : do not create a copy of the values pointed as belonging to the variance matrix*/
float * TreeEM::GetVariance(){
    if (this->GetNumberChildren()!=0) { /* If the considered element is not a leaf, there is no parameters since the distribution is of mixture type*/
        return NULL;
    }
    else {
        switch (this->GetDistributionType()) { /* For the moment only Gaussian distribution = default case might change in the future*/
                
            default:
                return this->GetParameters()->ValueParameters+this->GetNumberModalities();
                break;
        }
    }
}


int TreeEM::GetSizeParameters(){
    if (this->GetDistributionType()==1 && this->GetParameters()->SizeParameters==(int)((this->GetNumberModalities()+1)*this->GetNumberModalities())) {
        return this->GetParameters()->SizeParameters;
    }
    //cout<<"Not Gaussian or not appropriate size of parameters"<< this->GetDistributionType()<< endl;
    
        return this->GetParameters()->SizeParameters;
}

float * TreeEM::GetParametersValue(){
    return this->GetParameters()->ValueParameters;
}

// Returns in a vector of pointers to nifti_images the Priors of the children
vector<nifti_image *> TreeEM::GetPriorsVector(){
    vector<nifti_image *> PriorsVector;
    int numbchild=this->GetNumberChildren();
    for(int c=0;c<numbchild;c++){
        PriorsVector.push_back(this->GetChild(c)->GetPriors());
        //cout<<PriorsVector[c];
    }
    return PriorsVector;
}

nifti_image * TreeEM::GetDataImage(){
    //if (this->WhatTypeOfTreeComponent()==ROOT || this->WhatTypeOfTreeComponent()==INITIALNODE){
    if(this->GetParent()==NULL){
        //cout<< "The considered tree is a root or an initial node"<<endl;
        if( this->DataImage==NULL) {
        cout<<"This tree is not valid : no input image is given !"<<endl;
       
            return NULL;
    }
        else{
           return this->DataImage; 
        }
    }
    else{
        return this->GetParent()->GetDataImage();
    }
}

// Returns the max of modality m in data image after checking it can be found
float TreeEM::GetMaxDataModal(int m){
    float MaxResult=-1E32;
    if (this->GetDataImage()==NULL) { // no data image attributed in the tree
        cout<< "No data attributed to the tree, tree not valid"<< endl;
        return MaxResult;
    }
    else{
        // Check if modality wanted is compatible with Data
        if (m>=this->GetNumberModalities()){// modality wanted is out of bounds
            cout<< "Modality wanted is out of bounds"<<endl;
            return MaxResult;
        }
        // Everything has been checked so we look for the maximum
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        int numel=this->GetNumberElements();
        float * Data_PTRm=&Data_PTR[m*numel];

        for (int i=0; i<numel; i++,Data_PTRm++) {
            if((*Data_PTRm)>MaxResult){ // update of the maximum
                MaxResult=(*Data_PTRm);
            }
        }
        //cout << "the max is "<<MaxResult<<endl;
        return MaxResult;
    }
}

// Returns the min of modality m in data image after checking it can be found
float TreeEM::GetMinDataModal(int m){
    float MinResult=1E32;
    if (this->GetDataImage()==NULL) { // no data image attributed in the tree
        cout<< "No data attributed to the tree, tree not valid"<< endl;
        return MinResult;
    }
    else{
        // Check if modality wanted is compatible with Data
        if (m>=this->GetNumberModalities()){// modality wanted is out of bounds
            cout<< "Modality wanted is out of bounds"<<endl;
            return MinResult;
        }
        // Everything has been checked so we look for the maximum
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        int numel=this->GetNumberElements();
        float * Data_PTRm=&Data_PTR[m*numel];

        for (int i=0; i<numel; i++,Data_PTRm++) {
            if((*Data_PTRm)<MinResult){ // update of the maximum
                MinResult=(*Data_PTRm);
            }
        }
        //cout<<"the min is "<<MinResult<<endl;
        return MinResult;
    }
}

// Copy the data of modality m according to the mask in a float array and returns the pointer to the beginning
float * TreeEM::GetDataModalMasked(int m){
    if (this->GetDataImage()==NULL) { // no data image attributed in the tree
        cout<< "No data attributed to the tree, tree not valid"<< endl;
        return NULL;
    }
    else{
        // Check if modality wanted is compatible with Data
        if (m>=this->GetNumberModalities()){// modality wanted is out of bounds
            cout<< "Modality wanted is out of bounds"<<endl;
            return NULL;
        }
        int numel=this->GetNumberElements();
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        float * Data_PTRm=&Data_PTR[m*numel];
        float * DataMaskedModal=new float[numel];//{0};
        for (int i=0; i<numel; i++) {
            DataMaskedModal[i]=0;
        }
        float * DataMaskedModal_PTR=DataMaskedModal;
        int * L2S_PTR=this->GetL2S();
        for (int i=0; i<numel; i++,L2S_PTR++,Data_PTRm++,DataMaskedModal_PTR++) {
            if (*L2S_PTR>=0) {
                *DataMaskedModal_PTR=*Data_PTRm;
            }
        }
        return DataMaskedModal;
    }
}

int TreeEM::GetNumberElements(){
    if (this->GetDataImage()==NULL){
        cout<<"No image to look at"<<endl;
        return 0;
    }
    //cout<< this->GetDataImage()->nx<<endl;
    return this->GetDataImage()->nx*this->GetDataImage()->ny*this->GetDataImage()->nz;
}

int TreeEM::GetNumberModalities(){
    if (this->GetDataImage()==NULL){
        cout<<"No image to look at"<<endl;
        return 0;
    }
    return this->GetDataImage()->nu*this->GetDataImage()->nt;
}

int TreeEM::GetNumberMaskedElements(){
    int numel=this->GetDataImage()->nx*this->GetDataImage()->ny*this->GetDataImage()->nz;
    if(this->GetMask()==NULL || this->GetMask()->datatype!=DT_BINARY){
        //cout<<"Mask is NULL or not Binary "<<this->GetNumberElements()<<endl;
        if (this->GetMask()!=NULL && this->GetMask()->datatype!=DT_BINARY) {
            cout<<"Mask is not binary "<< this->GetMask()<<endl;
        }
        return this->GetNumberElements();
    }
    else{
        
        //cout<<"Mask is correct and binary"<<endl;
        int CountMaskedVoxels=0;
        bool * MaskPointer=static_cast<bool *>(this->GetMask()->data);
        for (int i=0; i<numel; i++) {
            if (MaskPointer[i]>0) {
                CountMaskedVoxels++;
            }
        }
        return CountMaskedVoxels;
    }
}

int * TreeEM::GetL2S(){
    if(this->GetMask()==this->GetMaskDirect()){
        if (this->L2S==NULL) {
            this->MakeL2S();
        }
        return this->L2S;
    }
    else{ // we are at a place where L2S should be NULL and we get it from the parent
        //cout<<"It must be a child"<<endl;
        if (this->L2S!=NULL) {
            delete [] this->L2S;
            this->L2S=NULL;
        }
        return this->GetParent()->GetL2S();
    }
}

void TreeEM::MakeL2S(){
    if(this->GetMask()==this->GetMaskDirect()){
        cout<<"we are at the basis level"<<endl;
        if (this->L2S!=NULL) {
            delete [] this->L2S;
        }
        int numel=this->GetNumberElements();
        this->L2S= new int [numel];
        if(this->GetMask()!=NULL && this->GetMask()->datatype==DT_BINARY){
            /* The compatibility between the dimensions of DataImage and MaskImage has to be checked when setting the MaskImage pointer beforehand*/
            bool * Maskptr = static_cast<bool *> (this->GetMask()->data);
            bool * Maskptrtmp = Maskptr;
            int * L2S_PTR = this->GetL2S();
            Maskptrtmp = Maskptr;
            int tempindex=0;


            for (int i=0; i<numel; i++,Maskptrtmp++,L2S_PTR++) {
                if ((*Maskptrtmp)>0) {
                    (*L2S_PTR)=tempindex;
                    tempindex++;
                    
                }
                else{
                    (*L2S_PTR)=-1;
                }
            }
        }
        else {
            cout<< "The datatype of the image used as mask is not binary or there is no mask to consider. All voxels are then considered as active" << endl;
            int * L2S_PTR = this->L2S;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++,L2S_PTR++) {
                
                (*L2S_PTR)=i;
            }
        }
    }
    else{
        if(this->L2S!=NULL){
            delete [] this->L2S;
            this->L2S=NULL;
        }
    }
}

int * TreeEM::GetS2L(){
    if (this->GetMask()==this->GetMaskDirect()) { // we are at the level where S2L must be calculated
        // if the pointer to the Mask is NULL, it means that the mask is the image itself and the full image must be considered
        if (this->S2L==NULL) {
            this->MakeS2L();
        }
        return this->S2L;

    }
        else{// S2L must be NULL and we get S2L from the parent
            if (this->S2L!=NULL) {
                delete [] S2L;
                this->S2L=NULL;
            }
            return this->GetParent()->GetS2L();
        }
}

void TreeEM::MakeS2L(){
    if (this->GetMask()==this->GetMaskDirect()) { // we are at the level where S2L must be calculated
        // if the pointer to the Mask is NULL, it means that the mask is the image itself and the full image must be considered
        if(this->GetMask()!=NULL && this->GetMask()->datatype==DT_BINARY){
            bool * Maskptr = static_cast<bool *> (this->GetMask()->data);
            bool * Maskptrtmp = Maskptr;
            int numelmasked=this->GetNumberMaskedElements();
            
            // reinitialise S2L
            if(this->S2L!=NULL){
                delete [] this->S2L;
            }
            this->S2L= new int[numelmasked];
            int * S2L_PTR = (int *)(this->S2L);
            Maskptrtmp = Maskptr;
            int tempindex=0;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++) {
                if ((*Maskptrtmp)>0) {
                    S2L_PTR[tempindex]=i;
                    tempindex++;
                }
                Maskptrtmp++;
            }
        }
        else {
            printf("the datatype is not binary, all the voxels are considered as active");
            // reinitialise S2L
            if(this->S2L!=NULL){
                delete [] this->S2L;
            }
            int numelmasked=this->GetNumberMaskedElements();
            this->S2L= new int[numelmasked];
            int * S2L_PTR=this->S2L;
            int numel=this->GetNumberElements();
            for(int i=0;i<numel;i++, S2L_PTR++){
                (*S2L_PTR)=i;
            }
        }
    }
        else{
            if(this->S2L!=NULL){
                delete [] this->S2L;
                this->S2L=NULL;
            }
        }
}

void TreeEM::GetLeaves(vector<TreeEM*> & LeavesVector){
    if (this->GetNumberChildren()==0) {
        LeavesVector.push_back(this);
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->GetLeaves(LeavesVector);
        }
    }
}

vector<TreeEM *> TreeEM::GetAllLeaves(){
    vector<TreeEM *> LeavesVector;
    if(this->GetNumberChildren()==0){
        return LeavesVector;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            if (this->GetChild(c)->GetNumberChildren()==0) {
                LeavesVector.push_back(this->GetChild(c));
            }
            else {
            vector<TreeEM *> ChildLeaveVector=this->GetChild(c)->GetAllLeaves();
            for (int i=0; i<ChildLeaveVector.size(); i++) {
                LeavesVector.push_back(ChildLeaveVector[i]);
            }
            }
        }
    }
    return LeavesVector;
}

int TreeEM::GetNumberAllLeaves(){
    int numbAllLeaves=0;
    if (this->GetNumberChildren()==0) {
        numbAllLeaves=0;
        return numbAllLeaves;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            if (this->GetChild(c)->GetNumberChildren()==0) {
                numbAllLeaves+=1;
            }
            else {
            numbAllLeaves=numbAllLeaves+this->GetChild(c)->GetNumberAllLeaves();
            }
        }
    }
    return numbAllLeaves;
}

// Return the number of free parameters in the model
int TreeEM::GetNumberFreeParameters(){
    int numbFreeParameters=0;
    int numbmodal=this->GetNumberModalities();
    if (this->GetNumberChildren()==0) {
        switch (this->GetDistributionType()) {
            case 0: // Should not be possible that a leaf is considered as a mixture but taken care of anyway
                return numbFreeParameters;
                break;
                
            default: // Gaussian case is default for leaves
                numbFreeParameters=(int)numbmodal+((numbmodal+1)*numbmodal)/2;
                break;
        }
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            numbFreeParameters=numbFreeParameters+this->GetChild(c)->GetNumberFreeParameters(); // Recursive part
        }
    }
    return numbFreeParameters;
}

//Returns the value of the loglikelihood for the considered model
float TreeEM::GetLogLikelihood(){
    int numelmasked=this->GetNumberMaskedElements();
    float LL=0;
    float * Distribution_PTR=this->GetDistribution();
    int CountNeedChangeValue=0;


    for (int i=0; i<numelmasked; i++,Distribution_PTR++) {
        if (*Distribution_PTR>=1E-6) {
            LL+=logf(*Distribution_PTR);
        }
        else {
            CountNeedChangeValue++;
            LL+=logf(1E-6);
        }
        }
    return LL;
}

vector<TreeEM *> TreeEM::GetAllDirectLeaves(){
    vector<TreeEM *> DirectLeavesVector;
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetNumberChildren()==0) {
            DirectLeavesVector.push_back(this->GetChild(c));
        }
    }
    return DirectLeavesVector;
}

int TreeEM::GetNumberDirectLeaves(){
    int numbDirectLeaves=0;
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetNumberChildren()==0) {
            numbDirectLeaves++;
        }
    }
    return numbDirectLeaves;
}

void TreeEM::GetNumberLeaves(int & numbleaves){
    
    if (this->GetNumberChildren()==0) {
        numbleaves++;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->GetNumberLeaves(numbleaves);
        }
    }
}

// Returns in a float array the shifted values according to the dimension of the shift and its direction
float * TreeEM::GetDataShifted(int dim,int dir){
    // First make sure that the input values are ok or change them accordingly
    dir=dir>0?1:-1;
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;

    // Initialisation of the needed values
    int * Dimensions=new int[3];
    Dimensions[0]=this->GetDataImage()->nx;
    Dimensions[1]=this->GetDataImage()->ny;
    Dimensions[2]=this->GetDataImage()->nz;
    int * Shift=new int[3];
    Shift[0]=1;
    Shift[1]=this->GetDataImage()->nx;
    Shift[2]=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int numbx=this->GetDataImage()->nx;
    int numby=this->GetDataImage()->ny;
    int numbz=this->GetDataImage()->nz;
    int indexShift=Shift[dim]*dir;
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    int nvox=this->GetDataImage()->nvox;
    int newIndex=0;
    int oldIndex=0;
    // Allocation of memory for shifted matrix
    float * DataShifted= new float[nvox];//{0};
    for (int i=0; i<nvox; i++) {
        DataShifted[i]=0;
    }
    float * DataPointed= static_cast<float*>(this->GetDataImage()->data);
    for (int m=0; m<numbmodal; m++) {
        for (int i=0; i<numbx; i++) {
            for (int j=0; j<numby; j++) {
                for (int k=0; k<numbz; k++) {
                    newIndex=i+j*Shift[1]+k*Shift[2]+indexShift;
                    oldIndex=i+j*Shift[1]+k*Shift[2];
                    if (newIndex>0 && newIndex<numel) {
                        DataShifted[m*numel+newIndex]=DataPointed[m*numel+oldIndex];
                    }
                }
            }
        }
    }
    delete Dimensions;
    delete Shift;
    return DataShifted;
}

float * TreeEM::GetDataCorrelation(int dim){
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    float * CorrelationResult=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        CorrelationResult[i]=0;
    }
    float * Mean=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        Mean[i]=0;
    }
    float * Std=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        Std[i]=0;
    }
    float * MeanShifted=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        MeanShifted[i]=0;
    }
    float * StdShifted=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        StdShifted[i]=0;
    }
    dim=dim<0?0:dim;
    dim=dim>2?2:dim;
    int numelmasked=this->GetNumberMaskedElements();
    
    // Handling the + part same afterwards for the - part
    float * DataShifted=this->GetDataShifted(dim, 1);
    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * DataShifted_PTR=DataShifted;
    float * Data_PTR=Data;
    int* L2S_PTR=this->GetL2S();
    
    // First calculate the mean and the standard deviation for shifted and not shifted
    for (int m=0; m<numbmodal; m++) {
        DataShifted_PTR=&DataShifted[m*numel];
        Data_PTR=&Data[m*numel];
        L2S_PTR=this->GetL2S();
        float tmpMeanShifted=0;
        float tmpMean=0;
        float tmpStd=0;
        float tmpStdShifted=0;
        for (int i=0; i<numel; i++,Data_PTR++,DataShifted_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpMeanShifted+=*DataShifted_PTR;
                tmpMean+=*Data_PTR;
                tmpStdShifted+=*DataShifted_PTR*(*DataShifted_PTR);
                tmpStd+=powf(*Data_PTR, 2);
            }
        }
        MeanShifted[m]=tmpMeanShifted/numelmasked;
        Mean[m]=tmpMean/numelmasked;
        Std[m]=powf(tmpStd/numelmasked-Mean[m]*Mean[m],0.5);
        StdShifted[m]=powf(tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m], 0.5);
    }
    
    // Then calculate the correlation

    for (int m=0; m<numbmodal; m++) {
        L2S_PTR=this->GetL2S();
        DataShifted_PTR=&DataShifted[m*numel];
        Data_PTR=&Data[m*numel];
        float tmpCorrelation=0;
        for (int i=0; i<numel; i++,Data_PTR++,DataShifted_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpCorrelation+=((*Data_PTR)-Mean[m])*((*DataShifted_PTR)-MeanShifted[m]);
            }
        }
        CorrelationResult[m]=tmpCorrelation/(StdShifted[m]*Std[m]);
    }
    
    delete [] DataShifted;
    DataShifted=NULL;
    DataShifted_PTR=NULL;
    
    float * DataShiftedMoins=this->GetDataShifted(dim, -1);
    L2S_PTR=this->GetL2S();
    nifti_image * ShiftMoins=SavePartialResult(DataShiftedMoins, this->GetDataImage(), (g_scratch_dir + "/TestShiftMoins").c_str());
    nifti_image_write(ShiftMoins);
    nifti_image_free(ShiftMoins);
    ShiftMoins=NULL;
    float * DataShiftedMoins_PTR=DataShiftedMoins;
    
    // Then calculate the mean and the standard deviation for shifted in other direction
    for (int m=0; m<numbmodal; m++) {
        DataShifted_PTR=&DataShifted[m*numel];
        L2S_PTR=this->GetL2S();
        float tmpMeanShifted=0;
        float tmpStdShifted=0;
        for (int i=0; i<numel; i++,DataShiftedMoins_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpMeanShifted+=*DataShiftedMoins_PTR;
                tmpStdShifted+=*DataShiftedMoins_PTR*(*DataShiftedMoins_PTR);
            }
        }
        MeanShifted[m]=tmpMeanShifted/numelmasked;
        StdShifted[m]=powf(tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m], 0.5);
    }
    
    // Then calculate the correlation
    for (int m=0; m<numbmodal; m++) {
        DataShiftedMoins_PTR=&DataShiftedMoins[m*numel];
        Data_PTR=&Data[m*numel];
        L2S_PTR=this->GetL2S();
        float tmpCorrelation=0;
        for (int i=0; i<numel; i++,Data_PTR++,DataShiftedMoins_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpCorrelation+=((*Data_PTR)-Mean[m])*((*DataShiftedMoins_PTR)-MeanShifted[m]);
            }

        }
        CorrelationResult[m]=(float)(CorrelationResult[m]+tmpCorrelation/(StdShifted[m]*Std[m]))/(2*(numelmasked-1));
    }
    delete [] DataShiftedMoins;
    DataShiftedMoins=NULL;
    DataShiftedMoins_PTR=NULL;
    delete [] Mean;
    Mean=NULL;
    delete [] Std;
    Std=NULL;
    delete [] MeanShifted;
    MeanShifted=NULL;
    delete [] StdShifted;
    StdShifted=NULL;
    return CorrelationResult;
}

// Returns the independence factor for dimension dim;
float TreeEM::GetIndFactor(int dim){
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;
    int numbmodal=this->GetNumberModalities();
    float * CorrelationDim=this->GetDataCorrelation(dim);
    float CorrMean=0;
    for (int m=0; m<numbmodal; m++) {
        CorrMean+=CorrelationDim[m];
    }
    CorrMean/=numbmodal;
    float IndFactorDim=0.9394/powf(-2*logf(2)/logf(CorrMean), 0.5);
    delete [] CorrelationDim;
    CorrelationDim=NULL;
    return IndFactorDim;
}

float TreeEM::GetIndFactorTot(){
    return this->GetIndFactor(0)*this->GetIndFactor(1)*GetIndFactor(2);
}

// Returns a float array of size the number of elements in the image (in order to save it as a nifti image afterwards) with the needed values
float * TreeEM::GetPartialResult(PartialResultType ResultType){
    if (!this->IsTreeValid()) {
        cout<<"Tree not valid : will not get any partial result"<<endl;
        return NULL;
    }
    if (ResultType>5 || ResultType<0) {
        cout<<"Not proper result type wanted"<<endl;
        return NULL;
    }
    int numel=this->GetNumberElements();

        float * PartialResult=new float[numel];//{0};
        for (int i=0; i<numel; i++) {
                    PartialResult[i]=0;
                }
    float * PartialResult_PTR=PartialResult;
    float * tmp_PTR;
    int * L2S_PTR=this->GetL2S();

    
    // Get the pointer to the data to put in the PartialResult array;
    switch (ResultType) {
        case DISTRIBUTION:
            tmp_PTR=this->GetDistribution();
            break;
        case NONNORMRESP:
            tmp_PTR=this->GetNonNormResp();
            break;
        case NORMRESP:
            tmp_PTR=this->GetNormResp();
            break;
        default:
            tmp_PTR=static_cast<float *>(this->GetDataImage()->data);
            break;
    }
    
    // Copying the corresponding data at the correct position in the PartialResult array
    for (int i=0; i<numel; i++,PartialResult_PTR++,L2S_PTR++) {
        if (*L2S_PTR>=0) {
            *PartialResult_PTR=*tmp_PTR;
            tmp_PTR++;
        }
    }
    return PartialResult;
}

float * TreeEM::GetBFCorrectionDirect(){
    return this->BFCorrection;
}

float * TreeEM::GetBFCoeffsDirect(){
    return this->BFCoeffs;
}

float * TreeEM::GetBFCorrection(){
    if (this->WhatTypeOfTreeComponent()!=ROOT && this->GetBFCorrectionDirect()==NULL) {
        return this->GetParent()->GetBFCorrection();
    }
    else{
        // if not yet calculated when needed, then create the right Basis Functions
        if (BFFlag && this->BFCorrection==NULL) {
            this->BFCorrection=this->MakeBFCorrection();
        }
        return this->BFCorrection;
    }
}

float * TreeEM::GetBFCoeffs(){
    if (this->WhatTypeOfTreeComponent()!=ROOT && !this->AreBFCoeffsDirectUseable()) {
        return this->GetParent()->GetBFCoeffs();
    }
    else {
//        if (BFFlag && !this->AreBFCoeffsDirectUseable()) {
////            if (this->GetBFCorrection()==NULL) {
////                this->BFCorrection=this->MakeBFCorrection();
////            }
//            this->BFCoeffs=this->MakeFinalBFCoeffsChildren();
//        }
        return this->BFCoeffs;
    }
}

/************** SET FUNCTIONS WITH NEEDED CHECKS *********************/

void TreeEM::SetParent(TreeEM *ParentInput){
    // First check if the ParentInput is a valid tree (Method to write, not done so far)
    // then change the pointer of parent
    this->Parent=ParentInput;
    // Reinitialise memory allocation and set everything back to zero
    this->ReinitialiseTree();
}

void TreeEM::ReinitialiseTree(){
    this->SetMask(this->GetMask());
    for (int c=0; c<this->GetNumberChildren(); c++) {
        this->GetChild(c)->ReinitialiseTree();
    }
}

void TreeEM::SetChildren(vector<TreeEM *> ChildrenInput){
    // First check if everything is right with each of the offered tree in the vector
    for (int c=0; c<ChildrenInput.size(); c++) { /*if valid, each in ChildrenInput is added to the Children vector of this 
                                                  WARNING : here the previous possible existing children are not suppressed see */
        if (ChildrenInput[c]->IsTreeValid() && ChildrenInput[c]->GetParent()==this) {
            this->Children.push_back(ChildrenInput[c]);
        }
    }
    if (ChildrenInput.size()==0) { // if there is no children as input, be sure that this is still in mode simple distribution
        if (this->GetNumberChildren()==0) {
            this->CreateAllocateAndInitializeParameters(1);
        }
    }
    else{
        this->MakeParametersMixture();
    }
}

// Add Child if valid to current vector of Children
void TreeEM::AddChild(TreeEM * ChildAdded){
    if (ChildAdded->GetParent()==this) { // check that the Child added is set for the right parent
        this->Children.push_back(ChildAdded);
        if (!this->IsNumbChildOKWithDistType()) {
            this->MakeParametersMixture();
        }
    }
    else{
        cout<<"The added child does not correspond to the right parent"<<endl;
    }
}

void TreeEM::SetData(nifti_image * DataInput){
    // Check if the considered Tree on which applied is of type Root otherwise data cannot be changed nor modified
    if (this->WhatTypeOfTreeComponent()!=ROOT&&this->WhatTypeOfTreeComponent()!=INITIALNODE) {
        cout<< "Data image can only be modified at root level and pointer is only stored there"<<endl;
        return;
    }
    else{ // the data can be modified since we are at Root level;
        this->DataImage=DataInput;
        if(!IsDataFloat()){
            this->ConvertDataImageToFloat();
        }
        cout<< "Re give data input"<<endl;
        if (!this->IsDataImageNormalised()) {
            cout<< "We have to normalise Image"<<endl;
            this->NormaliseDataImage();
        }
        cout<< "Image is already normalised" << endl;
        // if there was previously a mask, check if we can keep it
        if(this->GetMask()!=NULL){ // check on the dimensions
            if (this->GetMask()->nx!=DataInput->nx || this->GetMask()->ny!=DataInput->ny || this->GetMask()->nz!=DataInput->nz){ // Mask reput to NULL as default and so everything is reallocated as needed
                this->SetMask(NULL);
                this->ReinitialiseTree();
            }
            else { // in case of existing Mask, need to reinitialise the tree : needed mostly for the parameters structure in case the number of modalities has changed. Allows for keeping a similar structure but reinitialising and reallocating memory
                this->ReinitialiseTree();
            }
        }
    }
}

void TreeEM::NormaliseDataImage(){
    if (this->GetDataImage()==NULL){
        cout<<"WARNING : there is no data to segment !"<<endl;
        return;
    }
    if (!this->IsDataFloat()) {
        this->ConvertDataImageToFloat();
    }
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    cout<<"Data type image"<<this->GetDataImage()->datatype<< endl;
    float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
    float * Data_PTRtmp=Data_PTR;
    for (int m=0; m<numbmodal; m++) { 
        Data_PTRtmp=&Data_PTR[m*numel];
        if (this->GetMaxDataModal(m)==this->GetMinDataModal(m)) { // Case where all the values are equal in the image
            for (int i=0; i<numel; i++,Data_PTRtmp++) {
                (*Data_PTRtmp)=0;
            }
        }
        else{
            float MultiplicativeFactor=(exp(1.0f)-1.0f)/(this->GetMaxDataModal(m)-this->GetMinDataModal(m));
            float AdditiveFactor=1.0f-MultiplicativeFactor*this->GetMinDataModal(m);
            //cout<< MultiplicativeFactor << "and " <<AdditiveFactor << endl;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++,Data_PTRtmp++) {
                (*Data_PTRtmp)=MultiplicativeFactor*(*Data_PTRtmp)+AdditiveFactor;
                (*Data_PTRtmp)=logf((*Data_PTRtmp)); // Normally then between 0 and 1;
            }
        }
    }
    //cout<<this->GetMinDataModal(0)<<"and "<<this->GetMaxDataModal(0);
    return;
}

void TreeEM::ConvertDataImageToFloat(){
    if (this->GetDataImage()==NULL) {
        cout<<"WARNING : there is no data to convert !"<<endl;
        return;
    }
    if (this->IsDataFloat()) {
        cout<<"Data already of type float no need to change"<<endl;
        return;
    }
    
    
    cout<<"Image has to be converted"<<endl;
    // First changing the values to have only 0 and 1
    float * Data_PTR=static_cast<float *>(this->GetDataImage()->data); // consider that then data is float
    int numbvox=this->GetDataImage()->nvox;
    
    // Then changing datatype
    this->GetDataImage()->datatype=DT_FLOAT;
    
    // the initial array is saved and freeed
    float *initialValue = new float[this->GetDataImage()->nvox];
    float * initialValue_PTR=initialValue;
    for (int i=0; i<numbvox; i++,initialValue_PTR++,Data_PTR++) {
        *initialValue_PTR=(float)(*Data_PTR);
    }

    float *initialValue_tmp=initialValue;
    //cout<<initialValue<<endl;
    // the new array is allocated and then filled
    
    free(this->GetDataImage()->data);
    this->GetDataImage()->nbyper = sizeof(float);
    this->GetDataImage()->data = (void *)calloc(this->GetDataImage()->nvox,sizeof(float));
    float *dataPtr = static_cast<float *>(this->GetDataImage()->data);
    for (int i=0; i<numbvox; i++, dataPtr++,initialValue_tmp++) {
        (*dataPtr)=(float)(*initialValue_tmp);
    }
    delete [] initialValue;
    //cout<<this->GetDataImage()->datatype<<endl;
}



// After checking set the parameters to ParametersSet
void TreeEM::SetParameters(Parameters * ParametersToSet){
    
    // Check for compatibility between the Distribution Type and the Parameters and between the image and the Parameters
    if(!CheckForValidityOfParametersStructure(ParametersToSet)||!CheckForSizeParametersValidity( ParametersToSet)){
        cout<< "The structure of the parameters to set is not valid or the sizes are not compatible therefore creation and allocation of valid ones"<<endl;
        if(ParametersToSet !=NULL){
        delete ParametersToSet;
            ParametersToSet=NULL;
        }
        int TypeTreeComponent=this->WhatTypeOfTreeComponent();
        ParametersToSet=new Parameters();
        switch (TypeTreeComponent) {
            case 0: // Case of a root
                        this->CreateAllocateAndInitializeParameters(0);
                break;
            case 1:// Case of a branch
                        this->CreateAllocateAndInitializeParameters(0);
                break;
            default: // Case of a leaf;
                        this->CreateAllocateAndInitializeParameters( 1);
                break;
        }
    }
    else{
    // All checks are ok : Delete previous one and set new one instead;
        if(this->ParametersDistribution!=NULL){
            delete this->ParametersDistribution;
        }
    this->ParametersDistribution=new Parameters;
        this->ParametersDistribution->DistributionType=ParametersToSet->DistributionType;
        this->ParametersDistribution->SizeParameters=ParametersToSet->SizeParameters;
        this->ParametersDistribution->ValueParameters=new float [this->ParametersDistribution->SizeParameters];//{0};
        for (int i=0; i<this->ParametersDistribution->SizeParameters; i++) {
                    this->ParametersDistribution->ValueParameters[i]=0;
                }
        int SP=this->ParametersDistribution->SizeParameters;
        for (int sp=0; sp<SP; sp++) {
            this->ParametersDistribution->ValueParameters[sp]=ParametersToSet->ValueParameters[sp];
        }
    }
    return;
}

// Transforms the Parameters into a mixture (needed when creating child for a leaf becoming a branch)
void TreeEM::MakeParametersMixture(){
    if (CheckForValidityOfParametersStructure()&& this->GetDistributionType()==0) {
        return;
    }
    if(this->ParametersDistribution!=NULL){
        if (this->GetParametersValue()!=NULL) {
            delete [] this->GetParametersValue();
            this->ParametersDistribution->ValueParameters=NULL;
        }
        delete this->ParametersDistribution;
        this->ParametersDistribution=NULL;
    }
    this->CreateAllocateAndInitializeParameters(0);
}

// Set Mask to MaskImage after checking. For the moment : Mask can only be introduced at the root level and nowhere else
void TreeEM::SetMask(nifti_image *MaskImage){
    // First check if trying to set a mask at root or somewhere else
    if (this->GetParent()==NULL){
        cout<<"It is a root or an initial node"<<endl;
        // Then check for the dimension validity
        if (MaskImage!=NULL && (MaskImage->nx!=this->GetDataImage()->nx || MaskImage->ny!=this->GetDataImage()->ny ||MaskImage->nz!=this->GetDataImage()->nz) ) {
            cout<< "The dimensions between Mask and Image are not compatible"<<endl;
            return;
        }
       
        // Set the pointer to Mask to MaskImage
        this->Mask=MaskImage;
        if (MaskImage!=NULL && this->GetMask()->datatype!=DT_BINARY) {
            this->MakeMaskBinary();
            cout<<"Binarisation performed "<<endl;
        }
        // Calculate S2L and L2S
        if (this->L2S!=NULL) {
            delete [] this->L2S;
            this->L2S=NULL;
        }
        this->MakeL2S();
        if(this->S2L!=NULL){
            delete [] this->S2L;
            this->S2L=NULL;
        }
        this->MakeS2L();
        int numelmasked=this->GetNumberMaskedElements();
        //cout<< "the number of masked elements is "<<this->GetNumberMaskedElements()<<endl;
        // Allocate right amount of memory to Distribution, NonNormResp, NormResp
        if(this->GetDistribution()!=NULL){
            delete this->Distribution;
        }
        this->Distribution=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->Distribution[i]=0;
                }
        if(this->NonNormResp!=NULL){
            delete this->NonNormResp;
        }
        this->NonNormResp=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->NonNormResp[i]=0;
                }
        if(this->NormResp!=NULL){
            delete this->NormResp;
        }
        this->NormResp=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->NormResp[i]=0;
                }
        // Propagate Mask to children
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->SetMask(MaskImage);
        }
        return;
    }
    else {
        //cout<<"Parent Not NULL and address Mask is"<< this->GetMask()<<endl;
        int numelmasked=this->GetNumberMaskedElements();
        //cout<< "the number of masked elements is "<<this->GetNumberMaskedElements()<<endl;
        // Allocate right amount of memory to Distribution, NonNormResp, NormResp
        if(this->GetDistribution()!=NULL){
            delete this->Distribution;
        }
        this->Distribution=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->Distribution[i]=0;
                }
        if(this->NonNormResp!=NULL){
            delete this->NonNormResp;
        }
        this->NonNormResp=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->NonNormResp[i]=0;
                }
        if(this->NormResp!=NULL){
            delete this->NormResp;
        }
        this->NormResp=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->NormResp[i]=0;
                }
        
        // Ensuring that pointer to L2S and S2L are NULL if we are not at the node where the pointer to the Mask is stored
        if (this->L2S!=NULL) {
            delete [] this->L2S;
            this->L2S=NULL;
        }
        
        if(this->S2L!=NULL){
            delete [] this->S2L;
            this->S2L=NULL;
        }

        // Propagate memory allocation to children
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->SetMask(MaskImage);
        }
    }
    return;
}

void TreeEM::MakeMaskBinary(){
    if (this->GetMask()==NULL){
        cout<< "No mask so nothing to binarise"<<endl;
        return;
    }
    if(this->GetMask()->datatype==DT_BINARY){
        cout<< "Already binarised"<<endl;
        return;
    }
    else{
        cout<<"Mask has to be binarised"<<endl;
        // First changing the values to have only 0 and 1
        float * Mask_PTR=static_cast<float *>(this->GetMask()->data); // consider that then data is float
        float * Mask_PTRtmp=Mask_PTR;
        int numel=this->GetNumberElements();
        // if set as mask then it has been checked before that dimensions with dataimage are compatible
        bool *initialValue = new bool[this->GetMask()->nvox];
        bool *initialValue_PTR=initialValue;
        for (int i=0; i<numel; i++,Mask_PTRtmp++,initialValue_PTR++) {
            (*Mask_PTRtmp)=(*Mask_PTRtmp)>0?1:0;// Set to 1 if strictly positive, 0 otherwise
            (*initialValue_PTR)=(bool)(*Mask_PTRtmp);

        }
//        string FilenameTest=(g_scratch_dir + "/TryMask").c_str();
//        nifti_image * TestMask =SavePartialResult(Mask_PTR, this->GetDataImage(),(char *)FilenameTest.c_str() );
//        nifti_image_write(TestMask);
//        nifti_image_free(TestMask);
        // Then changing datatype
        this->GetMask()->datatype=DT_BINARY;

            // the initial array is saved and freeed

        initialValue_PTR=initialValue;
            free(this->GetMask()->data);
            this->GetMask()->nbyper = sizeof(bool);
            this->GetMask()->data = (void *)calloc(this->GetMask()->nvox,sizeof(bool));
            bool *dataPtr = static_cast<bool *>(this->GetMask()->data);
        for (int i=0; i<numel; i++, dataPtr++,initialValue_PTR++) {
            (*dataPtr)=(bool)(*initialValue_PTR);
        }
        delete [] initialValue;
        cout<<this->GetMask()->datatype<<endl;
        }
//    bool * test=static_cast<bool*>(this->GetMask()->data);
//    string FilenameTest=(g_scratch_dir + "/TryMask").c_str();
//    nifti_image * TestMask =SavePartialResult(test, this->GetDataImage(),(char *)FilenameTest.c_str() );
//    nifti_image_write(TestMask);
//    nifti_image_free(TestMask);
    
    }

// Setting the Priors to PriorsInput first checking the validity of the dimension and datatype, deleting any previous set of the priors and making sure then they are of probability type
void TreeEM::SetPriors(nifti_image *PriorsInput){
    if (PriorsInput == NULL) {
        this->Priors=NULL;
        return;
    }
    else{
    bool ValidityPriors=(PriorsInput->nx==this->GetDataImage()->nx)&&(PriorsInput->ny==this->GetDataImage()->ny)&&(PriorsInput->nz==this->GetDataImage()->nz);
    if(ValidityPriors){
        this->Priors=PriorsInput;
        if((PriorsInput->datatype!=DT_FLOAT32)){
            this->MakePriorsFloat();
        }
        this->MakePriorsProbabilityType();
    }
        return;
    }
}

// Change the priors referred so that we only handle float datatype
void TreeEM::MakePriorsFloat(){
    if (this->GetPriors()==NULL) {
        cout<<"WARNING : there is no data to convert !"<<endl;
        return;
    }
    if (this->IsPriorsFloat()) {
        cout<<"Priors already of type float no need to change"<<endl;
        return;
    }
    
    
    cout<<"Priors has to be converted"<<endl;
    // First changing the values to have only 0 and 1
    float * Priors_PTR=static_cast<float *>(this->GetPriors()->data); // consider that then data is float
    int numbvox=this->GetPriors()->nvox;
    
    // Then changing datatype
    this->GetPriors()->datatype=DT_FLOAT;
    
    // the initial array is saved and freeed
    float *initialValue = new float[numbvox];
    float * initialValue_PTR=initialValue;
    for (int i=0; i<numbvox; i++,initialValue_PTR++,Priors_PTR++) {
        *initialValue_PTR=(float)(*Priors_PTR);
    }
    
    float *initialValue_tmp=initialValue;
    //cout<<initialValue<<endl;
    // the new array is allocated and then filled
    
    free(this->GetPriors()->data);
    this->GetPriors()->nbyper = sizeof(float);
    this->GetPriors()->data = (void *)calloc(this->GetPriors()->nvox,sizeof(float));
    float *priorsPtr = static_cast<float *>(this->GetPriors()->data);
    for (int i=0; i<numbvox; i++, priorsPtr++,initialValue_tmp++) {
        (*priorsPtr)=(float)(*initialValue_tmp);
    }
    delete [] initialValue;
    //cout<<this->GetPriors()->datatype<<endl;
}




// Replaces the values in Priors so that they are between 0 and 1 and can be used as probabilities by rescaling the dispersion of the values
void TreeEM::MakePriorsProbabilityType(){
    if (this->GetPriors()==NULL) {
        return;
    }
    float * Priors_PTR=static_cast<float *>(this->GetPriors()->data);
    float * Priors_PTRtmp=Priors_PTR;
    // Initialisation of the maximal and minimal values
    float maxPriors=-1E32;
    float minPriors=1E32;
    int numel=this->GetNumberElements();
    // Determination of the value of the maximum and the minimum
    for (int i=0; i<numel; i++,Priors_PTRtmp++) {
        if ((*Priors_PTRtmp)>maxPriors) { // update of the maximum value
            maxPriors=(*Priors_PTRtmp);
        }
        if((*Priors_PTRtmp)<minPriors){ // update of the minimum value
            minPriors=(*Priors_PTRtmp);
        }
    }
    
    // reinitialisation of the pointer to the beginning of the Priors
    Priors_PTRtmp=Priors_PTR;
    // modification of the values
    if (minPriors==maxPriors) { // case (strange) where all the values are equal, then we put them all to 1
        for (int i=0; i<numel; i++,Priors_PTRtmp++) {
            (*Priors_PTRtmp)=1;
        }
    }
    else {// rescaling of the values
        for (int i=0;i<numel;i++,Priors_PTRtmp++){ // Putting all the values between 0 and 1
            (*Priors_PTRtmp)=(*Priors_PTRtmp-minPriors)/(maxPriors-minPriors);
        }
    }
}

// Copies the content of NonNormRespInput into the NonNormResp of the tree on which it is applied
void TreeEM::SetNonNormResp(float *NonNormRespInput){
    int numelmasked=this->GetNumberMaskedElements();
    if (NonNormRespInput==NULL) {
        if (this->GetNonNormResp()!=NULL) {
            delete []this->GetNonNormResp();
        }
        this->NonNormResp=NULL;
        return;
    }
    if (this->GetNonNormResp()==NULL) {
        this->NonNormResp=new float[this->GetNumberMaskedElements()];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->NonNormResp[i]=0;
                }
    }
    float * NonNormResp_PTR=this->GetNonNormResp();
    float * NonNormRespInput_PTR=NonNormRespInput;
    for (int i=0; i<numelmasked; i++,NonNormResp_PTR++,NonNormRespInput_PTR++) {
        *NonNormResp_PTR=*NonNormRespInput_PTR;
    }
    return;
}

// Copies the content of NormRespInput into NormResp of the tree on which it is applied
void TreeEM::SetNormResp(float * NormRespInput){
    int numelmasked=this->GetNumberMaskedElements();
    if (NormRespInput==NULL) {
        if (this->GetNormResp()!=NULL) {
            delete []this->GetNormResp();
        }
        this->NormResp=NULL;
        return;
    }
    if (this->GetNormResp()==NULL) {
        this->NormResp=new float[this->GetNumberMaskedElements()];//{0};
        for (int i=0; i<numelmasked; i++) {
                    this->NormResp[i]=0;
                }
    }
    float * NormResp_PTR=this->GetNormResp();
    float * NormRespInput_PTR=NormRespInput;
    
    for (int i=0; i<numelmasked; i++,NormResp_PTR++,NormRespInput_PTR++) {
        *NormResp_PTR=*NormRespInput_PTR;
    }
    return;
}

// At this specific level of the tree reinitialisation of the NormResp with the values of the priors after making them usable in a probabilistic manner
void TreeEM::InitialiseNormRespWithPriors(){
    float * NormResp_PTR=this->GetNormResp();
    int numelmasked=this->GetNumberMaskedElements();
    if (NormResp_PTR==NULL) {
        float * NormRespToSet=new float[numelmasked];
        for (int i=0; i<numelmasked; i++) {
                    NormRespToSet[i]=0;
                }
        this->SetNormResp( NormRespToSet);//{0});
    }
    if(!this->IsPriorsProbabilityType()){ // Conversion into probabilities if needed
        this->MakePriorsProbabilityType();
    }
    float * PriorsData_PTR=static_cast<float*>(this->GetPriors()->data);
    int * L2S_PTR=this->GetL2S();
    int numel=this->GetNumberElements();
    for (int i=0; i<numel; i++,PriorsData_PTR++,L2S_PTR++) {
        if ((*L2S_PTR)>=0) {
            (*NormResp_PTR)=(*PriorsData_PTR);
            NormResp_PTR++;
        }
    }
}

// Sets NonNormWeight of the tree on which it is applied to NonNormWeightInput
void TreeEM::SetNonNormWeight(float NonNormWeightInput){
        this->NonNormWeight=NonNormWeightInput;
}

// Sets NormWeight of the tree on which it is applied to NormWeightInput after checking that the value is between 0 and 1
void TreeEM::SetNormWeight(float NormWeightInput){
    if (NormWeightInput>=0&&NormWeightInput<=1) {
        this->NormWeight=NormWeightInput;
    }
    else{
        cout<<"the input value is not valid"<<endl;
    }
}

void TreeEM::SetDistribution(float * DistributionInput){
    int numelmasked=this->GetNumberMaskedElements();
    if (DistributionInput==NULL) {
        if (this->GetDistribution()!=NULL) {
            delete []this->GetDistribution();
        }
        this->Distribution=NULL;
        return;
    }
    if (this->GetDistribution()==NULL) {
        this->Distribution=new float[this->GetNumberMaskedElements()];//{0};
        for(int i=0;i<numelmasked;i++){
            this->Distribution[i]=0;
        }
    }
    float * Distribution_PTR=this->GetDistribution();
    float * DistributionInput_PTR=DistributionInput;

    //cout<<"Addresses distribution"<< Distribution_PTR << " "<< DistributionInput_PTR;
    for (int i=0; i<numelmasked; i++,Distribution_PTR++,DistributionInput_PTR++) {
        *Distribution_PTR=*DistributionInput_PTR;
    }
    return;
    
}

// After checking the validity of such an update change the BF coeffs
void TreeEM::SetBFCoeffs(float * BFCoeffsInput){
    int numbmodal=this->GetNumberModalities();

    // Clearing current space for BFCoeffs if allocated
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * BFCoeffsToChange=this->GetBFCoeffs();
    if(this->GetParent()!=NULL){ // if it is not a root we do not change the BF coeffs
        cout<<"BF coeffs only set at the root"<<endl;
        return;
    }

//    if (BFCoeffsToChange!=NULL) {
//        delete [] BFCoeffsToChange;
//        BFCoeffsToChange=NULL;
//    }
    // Fill with BF CoeffsInput
    if(BFCoeffsInput!=NULL) { // it is only in the case we are at the root (checked before)
        if(BFCoeffsToChange==NULL){
            cout<<"Memory allocated for BF"<<endl;
        this->BFCoeffs=new float[numbmodal*numbBF];
        }
        BFCoeffsToChange=this->BFCoeffs;
        for (int l=0; l<numbmodal*numbBF;l++) {
            BFCoeffsToChange[l]=BFCoeffsInput[l];
        }
    }

}

void TreeEM::SetBFCorrection(float * BFCorrectionInput){
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
//    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    // Clearing current space for BFCoeffs if allocated

    float * BFCorrectionToChange=this->BFCorrection;
    if(this->GetParent()!=NULL){ // if it is not a root we do not change the BF coeffs
        cout<<"BF correction only set at the root"<<endl;
        return;
    }

//    if (BFCorrectionToChange!=NULL) {
//        delete [] BFCorrectionToChange;
//        BFCorrectionToChange=NULL;
//    }
    // Fill with BF CoeffsInput
    if(BFCorrectionInput!=NULL) {
        if(BFCorrectionToChange==NULL){
        BFCorrectionToChange=new float[numelmasked*numbmodal];
        }
        for (int i=0; i<numelmasked*numbmodal;i++) {
            BFCorrectionToChange[i]=BFCorrectionInput[i];
        }
    }
}

/********************* EM RELATED METHODS ***************************/

// Update all the parameters of the model in a recursive way
void TreeEM::UpdateParameters(){
    int TypeTreeComponent=this->WhatTypeOfTreeComponent();
    if(TypeTreeComponent==LEAF){// Means it is a leaf
        /* For the moment only Gaussian distribution as simple distribution so no other case than the default one in the switch*/
        switch (this->GetDistributionType()) {
            default:
                this->UpdateGaussianParameters();
                //cout<<"Gaussian parameters updated";
                break;
        }
    }
    else{ // Recursive part
        for (int c=0; c<this->GetChildren().size(); c++) {
            //cout<<"Updating parameters of child "<<c<<endl;
            this->GetChild(c)->UpdateParameters();
        }
    }
    return;
}

void TreeEM::UpdateGaussianParameters(){
    if (this->GetDistributionType()!=1){
        cout<<"It is not a Gaussian distribution"<<endl;
        return;
    }
    this->UpdateGaussianMean();
    this->UpdateGaussianVariance();
    return;
}

void TreeEM::UpdateGaussianMean(){
    int numel=this->GetNumberElements(); // number of elements in the image
    int numbmodal=this->GetNumberModalities(); // number of modalities to consider
    //float * PointerToDataBegin=static_cast<float *>(this->GetDataImage()->data);
    float * PointerToDataBegin=this->MakeDataBFCorrected();
    for (int m=0;m<numbmodal;m++){
        float * PointerToImageBegin_PTR=&PointerToDataBegin[m*numel]; // Points the data to the beginning of each modality image
        float * PointerToNormRespBegin_PTR=this->GetNormResp(); // Points to the beginning of the NormalizedResponsabilities
        float meantmp=0; // Initialisation of the numerator in the update of the mean
        float sumResp=0; // Initialisation of the denominator in the update of the mean
        for(int i=0;i<numel;i++,PointerToImageBegin_PTR++){
            if (this->GetL2S()[i]>=0){
                /* Each time the considered voxel is active (<=> L2S contains an index), update the numerator and the denominator*/
                meantmp+=(*PointerToImageBegin_PTR)*(*PointerToNormRespBegin_PTR);
                if (*PointerToImageBegin_PTR>1) {
                    cout<<"Pb in normalisation image";
                }
                sumResp+=(*PointerToNormRespBegin_PTR);
                PointerToNormRespBegin_PTR++;
            }
        }
        if (sumResp!=0) {
            this->GetMean()[m]=meantmp/sumResp;
        }
        else{
            this->GetMean()[m]=meantmp/this->GetNumberMaskedElements();
        }
        //cout<<"Mean is now "<< this->GetMean()[m];
    }
    // Clearing memory
    if (PointerToDataBegin!=NULL) {
        delete [] PointerToDataBegin;
        PointerToDataBegin=NULL;
    }
    return;
}

void TreeEM::UpdateGaussianVariance(){
    
    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    float sumResp=0; // Initialisation of the denominator
    //float * PointerToDataBegin = static_cast<float *>(this->GetDataImage()->data);
    float * PointerToDataBegin=this->MakeDataBFCorrected();
    float * PointerToDataBegin_PTR1=PointerToDataBegin;
    float * PointerToDataBegin_PTR2=PointerToDataBegin;
    float * NormalisedResponsabilities_PTR=this->GetNormResp();
    int * L2S_PTR=this->GetL2S();
    
    // Calculation of the denominator (sum over the active voxels of the normalised responsabilities)
    for (int i=0; i<numel; i++) {
        if(L2S_PTR[i]>=0){
            sumResp+=(*NormalisedResponsabilities_PTR);
            NormalisedResponsabilities_PTR++;
        }
    }

    float * VarianceToUpdate=this->GetVariance();
    float * MeanToUse=this->GetMean();
    for(int m1=0;m1<numbmodal;m1++){
         // First data pointer to the beginning of the modality m1 considered
        for(int m2=0;m2<numbmodal;m2++){
            PointerToDataBegin_PTR1=&PointerToDataBegin[m1*numel];
            PointerToDataBegin_PTR2=&PointerToDataBegin[m2*numel]; // Second Data pointer to the beginning of the modality m2 considered
            NormalisedResponsabilities_PTR=this->GetNormResp(); // Reinitialisation of the responsabilities pointer to the beginning
            VarianceToUpdate[m1+m2*numbmodal]=0;

            for(int i=0;i<numel;i++,PointerToDataBegin_PTR1++,PointerToDataBegin_PTR2++){
                if (L2S_PTR[i]>=0) {
                    // Update of the numerator of the Variance Calculation only if in the case of an active voxel
                    VarianceToUpdate[m1+m2*numbmodal]+=(*NormalisedResponsabilities_PTR)*((*PointerToDataBegin_PTR1)-MeanToUse[m1])*((*PointerToDataBegin_PTR2)-MeanToUse[m2]);
                    NormalisedResponsabilities_PTR++;
                }
            }
            if (sumResp !=0) {
                VarianceToUpdate[m1+m2*numbmodal]/=sumResp;
                if (m1==m2) {
                    VarianceToUpdate[m1+m2*numbmodal]=VarianceToUpdate[m1+m2*numbmodal]<=1E-6?1E-6:VarianceToUpdate[m1+m2*numbmodal]; // in order to avoid going to 0 if too sharp distribution but not changing non diagonal of variance
                }
                

            }
            else{
                VarianceToUpdate[m1+m2*numbmodal]/=numelmasked;// this->GetNumberMaskedElements();
            }
            // Use of the symmetry property of the Variance matrix
            VarianceToUpdate[m2+m1*numbmodal]=VarianceToUpdate[m1+m2*numbmodal];
        }
    }
    // Clearing memory
    if (PointerToDataBegin!=NULL) {
        delete []  PointerToDataBegin;
        PointerToDataBegin=NULL;
    }
    return;
}

// Update the obtention of the bias field coefficients
void TreeEM::UpdateBFCoeffs(){
    float * BFCoeffsToUpdate=this->MakeFinalBFCoeffsChildren();
    this->SetBFCoeffs(BFCoeffsToUpdate);
    // Clearing memory not needed anymore.
        if (BFCoeffsToUpdate!=NULL) {
            delete [] BFCoeffsToUpdate;
            BFCoeffsToUpdate=NULL;
        }
}

// Update the Bias Field Correction
void TreeEM::UpdateBFCorrection(){
    float * BFCorrectionToUpdate=this->MakeBFCorrection();
    this->SetBFCorrection(BFCorrectionToUpdate);
    // Clearing memory not needed anymore
    if(BFCorrectionToUpdate!=NULL){
        delete [] BFCorrectionToUpdate;
        BFCorrectionToUpdate=NULL;
    }
}

// Recursively update the distribution
void TreeEM::UpdateDistribution(){
    this->ReputDistributionToZero();
    //cout<<"Distribution put back to 0"<<endl;
    int numbchild=this->GetNumberChildren();
    if (this->GetNumberChildren()==0) {/* If it is a leaf the distribution is then calculated according to the parameters stored. For the moment, their is only one choice of simple distribution : the Gaussian distribution*/
        switch (this->GetDistributionType()) {
            default:
                this->UpdateGaussianDistribution();
                break;
        }
    }
    else{/* It is not a leaf and the distribution is the weighted sum of the distributions below*/
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->UpdateDistribution();
            this->SumToDistributionWeightedDistributionOfChild(c);
        }
    }
//    string Filename=(g_scratch_dir + "/TryDistributionDataT1Nifty").c_str();
//    nifti_image * PartialRes=SavePartialResult(this->GetPartialResult(DISTRIBUTION),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);
    return;
}

void TreeEM::ReputDistributionToZero(){
    int numelmasked=this->GetNumberMaskedElements();
    float * Distribution_PTR=this->GetDistribution();
    for(int i=0;i<numelmasked;i++,Distribution_PTR++){
        *Distribution_PTR=0;
    }
    return;
}


void TreeEM::SumToDistributionWeightedDistributionOfChild(int c){
    // First need to check that child c exists
    if(c>=this->GetNumberChildren()||c<0){
        cout<<"Impossible to add weighted distribution of a child that does not exist !"<<endl;
        return;
    }
    float * Distribution_PTR=this->GetDistribution();
    float * DistributionChild_PTR=this->GetChild(c)->GetDistribution();
    int numelmasked=this->GetNumberMaskedElements();
    //cout<<"the multiplicative weight is "<<this->GetChild(c)->GetNormWeight()<<endl;
    for (int i=0; i<numelmasked; i++,Distribution_PTR++,DistributionChild_PTR++) {
        *Distribution_PTR=(*Distribution_PTR)+this->GetChild(c)->GetNormWeight()*(*DistributionChild_PTR);
    }
    return;
}

void TreeEM::UpdateGaussianDistribution(){
    //cout<<"Updating Gaussian Distribution"<<endl;
    if (this->GetDistributionType()!=1) {
        cout<<"This is not a leaf with Gaussian distribution : cannot be updated as wanted"<<endl;
        return;
    }
    // Distribution put back to zero every where
    this->ReputDistributionToZero();
    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
    float * DataPointer=this->MakeDataBFCorrected();
    //float * DataPointer=static_cast<float *>(this->GetDataImage()->data);

    
    matrix<float> VarianceMatrix=matrix<float>(numbmodal);
    for (int m1=0; m1<numbmodal; m1++) {
        for (int m2=0; m2<numbmodal; m2++) {
            VarianceMatrix.setvalue(m1, m2, this->GetVariance()[m1+m2*numbmodal]);
            
        }
    }
    
    // Calculation of the factor in front of the exponential
    float DeterminantVariance=VarianceMatrix.determinant();
    //cout<<"the Variance determinant is "<<DeterminantVariance<<endl;
    float NormalisationFactor=1.0/(float)(powf(2*M_PI , (float)((float)numbmodal/2.0))*powf(DeterminantVariance, 0.5));
    //cout <<"The normalisation factor is "<<NormalisationFactor<<endl;
    //Initialisation of the needed element to calculate the inside of the exponential
    float * Distribution_PTR=this->GetDistribution();
    float * InvertedVariance=new float[numbmodal*numbmodal];
    if (numbmodal>1) {
        VarianceMatrix.invert();
        
        bool success;
        for(int m1=0;m1<numbmodal;m1++){
            for (int m2=0; m2<numbmodal; m2++) {
                VarianceMatrix.getvalue(m1, m2, InvertedVariance[m1+m2*numbmodal], success);
            }
        }

    }
    else{
        bool success;
        VarianceMatrix.getvalue(0, 0, InvertedVariance[0], success);
        *InvertedVariance=1.0/(*InvertedVariance);
    }
        float temp;
    
    // Calculation of the inside of the exponential
    //cout<<InvertedVariance[0]<<endl;
    for (int i=0; i<numel; i++) {
        
        if (this->GetL2S()[i]>=0){
            for (int m1=0; m1<numbmodal; m1++) {
                temp=0;
                for (int m2=0; m2<numbmodal; m2++) {
                    temp+=(DataPointer[i+m2*numel]-this->GetMean()[m2])*InvertedVariance[m2+m1*numbmodal];
                }
                *Distribution_PTR=(*Distribution_PTR)+temp*(DataPointer[i+m1*numel]-this->GetMean()[m1]);
            }
            Distribution_PTR++;
        }
    }
    
    // Filling the distribution array with the value
    Distribution_PTR=this->GetDistribution();

    int numelmasked=this->GetNumberMaskedElements();
    for (int i=0;i<numelmasked;i++,Distribution_PTR++){
            *Distribution_PTR=NormalisationFactor*exp(-0.5*(*Distribution_PTR));
    }
    // Clearing space needing for inverse of variance
    if(InvertedVariance!=NULL){
        delete [] InvertedVariance;
        InvertedVariance=NULL;
    }
    if (DataPointer!=NULL) {
        delete [] DataPointer;
        DataPointer=NULL;
    }
//    string Filename=(g_scratch_dir + "/TryDistributionDataT1Nifty").c_str();
//    nifti_image * PartialRes=SavePartialResult(this->GetPartialResult(DISTRIBUTION),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);    return;
}


// Recursively update the non normalised responsabilities
void TreeEM::UpdateNonNormResp(){
    this->ReputNonNormRespToZero();
    int numbchild=this->GetNumberChildren();
    bool testNormPriors=this->ArePriorsNormalised();
    bool testNullPriors=(this->GetPriorsVector()[numbchild-1]==NULL);
    for (int c=0;c<numbchild;c++){
        if(this->GetChild(c)->GetNumberChildren()==0){ /* Case of a leaf : Need to calculate the multiplication of the distribution by the priors. If the priors are not normalised meaning we are at a subclass level, need also to be weighted*/
            switch (this->GetChild(c)->GetDistributionType()) {
                default:
                    //cout<<"Calculating NonNormResp for child "<<c<<endl;
                {this->GetChild(c)->CalculatingNonNormResp();
                    if (!testNormPriors || testNullPriors) {
                        //cout<<"Need to multiply by norm weight"<<endl;
                        this->GetChild(c)->MultiplyCurrentNonNormRespByNormWeight();
                    }
                    //cout<<" No need for multiplication by norm weight"<<endl;
                }
                    break;
            }
            this->SumToNonNormRespWeightedNonNormRespOfChild(c);
        }
        else{
            this->GetChild(c)->UpdateNonNormResp();
            this->SumToNonNormRespWeightedNonNormRespOfChild(c);
        }
    }
    if (!testNormPriors&&!testNullPriors) {/* if the priors are normalised, meaning that we are at first level or not handling subclasses with same priors, we do not need the weights (redundant with priors) */
        //this->MultiplyCurrentNonNormRespByNormWeight();
    }
//    string Filename=(g_scratch_dir + "/TryNonNormRespT1Nifty").c_str();
//    nifti_image * PartialRes=SavePartialResult(this->GetPartialResult(NONNORMRESP),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);
        return;
    }

void TreeEM::ReputNonNormRespToZero(){
    float * NonNormResp_PTR=this->GetNonNormResp();
    int numelmasked=this->GetNumberMaskedElements();
    for(int i=0;i<numelmasked;i++,NonNormResp_PTR++){
        *NonNormResp_PTR=0;
    }
}

void TreeEM::CalculatingNonNormResp(){
    int numel=this->GetNumberElements();
    float * Distribution_PTR=this->GetDistribution();
    float * NonNormResp_PTR=this->GetNonNormResp();
    if (this->GetPriors()==NULL) {
        for (int i=0; i<numel; i++) {
            if (this->GetL2S()[i]>=0){
                (*NonNormResp_PTR)=(*Distribution_PTR);
                NonNormResp_PTR++;
                Distribution_PTR++;
            }
        }
    }
    else{
        //cout<< "NonNormResp with priors used"<<endl;
        float * PriorsPointer=static_cast<float *>(this->GetPriors()->data);
        float * PriorsPointer_PTR=PriorsPointer;
        for (int i=0; i<numel; i++,PriorsPointer_PTR++) {
            if (this->GetL2S()[i]>=0){
                *NonNormResp_PTR=(*Distribution_PTR)*(*PriorsPointer_PTR);
                NonNormResp_PTR++;
                Distribution_PTR++;
            }
            // the data in Priors are of size numel and not numelmasked as NonNormResp and Distribution
        }
    }
}

void TreeEM::MultiplyCurrentNonNormRespByNormWeight(){
    int numelmasked=this->GetNumberMaskedElements();
    float * NonNormResp_PTR=this->GetNonNormResp();
    for (int i=0; i<numelmasked; i++,NonNormResp_PTR++) {
        (*NonNormResp_PTR)=(*NonNormResp_PTR)*this->GetPartNormWeight();   
    }
    //cout<<"PartNormWeight is "<<this->GetPartNormWeight()<<endl;
    return;
}

void TreeEM::SumToNonNormRespWeightedNonNormRespOfChild(int c){
    // First need to check that child c exists
    if(c>=this->GetNumberChildren()||c<0){
        cout<<"Impossible to add weighted non normalised responsabilities of a child that does not exist !"<<endl;
        return;
    }
    float * NonNormResp_PTR=this->GetNonNormResp();
    float * NonNormRespChild_PTR=this->GetChild(c)->GetNonNormResp();
    int numelmasked=this->GetNumberMaskedElements();
    for (int i=0; i<numelmasked; i++,NonNormResp_PTR++,NonNormRespChild_PTR++) {
        *NonNormResp_PTR=(*NonNormResp_PTR)+(*NonNormRespChild_PTR);
    }
    return;
}

/* Using the NonNormResp calculated before, goes back to the root to calculate the normalisation factor and divide elementwise to obtain the normalised responsabilities*/
void TreeEM::UpdateNormResp(){
    // First find the root of the looked element
    int numbchildRoot=this->FindRoot()->GetNumberChildren();
    //cout<<"the number of children in Root is "<<numbchildRoot<<endl;
    // Then Calculate the overall sum of the non normalised responsabilities to obtain the normalisation factor
    this->FindRoot()->ReputNonNormRespToZero();
    for (int c=0; c<numbchildRoot; c++) {
        (this->FindRoot())->SumToNonNormRespWeightedNonNormRespOfChild(c);
    }
//    string Filename=(g_scratch_dir + "/TryRootNonNormRespT1Nifty").c_str();
//    nifti_image * PartialRes=SavePartialResult(this->FindRoot()->GetPartialResult(NONNORMRESP),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        this->GetChild(c)->UpdateNormResp();
    }
    // Finally Divide the NonNormResp when going through this by the SumNonNormResp found just before;
    this->DivideNonNormRespByRootNonNormResp();
//    Filename=(g_scratch_dir + "/TryNormRespT1Nifty3").c_str();
//    PartialRes=SavePartialResult(this->GetPartialResult(NORMRESP),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);
//    //this->NormaliseNonNormResp();
//     Filename=(g_scratch_dir + "/TryNormRespT1Nifty2").c_str();
//     PartialRes=SavePartialResult(this->FindRoot()->GetPartialResult(NORMRESP),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);
    return;
    
}

// Recursively divide NonNormResp of the Children
void TreeEM::NormaliseNonNormResp(){
    
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        //cout<<"Normalising NonNormResp of child "<<c<<endl;
        this->GetChild(c)->NormaliseNonNormResp(); // Recursive part
    }
    this->DivideNonNormRespByRootNonNormResp();
    return;
}

void TreeEM::DivideNonNormRespByRootNonNormResp(){
    // No need to reput NormResp to zero, only replacement of the contained values
    float * NormResp_PTR=this->GetNormResp();
    float * NonNormResp_PTR=this->GetNonNormResp();
    float * RootNonNormResp_PTR=this->FindRoot()->GetNonNormResp();
//    string Filename=(g_scratch_dir + "/TryRootNonNormRespT1Nifty2").c_str();
//    nifti_image * PartialRes=SavePartialResult(this->FindRoot()->GetPartialResult(NONNORMRESP),this->GetDataImage(),(char *)Filename.c_str());
//    nifti_image_write(PartialRes);
//    nifti_image_free(PartialRes);
    int numelmasked=this->GetNumberMaskedElements();
    for (int i=0; i<numelmasked; i++,NormResp_PTR++,NonNormResp_PTR++,RootNonNormResp_PTR++) {
        if(*RootNonNormResp_PTR>0){
            (*NormResp_PTR)=(float)(*NonNormResp_PTR)/(*RootNonNormResp_PTR);
        }
        else{
            cout<<"Root non norm resp is 0, there is a problem somewhere..."<<endl;
            (*NormResp_PTR)=1.0;
        }
    }
}

// Recursively obtain the NonNormWeight (sum over the voxels for one class of the normalised responsabilities)
void TreeEM::UpdateNonNormWeights(){
    int numelmasked=this->GetNumberMaskedElements();
    int numbchild=this->GetNumberChildren();
    this->SetNonNormWeight(0);
    float * NormResp_PTR=this->GetNormResp();
    for(int i=0;i<numelmasked;i++,NormResp_PTR++){
        this->SetNonNormWeight(this->GetNonNormWeight()+(*NormResp_PTR));
    }
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->UpdateNonNormWeights();
    }
    return;
}

// Recursively Update the normalised weights when the non normalised ones are obtained
void TreeEM::UpdateNormWeights(){
    // First determination of the normalisation factor for the children
    float SumNonNormWeights=0;
    int numbchild=this->GetNumberChildren();
    for (int c=0;c<numbchild; c++) {
        SumNonNormWeights+=this->GetChild(c)->GetNonNormWeight();
    }
    // Then normalisation of the not normalised weights for the children
    for (int c=0; c<numbchild; c++) {
        if (SumNonNormWeights>0) {
            this->GetChild(c)->SetNormWeight(this->GetChild(c)->GetNonNormWeight()/SumNonNormWeights);
        }
        else{// in case sum of the weights is 0
            this->GetChild(c)->SetNormWeight(1/this->GetNumberChildren());
        }
        //cout<< "Norm weight of child "<<c<<" is of value "<<this->GetChild(c)->GetNormWeight()<<endl;
        this->GetChild(c)->UpdateNormWeights(); // Recursive part
    }
    
    return;
}

//void TreeEM::MakeNonNormWeightedSum(float & NonNormSum){
//    int numbchild=this->GetNumberChildren();
//    int numelmasked=this->GetNumberMaskedElements();
//    int numel=this->GetNumberElements();
//    if(numbchild==0){// case where it is a leaf
//        switch (this->GetDistributionType()){
//                default : {
//                float * DistSum_tmp=this->MakeGaussianDistribution();
//                float * DistSum_tmp_PTR=DistSum_tmp;
//                if(this->GetPriorsDirect()==NULL){// case we have to multiply by the normalised weight
//                    float NormWeightUsed=this->GetNormWeight();
//                    for(int i=0;i<numelmasked;i++,DistSum_tmp_PTR++){
//                        *DistSum_tmp_PTR*=NormWeightUsed;
//                    }
//                    if(this->GetPriors()!=NULL){// case where we are using atlases
//                        DistSum_tmp_PTR=DistSum_tmp;
//                        L2S_PTR=this->GetL2S();
//                        float * Priors_PTR=static_cast<float*>(this->GetPriors()->data);
//                        for(int i=0;i<numel;i++,L2S_PTR++,Priors_PTR++){
//                            if(*L2S_PTR>=0){
//                               (*DistSum_tmp_PTR)*=*Priors_PTR;
//                                DistSum_tmp_PTR++;
//                            }
//                        }
//                    }

//                }
//        }

//        }
//        // Copying the final result in the destined float array
//        float * NonNormSum_PTR=NonNormSum;
//        DistSum_tmp_PTR=DistSum_tmp;
//        for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
//            *NonNormSum_PTR=*DistSum_tmp_PTR;
//        }
//        // Clearing memory not needed anymore
//        delete[] DistSum_tmp;
//        DistSum_tmp=NULL;
//    }
//    else {
//        for(int c=0;c<numbchild;c++){
//            // Initialisation of temporary result
//            float * DistSum_tmp=new float[numelmasked];
//            for(int i=0;i<numelmasked;i++){
//                DistSum_tmp[i]=0;
//            }
//            this->GetChild(c)->MakeNonNormWeightedSum(DistSum_tmp); // recursive part
//            if(this->GetChild(c)->GetPriorsDirect()==NULL){ // need to multiply by the normalised weight
//                float * NonNormSum_PTR=NonNormSum;
//                float * DistSum_tmp_PTR=DistSum_tmp;
//                float NormWeightUsed=this->GetChild(c)->GetNormWeight();
//                for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
//                    (*NonNormSum_PTR)+=NormWeightUsed*(*DistSum_tmp_PTR);
//                }
//            }
//            else{// case we do not have to multiply by the normalised weight before adding
//                float * NonNormSum_PTR=NonNormSum;
//                float * DistSum_tmp_PTR=DistSum_tmp;
//                for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
//                    (*NonNormSum_PTR)+=(*DistSum_tmp_PTR);
//                }
//            }
//            delete [] DistSum_tmp;
//            DistSum_tmp=NULL;
//        }
//    }
//  return;
//}

//void TreeEM::UpdateNewNormResp(){

//}

// REMARK : Be really careful when Mask is changed from one layer to another. Always check that the right amount of memory is then allocated
float TreeEM::CompleteLogLikelihood(){
    float CompleteLogLikelihood=0;
    int numelmasked=this->GetNumberMaskedElements();
    int numel=this->GetNumberElements();
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        int * L2S_PTR=this->GetL2S();
        float * NormResp_PTR=this->GetChild(c)->GetNormResp();
        float * Priors_PTR=NULL;
        if(this->GetChild(c)->GetPriors()!=NULL){
         Priors_PTR=static_cast<float *>(this->GetChild(c)->GetPriors()->data);
        }
        float * Distribution_PTR=this->GetChild(c)->GetDistribution();
        
        if (Priors_PTR!=NULL) {
            for(int i=0;i<numel;i++){
                if ((*L2S_PTR)>=0){
                    if((*Distribution_PTR)>0.00001){
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(*Distribution_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(0.00001);
                    }
                    
                    if((*Priors_PTR)>0.00001){
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(*Priors_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(0.00001);
                    }
                    
                    NormResp_PTR++;
                    Distribution_PTR++;
                }
                Priors_PTR++;
                L2S_PTR++;
            }
        }
        
        else {
            for(int i=0;i<numelmasked;i++,Distribution_PTR++,NormResp_PTR++){
                
                    if((*Distribution_PTR)>0){
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(*Distribution_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(0.00001);
                    }
                    
                    if(this->GetChild(c)->GetNormWeight()>0){
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(this->GetChild(c)->GetNormWeight());
                    }
                    else{
                        CompleteLogLikelihood+=(*NormResp_PTR)*logf(0.00001);
                    }

            }

        }
    }
    return CompleteLogLikelihood;
}

void TreeEM::EMCompleteIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration){
    
    OldCompleteLogLikelihood=CompleteLogLikelihood;
    this->UpdateParameters();
    //cout<<"Update Parameters done"<<endl;

    this->UpdateDistribution();
    //cout<<"Update distribution done"<<endl;
    this->UpdateNonNormResp();
    //cout<<"Update NonNormResp done"<<endl;
    this->UpdateNormResp();
    //cout<<"Update NormResp done"<<endl;
    this->UpdateNonNormWeights();
    //cout<<"Update NonNormWeights done"<<endl;
    this->UpdateNormWeights();
    if (BFFlag) {
        this->UpdateBFCoeffs();
        this->UpdateBFCorrection();
    }
    CompleteLogLikelihood=this->CompleteLogLikelihood();
    Iteration++;
    //this->SaveAllClasses((g_scratch_dir + "/TestBRATSResult").c_str());
}

void TreeEM::RunFullEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration){
    bool ValidityInitialisationTree = 1;// this->IsTreeWellInitialised();
    if (!ValidityInitialisationTree){
        cout<< "Tree not well initialised, EM cannot be performed" << endl;
        return;
    }
    cout<<"Tree well initialised we can run the EM"<<endl;
//    do {
//        this->EMCompleteIteration(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration);
//    }
    while (((CompleteLogLikelihood-OldCompleteLogLikelihood)/OldCompleteLogLikelihood>0.0001 || Iteration<6) && Iteration<35) {
        cout<<"we are at iteration "<<Iteration<<endl;
        cout<< CompleteLogLikelihood <<" and the old CLL "<< OldCompleteLogLikelihood<<endl;
        this->EMCompleteIteration(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration);
    }
}

TreeEM * TreeEM::RunFullBiASM(){
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration);
    if(BFFlag){
        this->SaveBFCorrectedData((g_scratch_dir + "/BRATS_DataCorrected").c_str());
        this->SaveBFCorrection((g_scratch_dir + "/BRATS_BFCorrection").c_str());
    }
    BFFlag=0;
    TreeEM * TreeSplit=this;
    TreeEM * TreeMerge=this;
    MergeKLD * MergeTest=TreeMerge->GetToMerge();
    SplitKLD * SplitTest=TreeSplit->GetToSplit();
    int numberAllLeaves=this->GetNumberAllLeaves();
    while ((SplitTest!=NULL || MergeTest!=NULL) && numberAllLeaves<=15) {
        if(SplitTest!=NULL){
            delete SplitTest;
            SplitTest=NULL;
        }
        if(MergeTest!=NULL){
            delete MergeTest;
            MergeTest=NULL;
        }
        TreeSplit=TreeSplit->RunSplitOperation();
        int numbchild=TreeSplit->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            (TreeSplit->GetChild(c))->PutAllLeavesToChildrenLevel();
        }
        TreeMerge=TreeSplit->RunMergeOperation();
        TreeSplit=TreeMerge;
        TreeSplit->SaveTreeInTextFile((g_scratch_dir + "/DataTree.txt").c_str());

        SplitTest=TreeSplit->GetToSplit();
        MergeTest=TreeMerge->GetToMerge();
        numberAllLeaves=TreeSplit->GetNumberAllLeaves();
    }
    if(MergeTest!=NULL){
        delete MergeTest;
        MergeTest=NULL;
    }
    if(SplitTest!=NULL){
        delete SplitTest;
        SplitTest=NULL;
    }
    return TreeSplit;
}


/**************  TREE RELATED METHODS ******************/


/*A root element is defined by the fact that the pointer to Parent is NULL. FindRoot recursively finds the element whose Parent is NULL. Should be the only one to use links upwards*/
TreeEM * TreeEM::FindRoot(){
    if(this->GetParent()==NULL){
        return this;
    }
    else{
        return this->GetParent()->FindRoot();
    }
}

// If the current node is considered as a mixture but has only one child, the node is replaced by the child itself and transformed as a simple distribution.
void TreeEM::CollapseOnlyChild(){
    if (this->GetNumberChildren()!=1) {
        cout<<"Collapsing only when number of children is 1"<<endl;
        return;
    }
    //this->SetDistribution(this->GetChild(0)->GetDistribution());
    //this->SetNormResp(this->GetChild(0)->GetNormResp());
    //this->SetNonNormResp(this->GetChild(0)->GetNonNormResp());
    //this->SetNormWeight(this->GetChild(0)->GetNormWeight());
    //this->SetNonNormWeight(this->GetChild(0)->GetNonNormWeight());
    this->SetParameters(this->GetChild(0)->GetParameters());
    delete this->GetChild(0);
    this->Children.erase(this->Children.begin());

}

// Recursively collapsing all only children in the tree to avoid to have only one child
void TreeEM::CollapseOnlyChildTot(){
    if (this->GetNumberChildren()==1) {
        this->CollapseOnlyChild();
        return;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->CollapseOnlyChildTot();
        }
    }
    return;
}

void TreeEM::CreateAndAddChild(nifti_image * PriorsInput=NULL,int DistributionType=1){
    
    // First taking care of all needed changes on the created child to make a coherent tree
    TreeEM * ChildCreated=new TreeEM();
    //cout<<ChildCreated->GetDataDirect()<<endl;
    ChildCreated->SetParent(this); // set the parent and allocate the proper amount of memory
    //cout<<ChildCreated->GetDataDirect()<<endl;
    cout<<"Allocate memory when creating child"<<endl;
    // Allocation of good amount of memory for Distribution, NonNormResp, NormResp
    int numelmasked=ChildCreated->GetNumberMaskedElements();
    //cout<<"Setting distribution";
    if (ChildCreated->GetDistribution()!=NULL) {
        delete [] ChildCreated->Distribution;
        ChildCreated->Distribution=NULL;
    }
    float * DistributionToSet=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        DistributionToSet[i]=0;
    }
    ChildCreated->SetDistribution(DistributionToSet);
    if(DistributionToSet!=NULL){
        delete [] DistributionToSet;
    }
    //cout<<"Distribution Set";
    if (ChildCreated->GetNonNormResp()!=NULL) {
        delete []ChildCreated->NonNormResp;
        ChildCreated->NonNormResp=NULL;
    }
    float * NonNormRespToSet=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        NonNormRespToSet[i]=0;
    }
    ChildCreated->SetNonNormResp(NonNormRespToSet);
    if(NonNormRespToSet!=NULL){
        delete [] NonNormRespToSet;
        NonNormRespToSet=NULL;
    }
    //cout<<"NonNormResp set";
    if (ChildCreated->NormResp!=NULL) {
        delete [] ChildCreated->NormResp;
        ChildCreated->NormResp=NULL;
    }
    float * NormRespToSet=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        NormRespToSet[i]=0;
    }
    ChildCreated->SetNormResp(NormRespToSet);
    if(NormRespToSet!=NULL){
        delete [] NormRespToSet;
        NormRespToSet=NULL;
    }
    //cout<< "NormResp set";
    
    // Handling the Priors input
    ChildCreated->SetPriors(PriorsInput); // Takes care of the validity checking of the priors and of making usable as probabilities
    
    // Handling the distribution type input. Must necessary be a simple distribution. By default Gaussian distribution chosen
    if (DistributionType==0){ // Cannot be a mixture if leaf
        DistributionType=1;
    }
    //cout<<"Normally the distribution type is now "<<DistributionType;
    ChildCreated->CreateAllocateAndInitializeParameters(DistributionType);
    cout<<"Parameters allocated for child"<<endl;
    // Then taking care of changes needed on parent on which ChildCreated is added
    this->MakeParametersMixture(); // Modify parameters structure : check if it is already a mixture otherwise make it a mixture parameters
    this->Children.push_back(ChildCreated);
}

void TreeEM::InitialiseBasicTree(){
    if (this->IsBasicTree()) {
        srand (time(NULL));
        cout<<"This is a basic Tree !!";
        //vector<nifti_image*> PriorsVector=this->GetPriorsVector();
        //cout<<this->GetPriorsVector().size();
        if(ArePriorsNormalised()){
            if (this->GetPriorsVector()[this->GetNumberChildren()-1]==NULL) {// <=> no priors considered
                /* If there is no prior considered, the same value is put everywhere, equal 1/number of initial children*/
                cout<<"No priors set beforehand"<<endl;
                float * NormResp_PTR=NULL;
                int numbchild=this->GetNumberChildren();
                int numbmodal=this->GetNumberModalities();
                int numelmasked=this->GetNumberMaskedElements();
                for (int c=0;c<numbchild; c++) {
                    NormResp_PTR=this->GetChild(c)->GetNormResp();
                    for (int i=0; i<numelmasked; i++,NormResp_PTR++) {
                        (*NormResp_PTR)=(float)1.0/numbchild;
                    }
                    this->GetChild(c)->SetNormWeight((float)(1.0/(float)this->GetNumberChildren()));

                    // Random initialisation of the mean, variance set to the identity matrix
                    if(!this->GetChild(c)->CheckForValidityOfParametersStructure()){
                        this->GetChild(c)->CreateAllocateAndInitializeParameters(1);
                    }
                    for (int m=0; m<numbmodal; m++) {
                        
                        this->GetChild(c)->GetMean()[m]=(float)((float)rand()/(float)RAND_MAX);
                        //cout<<"Address mean "<<this->GetChild(c)->GetMean()<<" and "<<this->GetChild(c)->ParametersDistribution->ValueParameters;
                        //cout<<"Value mean = "<<this->GetChild(c)->GetMean()[m]<<endl;
                        this->GetChild(c)->GetVariance()[numbmodal*m]=0.1;
                    }
                }
               this->UpdateDistribution();
                this->UpdateNonNormResp();
                this->UpdateNonNormWeights();
                this->UpdateNormResp();
                this->UpdateNormWeights();
                if (BFFlag) {
                    this->UpdateBFCoeffs();
                    this->UpdateBFCorrection();
                }
            }
            else {// Priors are normalised and point to valid nifti_images
                int numel=this->GetNumberElements();
                int numbchild=this->GetNumberChildren();
                int numelmasked=this->GetNumberMaskedElements();
                float * NormResp_PTR=NULL;
                float * Priors_PTR=NULL;
                float tmpNormWeight;
                for (int c=0; c<numbchild; c++) {
                    //cout<<"Initialisation of the pointers"<<endl;
                    this->GetChild(c)->SetNormWeight(0);
                    tmpNormWeight=0;
                    NormResp_PTR=this->GetChild(c)->GetNormResp();
                    //cout<<"the first value of NormResp is "<<(*NormResp_PTR)<<endl;
                    Priors_PTR=static_cast<float *>(this->GetPriorsVector()[c]->data);
                    int* L2S_PTR=this->GetChild(c)->GetL2S();
                    for (int i=0; i<numel; i++, L2S_PTR++,Priors_PTR++) {
                        if ((*L2S_PTR)>=0) {
                            //cout<< "the index is "<<(*L2S_PTR)<< " and Priors are "<<(*Priors_PTR)<<endl;
                            (*NormResp_PTR)=(*Priors_PTR);
                            tmpNormWeight+=(*Priors_PTR);
                            NormResp_PTR++;
                        }
                    }
                    this->GetChild(c)->SetNormWeight(tmpNormWeight/numelmasked);
                    //cout<<"Initialisation done for this prior"<<endl;
// In basic tree there is always more than one child so no division by 0;
                }
            }
        }
        else{ // In case the priors are not normalised, then need to normalise them
            this->NormalisePriors(); // We normalise and call back InitialiseBasicTree
            cout<<this->ArePriorsNormalised();
            cout<<"Normalised Priors"<<endl;
            this->InitialiseBasicTree();
        }
    }
    else{
        cout<<"This tree is not to initialise"<<endl;
    }
    return;
}

/* Reinitialise the checks for the SM operations.*/
void TreeEM::ClearSMChecks(){
    int numbDirectLeaves=this->GetNumberDirectLeaves();
    int numbchild=this->GetNumberChildren();

    if (numbDirectLeaves!=0) {
        if (this->GetSplitCheck()!=NULL) {
            delete [] this->GetSplitCheck();
        }
        this->SplitCheck=NULL;
        this->SplitCheck=new bool[numbDirectLeaves];//{0};
        for(int i=0;i<numbDirectLeaves;i++){
            this->SplitCheck[i]=0;
        }
        if (this->GetMergeCheck()!=NULL) {
            delete [] this->GetMergeCheck();
        }
        this->MergeCheck=NULL;
        this->MergeCheck=new bool[numbDirectLeaves*numbDirectLeaves];//{0};
        for(int i=0;i<numbDirectLeaves*numbDirectLeaves;i++){
            this->MergeCheck[i]=0;
        }
        for (int dl=0; dl<numbDirectLeaves; dl++) {
            this->GetMergeCheck()[dl+dl*numbDirectLeaves]=1;
        }
    }
    else{
        if (this->GetMergeCheck()!=NULL) {
            delete [] this->GetMergeCheck();
        }
        this->MergeCheck=NULL;
        if (this->GetSplitCheck()!=NULL) {
            delete [] this->GetSplitCheck();
        }
        this->SplitCheck=NULL;
    }
    for (int c=0; c<numbchild; c++) {
        this->GetChild(c)->ClearSMChecks();
    }
    return;
}

// Modify NormWeights in the tree in order to just afterwards collapse everything to the children level
void TreeEM::ModifyNormWeightsForChildrenLevel(){
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetPriorsDirect()==NULL) {
            this->GetChild(c)->SetNormWeight(this->GetNormWeight()*this->GetChild(c)->GetNormWeight());
        }
        else{
            this->GetChild(c)->SetPriors(this->GetPriorsDirect());
        }
        this->GetChild(c)->ModifyNormWeightsForChildrenLevel(); // Recursive part
    }
}

// Using the preliminary function ModifyNormWeights for ChildrenLevel then collect all the leaves, copy them, delete the children and put the copied leaves instead. Allows for changing the value of parent for the copied leaves right away before adding them as children
void TreeEM::PutAllLeavesToChildrenLevel(){
    vector<TreeEM *> Leaves=this->GetAllLeaves();
    int numbLeaves=this->GetNumberAllLeaves();
    vector<TreeEM *> CopiedLeaves;
    // Copy the leaves with the right weights
    for (int l=0; l<numbLeaves; l++) {
        CopiedLeaves.push_back(Leaves[l]->CopyTree(this));
    }
    // Delete the children (not useful anymore)
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        delete this->GetChild(0);
        this->Children.erase(this->Children.begin());
    }
    // Set the copied leaves as the children
    for (int l=0; l<numbLeaves; l++) {
        this->AddChild(CopiedLeaves[l]);
    }
}

void TreeEM::RecPutAllLeavesToChildrenLevel(){
    if (this->GetNumberChildren()==this->GetNumberDirectLeaves()) {
        return;
    }
    else{
        for (int c=0; c<this->GetNumberChildren(); c++) {
            if (this->GetChild(c)->GetNumberChildren()!=0) {
                int numbchild2=this->GetChild(c)->GetNumberChildren();
                for (int c2=0; c2<numbchild2; c2++) {
                    this->GetChild(c)->GetChild(c2)->SetNormWeight(this->GetChild(c)->GetNormWeight()*this->GetChild(c)->GetChild(c2)->GetNormWeight());
                    TreeEM * CopiedChild=this->GetChild(c)->GetChild(c2)->CopyTree(this);
                    this->AddChild(CopiedChild);
                }
                delete this->GetChild(c);
                this->Children.erase(this->Children.begin()+c);
            }
        }
    }
}


/**************** TESTS ON TREE *************************/

/* Test if the tree on which applied is a leaf. A leaf is defined as a tree that does not have any children*/
bool TreeEM::IsLeaf(){
    if(this->GetNumberChildren()==0){
        return 1;
    }
    return 0;
}

/* Test if the tree on which applied is a branch. A branch is defined as a tree that has both a defined Parent and Children*/
bool TreeEM::IsBranch(){
    if(this->GetParent()!=NULL && this->GetNumberChildren()>0){
        return 1;
    }
    return 0;
}

/* Test if the tree on which applied is a root. A root is defined as a tree that has no Parent*/
bool TreeEM::IsRoot(){
    if(this->GetParent()==NULL){
        return 1;
    }
    return 0;
}

/*Returns as an integer the type of component (root,branch,leaf,initialNode) considered. Could be coded afterwards as enum*/
int TreeEM::WhatTypeOfTreeComponent(){
    if (this->IsBranch()) {
        return BRANCH;
    }
    if (this->IsRoot()&&this->IsLeaf()) {
        return INITIALNODE; // is an initial node with no parent and no children
    }
    if (this->IsLeaf()) {
        return LEAF;
    }

    return ROOT; // Is a root
}



/* Tests if the tree considered is a basic tree or not. A basic tree is defined as a tree with only a root and one level of children*/
bool TreeEM::IsBasicTree(){
    int numbchild=this->GetNumberChildren();
    if (this->WhatTypeOfTreeComponent()!=ROOT) {
        return 0;
    }
    for (int c=0; c<numbchild; c++) {
        if (this->Children[c]->GetNumberChildren()!=0) {
            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsTreeValid(){
    if (this->GetDataImage()==NULL){ // Force the tree to make reference to an image to segment
        cout<< "Not valid because no data image referenced"<<endl;
        return 0;
    }
    if (!this->IsDataImageNormalised()) {
        cout<<"Not valid because data not normalised"<<endl;
        return 0;
    }
    if (this->WhatTypeOfTreeComponent()==LEAF||this->WhatTypeOfTreeComponent()==INITIALNODE){ // Mask is only set at the root level for allocation problems and sizes incompatibilities if changed through the tree
        if(this->GetMaskDirect()!=NULL){
            cout<<"Not valid because Mask defined at wrong place"<<endl;
            return 0;
        }
    }
    if (!this->IsPriorsProbabilityType()) { // Not valid if priors not of probability type
        cout<<"Not valid because Priors not set as probabilities"<<endl;
        return 0;
    }
    if(!this->IsNumbChildOKWithDistType()){ // not valid if number of children not valid with parameters
        cout<<"Not valid of incompatibility between number of children and distribution type"<<endl;
        return 0;
    }
    if(!this->CheckForValidityOfParametersStructure()){ // not valid if parameters structure not valid
        cout<<"Not valid because of unvalidity of Parameters structure"<<endl;
        return 0;
    }
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (!this->GetChild(c)->IsTreeValid()) {
            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsNumbChildOKWithDistType(){
    if (this->GetNumberChildren()==0&&this->GetDistributionType()==0){ // case considered as mixture whereas does not have any child
        return 0;
    }
    if (this->GetNumberChildren()>0 && this->GetDistributionType()>0){ // case where considered as simple distribution whereas has children
        return 0;
    }
    return 1;
}

bool TreeEM::IsDataImageNormalised(){
    if (this->GetDataImage()==NULL) {// no data image in the tree to look at
        cout<<"No data image as reference"<<endl;
        return 0;
    }
    int numbmodal=this->GetNumberModalities();
    for (int m=0; m<numbmodal; m++) {
        //cout << "The number of modalities is "<<this->GetNumberModalities()<<endl;
        float minRes=this->GetMinDataModal(m);
        float maxRes=this->GetMaxDataModal(m);
        //cout<<(minRes-0.0)<<"and "<<(maxRes-1.0)<<endl;
        if (minRes>1E-6||minRes<-1E-6){
            
            return 0;
        }
        if (maxRes>1+1E-6||maxRes<1-1E-6){
            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsDataFloat(){
    if (this->GetDataImage()==NULL) {
        cout<<"No data to check"<<endl;
        return 0;
    }
    if (this->GetDataImage()->datatype==DT_FLOAT) {
        return 1;
    }
    return 0;
}

bool TreeEM::IsPriorsFloat(){
    if (this->GetPriors()==NULL) {
        cout<<"No priors to check"<<endl;
        return 0;
    }
    if (this->GetPriors()->datatype==DT_FLOAT) {
        return 1;
    }
    return 0;
}

// Check if the direct BF coefficients obtained are directly useable or not
bool TreeEM::AreBFCoeffsDirectUseable(){
//    vector<float*> BFCoeffsToTest=this->GetBFCoeffsDirect();
//    // 1. Test on the size of the vector
//    if (BFCoeffsToTest.size()==0) {
//        return 0;
//    }
//    int numbmodal=this->GetNumberModalities();
//    // 2. Test on the compatibility of the dimensions
//    if (BFCoeffsToTest.size()!=numbmodal) {
//        return 0;
//    }
//    // 3. Test on the containts of the vector if of the right size.
//    for (int m=0; m<numbmodal; m++) {
//        if (BFCoeffsToTest[m]==NULL) {
//            return 0;
//        }
//    }
    if(this->GetBFCoeffsDirect()==NULL){
        return 0;
    }
    return 1;
}
/********************** METHODS ON PARAMETERS STRUCTURE ***********************/

int TreeEM::CalculateSizeParameters(){
    int SizeParameters;
    switch (this->GetDistributionType()) {
        case 0:{
            //cout << "We are in the case of a mixture" << endl;
            SizeParameters=0;
        }
            break;
            
        default:
        {
            //cout<< "We are looking at a Gaussian distribution" << endl;
            int numbmodal=this->GetNumberModalities();
            SizeParameters=(int)(numbmodal*(numbmodal+1)); // Size of the mean vector + size of the variance matrix (numbmodal)
            //cout<< "Size parameters is "<<SizeParameters;
            /* REMARK :
             this might be changed if decided then to consider only number of free parameters, then the Gaussian Parameters calculation function must be changed consequently */
        }
            
    }
    return SizeParameters;
    
}

int TreeEM::CalculateSizeParameters(int DistributionTypeInput){
    int SizeParameters;
    
    switch (DistributionTypeInput) {
        case 0:{
            //cout << "We are in the case of a mixture" << endl;
            SizeParameters=0;
        }
            break;
            
        default:
        {
            //cout<< "We are looking at a Gaussian distribution" << endl;
            int numbmodal=this->GetNumberModalities();
            SizeParameters=(int)(numbmodal*(numbmodal+1)); // Size of the mean vector + size of the variance matrix (numbmodal)
            //cout<< "Size parameters is "<<SizeParameters;
            /* REMARK :
             this might be changed if decided then to consider only number of free parameters, then the Gaussian Parameters calculation function must be changed consequently */
        }
            
    }
    return SizeParameters;
}

void TreeEM::CreateAllocateAndInitializeParameters(int DistributionTypeInput){
    if (this->GetNumberChildren()!=0) { // if there are children, it is necessarily a mixture
        cout<<"We have to modify Distribution type"<<endl;
        DistributionTypeInput=0;
    }
    if (this->GetParameters()==NULL) {
        //cout<<"Parameters are NULL at the momemt"<<endl;
        this->ParametersDistribution=new Parameters();
        this->ParametersDistribution->DistributionType=DistributionTypeInput;
        //cout<<this->GetDistributionType()<< "and " << this->ParametersDistribution->DistributionType;
        this->ParametersDistribution->SizeParameters=this->CalculateSizeParameters();
        if (this->GetSizeParameters()>0) {
            this->GetParameters()->ValueParameters=new float[this->GetSizeParameters()];//{0};
            for(int i=0;i<this->GetSizeParameters();i++){
                this->GetParameters()->ValueParameters[i]=0;
            }
        }
        else{
            this->GetParameters()->ValueParameters=NULL;
        }
        return;
    }
    this->ParametersDistribution->DistributionType=DistributionTypeInput;
    this->ParametersDistribution->SizeParameters=this->CalculateSizeParameters();
    //cout<<"the size of parameters is "<<this->GetSizeParameters()<<endl;
    if (this->GetSizeParameters()>0) {
        if(this->GetParametersValue()!=NULL){
            delete [] this->GetParametersValue();
            this->GetParameters()->ValueParameters=NULL;
        }
        int numbsp=this->GetSizeParameters();
        this->GetParameters()->ValueParameters=new float[this->GetSizeParameters()];//{0};
        for(int i=0;i<numbsp;i++){
            this->GetParameters()->ValueParameters[i]=0;
        }
    }
    else {
        if(this->GetParametersValue()!=NULL){
            delete [] this->GetParametersValue();
        }
        this->GetParameters()->ValueParameters=NULL;
    }
    return;
}

bool TreeEM::CheckForValidityOfParametersStructure(){
    // Incompatibility between type mixture and existence of values for the parameters
    
    if(this->GetDistributionType()==0 && this->GetParametersValue()!=NULL){
        return 0;
    }
    // Incompatibility between mixture DistributionType and SizeParameters > 0
    if (this->GetDistributionType()==0 && this->GetSizeParameters()>0) {
        return 0;
    }
    // Incompatibility between simple distribution and SizeParameters=0
    if (this->GetDistributionType()>0 && this->GetSizeParameters()==0) {
        return 0;
    }
    // Incompatibility between simple distribution and non allocation of values for parameters
    if(this->GetDistributionType()>0 && this->GetParametersValue()==NULL){
        return 0;
    }
    // Incompatibility if SizeParameters is 0 and ValueParameters is allocated
    if (this->GetSizeParameters()==0 && this->GetParametersValue()!=NULL) {
        return 0;
    }
    // Incompatibility if SizeParameters is >0 and ValueParameters is not allocated
    if (this->GetSizeParameters()>0 && this->GetParametersValue()==NULL) {
        return 0;
    }
    else {
        return 1;
    }
}

// Check for the compatibility of the data and the size of the parameters to set
bool TreeEM::CheckForSizeParametersValidity(){
    // First check for the validity of the structure of ParametersToSet
    if(!CheckForValidityOfParametersStructure()){
        return 0;
    }
    // Once checked for the validity of the structure check if mixture type or something else (if mixture, ok by definition)
    if(this->GetDistributionType()==0){
        return 1;
    }
    if (this->GetSizeParameters()!=this->CalculateSizeParameters()){
        return 0;
    }
    return 1;
}

bool TreeEM::ArePriorsNormalised(){
    int numel=this->GetNumberElements();
    if (this->IsPriorsVectorFilledWithNULL()) {
        return 1;
    }
    int numbclasses=this->GetPriorsVector().size();
    /*Check if some or all of them are NULL pointers in the vector if not all of them, must be normalised (all must be put to NULL)*/
    int indexFirstNull=numbclasses;
    for(int c=numbclasses-1;c>=0;c--){
        if (this->GetPriorsVector()[c]==NULL){
            indexFirstNull=c;
        }
    }
    if(indexFirstNull<numbclasses){// Means one of the pointer is NULL but one at least the first is not NULL
        cout<<"One of the Pointer to priors is NULL and is not the first one"<<endl;
        return 0;
    }
    if(indexFirstNull==0){// the first pointer is NULL, we must check for all the others in the vector
        for (int c=0;c<numbclasses;c++){
            if (this->GetPriorsVector()[c]!=NULL){
                cout<<"One of the Priors is not NULL whereas the first one is NULL"<<endl;
                return 0;
            }
        }
        return 1; // All the pointers in the vector are NULL
    }
    /* if get there : all Priors pointers in the vector correspond to valid pointer to nifti images*/
    vector<float *>PriorsVectorData_PTR;
    //cout << "All pointers correspond to a valid nifti_image"<<endl;
    for (int c=0;c<numbclasses;c++){
        PriorsVectorData_PTR.push_back(static_cast<float *>(this->GetPriorsVector()[c]->data));
    }
    float sumCheck=0;
    for (int i=0;i<numel;i++){
        sumCheck=0;
        for(int c=0;c<numbclasses;c++){
            sumCheck+=*PriorsVectorData_PTR[c];
            PriorsVectorData_PTR[c]++;
        }
        if(fabs(sumCheck-1)>1E-6){
            //cout<<"The priors are not well normalised, the difference to 1 is "<<sumCheck-1<<endl;
            
            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsPriorsVectorFilledWithNULL(){
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetPriors()!=NULL) {
            return 0;
        }
    }
    return 1; // By default return 1 if it is a leaf or an initial Node we are looking at
}

bool TreeEM::IsOneOfThePriorsNULL(){
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetPriors()==NULL) {
            return 1;
        }
    }
    return 0;
}

void TreeEM::NormalisePriors(){
    int numbchild=this->GetNumberChildren();
    int numel=this->GetNumberElements();
    
    if(this->ArePriorsNormalised()){
        cout<<"Priors are already normalised";
        return;
    }
    else{
        
        // First if one of them is a pointer to NULL, all of the priors are deleted and set to NULL
        if( this->IsOneOfThePriorsNULL()){
            int numbchild=this->GetNumberChildren();
            for (int c=0; c<numbchild; c++) {
                if(this->GetChild(c)->Priors!=NULL){
                    delete this->GetChild(c)->Priors;
                    this->GetChild(c)->Priors=NULL;
                }
            }
        }
        else{ // Priors are not normalised but all pointing to valid nifti images
        vector<float *> PriorsData_PTR;
        for (int c=0; c<numbchild; c++) { // Pushing pointers to begin of each of the Priors in the vector PriorsData_PTR
            if (!this->GetChild(c)->IsPriorsProbabilityType()) {
                cout<<"Priors is not of probability type";
                this->GetChild(c)->MakePriorsProbabilityType();
                cout << "and now is of probability type "<<this->GetChild(c)->IsPriorsProbabilityType();
            }
            //cout<< "Priors already of probability type"<<endl;
            PriorsData_PTR.push_back(static_cast<float*>(this->GetChild(c)->GetPriors()->data));
        }
            int * L2S_PTR=this->GetL2S(); // Pointer to beginning of L2S
            float tmpSum=0;
            for (int i=0; i<numel; i++,L2S_PTR++) { // for each of the voxels
                
                    tmpSum=0;
                    for (int c=0; c<numbchild; c++) { // calculation of the normalising factor
                        tmpSum+=(*PriorsData_PTR[c]);
                    }
                    for (int c=0; c<numbchild; c++) { // Normalisation of the priors data
                        if(tmpSum>0){
                            (*PriorsData_PTR[c])/=tmpSum;
                        }
                        else{ // taking care of breaking case where normalisation factor is 0
                            //cout<<"Normalisation factor is 0"<<endl;
                            (*PriorsData_PTR[c])=(float)1.0/numbchild;
                        }
                    }
                for (int c=0; c<numbchild; c++) { // incrementation of the pointer
                    PriorsData_PTR[c]++;
                }
            }
        }
        return;
    }
}

// Check if the priors is a probability type Prior (everything float between 0 and 1)
bool TreeEM::IsPriorsProbabilityType(){
    if (this->GetPriors()==NULL) {
        return 1;
    }
    else{
        int numel=this->GetNumberElements();
        float * Priors_PTR=static_cast<float*>(this->GetPriors()->data);
        float maxPriors=-1E32;
        float minPriors=1E32;
        for (int i=0; i<numel; i++,Priors_PTR++) {
            if ((*Priors_PTR)>maxPriors) {
                maxPriors=(*Priors_PTR);
            }
            if ((*Priors_PTR)<minPriors) {
                minPriors=(*Priors_PTR);
            }
        }
        bool testminPriors=(minPriors>-1E-6);
        bool testmaxPriors=(maxPriors<1+1E-6);
        if (testmaxPriors&&testminPriors) {
            //cout<<"Truly of probability type"<<endl;
            return 1;
        }
    }
    return 0;
}

bool TreeEM::CheckForValidityOfParametersStructure(Parameters * ParametersToCheck){
    // Incompatibility between type mixture and existence of values for the parameters
    if(ParametersToCheck->DistributionType==0 && ParametersToCheck->ValueParameters!=NULL){
        return 0;
    }
    // Incompatibility between mixture DistributionType and SizeParameters > 0
    if (ParametersToCheck->DistributionType==0 && ParametersToCheck->SizeParameters>0) {
        return 0;
    }
    // Incompatibility between simple distribution and SizeParameters=0
    if (ParametersToCheck->DistributionType>0 && ParametersToCheck->SizeParameters==0) {
        return 0;
    }
    // Incompatibility between simple distribution and non allocation of values for parameters
    if(ParametersToCheck->DistributionType>0 && ParametersToCheck->ValueParameters==NULL){
        return 0;
    }
    // Incompatibility if SizeParameters is 0 and ValueParameters is allocated
    if (ParametersToCheck->SizeParameters==0 && ParametersToCheck->ValueParameters!=NULL) {
        return 0;
    }
    // Incompatibility if SizeParameters is >0 and ValueParameters is not allocated
    if (ParametersToCheck->SizeParameters>0 && ParametersToCheck->ValueParameters==NULL) {
        return 0;
    }
    else {
        return 1;
    }
}

// Check for the compatibility of the data and the size of the parameters to set
bool TreeEM::CheckForSizeParametersValidity( Parameters * ParametersToSet){
    // First check for the validity of the structure of ParametersToSet
    if(!CheckForValidityOfParametersStructure(ParametersToSet)){
        return 0;
    }
    // Once checked for the validity of the structure check if mixture type or something else (if mixture, ok by definition)
    if(ParametersToSet->DistributionType==0){
        return 1;
    }
    if (ParametersToSet->SizeParameters!=this->CalculateSizeParameters(ParametersToSet->DistributionType)){
        return 0;
    }
    return 1;
}


/***************************** CALCULATION OF KLD ********************************/

// Update the histogram of the data weighted with NormResp for each of the modality in a vector of float pointer to each of the monomodal histogram
vector<float *> TreeEM::GetDataHistogram(){
    vector<float *> DataHistogram_VEC;
    // First check if there is Data to make histogram
    if (this->GetDataImage()==NULL) {
        cout<<"No data to make histogram from"<<endl;
        return DataHistogram_VEC;
    }
    
    // Then if modality wanted in histogram is available
    int numbbins=this->GetNumbbins();
    if (numbbins<=0) {
        cout<<"No proper number of bins asked for"<<endl;
        return DataHistogram_VEC;
    }
    float sizeBin=(float)1.0/numbbins;
    int *L2S_PTR=this->GetL2S();
    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
//    int numbchild=this->GetNumberChildren();
    float tmpValue=0;
    float * Data=this->MakeDataBFCorrected();
    //float * Data_PTR =static_cast<float*>(this->GetDataImage()->data);
    //vector<float *> DataHistogram_VEC;
    for (int m=0; m<numbmodal; m++) {
        float * InitDataHistToPushBack=new float[numbbins];
        for(int i=0;i<numbbins;i++){
            InitDataHistToPushBack[i]=0;
        }
        DataHistogram_VEC.push_back(InitDataHistToPushBack);
    }
    float * NormResp_PTR=this->GetNormResp();
    for (int m=0; m<numbmodal; m++) {
        float * DataHistogram_PTR=DataHistogram_VEC[m];
        //float * DataHistogram=new float[numbbins]{0};
        
        L2S_PTR=this->GetL2S();
        
        //Data_PTR=&static_cast<float*>(this->GetDataImage()->data)[m*numel];
        float * Data_PTR=&Data[m*numel];
        NormResp_PTR=this->GetNormResp();
        float SumNormResp=0;
        tmpValue=0;
        // Case when NormResp is NULL
        if (this->GetNormResp()==NULL) { // then consider that there are 1 everywhere
            for (int i=0; i<numel; i++, Data_PTR++,L2S_PTR++) {
                if ((*L2S_PTR)>=0) {
                    tmpValue=(float)*Data_PTR/sizeBin-0.5;
                    if (tmpValue<0) {
                        DataHistogram_PTR[0]+=1;
                    }
                    else if (tmpValue>numbbins-1){
                        DataHistogram_PTR[numbbins-1]+=1;
                    }
                    else{
                        DataHistogram_PTR[(int)floorf(tmpValue)]+=(1-(tmpValue-floorf(tmpValue)));
                        DataHistogram_PTR[(int)floorf(tmpValue)+1]+=(tmpValue-floorf(tmpValue));
                    }
                    SumNormResp++;
                    //NormResp_PTR++;
                }
            }
            
        }
        else{
            for (int i=0; i<numel; i++, Data_PTR++,L2S_PTR++) {
                if ((*L2S_PTR)>=0) {
                    tmpValue=(float)*Data_PTR/sizeBin-0.5;
                    if (tmpValue<0) {
                        DataHistogram_PTR[0]+=*NormResp_PTR;
                    }
                    else if (tmpValue>numbbins-1){
                        DataHistogram_PTR[numbbins-1]+=*NormResp_PTR;
                    }
                    else{
                        DataHistogram_PTR[(int)floorf(tmpValue)]+=*NormResp_PTR*(1-(tmpValue-floorf(tmpValue)));
                        DataHistogram_PTR[(int)floorf(tmpValue)+1]+=*NormResp_PTR*(tmpValue-floorf(tmpValue));
                    }
                    SumNormResp+=*NormResp_PTR;
                    NormResp_PTR++;
                }
            }

        }
        DataHistogram_PTR=DataHistogram_VEC[m];
        //int numelmasked=this->GetNumberMaskedElements();
        float SumData=0;
        for (int i=0; i<numbbins; i++,DataHistogram_PTR++) {
            *DataHistogram_PTR/=SumNormResp;
            SumData+=*DataHistogram_PTR;
        }
        //cout<<"SumData in GetDataHistogram is"<<SumData;
    }
    if (Data!=NULL) {
        delete [] Data;
        Data=NULL;
    }
    return DataHistogram_VEC;
}

/* Update the DataHistogram with all modalities weighted by NormResp.
 At the root, it will be the true histogram of the masked image*/
float * TreeEM::GetDataHistogramTotal(){
    // First check if there is Data to make histogram
    if (this->GetDataImage()==NULL) {
        cout<<"No data to make histogram from"<<endl;
        return NULL;
    }
    // then check that the number of bins set for the tree is compatible with 
    int numbbins=this->GetNumbbins();
    if (numbbins<=0) {
        return NULL;
    }
    
    else{
        int numbmodal=this->GetNumberModalities();
//        int numbchild=this->GetNumberChildren();
        int numel=this->GetNumberElements();
        //float * DataHistogram=new float[(int)powf(numbbins,numbmodal)]{0};
        float * DataHistogramTotal_PTR=new float[(int)powf(numbbins, numbmodal)];//{0};
        for(int i=0;i<(int)powf(numbbins,numbmodal);i++){
            DataHistogramTotal_PTR[i]=0;
        }
        float sizeBin=(float)1.0/numbbins;
        int *L2S_PTR=this->GetL2S();

        vector<float *> Data_PTRModalVector;
        float * Data=this->MakeDataBFCorrected();
        for (int m=0; m<numbmodal; m++) {
            //Data_PTRModalVector.push_back(&static_cast<float*>(this->GetDataImage()->data)[m*numel]);
            Data_PTRModalVector.push_back(&Data[m*numel]);
        }
        float * NormResp_PTR=this->GetNormResp();
        float * tmpValue=new float[numbmodal];//{0};
        for(int i=0;i<numbmodal;i++){
            tmpValue[i]=0;
        }
        int * indexTab=new int[(int)powf(2, numbmodal)];//{0};
        for(int i=0;i<(int)powf(2,numbmodal);i++){
            indexTab[i]=0;
        }
        float * valuePerc=new float[(int)powf(2, numbmodal)];//{1};
        for(int i=0;i<(int)powf(2,numbmodal);i++){
            valuePerc[i]=1;
        }
        int CountNormRespZero=0;
        float SumNormResp=0;
        for (int i=0; i<numel; i++,L2S_PTR++) {
            if ((*L2S_PTR>=0)) {
                for (int m=0; m<numbmodal; m++) {
                    tmpValue[m]=*Data_PTRModalVector[m]/sizeBin-0.5;
                    tmpValue[m]=tmpValue[m]>(numbbins-1)?numbbins-1:tmpValue[m];
                    tmpValue[m]=tmpValue[m]<0?0:tmpValue[m];
                    for (int c=0; c<(int)powf(2, m); c++) {
                        indexTab[c+(int)(powf(2, m))]=indexTab[c];
                        valuePerc[c+(int)(powf(2, m))]=valuePerc[c];
                    }
                    for (int c=0; c<(int)powf(2, m); c++) {
                        indexTab[c]+=(int)powf(numbbins, m)*floorf(tmpValue[m]);
                        indexTab[c+(int)powf(2, m)]+=(int)powf(numbbins, m)*floorf(tmpValue[m])+1;
                        valuePerc[c]*=(1-tmpValue[m]+floorf(tmpValue[m]));
                        valuePerc[c+(int)powf(2, m)]*=(tmpValue[m]-floorf(tmpValue[m]));
                    }
                }
                float sumValuePerc=0;
                for (int c=0; c<(int)powf(2, numbmodal); c++) {
                    sumValuePerc+=valuePerc[c];
                }
                //cout<<"sumValuePerc "<< sumValuePerc<<" for index "<<i<<endl;
                if(NormResp_PTR!=NULL){
                    
                for (int p=0; p<(int)powf(2, numbmodal); p++) {
                        DataHistogramTotal_PTR[indexTab[p]]+=valuePerc[p]*(*NormResp_PTR);
                        
                    }
                    if (*NormResp_PTR<1-1E-6) {
                        CountNormRespZero++;
                    }
                    SumNormResp+=*NormResp_PTR;
                    NormResp_PTR++;
                }
                    else{
                        for(int p=0; p<(int)powf(2, numbmodal); p++) {
                        DataHistogramTotal_PTR[indexTab[p]]+=valuePerc[p];
                            
                    }
                   SumNormResp++; 
                }
            }
            for (int p=0; p<(int)powf(2, numbmodal); p++) {
                indexTab[p]=0;
                valuePerc[p]=1;
            }
            for (int m=0; m<numbmodal; m++) {
                Data_PTRModalVector[m]++;
            }
        }
        //cout<<"CountNormRespZero is"<<CountNormRespZero<<endl;
        delete [] indexTab;
        indexTab=NULL;
        delete [] valuePerc;
        valuePerc=NULL;
        delete [] tmpValue;
        tmpValue=NULL;
        float SumDataHistogram=0;
//        float NumberDataHistogramZero=0;
//        int numelmasked=this->GetNumberMaskedElements();
        for (int i=0; i<(int)powf(numbbins,numbmodal); i++) {
            DataHistogramTotal_PTR[i]/=SumNormResp;
            SumDataHistogram+=DataHistogramTotal_PTR[i];
        }
        
        //cout<<"Sum DataHistogram is "<<SumDataHistogram<<endl;
        if (Data!=NULL) {
            delete [] Data;
            Data=NULL;
        }
        return DataHistogramTotal_PTR;
    }
}

// Returns the histogram of the distribution for numbbins
vector<float *> TreeEM::GetDistHistogram(){
    vector<float *> DistHistogramVector;
    // First check if there is Data to make histogram
    if (this->GetDataImage()==NULL) {
        cout<<"No data to make histogram from"<<endl;
        return DistHistogramVector;
    }
    // Then if modality wanted in histogram is available
    //vector<float *> DistHistogramVector;
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    int numbbins=this->GetNumbbins();
    if (this->GetNumberChildren()==0) {/* If it is a leaf the distribution is then calculated according to the parameters stored. For the moment, their is only one choice of simple distribution : the Gaussian distribution*/
        switch (this->GetDistributionType()) {
            default:
                vector<float*>ChildHistogram_VEC=this->GetGaussianDistributionHist();
                for (int m=0; m<numbmodal; m++) {
                    DistHistogramVector.push_back(ChildHistogram_VEC[m]);
                }
                break;
        }
    }
    else{/* It is not a leaf and the distribution is the weighted sum of the distributions below*/
        float * ChildHistogram_PTR;
        float * DistHistogram_PTR;
     
        float NormWeightChild=1;
        for (int m=0; m<numbmodal; m++) {
            float * InitDistHistToPushBack=new float[numbbins];
            for(int i=0;i<numbbins;i++){
                InitDistHistToPushBack[i]=0;
            }
            DistHistogramVector.push_back(InitDistHistToPushBack);
        }
        for (int c=0; c<numbchild; c++) {
            vector<float*>ChildHistogram_VEC=this->GetChild(c)->GetDistHistogram();
            for (int m=0; m<numbmodal; m++) {
                DistHistogram_PTR=DistHistogramVector[m];
                ChildHistogram_PTR=ChildHistogram_VEC[m];
                NormWeightChild=this->GetChild(c)->GetNormWeight();
                for (int i=0; i<numbbins; i++,DistHistogram_PTR++,ChildHistogram_PTR++) {
                    *DistHistogram_PTR+=NormWeightChild*(*ChildHistogram_PTR);
                }
                // ClearingMemory
                if(ChildHistogram_VEC[m]!=NULL){
                    delete ChildHistogram_VEC[m];
                    ChildHistogram_VEC[m]=NULL;
                    cout<< " Delete done "<<endl;
                }
//                if (ChildHistogram_VEC[m]!=NULL && ((*ChildHistogram_VEC[m])!= NULL) ) {
//                                    delete [] ChildHistogram_VEC[m];
//                                    cout<< "Delete ChildHistogram_VEC partially done for modality "<<m<<" and child "<<c<<endl;
//                                   ChildHistogram_VEC[m]=NULL;
//                }
//                if (ChildHistogram_VEC[m]==NULL) {
//                    cout<<"ChildHistogram_VEC pointer already NULL"<<endl;
//                }
//                if ((*ChildHistogram_VEC[m])==NULL) {
//                    cout<<"Pb with memory allocation"<<endl;
//                }

            }
        }
    }
    return DistHistogramVector;
}

// Returns the Gaussian Distribution for each modality independently in a vector of float pointers
vector<float *> TreeEM::GetGaussianDistributionHist(){
    vector<float *> GaussianDistributionHist;
    if (this->GetNumberChildren()!=0) {
        cout<< "not leaf so cannot return gaussian distribution"<<endl;
        return GaussianDistributionHist;
    }
    if(this->GetDistributionType()!=1){
        cout<<"no gaussian distribution type"<<endl;
        return GaussianDistributionHist;
    }
    //vector<float *> GaussianDistributionHist;
    int numbmodal=this->GetNumberModalities();
    
    float * Variance=this->GetVariance();
    float * Mean=this->GetMean();
    float Value=0;
    float sizeBin=1.0/(float)numbbins;
    for (int m=0; m<numbmodal; m++) {
        float * GaussianDist_tmp=new float[numbbins];
        for(int i=0;i<numbbins;i++){
            GaussianDist_tmp[i]=0;
        }
        float * GaussianDist_tmpPTR=GaussianDist_tmp;
        float NormalisationFactor=1.0/powf(2*M_PI*Variance[m+m*numbmodal], 0.5);
        GaussianDist_tmpPTR=GaussianDist_tmp;
        float sumDist=0;
        for (int i=0; i<numbbins; i++,GaussianDist_tmpPTR++) {
            Value=i*sizeBin+0.5*sizeBin;
            *GaussianDist_tmpPTR=NormalisationFactor*expf(-(float)0.5/Variance[m+m*numbmodal]*(Value-Mean[m])*(Value-Mean[m]))*sizeBin;
            sumDist+=*GaussianDist_tmpPTR;
        }
        //cout<<"sumDist in vectorial Gaussian hist is "<<sumDist<<endl;
        GaussianDistributionHist.push_back(GaussianDist_tmp);
    }
//    delete [] GaussianDist_tmp;
//    GaussianDist_tmp=NULL;
    return GaussianDistributionHist;
}

// Returns the complete histogram corresponding to the distribution with the given number of bins. To be used for the derivation of the KLD
float * TreeEM::GetGaussianDistributionHistTotal(){
    if (this->GetNumberChildren()!=0) {
        cout<< "not leaf so cannot return gaussian distribution"<<endl;
        return NULL;
    }
    if(this->GetDistributionType()!=1){
        cout<<"no gaussian distribution type"<<endl;
        return NULL;
    }
    
    int numbmodal=this->GetNumberModalities();
    
    // First initialise the Variance Matrix and inverted one
    matrix<float> VarianceMatrix=matrix<float>(numbmodal);
    for (int m1=0; m1<numbmodal; m1++) {
        for (int m2=0; m2<numbmodal; m2++) {
            VarianceMatrix.setvalue(m1, m2, this->GetVariance()[m1+m2*numbmodal]);
        }
    }
    
    // Calculation of the factor in front of the exponential
    float DeterminantVariance=VarianceMatrix.determinant();
    float * Mean=this->GetMean();
    float sizeBin=1.0/(float)numbbins;
    //cout<<"the Variance determinant is "<<DeterminantVariance<<endl;
    float NormalisationFactor=1.0/(float)(powf(2*M_PI , (float)((float)numbmodal/2.0))*powf(DeterminantVariance, 0.5));
    //cout <<"The normalisation factor is "<<NormalisationFactor<<endl;
    //Initialisation of the needed element to calculate the inside of the exponential
    float * GaussianDistributionTotal = new float[(int)powf(numbbins, numbmodal)];//{0};
    for(int i=0;i<(int)powf(numbbins,numbmodal);i++){
        GaussianDistributionTotal[i]=0;
    }
    float * GaussianDistribution_PTR=GaussianDistributionTotal;
    float * InvertedVariance=new float[numbmodal*numbmodal];
    if (numbmodal>1) {
        VarianceMatrix.invert();
        bool success;
        for(int m1=0;m1<numbmodal;m1++){
            for (int m2=0; m2<numbmodal; m2++) {
                VarianceMatrix.getvalue(m1, m2, InvertedVariance[m1+m2*numbmodal], success);
            }
        }
    }
    else{
        bool success;
        VarianceMatrix.getvalue(0, 0, InvertedVariance[0], success);
        *InvertedVariance=(float)1.0/DeterminantVariance;
        //*InvertedVariance=1/(*InvertedVariance);
    }
    float temp;
    int * IndexConversion=new int [numbmodal];//{0};
    for(int i=0;i<numbmodal;i++){
        IndexConversion[i]=0;
    }
    int Rem;
    // Calculation of the inside of the exponential
    //cout<<InvertedVariance[0]<<endl;
    int CountNegativeValues=0;
    for (int i=0; i<(int)powf(numbbins,numbmodal); i++) {
        Rem=i;
        
        // Conversion of the index i into the different index for each modality giving then the value to consider
        for (int m=numbmodal-1; m>=0; m--) {
            IndexConversion[m]=Rem/(int)powf(numbbins, m);
            Rem-=IndexConversion[m]*(int)powf(numbbins, m);
        }
            for (int m1=0; m1<numbmodal; m1++) {
                temp=0;
                for (int m2=0; m2<numbmodal; m2++) {
                    temp+=((sizeBin*IndexConversion[m2]+0.5*sizeBin)-Mean[m2])*InvertedVariance[m2+m1*numbmodal];
                }
                *GaussianDistribution_PTR=(*GaussianDistribution_PTR)+temp*((sizeBin*IndexConversion[m1]+0.5*sizeBin)-Mean[m1]);
            }
        if (*GaussianDistribution_PTR<0) {
            CountNegativeValues++;
        }
            GaussianDistribution_PTR++;
        }
    //cout<<"The number of negative values is "<<CountNegativeValues<<endl;
    // Filling the distribution array with the value
    GaussianDistribution_PTR=GaussianDistributionTotal;
    //cout << "Normalisation factor in GetGaussianDistHistTotal is "<<NormalisationFactor<<endl;
    float sumDist=0;
//    int numelmasked=this->GetNumberMaskedElements();
    for (int i=0;i<(int)powf(numbbins,numbmodal);i++,GaussianDistribution_PTR++){
        *GaussianDistribution_PTR=NormalisationFactor*exp(-0.5*(*GaussianDistribution_PTR))*powf(sizeBin,numbmodal);
        sumDist+=*GaussianDistribution_PTR;
    }
    //cout<<"sumDist in GetGaussianTotal is "<<sumDist<<endl;
    // Clearing space needing for inverse of variance
    if(InvertedVariance!=NULL){
        delete [] InvertedVariance;
        InvertedVariance=NULL;
    }

    if (IndexConversion!=NULL) {
        delete [] IndexConversion;
        IndexConversion=NULL;
    }
    
    // Return result
    return GaussianDistributionTotal;
}

float * TreeEM::GetDistHistogramTotal(){

// Then if modality wanted in histogram is available
int numbchild=this->GetNumberChildren();
int numbmodal=this->GetNumberModalities();
int numbbins=this->GetNumbbins();
float * DistHistogramTotal=new float[(int)powf(numbbins, numbmodal)];//{0};
for(int i=0;i<(int)powf(numbbins,numbmodal);i++){
    DistHistogramTotal[i]=0;
}

if (this->GetNumberChildren()==0) {/* If it is a leaf the distribution is then calculated according to the parameters stored. For the moment, their is only one choice of simple distribution : the Gaussian distribution*/
    switch (this->GetDistributionType()) {
        default:
            float * LeafHistogramTotal=this->GetGaussianDistributionHistTotal();
//            vector<float *> LeafVectorHistogram=this->GetGaussianDistributionHist();
//            float * LeafVectorHistogram_PTR=LeafVectorHistogram[0];
            float * LeafHistogramTotal_PTR=LeafHistogramTotal;
//            int CountDifference=0;
//            for (int i=0; i<numbbins; i++,LeafHistogramTotal_PTR++,LeafVectorHistogram_PTR++) {
//                if (fabsf((*LeafHistogramTotal_PTR)-(*LeafVectorHistogram_PTR))>1E-6) {
//                    CountDifference++;
//                }
//            }
//            //cout<<"Count difference direct gaussian is "<<CountDifference<<endl;
//            LeafHistogramTotal_PTR=LeafHistogramTotal;
            float SumLeaf=0;
            for (int i=0; i<(int)powf(numbbins,numbmodal); i++,LeafHistogramTotal_PTR++) {
                DistHistogramTotal[i]=*LeafHistogramTotal_PTR;
                //SumLeaf+=*LeafHistogramTotal_PTR;
            }
            //cout << "Sum leaf is "<<SumLeaf<<endl;
            delete [] LeafHistogramTotal;
            LeafHistogramTotal=NULL;
            break;
    }
}
else{/* It is not a leaf and the distribution is the weighted sum of the distributions below*/
    
    float * DistHistogram_PTR=DistHistogramTotal;
    float * ChildDistHistogram;
    float * ChildDistHistogram_PTR;
    float NormWeightChild;
    for (int c=0; c<numbchild; c++) {
        
    ChildDistHistogram=(this->GetChild(c))->GetDistHistogramTotal();
         ChildDistHistogram_PTR=ChildDistHistogram;
        DistHistogram_PTR=DistHistogramTotal;
        NormWeightChild=this->GetChild(c)->GetNormWeight();
        for (int i=0; i<(int)powf(numbbins, numbmodal); i++) {
            *DistHistogram_PTR+=NormWeightChild*(*ChildDistHistogram_PTR);
            DistHistogram_PTR++;
            ChildDistHistogram_PTR++;
        }
        // Clearing
        ChildDistHistogram_PTR=NULL;
        if (ChildDistHistogram!=NULL) {
            delete [] ChildDistHistogram;
            ChildDistHistogram=NULL;
            cout<< "Delete done for child "<<c<<endl;
        }
//        if ((ChildDistHistogram!=NULL) && ((*ChildDistHistogram)!= NULL)) {
//        
//            delete [] ChildDistHistogram;
//            ChildDistHistogram=NULL;
//            cout<<"Delete done for child "<<c<<endl;
//        }
//        if (ChildDistHistogram==NULL) {
//            cout<<"ChildDistHistogram is NULL"<<endl;
//        }
//        cout<<ChildDistHistogram<<" and "<< *ChildDistHistogram<<endl;
        
        
       }
    }
    return DistHistogramTotal;
}

// Returns the KLD when comparing the total histogram of the distribution modelled and the data
float TreeEM::GetDistCompKLDTotal(){
    float * DistHistTotal=this->GetDistHistogramTotal();
    float * DistHistTotal_PTR=DistHistTotal;
    float * DataHistTotal=this->GetDistHistogramTotal();
    float * DataHistTotal_PTR=DataHistTotal;
    float KLD=0;
    int numbmodal=this->GetNumberModalities();
    for (int i=0; i<(int)powf(numbbins, numbmodal); i++,DataHistTotal_PTR++,DistHistTotal_PTR++) {
        KLD+=logf(((*DataHistTotal_PTR)/(*DistHistTotal_PTR)))*(*DataHistTotal_PTR);
    }
    delete [] DistHistTotal;
    DistHistTotal=NULL;
    delete [] DataHistTotal;
    DataHistTotal = NULL;
    return KLD;
}

float * TreeEM::GetDistCompKLD(){
    int numbmodal=this->GetNumberModalities();
    float * KLD=new float[numbmodal];//{0};
    for(int i=0;i<numbmodal;i++){
        KLD[i]=0;
    }
    vector<float *> DistHist_VEC=this->GetDistHistogram();
    vector<float *> DataHist_VEC=this->GetDataHistogram();
    for (int m=0; m<numbmodal; m++) {
        float * DistHist_PTR=DistHist_VEC[m];
        float * DataHist_PTR=DataHist_VEC[m];
        for (int i=0; i<(int)powf(numbbins, numbmodal); i++, DataHist_PTR++,DistHist_PTR++) {
            KLD[m]+=logf(((*DataHist_PTR)/(*DistHist_PTR)))*(*DataHist_PTR);
        }
        delete [] DistHist_VEC[m];
        DistHist_VEC[m]=NULL;
        delete [] DataHist_VEC[m];
        DataHist_VEC[m]=NULL;
    }
    return KLD;
}

float TreeEM::CriterionCalculation(){
    int numbFreeParams=this->GetNumberFreeParameters();
    int numel=this->GetNumberElements();
    float IndFactor=this->GetIndFactorTot();
    float LL=this->GetLogLikelihood();
    float BIC=LL*IndFactor-(float)numbFreeParams/2.0*logf(numel*IndFactor);
    return BIC;
}

/************************** METHODS FOR MERGING OPERATIONS ******************************/

// Returns the comparison between the two modelled distributions if comparable meaning that the two trees compared have the same Parent
MergeKLD * TreeEM::GetKLDMerge(int Child1Input,int Child2Input){
    // Check that the given Children exist
    int numbchild=this->GetNumberChildren();
    if (Child1Input>=numbchild||Child2Input>=numbchild||Child1Input<0||Child2Input<0) {
        cout<<"Out of bound for the children to merge"<<endl;
        return NULL;
    }
    // Check if the merging try is made on leaves
    if (this->GetChild(Child1Input)->GetNumberChildren()!=0||this->GetChild(Child2Input)->GetNumberChildren()!=0) {
        cout<< "One of the children is not a leaf, merging cannot be done"<<endl;
        return NULL;
    }
    // Check if the merging is tried on leaves with same priors
//    cout<< this->GetChild(Child1Input)->GetPriors()<<" and "<<this->GetChild(Child2Input)->GetPriors();
    if (this->GetChild(Child1Input)->GetPriors()!=this->GetChild(Child2Input)->GetPriors()) {
        cout<<"Impossible to merge leaves with not same priors"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    // Initialisation of ResultMergeKLD
    MergeKLD * ResultMergeKLD=new MergeKLD;
    ResultMergeKLD->Parent=this;
    ResultMergeKLD->Child1=Child1Input;
    ResultMergeKLD->Child2=Child2Input;
    ResultMergeKLD->KLDTot=0;
    ResultMergeKLD->KLD=new float[numbmodal];//{0};
    for(int i=0;i<numbmodal;i++){
        ResultMergeKLD->KLD[i]=0;
    }
    
    // Calculation of the needed histograms for the calculation of the KLDcompTot
    float * DistHistTotal1=this->GetChild(Child1Input)->GetDistHistogramTotal();
    float * DistHistTotal1_PTR=DistHistTotal1;
    float * DistHistTotal2=this->GetChild(Child2Input)->GetDistHistogramTotal();
    float * DistHistTotal2_PTR=DistHistTotal2;
//    int CountNan1=0;
//    int CountNan2=0;
//    int CountNanAdd=0;
    for (int i=0; i<(int)powf(numbbins, numbmodal); i++,DistHistTotal1_PTR++,DistHistTotal2_PTR++) {
        float DistValue1=(*DistHistTotal1_PTR)<=0.00001?0.00001:(*DistHistTotal1_PTR);
        float DistValue2=(*DistHistTotal2_PTR)<=0.00001?0.00001:(*DistHistTotal2_PTR);
//        if (logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)!=logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)) {
//            CountNan1++;
//        }
//        if (logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR)!=logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR)) {
//            CountNan2++;
//        }
//        if (0.5*(logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)+logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR))!=0.5*(logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)+logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR))) {
//            CountNanAdd++;
//        }
        ResultMergeKLD->KLDTot+=0.5*(logf(((DistValue1)/(DistValue2)))*(*DistHistTotal1_PTR)+logf(((DistValue2)/(DistValue1)))*(*DistHistTotal2_PTR)); // Use of the symmetric KLD between the two distributions
        //cout<<ResultMergeKLD->KLDTot<< "   -   ";
    }
    cout<<endl;
    delete [] DistHistTotal1;
    DistHistTotal1=NULL;
    delete [] DistHistTotal2;
    DistHistTotal2=NULL;
    
    // Calculation of KLDcomp
    vector<float *> DistHist1_VEC=this->GetChild(Child1Input)->GetDistHistogram();
    vector<float *> DistHist2_VEC=this->GetChild(Child2Input)->GetDistHistogram();
    for (int m=0; m<numbmodal; m++) {
        float * DistHist1_PTR=DistHist1_VEC[m];
        float * DistHist2_PTR=DistHist2_VEC[m];
        for (int i=0; i<numbbins; i++,DistHist1_PTR++,DistHist2_PTR++) {
            float DistValue1=(*DistHist1_PTR)<=0.00001?0.00001:(*DistHist1_PTR);
            float DistValue2=(*DistHist2_PTR)<=0.00001?0.00001:(*DistHist2_PTR);
            ResultMergeKLD->KLD[m]+=0.5*(logf(((DistValue1/DistValue2)))*(*DistHist1_PTR)+logf(((DistValue2)/(DistValue1)))*(*DistHist2_PTR)); // Use of the symmetric KLD between the two distributions
        }
        delete [] DistHist1_VEC[m];
        DistHist1_VEC[m]=NULL;
        delete [] DistHist2_VEC[m];
        DistHist2_VEC[m]=NULL;
    }
    return ResultMergeKLD;
}

// Returns for all children of a TreeEM the pairwise MergeKLD in a vector;
vector<MergeKLD *> TreeEM::GetVectorKLDMergeChildren(){
    int numbchild=this->GetNumberChildren();
    vector<MergeKLD *> ResultCompChildrenKLD;
    for (int c1=0; c1<numbchild; c1++) {
        for (int c2=c1+1; c2<numbchild; c2++) {
            MergeKLD * MergeKLDToAdd=this->GetKLDMerge(c1, c2);
            if (MergeKLDToAdd!=NULL) {
                ResultCompChildrenKLD.push_back(MergeKLDToAdd);
            }
            
        }
    }
    return ResultCompChildrenKLD;
}

vector<MergeKLD*> TreeEM::GetVectorKLDMergeLeaves(){
    vector<MergeKLD *> VectorKLDMergeLeaves;
    vector<MergeKLD *> VectorKLDMergeChildren=this->GetVectorKLDMergeChildren();
    for (int i=0; i<VectorKLDMergeChildren.size(); i++) {
        VectorKLDMergeLeaves.push_back(VectorKLDMergeChildren[i]);
    }
    //VectorKLDMergeChildren.clear();
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        VectorKLDMergeChildren=(this->GetChild(c))->GetVectorKLDMergeLeaves();
        for (int i=0; i<VectorKLDMergeChildren.size(); i++) {
            VectorKLDMergeLeaves.push_back(VectorKLDMergeChildren[i]);
        }
        //VectorKLDMergeChildren.clear();
    }
    return VectorKLDMergeLeaves;
}

MergeKLD * TreeEM::GetToMerge(){
    vector<MergeKLD *> VectorKLDMergeLeaves=this->GetVectorKLDMergeLeaves();
    if (VectorKLDMergeLeaves.size()==0) {
        return NULL;
    }
    else{
        int numbmodal=this->GetNumberModalities();
        MergeKLD * MergeMin=new MergeKLD();
        MergeMin->Parent=NULL;
        MergeMin->Child1=-1;
        MergeMin->Child2=-1;
        MergeMin->KLDTot=1E32;
        MergeMin->KLD=new float[numbmodal];//{0};
        for(int i=0;i<numbmodal;i++){
            MergeMin->KLD[i]=0;
        }
        
        for (int i=0; i<VectorKLDMergeLeaves.size(); i++) {
            int Child1=VectorKLDMergeLeaves[i]->Child1;
            int Child2=VectorKLDMergeLeaves[i]->Child2;
            int numbDirectLeaves=VectorKLDMergeLeaves[i]->Parent->GetNumberDirectLeaves();
            float KLDTotTry=VectorKLDMergeLeaves[i]->KLDTot;
            float * KLDTry=VectorKLDMergeLeaves[i]->KLD;
            // Update if not already checked and is lowest KLDTot
            bool MergeCheckTest=VectorKLDMergeLeaves[i]->Parent->GetMergeCheck()[Child1+Child2*numbDirectLeaves];
            bool MergingPossibility=VectorKLDMergeLeaves[i]->Parent->GetChild(Child1)->GetPriors()==VectorKLDMergeLeaves[i]->Parent->GetChild(Child2)->GetPriors();
            if (MergingPossibility) {
//                cout<<"Possible merging"<<endl;
                if (!MergeCheckTest&& KLDTotTry<MergeMin->KLDTot) {
                    MergeMin->Parent=VectorKLDMergeLeaves[i]->Parent;
                    MergeMin->Child1=Child1;
                    MergeMin->Child2=Child2;
                    MergeMin->KLDTot=KLDTotTry;
                    for (int m=0; m<numbmodal; m++) {
                        MergeMin->KLD[m]=KLDTry[m];
                    }
                }
            }
            
            // clear space needed by ith element of VectorKLDMergeLeaves once it has been looked at
            delete VectorKLDMergeLeaves[i];
            VectorKLDMergeLeaves[i]=NULL;
        }
        // if MergeMin not modified <=> no possible merging return NULL
        if (MergeMin->Child1==-1) {
            delete MergeMin;
            return NULL;
        }
        else{
            return MergeMin;
        }
    }
}

void TreeEM::MergeOperation(int Child1, int Child2){
    int numbchild=this->GetNumberChildren();
    if (Child1>=numbchild||Child1<0||Child2<0||Child2>=numbchild) {
        cout<<"Index children out of bound, no merging possible"<<endl;
        return;
    }
    if ((this->GetChild(Child1))->GetNumberChildren()!=0|| (this->GetChild(Child2))->GetNumberChildren()!=0) {
        cout << "Nothing other than leaf can be merged"<<endl;
        return;
    }
    if (Child1==Child2) {
        cout<<"Impossible to merge the same child"<<endl;
        return;
    }
    if (Child1>Child2) {// Reorder Child index in order to properly erase them afterwards
        int tmp=Child1;
        Child1=Child2;
        Child2=tmp;
    }
    this->CreateAndAddChild(this->GetChild(Child1)->GetPriorsDirect(),1);
    int newnumbchild=this->GetNumberChildren();
    this->GetChild(newnumbchild-1)->NormWeight=(this->GetChild(Child1))->GetNormWeight()+(this->GetChild(Child2))->GetNormWeight();
    Parameters * ParametersMerge=this->ParametersForMerging(Child1,Child2);
    this->GetChild(newnumbchild-1)->SetParameters(ParametersMerge);
//    delete[] ParametersMerge->ValueParameters;
//    ParametersMerge->ValueParameters=NULL;
    delete ParametersMerge;
    ParametersMerge=NULL;
    delete this->GetChild(Child1);
    //this->GetChild(Child1)=NULL;
    delete this->GetChild(Child2);
    //this->GetChild(Child2)=NULL;
    //cout<<this->GetChildren().begin()+Child1<<" is the place to look at in the vector"<<endl;
    (this->Children).erase(this->Children.begin()+Child1);
    this->Children.erase(this->Children.begin()+Child2-1);
    newnumbchild=this->GetNumberChildren();
    if (newnumbchild==1) {
        this->CollapseOnlyChild();
    }

}

Parameters * TreeEM::ParametersForMerging(int Child1, int Child2){
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    if (Child1>=numbchild||Child1<0||Child2<0||Child2>=numbchild) {
        cout<<"Index children out of bound, no merging possible"<<endl;
        return NULL;
    }
    if ((this->GetChild(Child1))->GetNumberChildren()!=0|| (this->GetChild(Child2))->GetNumberChildren()!=0) {
        cout << "Nothing other than leaf can be merged"<<endl;
        return NULL;
    }
    if (Child1==Child2) {
        cout<<"Impossible to merge the same child"<<endl;
        return NULL;
    }
    if(this->GetChild(Child1)->GetDistributionType()!=this->GetChild(Child2)->GetDistributionType()){
        cout<< "Impossible to merge leaves of different distribution type"<<endl;
        return NULL;
    }
    Parameters * ParametersMerge=new Parameters;
    ParametersMerge->DistributionType=(this->GetChild(Child1))->GetDistributionType();
    ParametersMerge->SizeParameters=(this->GetChild(Child1))->GetSizeParameters();
    ParametersMerge->ValueParameters=new float[ParametersMerge->SizeParameters];//{0};
    for(int i=0;i<ParametersMerge->SizeParameters;i++){
        ParametersMerge->ValueParameters[i]=0;
    }
    // As checked before it is a leaf, necessary of simple distribution type so there is no problem with SizeParameters that is necessarily strictly positive
    switch (ParametersMerge->DistributionType) {
 
            
        default: // Gaussian case as default
            float * Mean1_PTR=(this->GetChild(Child1))->GetMean();
            float * Mean2_PTR=(this->GetChild(Child2))->GetMean();
            float * Variance1_PTR=(this->GetChild(Child1))->GetVariance();
            float * Variance2_PTR=(this->GetChild(Child2))->GetVariance();
            float Weight1=this->GetChild(Child1)->GetNormWeight();
            float Weight2=this->GetChild(Child2)->GetNormWeight();
            // Mean for ParametersMerge
            for (int m=0; m<numbmodal; m++) {
                ParametersMerge->ValueParameters[m]=(Weight1*Mean1_PTR[m]+Weight2*Mean2_PTR[m])/(Weight1+Weight2);
            }
            // Variance for ParametersMerge
            for (int m1=0; m1<numbmodal; m1++) {
                for (int m2=0; m2<numbmodal; m2++) {
                    ParametersMerge->ValueParameters[numbmodal+m1+m2*numbmodal]=(Weight1*(Variance1_PTR[m1+m2*numbmodal]+(Mean1_PTR[m1]-ParametersMerge->ValueParameters[m1])*(Mean1_PTR[m2]-ParametersMerge->ValueParameters[m2]))+Weight2*(Variance2_PTR[m1+m2*numbmodal]+(Mean2_PTR[m1]-ParametersMerge->ValueParameters[m1])*(Mean2_PTR[m2]-ParametersMerge->ValueParameters[m2])))/(Weight1+Weight2);
                }
            }
            break;
    }
    return ParametersMerge;
}

// When testing a S operation returns the best of the two models obtained.
TreeEM* TreeEM::RunMergeOperation(){
    MergeKLD * MergeTry=this->GetToMerge();
    if (MergeTry==NULL) {
        cout<<"No merging to try"<<endl;
        return this;
    }
    // Updating of the checking
    cout<< "Trying merging"<<endl;
    int numbDirectLeaves=(MergeTry->Parent)->GetNumberDirectLeaves();
    (MergeTry->Parent)->GetMergeCheck()[MergeTry->Child1+numbDirectLeaves*MergeTry->Child2]=1;
    (MergeTry->Parent)->GetMergeCheck()[MergeTry->Child2+numbDirectLeaves*MergeTry->Child1]=1;
    // Copy of the current tree
    TreeEM * CopiedTree=this->CopyTree(NULL);
    (MergeTry->Parent)->MergeOperation(MergeTry->Child1,MergeTry->Child2);
    this->CollapseOnlyChildTot();
    this->ClearSMChecks();
    
    // Clear memory for MergeTry
//    delete [] MergeTry->KLD;
//    MergeTry->KLD=NULL;
    delete MergeTry;
    MergeTry=NULL;
    
    // Partial EM first
//    float PartCLL=0;
//    float PartOldCLL=0;
//    int PartIteration=0;
    //(MergeTry->Parent)->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    // Run EM
    this->UpdateDistribution();
    this->UpdateNonNormResp();
    this->UpdateNormResp();
    this->SaveAllClasses((g_scratch_dir + "/TestMergeResult").c_str());
    this->UpdateNonNormWeights();
    this->UpdateNormWeights();
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration);
    this->SaveAllClasses((g_scratch_dir + "/TestMergeResult").c_str());
    // Test of the operation
    float OldCriterion=CopiedTree->CriterionCalculation();
    float NewCriterion=this->CriterionCalculation();
    if (OldCriterion>=NewCriterion) { // Criterion of change not met, delete the modified tree and return the copied one
        delete this;
        cout<<"No merging accepted"<<endl;
        return CopiedTree;
    }
    else{ // Criterion of change met delete the not modified tree stored in CopiedTree and return the new model.
        delete CopiedTree;
        cout<<"Merging accepted"<<endl;
        this->CollapseOnlyChildTot();
        this->ClearSMChecks();
        return this;
    }
}


/******** METHODS FOR SPLITTING OPERATIONS *************************/
SplitKLD * TreeEM::GetKLDSplit(int ChildInput){
    // First check if input integer is correct
    int numbchild=this->GetNumberChildren();
    if (ChildInput<0||ChildInput>=numbchild) {
        cout<<"Out of bound for child to split"<<endl;
        return NULL;
    }
    // then check if can be split <=> check if split a leaf
    if (this->GetChild(ChildInput)->GetNumberChildren()!=0) {
        cout<<"It is not a leaf so cannot be split"<<endl;
        return NULL;
    }
    SplitKLD * SplitKLDResult=new SplitKLD();
    SplitKLDResult->Parent=this;
    SplitKLDResult->ChildToSplit=ChildInput;
    SplitKLDResult->KLDTot=0;
    int numbmodal=this->GetNumberModalities();
    SplitKLDResult->KLD=new float[numbmodal];//{0};
    for(int i=0;i<numbmodal;i++){
        SplitKLDResult->KLD[i]=0;
    }
    
    // Calculating KLDTot and deleting the corresponding float arrays needed after use
    float * DistHistTotal=(this->GetChild(ChildInput))->GetDistHistogramTotal();
    float * DistHistTotal_PTR=DistHistTotal;
    float * DataHistTotal=(this->GetChild(ChildInput))->GetDataHistogramTotal();
    float * DataHistTotal_PTR=DataHistTotal;
    float SumDist=0;
    float SumData=0;
    

    
    for (int i=0; i<(int)powf(numbbins, numbmodal); i++,DistHistTotal_PTR++,DataHistTotal_PTR++) {
        SumDist+=*DistHistTotal_PTR;
        SumData+=*DataHistTotal_PTR;
        float DistValue=(*DistHistTotal_PTR)<=0.00001?0.00001:(*DistHistTotal_PTR);
        float DataValue=(*DataHistTotal_PTR)<=0.00001?0.00001:(*DataHistTotal_PTR);
        SplitKLDResult->KLDTot+=0.5*(logf(DataValue/DistValue)*(*DataHistTotal_PTR)+logf(DistValue/DataValue)*(*DistHistTotal_PTR)); // Use of the symmetric KLD between the data and the distribution
        if (SplitKLDResult->KLDTot!=SplitKLDResult->KLDTot) {
            cout<<"Pb at index "<<i<<" : "<<DataValue<<" and "<<DistValue<<endl;
        }
    }
    //cout<<SplitKLDResult->KLDTot<< " and SumData "<<SumData<< "and SumDist are "<<SumDist<<endl;

    
    
    // Calculation of KLDcomp
    vector<float *> DistHist_VEC=this->GetChild(ChildInput)->GetDistHistogram();
    vector<float *> DataHist_VEC=this->GetChild(ChildInput)->GetDataHistogram();
    DistHistTotal_PTR=DistHistTotal;
    DataHistTotal_PTR=DataHistTotal;
    
//    float * DiffDistHist=new float[numbbins];//{0};
//    for(int i=0;i<numbbins;i++){
//        DiffDistHist[i]=0;
//    }
//    float * DiffDataHist=new float[numbbins];//{0};
//    for(int i=0;i<numbbins;i++){
//        DiffDataHist[i]=0;
//    }

//    float * DiffDistHist_PTR=DiffDistHist;
//    float * DiffDataHist_PTR=DiffDataHist;
//    float * DistHist_PTR=DistHist_VEC[0];
//    float * DataHist_PTR=DataHist_VEC[0];
//    int CountNonZeroDiffDist=0;
//    int CountNonZeroDiffData=0;
//      SumData=0;
//    for (int i=0; i<numbbins; i++,DistHist_PTR++,DistHistTotal_PTR++,DataHistTotal_PTR++,DataHist_PTR++,DiffDataHist_PTR++,DiffDistHist_PTR++) {
        
//        *DiffDistHist_PTR=(*DistHist_PTR)-(*DistHistTotal_PTR);
//        *DiffDataHist_PTR=(*DataHist_PTR)-(*DataHistTotal_PTR);
//        SumData+=*DistHistTotal_PTR;
//        if (fabsf(*DiffDistHist_PTR)>1E-6) {
//            CountNonZeroDiffDist++;
//            //cout<<*DiffDistHist_PTR<<endl;
//        }
//        if (fabsf(*DiffDataHist_PTR)>1E-6) {
//            CountNonZeroDiffData++;
//        }
//    }
//    delete [] DiffDataHist;
//    DiffDataHist=NULL;
//    delete [] DiffDistHist;
//    DiffDistHist=NULL;
    
    for (int m=0; m<numbmodal; m++) {
        float * DistHist_PTR=DistHist_VEC[m];
        float * DataHist_PTR=DataHist_VEC[m];
        for (int i=0; i<numbbins; i++,DistHist_PTR++,DataHist_PTR++) {
            float DistValue=(*DistHist_PTR)<=0.00001?0.00001:(*DistHist_PTR);
            
            float DataValue=(*DataHist_PTR)<=0.00001?0.00001:(*DataHist_PTR);
            
            SplitKLDResult->KLD[m]+=0.5*(logf(DataValue/DistValue)*(*DataHist_PTR)+logf(DistValue/DataValue)*(*DistHist_PTR)); // Use of the symmetric KLD between the two distributions
        }
        delete [] DistHist_VEC[m];
        DistHist_VEC[m]=NULL;
        delete [] DataHist_VEC[m];
        DataHist_VEC[m]=NULL;
    }
    delete [] DistHistTotal;
    DistHistTotal=NULL;
    delete [] DataHistTotal;
    DataHistTotal=NULL;
    
    return SplitKLDResult;
}

// Returns the vector containing the pointers to all the SplitKLD of the children of the considered tree. Used to obtain ordered list of Splitting to try.
vector<SplitKLD *> TreeEM::GetVectorKLDSplitChildren(){
    int numbchild=this->GetNumberChildren();
    vector<SplitKLD *> KLDSplitChildren;
    for (int c=0; c<numbchild; c++) {
        KLDSplitChildren.push_back(this->GetKLDSplit(c));
    }
    return KLDSplitChildren;
}

// Returns the vector of pointers to all SplitKLD structures for the leaves of a given tree.
vector<SplitKLD *> TreeEM::GetVectorKLDSplitLeaves(){
    vector<SplitKLD *> KLDSplitLeaves;
    //cout<<KLDSplitLeaves.size()<<endl;
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetNumberChildren()==0) {
            //cout<<"Leaf for KLD calculation"<<endl;
            KLDSplitLeaves.push_back(this->GetKLDSplit(c));
        }
        else{
            vector<SplitKLD*> KLDSplitLeavesChild=this->GetChild(c)->GetVectorKLDSplitLeaves();
            for (int i=0; i<KLDSplitLeavesChild.size(); i++) {
                KLDSplitLeaves.push_back(KLDSplitLeavesChild[i]);
            }
        }
    }
    //cout<<"Size of KLDVector before return "<<KLDSplitLeaves.size()<<endl;
    return KLDSplitLeaves;
}

SplitKLD * TreeEM::GetToSplit(){
    vector<SplitKLD *> VectorSplit( this->GetVectorKLDSplitLeaves());
    int numbLeaves=VectorSplit.size();
    //cout<<VectorSplit.size()<<"is the size of Vector Split"<<endl;
    if (numbLeaves<=0) {
        return NULL;
    }
    // Initialisation of the result
    SplitKLD * SplitMax=new SplitKLD;
    SplitMax->KLDTot=-1E32;
    int numbmodal=this->GetNumberModalities();
    SplitMax->KLD=new float[numbmodal];//{-1E32};
    for(int i=0;i<numbmodal;i++){
        SplitMax->KLD[i]=-1E32;
    }
    for (int i=0; i<numbLeaves; i++) {
        bool SplitCheckTest=(VectorSplit[i]->Parent)->GetSplitCheck()[VectorSplit[i]->ChildToSplit];
        if (!SplitCheckTest&& SplitMax->KLDTot<VectorSplit[i]->KLDTot) { // if the leave has not been checked yet for the split and the KLD is more than the previous one, we update SplitMax by copying the structure of VectorSplit[i]
            SplitMax->Parent=VectorSplit[i]->Parent;
            SplitMax->ChildToSplit=VectorSplit[i]->ChildToSplit;
            SplitMax->KLDTot=VectorSplit[i]->KLDTot;
            for (int m=0; m<numbmodal; m++) {
                SplitMax->KLD[m]=VectorSplit[i]->KLD[m];
            }
        }
        // Delete and free the considered element of VectorSplit.
        //delete [] VectorSplit[i]->KLD;
        delete VectorSplit[i];
        VectorSplit[i]=NULL;
    }
    if (SplitMax->KLDTot<0) {
        delete SplitMax;
        SplitMax=NULL;
    }
    return SplitMax;
}

// Performs the Split operation and modify accordingly the tree
void TreeEM::SplitOperation(int ChildToSplit){
    // Check if it is effectively a child that can be split
    int numbchild=this->GetNumberChildren();
    if (ChildToSplit>=numbchild||ChildToSplit<0) {
        cout<<"Child to split out of bounds"<<endl;
        return;
    }
    if (this->GetChild(ChildToSplit)->GetNumberChildren()!=0) {
        cout<<"Impossible to split something else than leaf"<<endl;
        return;
    }
    Parameters ** ParametersSplit=(this->GetChild(ChildToSplit))->ParametersForSplitting();
    (this->GetChild(ChildToSplit))->CreateAndAddChild(NULL,1);
    (this->GetChild(ChildToSplit))->CreateAndAddChild(NULL,1);
    (this->GetChild(ChildToSplit)->GetChild(0))->SetParameters(ParametersSplit[0]);
    (this->GetChild(ChildToSplit)->GetChild(1))->SetParameters(ParametersSplit[1]);
    
    // Clearing memory dedicated to Parameters split
    delete ParametersSplit[0];
    ParametersSplit[0]=NULL;
    delete ParametersSplit[1];
    ParametersSplit[1]=NULL;
    delete ParametersSplit;
    ParametersSplit=NULL;
    this->GetChild(ChildToSplit)->GetChild(0)->SetNormWeight(0.5);
    this->GetChild(ChildToSplit)->GetChild(1)->SetNormWeight(0.5);
//    this->PutAllLeavesToChildrenLevel();
    return;
}

// Returns the pointer to the array of Parameters pointer initialising the classes formed by the split
Parameters ** TreeEM::ParametersForSplitting(){
    
    // Check if it is really a leaf to split
    if (this->GetNumberChildren()!=0) {
        return NULL;
    }
    Parameters ** ParametersSplit=new Parameters* [2];
    ParametersSplit[0]=new Parameters();
    ParametersSplit[1]=new Parameters();
    ParametersSplit[0]->DistributionType=this->GetDistributionType();
    ParametersSplit[1]->DistributionType=this->GetDistributionType();
    ParametersSplit[0]->SizeParameters=this->GetSizeParameters();
    ParametersSplit[1]->SizeParameters=this->GetSizeParameters();
    ParametersSplit[0]->ValueParameters=new float[ParametersSplit[0]->SizeParameters];//{0};
    for(int i=0;i<ParametersSplit[0]->SizeParameters;i++){
        ParametersSplit[0]->ValueParameters[i]=0;
    }
    ParametersSplit[1]->ValueParameters=new float[ParametersSplit[1]->SizeParameters];//{0};
    for(int i=0;i<ParametersSplit[1]->SizeParameters;i++){
        ParametersSplit[1]->ValueParameters[i]=0;
    }
    float * ParametersMean=this->GetMean();
    float * ParametersVariance=this->GetVariance();
    
    // Temporary values for splitting initialisation :
    int numbmodal=this->GetNumberModalities();
//    for (int m=0; m<numbmodal; m++) {
//        
//        //Mean Initialisation
//        ParametersSplit[0]->ValueParameters[m]=ParametersMean[m]-0.5*powf(ParametersVariance[m+numbmodal*m], 0.5);
//        ParametersSplit[1]->ValueParameters[m]=ParametersMean[m]+0.5*powf(ParametersVariance[m+numbmodal*m], 0.5);
//        
//        //Variance Initialisation
//        ParametersSplit[0]->ValueParameters[numbmodal+m+numbmodal*m]=ParametersVariance[m+numbmodal*m]/2.0;
//        ParametersSplit[1]->ValueParameters[numbmodal+m+numbmodal*m]=ParametersVariance[m+numbmodal*m]/2.0;
//        
//    }
    
    // Splitting initialisation with Richardson initialisation
    float * MeanInitSplit=this->MeanForSplitInitialisation();
    float * VarianceInitSplit=this->VarianceForSplitInitialisation();
    
    // Filling in the parameters
    for (int m1=0; m1<numbmodal; m1++) {
        ParametersSplit[0]->ValueParameters[m1]=MeanInitSplit[m1];
        ParametersSplit[1]->ValueParameters[m1]=MeanInitSplit[m1+numbmodal];
        for (int m2=0; m2<numbmodal; m2++) {
            ParametersSplit[0]->ValueParameters[numbmodal+m1+m2*numbmodal]=VarianceInitSplit[m1+m2*numbmodal];
            ParametersSplit[1]->ValueParameters[numbmodal+m1+m2*numbmodal]=VarianceInitSplit[numbmodal*numbmodal+m1+m2*numbmodal];
        }
    }
    // Memory clearing
    delete [] MeanInitSplit;
    MeanInitSplit=NULL;
    delete [] VarianceInitSplit;
    VarianceInitSplit=NULL;
    
    return ParametersSplit;
}

// Applied on the leaf to split returns a pointer to the float array containing the initialisation to use for the mean
float * TreeEM::MeanForSplitInitialisation(){
    if (this->GetNumberChildren()!=0) {
        cout<<"Nothing to split since not leaf"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * MeanSplitInit=new float[2*numbmodal];
    float * MeanToSplit=this->GetMean();
    SVD SVDForSplit=SVD(this->GetVariance(), numbmodal,numbmodal);
    float * SVDEigen=SVDForSplit.getSingularValues();
    float * SVDVec=SVDForSplit.getV();
    int IndMaxEigen=0;
    float MaxEigen=1E-6;
    float * MaxVec=SVDVec;
    for (int i=0; i<numbmodal; i++) {
        if (SVDEigen[i]>MaxEigen) {
            IndMaxEigen=i;
            MaxEigen=SVDEigen[i];
            MaxVec=&SVDVec[i*numbmodal];
        }
    }
    for (int m=0; m<numbmodal; m++) {
        MeanSplitInit[m]=MeanToSplit[m]-0.5*MaxVec[m]*sqrtf(MaxEigen);
        MeanSplitInit[m+numbmodal]= MeanToSplit[m]+0.5*MaxVec[m]*sqrtf(MaxEigen);
    }
//    delete [] SVDForSplit.getU();
//    delete [] SVDForSplit.getV();
//    delete [] SVDForSplit.getSingularValues();
    return MeanSplitInit;
}

float * TreeEM::VarianceForSplitInitialisation(){
    if (this->GetNumberChildren()!=0) {
        cout<<"Nothing to split since not leaf"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * VarianceSplitInit=new float[2*numbmodal*numbmodal];
    float * VarianceToSplit=this->GetVariance();
    SVD SVDForSplit=SVD(this->GetVariance(), numbmodal,numbmodal);
    float * SVDEigen=SVDForSplit.getSingularValues();
    float * SVDVec=SVDForSplit.getV();
    int IndMaxEigen=0;
    float MaxEigen=1E-6;
    float * MaxVec=SVDVec;
    for (int i=0; i<numbmodal; i++) {
        if (SVDEigen[i]>MaxEigen) {
            IndMaxEigen=i;
            MaxEigen=SVDEigen[i];
            MaxVec=&SVDVec[i*numbmodal];
        }
    }
    for (int m1=0; m1<numbmodal; m1++) {
        for (int m2=0; m2<numbmodal; m2++) {
            VarianceSplitInit[m1+m2*numbmodal]=VarianceToSplit[m1+m2*numbmodal]+2*(0.5-1-0.5*0.5*0.5)*MaxVec[m1]*MaxVec[m2]*MaxEigen+MaxVec[m1]*MaxVec[m2]*MaxEigen;
            VarianceSplitInit[m1+m2*numbmodal+numbmodal*numbmodal]=VarianceToSplit[m1+m2*numbmodal]+2*(0.5*0.5*0.5-0.5-0.5*0.5)*MaxVec[m1]*MaxVec[m2]*MaxEigen+MaxVec[m1]*MaxVec[m2]*MaxEigen;
        }
    }
    //delete &SVDForSplit;
//    delete [] SVDForSplit.getU();
//    delete [] SVDForSplit.getV();
//    delete [] SVDForSplit.getSingularValues();
    return VarianceSplitInit;
}



// When testing a S operation returns the best of the two models obtained.
TreeEM* TreeEM::RunSplitOperation(){
    SplitKLD * SplitTry=this->GetToSplit();
    if (SplitTry==NULL) {
        return this;
    }
    // Updating of the checking
    cout<< "Trying split operation"<<endl;
    SplitTry->Parent->GetSplitCheck()[SplitTry->ChildToSplit]=1;
    // Copy of the current tree
    TreeEM * CopiedTree=this->CopyTree(NULL);
    (SplitTry->Parent)->SplitOperation(SplitTry->ChildToSplit);
    this->ClearSMChecks();
    


    
    // Partial EM first
    float PartCLL=0;
    float PartOldCLL=0;
    int PartIteration=0;
    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateDistribution();
    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNonNormResp(); 
    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormResp();
    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNonNormWeights();
    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormWeights();
    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    // Run EM
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    this->UpdateDistribution();
    this->UpdateNonNormResp();
    this->UpdateNormResp();
    this->SaveAllClasses((g_scratch_dir + "/TestSplitResult").c_str());
    this->UpdateNonNormWeights();
    this->UpdateNormWeights();
    //this->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration);
    this->SaveAllClasses((g_scratch_dir + "/TestSplitResult").c_str());
    // Test of the operation
    float OldCriterion=CopiedTree->CriterionCalculation();
    float NewCriterion=this->CriterionCalculation();
    
    // Clear SplitTry
    if (SplitTry!=NULL) {
        delete SplitTry;
        SplitTry=NULL;
    }
    
    if (OldCriterion>=NewCriterion) { // Criterion of change not met, delete the modified tree and return the copied one
        delete this;
        cout<<"Split not accepted"<<endl;
        return CopiedTree;
    }
    else{ // Criterion of change met delete the not modified tree stored in CopiedTree and return the new model.
        delete CopiedTree;
        cout<<"Split accepted"<<endl;
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->PutAllLeavesToChildrenLevel();
        }
        this->ClearSMChecks();
        return this;
    }
}


/***************** METHODS FOR BIAS FIELD CORRECTION *************/

// Returns a float array with all the basis functions to be used in the bias field correction
/** REMARK : To be thought again to see if we want to get it at each iteration and delete it afterwards or store it at the root as L2S or S2L **/
float * TreeEM::MakeBasisFunctions(){
    int numbFunctions=(int)((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numelmasked=this->GetNumberMaskedElements();
    float * BasisFunctions=new float[numelmasked*numbFunctions];//{0};
    for(int i=0;i<numelmasked*numbFunctions;i++){
        BasisFunctions[i]=0;
    }
    // Find the factors for the affine transformation of the indices so that the power functions are applied between 0 and 1
    float FactorX=(float)(2.0/(this->GetDataImage()->nx-1));
    float FactorY=(float)(2.0/(this->GetDataImage()->ny-1));
    float FactorZ=(float)(2.0/(this->GetDataImage()->nz-1));
    
    int SizePlane=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int SizeColumn=this->GetDataImage()->nx;
    int zValue=0;
    int yValue=0;
    int xValue=0;
    
    int l=0;//Represent the index of the considered basis function
    for (int order=0; order<=BForder; order++) {
        for (int xorder=0; xorder<=order; xorder++) {
            for (int yorder=0; yorder<=order-xorder; yorder++) {
                int zorder=order-xorder-yorder; // the sum of the orders for the different directions must be at most BForder.
                // Reinitialisation of the pointers
                int * S2L_PTR=this->GetS2L();
                float * BasisFunctions_PTR=&BasisFunctions[numelmasked*l];
                for (int i=0; i<numelmasked; i++,S2L_PTR++,BasisFunctions_PTR++) {
                    // conversion of S2L value into x,y and z component
                    zValue=*S2L_PTR/SizePlane;
                    yValue=(*S2L_PTR-zValue*SizePlane)/SizeColumn;
                    xValue=(*S2L_PTR)-zValue*SizePlane-yValue*SizeColumn;
                    
                    // Filling the BasisFunctions matrix with the proper values modified to be symmetric and between -1 and 1.
                    *BasisFunctions_PTR=powf(xValue*FactorX-1, xorder)*powf(yValue*FactorY-1, yorder)*powf(zValue*FactorZ-1, zorder);
                }
                l++;
            }
        }
    }
    return BasisFunctions;
}

// Returns the vector of numbmodal pointers to float array containing the values for the W matrix used in VL BF correction.
float * TreeEM::MakeWMatrixChildren(){
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    int numbchild=this->GetNumberChildren();

    float * WMatrix=new float[numelmasked*numbmodal];//{0};
    if(WMatrix==NULL){
        cout<<"Did not allocate Wmatrix properly"<<endl;
    }
        for(int i=0;i<numelmasked*numbmodal;i++){
            WMatrix[i]=0;
        }
        for(int m=0;m<numbmodal;m++){
        for (int c=0; c<numbchild; c++) {
            //Reinitialisation of the pointers
            float * WMatrix_PTR=&WMatrix[m*numelmasked];
            float * NormResp_PTR=this->GetChild(c)->GetNormResp();
            float * DiagVarianceValue=this->GetChild(c)->GetDiagVarianceDirect();
            float invVariance=DiagVarianceValue[m]>1E-6?1.0/DiagVarianceValue[m]:1E6;
            for (int i=0; i<numelmasked; i++,WMatrix_PTR++,NormResp_PTR++) {
                // Filling of the matrix as sum of normresp for the classes divided by the value of the diagonal in the covariance matrix.
                *WMatrix_PTR+=*NormResp_PTR*invVariance;
            }
            delete [] DiagVarianceValue;
        }
    }
    return WMatrix;
}

// Returns in a vector of pointers  the RMatrix needed for VL BF correction according to the different modalities.
float * TreeEM::MakeRMatrixChildren(){
    

    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    float * RMatrix=new float[numelmasked*numbmodal];
    int numbchild=this->GetNumberChildren();
    int numel=this->GetNumberElements();
    float * WMatrix=this->MakeWMatrixChildren();
    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * Data_PTR=Data;
    int * L2S_PTR=this->GetL2S();
    


        for(int i=0;i<numelmasked*numbmodal;i++){
            RMatrix[i]=0;
        }
        for(int m=0;m<numbmodal;m++){
        L2S_PTR=this->GetL2S();
        Data_PTR=&Data[m*numel];
        float * RMatrix_PTR=&RMatrix[m*numelmasked];
        for (int c=0; c<numbchild; c++) {
            //Reinitialisation of the pointers
            RMatrix_PTR=&RMatrix[m*numelmasked];
            float * NormResp_PTR=this->GetChild(c)->GetNormResp();
            float * WMatrix_PTR=&WMatrix[m*numelmasked];
            float *  Mean=this->GetChild(c)->GetMeanDirect();
            float * DiagVarianceValue=this->GetChild(c)->GetDiagVarianceDirect();
            float invVariance=DiagVarianceValue[m]>1E-6?1.0/DiagVarianceValue[m]:1E6;
            for (int i=0; i<numelmasked; i++,WMatrix_PTR++,NormResp_PTR++,RMatrix_PTR++) {
                // Filling of the matrix as normalised weighted sum of the means 
                *RMatrix_PTR+=(*NormResp_PTR*Mean[m]/invVariance)/(*WMatrix_PTR);
            }
            delete [] Mean;
            delete [] DiagVarianceValue;
        }
        RMatrix_PTR=&RMatrix[m*numelmasked];
        for (int i=0; i<numel; i++, Data_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                *RMatrix_PTR=*Data_PTR-*RMatrix_PTR;
                RMatrix_PTR++;
            }
        }
    }
        if(WMatrix!=NULL){
            delete [] WMatrix;
            WMatrix=NULL;
        }
    return RMatrix;
}

// Returns a vector of pointers to the AtWA matrices used for the BF correction. The obtained matrix will need to be inverted afterwards.
float * TreeEM::MakeAtWAMatrixChildren(){
    float * WMatrix=this->MakeWMatrixChildren();
    float * BFFunctions=this->MakeBasisFunctions();
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float* AtWAMatrix=new float[numbBF*numbBF*numbmodal];

        for(int i=0;i<numbBF*numbBF*numbmodal;i++){
            AtWAMatrix[i]=0;
        }
        
        // Not optimal since all calculated whereas symmetric matrix obtained
        for(int m=0;m<numbmodal;m++){
        for (int j1=0; j1<numbBF; j1++) {

            for (int j2=0; j2<numbBF; j2++) {
                float * WMatrix_PTR=&WMatrix[m*numelmasked];
                float * BFFunctions_PTR1=&BFFunctions[j1*numelmasked];
                float * BFFunctions_PTR2=&BFFunctions[j2*numelmasked];
                for (int i=0; i<numelmasked; i++,BFFunctions_PTR1++,WMatrix_PTR++,BFFunctions_PTR2++) {
                    AtWAMatrix[j1+numbBF*j2+m*numbBF*numbBF]+=(*BFFunctions_PTR1)*(*BFFunctions_PTR2)*(*WMatrix_PTR);
                }
            }
        }
    }
        //clearing WMatrix
        if(WMatrix!=NULL){
            delete [] WMatrix;
            WMatrix=NULL;
        }
    // clearing BFFunctions
    delete [] BFFunctions;
    BFFunctions=NULL;
    return AtWAMatrix;
}

// Returns a vector of float pointers to the AtWR vectors needed for the VL BF correction
float * TreeEM::MakeAtWRVectorChildren(){
    // Initialisation and getting all needed presteps with other matrices (A, W and R)

    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities();
    float * BFFunctions=this->MakeBasisFunctions();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * WMatrix=this->MakeWMatrixChildren();
    float * RMatrix=this->MakeRMatrixChildren();
    

        float * AtWRVector=new float[numbBF*numbmodal];//{0};
        for(int i=0;i<numbBF*numbmodal;i++){
            AtWRVector[i]=0;
        }
        for (int m=0; m<numbmodal; m++) {
        for (int j=0; j<numbBF; j++) {
            float * BFFunctions_PTR=&BFFunctions[j*numelmasked];
            float * WMatrix_PTR=&WMatrix[m*numelmasked];
            float * RMatrix_PTR=&RMatrix[m*numelmasked];
            for (int i=0; i<numelmasked; i++,WMatrix_PTR++,BFFunctions_PTR++,RMatrix_PTR++) {
                AtWRVector[j+m*numbBF]+=(*BFFunctions_PTR)*(*WMatrix_PTR)*(*RMatrix_PTR);
            }


    }
        }
        // Clearing memory allocation
        delete [] WMatrix;
        WMatrix=NULL;
        delete [] RMatrix;
        RMatrix=NULL;
        delete [] BFFunctions;
        BFFunctions=NULL;

    return AtWRVector;
}

// Returns in a vector of pointers the inverted matrix AtWA for each modality.
float* TreeEM::MakeInvAtWAMatrixChildren(){
    float * AtWAMatrix=this->MakeAtWAMatrixChildren();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
   // float * invAtWAMatrixResult;
    float * invAtWAMatrix=new float[numbBF*numbBF*numbmodal];//{0};
    for(int i=0;i<numbBF*numbBF*numbmodal;i++){
        invAtWAMatrix[i]=0;
    }
    // Allocation of the matrix
    for (int m=0; m<numbmodal; m++) {
        matrix<float> AtWAMatrix_mat=matrix<float>(numbBF);

        float * AtWAMatrix_tmp=&AtWAMatrix[m*numbBF*numbBF];
        bool success=0;
        for (int j1=0; j1<numbBF; j1++) {
            for (int j2=0; j2<numbBF; j2++) {
                AtWAMatrix_mat.setvalue(j1, j2, AtWAMatrix_tmp[j1+j2*numbBF]);
            }
        }
        AtWAMatrix_mat.invert();
        for (int j1=0; j1<numbBF; j1++) {
            for (int j2=0; j2<numbBF; j2++) {
                AtWAMatrix_mat.getvalue(j1, j2, invAtWAMatrix[j1+j2*numbBF+m*numbBF*numbBF], success);
            }
        }
    }
    // Clearing memory
    delete [] AtWAMatrix;
    AtWAMatrix=NULL;
    return invAtWAMatrix;
}

// Returns in a vector of pointers the coefficients for the bias field
float * TreeEM::MakeFinalBFCoeffsChildren(){
    // Initialisation and obtention of the needed values

    int numbmodal=this->GetNumberModalities();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * FinalBFCoeffsResult = new float[numbBF*numbmodal];
    for(int i=0;i<numbBF*numbmodal;i++){
        FinalBFCoeffsResult[i]=0;
    }
    float * invAtWAMatrixVector=this->MakeInvAtWAMatrixChildren();
    float * AtWRVector=this->MakeAtWRVectorChildren();
    for (int m=0; m<numbmodal; m++) {
        float * invAtWA=&invAtWAMatrixVector[m*numbBF*numbBF];
        float * AtWRV_PTR=&AtWRVector[m*numbBF];
        float * invAtWA_PTR=invAtWA;
        for (int j1=0; j1<numbBF; j1++) {
            invAtWA_PTR=&invAtWA[j1*numbBF];
            AtWRV_PTR=&AtWRVector[m*numbBF];
            for (int j2=0; j2<numbBF; j2++,AtWRV_PTR++,invAtWA_PTR++) {
                FinalBFCoeffsResult[j1+m*numbBF]+=(*AtWRV_PTR)*(*invAtWA_PTR);
            }
        }
    }
    // Clearing memory
    delete [] invAtWAMatrixVector;
    invAtWAMatrixVector=NULL;
    delete [] AtWRVector;
    AtWRVector=NULL;
    return FinalBFCoeffsResult;
}

float * TreeEM::MakeDataBFCorrected(){
    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
    float * DataBFCorrected=new float[numel*numbmodal];//{0};
    for(int i=0;i<numel*numbmodal;i++){
        DataBFCorrected[i]=0;
    }
    float * DataBFCorrected_PTR=DataBFCorrected;
    float * DataInit=static_cast<float *>(this->GetDataImage()->data);
    float * DataInit_PTR=DataInit;
    for (int m=0; m<numbmodal; m++) {
        DataBFCorrected_PTR=&DataBFCorrected[m*numel];
        DataInit_PTR=&DataInit[m*numel];
        for (int i=0; i<numel; i++,DataBFCorrected_PTR++,DataInit_PTR++) {
            *DataBFCorrected_PTR=*DataInit_PTR;
        }
    }
    // If there is a correction of the bias field to be done
    float * BFCoeffsUsed=this->GetBFCoeffs();
    

        //float * BFCoeffs=this->GetBFCoeffs();
        int numelmasked=this->GetNumberMaskedElements();
        //int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
        if (BFCoeffsUsed!=NULL) { // Check first if the correction is effectively possible
            float * BFCorrection_new=this->GetBFCorrection();
            for (int m=0; m<numbmodal; m++) {
                float * BFCorrection_PTR=&BFCorrection_new[m*numelmasked];
                    int * L2S_PTR=this->GetL2S();
                    for (int i=0; i<numel; i++,DataBFCorrected_PTR++) {
                        if (*L2S_PTR>=0) {
                            *DataBFCorrected_PTR-=(*BFCorrection_PTR);
                            BFCorrection_PTR++;
                        }
                    }
                }
        }

                return DataBFCorrected;
        }


float * TreeEM::MakeBFCorrection(){

    float * BFCoeffs_new=this->GetBFCoeffs();
    if(BFCoeffs_new==NULL){
        return NULL;
    }
    float * BFBasisFunctions=this->MakeBasisFunctions();
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
    float * BFCorrectionResult=new float[numbmodal*numelmasked];
    for(int i=0;i<numbmodal*numelmasked;i++){
        BFCorrectionResult[i]=0;
    }
    if(BFCoeffs_new!=NULL){
    for(int m=0;m<numbmodal;m++){
        for(int j=0;j<numbBF;j++){
            float * BasisFunctions_PTR=&BFBasisFunctions[j*numelmasked];
            float * BFCoeffs_PTR=&BFCoeffs_new[m*numbBF];
            float * BFCorrectionResult_PTR=&BFCorrectionResult[m*numelmasked];
            for(int i=0;i<numelmasked;i++,BFCorrectionResult_PTR++,BasisFunctions_PTR++){
                *BFCorrectionResult_PTR+=BFCoeffs_PTR[j]*(*BasisFunctions_PTR);
            }
        }
    }
    }
    delete[] BFBasisFunctions;
    BFBasisFunctions=NULL;
    return BFCorrectionResult;
}





/******************* SAVING RESULTS ******************************/
void TreeEM::SaveAllClasses(char * filenameOut){
    int numbTotalClasses=0;
    //this->GetNumberLeaves(numbTotalClasses);
    numbTotalClasses=this->GetNumberAllLeaves();
    vector<TreeEM *> LeavesVector;
    //this->GetLeaves(LeavesVector);
    LeavesVector=this->GetAllLeaves();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbTotalClasses;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    //float * Input_PTR=LeavesVector[0]->GetPartialResult(NORMRESP);
    for (int c=0; c<numbTotalClasses; c++) {
        Result_PTRtmp=&Result_PTR[c*numel];
        float * Input=LeavesVector[c]->GetPartialResult(NORMRESP);
        float * Input_PTR=Input;
        for (int i=0; i<numel; i++, Result_PTRtmp++,Input_PTR++) {
            *Result_PTRtmp=*Input_PTR;
        }
        if(Input!=NULL){
            delete[]Input;
            Input=NULL;
        }
    }
    nifti_set_filenames(Result, filenameOut, 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
    LeavesVector.clear();
}

void TreeEM::SaveGeneralClasses(char *filenameOut){
    int numbGeneralClasses=this->GetNumberChildren();
    vector<TreeEM *> ChildrenVector=this->GetChildren();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbGeneralClasses;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    //float * Input_PTR=ChildrenVector[0]->GetPartialResult(NORMRESP);
    for (int c=0; c<numbGeneralClasses; c++) {
        Result_PTRtmp=&Result_PTR[c*numel];
        float *Input=ChildrenVector[c]->GetPartialResult(NORMRESP);
        float * Input_PTR=Input;
        for (int i=0; i<numel; i++, Result_PTRtmp++,Input_PTR++) {
            *Result_PTRtmp=*Input_PTR;
        }
        if(Input!=NULL){
            delete [] Input;
            Input=NULL;
        }
    }
    nifti_set_filenames(Result, filenameOut, 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
    ChildrenVector.clear();
}

void TreeEM::SaveBFBasisFunctions(char *filenameBF){
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * BasisFunctionsToSave=this->MakeBasisFunctions();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbBF;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    float * Input_PTR=BasisFunctionsToSave;
    int * L2S_PTR=this->GetL2S();
    for (int j=0; j<numbBF; j++) {
        Result_PTRtmp=&Result_PTR[j*numel];
        Input_PTR=&BasisFunctionsToSave[j*numelmasked];
        L2S_PTR=this->GetL2S();
        for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                *Result_PTRtmp=*Input_PTR;
                Input_PTR++;
            }
            else{
                *Result_PTRtmp=0;
            }
            
        }
    }
    if(BasisFunctionsToSave!=NULL){
        delete[] BasisFunctionsToSave;
        BasisFunctionsToSave=NULL;
    }
    nifti_set_filenames(Result, filenameBF, 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
}

void TreeEM::SaveBFBasisFunctions(int BFnumb,char *filenameBF){
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    if (BFnumb>=numbBF) {
        return;
    }
    float * BasisFunctionsToSave=this->MakeBasisFunctions();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=3;
    Result->dim[4]=1;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    float * Input_PTR=&BasisFunctionsToSave[BFnumb*numelmasked];
    int * L2S_PTR=this->GetL2S();


        for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                *Result_PTRtmp=*Input_PTR;
                Input_PTR++;
            }
            else{
                *Result_PTRtmp=0;
            }
            
        }
        if(BasisFunctionsToSave!=NULL){
            delete[]BasisFunctionsToSave;
            BasisFunctionsToSave=NULL;
        }
    nifti_set_filenames(Result, filenameBF, 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
}

void TreeEM::SaveBFCorrection(char *filenameBF){
    int numelmasked=this->GetNumberMaskedElements();
//    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
    if (BFFlag) {

        
            float * BasisFunctionsToSave=this->MakeBasisFunctions();
            float * BFCoeffsToSave=this->GetBFCoeffs();
            if (this->AreBFCoeffsDirectUseable()) {
                float * BFCoeffs_tmp=BFCoeffsToSave;
                nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
                Result->dim[0]=4;
                Result->dim[4]=numbmodal;
                Result->dim[5]=1;
                nifti_update_dims_from_array(Result);
                Result->data = (void *) calloc(Result->nvox, sizeof(float));
                int numel=this->GetNumberElements();
                int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
                float * Result_PTR=static_cast<float *>(Result->data);
                float * Result_PTRtmp=Result_PTR;
                for (int i=0; i<numel; i++, Result_PTRtmp++) {
                    *Result_PTRtmp=0;
                }
                float * Input_PTR=BasisFunctionsToSave;
                int * L2S_PTR=this->GetL2S();
                
                for (int m=0; m<numbmodal; m++) {
                    for (int j=0; j<numbBF; j++) {
                        BFCoeffs_tmp=&BFCoeffsToSave[m*numbBF];
                        Result_PTRtmp=&Result_PTR[m*numel];
                        Input_PTR=&BasisFunctionsToSave[j*numelmasked];
                        L2S_PTR=this->GetL2S();
                        for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
                            if (*L2S_PTR>=0) {
                                *Result_PTRtmp+=*Input_PTR*BFCoeffs_tmp[j];
                                Input_PTR++;
                            }
                            else{
                                *Result_PTRtmp=0;
                            }
                            
                        }
                        
                    }
                }
                nifti_set_filenames(Result, filenameBF, 0, 0);
                nifti_image_write(Result);
                nifti_image_free(Result);
            }
        }


        }

void TreeEM::SaveBFCorrectedData(char * filenameBF){
//    int numelmasked=this->GetNumberMaskedElements();
//    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
        if (this->AreBFCoeffsDirectUseable()) {
            nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
            Result->dim[0]=4;
            Result->dim[4]=numbmodal;
            Result->dim[5]=1;
            nifti_update_dims_from_array(Result);
            Result->data = (void *) calloc(Result->nvox, sizeof(float));
            int numel=this->GetNumberElements();
            float * Result_PTR=static_cast<float *>(Result->data);
            float * Result_PTRtmp=Result_PTR;
            for (int i=0; i<numel; i++, Result_PTRtmp++) {
                *Result_PTRtmp=0;
            }
            float * Input=this->MakeDataBFCorrected();
            float * Input_PTR=Input;
//            int * L2S_PTR=this->GetL2S();
            
            for (int m=0; m<numbmodal; m++) {
                
                    
                    
                    Result_PTRtmp=&Result_PTR[m*numel];
                    Input_PTR=&Input[m*numel];
                    
                    for (int i=0; i<numel; i++, Result_PTRtmp++,Input_PTR++) {
                        
                            *Result_PTRtmp=*Input_PTR;
       
                        
                    }
                    
                }
            
            nifti_set_filenames(Result, filenameBF, 0, 0);
            nifti_image_write(Result);
            nifti_image_free(Result);
        }
    
}

void TreeEM::SaveTreeInTextFile(char *filenameOut){
    std::ofstream fs(filenameOut);
    
    if(!fs)
    {
        std::cout<<"Cannot open the output file."<<std::endl;
        return;
    }
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    fs<<"There are "<<numbchild<<" general classes \n";
    fs<<"The studied data contains "<<numbmodal<<" modalities\n";
    // Writing BF coefficients if exist
    if (this->AreBFCoeffsDirectUseable()) {
        int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
        float * BFCoeffs_new=this->GetBFCoeffs();
        fs<<" The BF coefficients are : \n";
        for (int m=0; m<numbmodal; m++) {
            float * BFCoeffs_PTR=&BFCoeffs_new[m*numbBF];
            fs << "For modality "<<m<<" \n";
            for (int l=0; l<numbBF; l++,BFCoeffs_PTR++) {
                fs<<*BFCoeffs_PTR<<"   -    ";
            }
            fs<<"\n";
        }
    }
    // Writing parameters and weights of the different classes
    for (int c=0; c<numbchild; c++) {
        int numbLeaves=this->GetChild(c)->GetNumberAllLeaves();
        vector<TreeEM *> LeavesVector=this->GetChild(c)->GetAllLeaves();
        fs<<"General class "<<c<<" contains "<<numbLeaves<<" subclasses and its weight is "<<this->GetChild(c)->GetNormWeight()<<" \n";
        if (numbLeaves==0) {
            fs<<"The general class "<<c<<" is not a mixture and its parameters are : \n Mean :\n";
            float * Mean=this->GetChild(c)->GetMean();
            float * Variance=this->GetChild(c)->GetVariance();
            for (int m=0; m<numbmodal; m++) {
                fs<<Mean[m]<<"    ";
            }
            fs<<"\n Variance \n";
            for (int m1=0; m1<numbmodal; m1++) {
                for (int m2=0; m2<numbmodal; m2++) {
                    fs<<Variance[m1+m2*numbmodal]<<"     ";
                }
                fs<<"\n";
            }

        }
        else{
        for (int l=0; l<numbLeaves; l++) {
            fs<<"The subclass "<<l<<" of general class "<<c<<" is of weight "<<LeavesVector[l]->GetNormWeight()<<" \n";
            fs<<"The parameters of the distribution are : \n Mean \n";
            float * Mean=LeavesVector[l]->GetMean();
            float * Variance=LeavesVector[l]->GetVariance();
            for (int m=0; m<numbmodal; m++) {
                fs<<Mean[m]<<"    ";
            }
            fs<<"\n Variance \n";
            for (int m1=0; m1<numbmodal; m1++) {
                for (int m2=0; m2<numbmodal; m2++) {
                    fs<<Variance[m1+m2*numbmodal]<<"     ";
                }
                fs<<"\n";
            }
        }
    }
    }
    //fs<<"ghgh";
    fs.close();
    return;
}










