#ifndef _SEG_IR_CPP
#define _SEG_IR_CPP
//
//  Seg_InitAndRun.cpp
//  TreeSeg
//
//  Created by Carole Sudre on 15/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 08/05/2026 - Fixed MRFImage, GMatrixIn, and GMatrixPost path resolution with ChangePath parameter.
//

#include "_Seg_InitAndRun.h"




vector<nifti_image *> ReadFromFilenamesVector(vector<string> FilenamesVector){
    int sizeVector=FilenamesVector.size();
    std::cout<< "Number of Filenames is "<<sizeVector<<endl;
    vector<nifti_image*> ImagesVector;
    for (int f=0; f<sizeVector; f++) {
        nifti_image * ReadImage=nifti_image_read(FilenamesVector[f].c_str(), true);
        ImagesVector.push_back(ReadImage);
//        if (ReadImage->dim[0]>=4) {
//            break;
//        }
        //cout<< "Size of dimension x in image "<< nifti_image_read(FilenamesVector[f], true)->nx<<endl;
    }
    return ImagesVector;
}

nifti_image * ReadFromFilename(string Filename){
//    std::cout<<Filename<<endl;
    return nifti_image_read(Filename.c_str(),true);
}

int * ReadCountfromFile(string DPFile){
    ifstream f ;
    std::cout<< DPFile<<endl;
    f.open( DPFile.c_str()) ;
    if(!f.is_open()){
        cout<<"Pb in opening DP file"<<endl;
        return NULL;
    }
    else{
        vector<int> DPValues;
        int num;
        while(f >> num){
            DPValues.push_back(num);
        }
        int DPsize=DPValues.size();
        int * DPRead=new int[DPsize];
        for(int i=0;i<DPsize;i++){
            DPRead[i]=DPValues[i];
        }
        DPValues.clear();
        return DPRead;
    }
    
}

vector<int *> ReadCountfromFiles(vector<string> DPFiles){
    int numbFiles=DPFiles.size();
    vector<int*> DPChildren;
    for(int c=0;c<numbFiles;c++){
        int * DPRead=ReadCountfromFile( DPFiles[c]);
        DPChildren.push_back(DPRead);
    }
    return DPChildren;
}

// Consider only the first numbImagesToConsider to build the DataImage and concatenate them on the 4th dimension
nifti_image * CreateDataImageProgressive(vector<nifti_image *> ImagesToSegment, int numbImagesToConsider){
    // Different treatment according to numbImagesToConsider
    if (numbImagesToConsider<=0) { // No images to consider
        return NULL;
    }
    int numbSeg=ImagesToSegment.size();
    if (numbImagesToConsider>numbSeg) { // All images contained in vector
        return CreateDataImage(ImagesToSegment);
    }
    else{
        vector<nifti_image *> PartialImageVector;
        for (int i=0; i<numbImagesToConsider; i++) { // Create the vector of images to consider to build the new data image.
            PartialImageVector.push_back(ImagesToSegment[i]);
        }
        return CreateDataImage(PartialImageVector);
    }
}


nifti_image * CreateDataImage(vector<nifti_image*> ImagesToSegment){
    int finn=FirstNotNULL(ImagesToSegment);
    int linn=LastNotNULL(ImagesToSegment);
    std::cout<< "First not Null is "<<finn<< " and last is "<< linn<<endl;
    if (linn-finn<0) {
        std::cout<<"No Images To Segment"<<endl;
        return NULL;
    }
    int nx=ImagesToSegment[finn]->nx;
    int ny=ImagesToSegment[finn]->ny;
    int nz=ImagesToSegment[finn]->nz;
    nifti_image * CreatedImage=NULL;
    if ((linn-finn)<0) { // no image in the vector
        std::cout<<"no image in the vector"<<endl;
        return CreatedImage;
    }
    else if((finn-linn)==0){ // only one image in the vector but can already be multimodal
        std::cout<<"only one image in the vector"<<endl;
        if (ImagesToSegment[finn]->datatype!=NIFTI_TYPE_FLOAT32){
            
            seg_changeDatatype<float>(ImagesToSegment[finn]);
        }
        CreatedImage =nifti_copy_nim_info(ImagesToSegment[finn]);
        float * ImageDataPtr=static_cast<float *>(ImagesToSegment[finn]->data);
        CreatedImage->data = (void *) calloc(CreatedImage->nvox, sizeof(float));
        
        float * CreatedImage_PTR_start=static_cast<float *>(CreatedImage->data);
        
        float * CreatedImage_PTR=CreatedImage_PTR_start;
        int numbvoxCreated=CreatedImage->nvox;
        for (int i=0; i<numbvoxCreated; i++ ,CreatedImage_PTR++, ImageDataPtr++) {
            if(*ImageDataPtr!=*ImageDataPtr){
                *CreatedImage_PTR=0;
            }
            else{
            (*CreatedImage_PTR)=(*ImageDataPtr);
            }
        }
        return CreatedImage;
    }
    else{ // for the moment only handle images with only 3 D and not 4D
        std::cout<<"There are "<<linn-finn+1<<" images to consider"<<endl;
        
        int numbmulti=0;
        for (int i=finn; i<=linn; i++) {
            if(ImagesToSegment[i]!=NULL){
                
                if (ImagesToSegment[finn]->datatype!=NIFTI_TYPE_FLOAT32){
                    seg_changeDatatype<float>(ImagesToSegment[i]);
                }
                //                numbmulti+=ImagesToSegment[i]->nu*ImagesToSegment[i]->nt;
                numbmulti++;
            }
        }
        CreatedImage =nifti_copy_nim_info(ImagesToSegment[finn]);
        CreatedImage->dim[0]=4;
        CreatedImage->dim[4]=numbmulti;
        CreatedImage->dim[5]=1;
        nifti_update_dims_from_array(CreatedImage);
        CreatedImage->data = (void *) calloc(CreatedImage->nvox, sizeof(float));
        float * CreatedImage_PTR_start=static_cast<float *>(CreatedImage->data);
        float * ImageDataPtr;
        float * CreatedImage_PTR=CreatedImage_PTR_start;
        for (int m=finn; m<=linn; m++) {
            if (ImagesToSegment[m]!=NULL) {
                if( ImagesToSegment[m]->nx!=nx || ImagesToSegment[m]->ny!=ny || ImagesToSegment[m]->nz!=nz){
                    std::cout<<"Not compatible dimensions"<<endl;
                }
                else{
                    int numel=ImagesToSegment[m]->nx*ImagesToSegment[m]->ny*ImagesToSegment[m]->nz;
                    std::cout<<"Compatible dimensions"<<endl;
                    seg_changeDatatype<float>(ImagesToSegment[m]);
                    ImageDataPtr=static_cast<float *>(ImagesToSegment[m]->data);
                    //                    for (int i=0; i<ImagesToSegment[m]->nvox; i++) {
                    for(int i=0;i<numel;i++){
                        (*CreatedImage_PTR)=(*ImageDataPtr);
                        CreatedImage_PTR++;
                        ImageDataPtr++;
                    }
                }
            }
        }
        //        for (int i=finn; i<=linn; i++) {
        //            nifti_image_free(ImagesToSegment[i]);
        //        }
        //        ImagesToSegment.clear();
        //        nifti_set_filenames(CreatedImage, (g_scratch_dir + "/CreatedImage.nii.gz").c_str(),0,0);
        //        nifti_image_write(CreatedImage);
        return CreatedImage;
    }
}

nifti_image * NoiseMaskOtsu(nifti_image * CreatedImage ){
    int numbmodal=CreatedImage->nu*CreatedImage->nt;
    int numel=CreatedImage->nx*CreatedImage->ny*CreatedImage->nz;
    nifti_image * TmpNoiseMask=nifti_copy_nim_info(CreatedImage);
    TmpNoiseMask->dim[0]=3;
    TmpNoiseMask->dim[4]=1;
    TmpNoiseMask->dim[5]=1;
    TmpNoiseMask->data=(void *)calloc(numel, sizeof(float));
    nifti_update_dims_from_array(TmpNoiseMask);
    float * TmpNoiseMask_PTR=static_cast<float *>(TmpNoiseMask->data);
    //Initialisation of noise mask
    for (int i=0; i<numel; i++) {
        TmpNoiseMask_PTR[i]=1;
    }
    
    float * Image_PTR=static_cast<float *>(CreatedImage->data);
    for (int m=0; m<numbmodal; m++) {
        //find image max and min
        float * Image=&Image_PTR[m*numel];
        float tempmax=-(1e32);
        float tempmin=1e32;
        for (int i=0; i<numel; i++) {
            if (Image[i]<tempmin) {
                tempmin=Image[i];
            }
            if (Image[i]>tempmax) {
                tempmax=Image[i];
            }
            
        }
        // Fill histogram
        float histsize=1002.0f;
        float histo[1002]={0};
        for(int i=0; i<(int)histsize; i++) histo[i]=0;
        
        for(int i=0; i<numel; i++){
            float location=(histsize-2)*(Image[i]-tempmin)/(tempmax-tempmin)+1;
            float weight=location-floorf(location);
            histo[(int)floor(location)]+=(1-weight);
            histo[(int)ceil(location)]+=(weight);
        }
        
        // Normalise histogram
        float sumhisto=0;
        for(int i=0; i<(int)histsize; i++)
            sumhisto+=histo[i];
        for(int i=0; i<(int)histsize; i++)
            histo[i]=histo[i]/sumhisto;
        
        
        float  w = 0;                // first order cumulative
        float  u = 0;                // second order cumulative
        float  uT = 0;               // total mean level
        float  work1, work2;		// working variables
        double work3 = 0.0;
        
        
        
        // Calculate total mean level
        for (int i=1; i<=histsize; i++)
            uT+=(i*histo[i-1]);
        
        
        // Find optimal threshold value
        for (int i=1; i<histsize; i++) {
            w+=histo[i-1];
            u+=(i*histo[i-1]);
            work1 = (uT * w - u);
            if(w==0||w==1){
                work2=0;
            }
            else{
                work2 = (work1 * work1) / ( w * (1.0f-w) );
            }
            if (work2>work3) work3=work2;
            //cout<<work3<<endl;
        }
        
        float threshold=(sqrt(work3)-1)/(histsize-2)*(tempmax-tempmin)+tempmin;
        
        
        cout<< "threshold = "<<threshold<<endl;
        // Convert the final value to an integer
        for(int i=0; i<numel; i++){
            TmpNoiseMask_PTR[i]=Image[i]>threshold?0:TmpNoiseMask_PTR[i];
        }
    }
    return TmpNoiseMask;
}

float * NoiseCovPriorsOtsu(nifti_image * CreatedImage){
    nifti_image * NoiseMask=NoiseMaskOtsu( CreatedImage);
    nifti_set_filenames(NoiseMask, (g_scratch_dir + "/NoiseMask.nii.gz").c_str(), 0, 0);
    nifti_image_write(NoiseMask);
    int numbmodal=CreatedImage->nu*CreatedImage->nt;
    int numel=CreatedImage->nx*CreatedImage->ny*CreatedImage->nz;
    float * NoiseMask_PTR=static_cast<float *>(NoiseMask->data);
    float * NoiseMask_PTRUsed=NoiseMask_PTR;
    float * VarianceNoise=new float[numbmodal*numbmodal];
    float * MeanNoise=new float[numbmodal];
    float * CreatedImage_PTR=static_cast<float *>(CreatedImage->data);
    float * CreatedImage_PTRUsed=CreatedImage_PTR;
    
    // Determination of mean for noise part of image
    int numbNoise=0;
    for (int m=0; m<numbmodal; m++) {
        MeanNoise[m]=0;
        numbNoise=0;
        CreatedImage_PTRUsed=&CreatedImage_PTR[m*numel];
        for (int i=0; i<numel; i++,CreatedImage_PTRUsed++, NoiseMask_PTRUsed++) {
            if (*NoiseMask_PTRUsed) {
                MeanNoise[m]+=*CreatedImage_PTRUsed;
                numbNoise++;
            }
        }
        MeanNoise[m]/=numbNoise;
        cout<<MeanNoise[m]<<" ";
    }
    cout<<endl;
    cout<<numbNoise<<" numbNoise"<<endl;
    
    // Determination of covariance matrix for noise part of image
    for(int m1=0;m1<numbmodal;m1++){
        // First data pointer to the beginning of the modality m1 considered
        for(int m2=0;m2<numbmodal;m2++){
            float * CreatedImage_PTRUsed1=&CreatedImage_PTR[m1*numel];
            float * CreatedImage_PTRUsed2=&CreatedImage_PTR[m2*numel];
            NoiseMask_PTRUsed=NoiseMask_PTR;
            PrecisionTYPE VarianceToUpdate_tmp=0;
            for(int i=0;i<numel;i++,CreatedImage_PTRUsed1++,CreatedImage_PTRUsed2++,NoiseMask_PTRUsed++){
                // Update of the numerator of the Variance Calculation only if in the case of an active voxel
                if (*NoiseMask_PTRUsed) {
                    VarianceToUpdate_tmp+=(PrecisionTYPE)((*CreatedImage_PTRUsed1)-MeanNoise[m1])*((*CreatedImage_PTRUsed2)-MeanNoise[m2]);
                }
            }
            VarianceNoise[m1+m2*numbmodal]=VarianceToUpdate_tmp/numbNoise;
        }
    }
    nifti_image_free(NoiseMask);
    delete [] MeanNoise;
    MeanNoise=NULL;
    return VarianceNoise;
}





vector<nifti_image *> ReadFromFilenamesVectorProgressive(vector<string> Filenames, int numbElementsToRead){
    vector<nifti_image *> VectorImages;
    int numbMaxFilenames=Filenames.size();
    if (numbElementsToRead<=0) {
        return VectorImages;
    }
    if (numbMaxFilenames<=numbElementsToRead) {
        VectorImages=ReadFromFilenamesVector(Filenames);
        return VectorImages;
    }
    vector<string> PartVectorFilenames;
    for (int i=0; i<numbElementsToRead; i++) {
        PartVectorFilenames.push_back(Filenames[i]);
    }
    VectorImages=ReadFromFilenamesVector(PartVectorFilenames);
    return VectorImages;
    
}

// Creates the pointer to the tree that will have the ImagesToSegment as input data, the PriorsVector as population priors and the Mask as mask.
TreeEM * CreateTreeSegmentationInit(SEG_PARAMETERS * segment_param){
    
//    vector<nifti_image* > VectorImageToSegment=ReadFromFilenamesVector(segment_param->filename_input);
    vector<nifti_image *> VectorImageToSegment=ReadFromFilenamesVectorProgressive(segment_param->filename_input,segment_param->NumberAddedModalities[0]);
    nifti_image * DataImage=CreateDataImage(VectorImageToSegment);
    int numbImages=VectorImageToSegment.size();
    for(int i=0;i<numbImages;i++){
        nifti_image_free(VectorImageToSegment[i]);
        VectorImageToSegment[i]=NULL;
    }
    VectorImageToSegment.clear();
    vector<nifti_image *> VectorPriors;
    if(segment_param->flag_manual_priors==1){
        VectorPriors=ReadFromFilenamesVector(segment_param->filename_priors);
    }
    vector<nifti_image *> VectorPriorsOut;
    if(segment_param->flag_manual_priors==1){
        VectorPriorsOut=ReadFromFilenamesVector(segment_param->filename_priors_out);
    }
    if (segment_param->OutliersMod==6 && VectorPriorsOut.size()==0) { // Correspond to incompatibility between OutliersMod and options
        segment_param->OutliersMod=3;
    }
    nifti_image * MaskToUse;
    if(segment_param->flag_mask==1){
        MaskToUse=ReadFromFilename(segment_param->filename_mask);
        int numbmodal=DataImage->dim[4];
        if(numbmodal>0){
            int numel=DataImage->nx*DataImage->ny*DataImage->nz;
            float * PointerToData=static_cast<float*>(DataImage->data);
            bool * PointerToMask=static_cast<bool*>(MaskToUse->data);
            for(int i=0;i<numel;i++){
                for(int m=0;m<numbmodal;m++){
                    if(PointerToData[m*numel+i]==0){
                        PointerToMask[i]=0;
                    }
                }
            }
            
        }
        std::cout<<"Mask number of elements is "<<MaskToUse->nx*MaskToUse->ny*MaskToUse->nz<<endl;
    }
    else{
        MaskToUse=NULL;
    }
    TreeEM * TreeToSegment=NULL;
    if (DataImage==NULL) {
        std::cout<<"Nothing to segment"<<endl;
        return NULL;
    }
    // creation of the node root
    TreeToSegment=new TreeEM();
    if(segment_param->flag_Outliers){
        TreeToSegment->SetFlagOutliers(segment_param->OutliersMod);
    }
    TreeToSegment->SetFlagUTC(segment_param->uniformTypeChange);
    TreeToSegment->SetData(DataImage);
    TreeToSegment->SetMask(MaskToUse);
    TreeToSegment->QuantilizeDataImage( segment_param);
    if(segment_param->flag_NormMask){
        if (!TreeToSegment->IsDataImageNormalised() && ! segment_param->flag_inDC) {
            cout<<"We have to normalise Image masked"<<endl;
            TreeToSegment->NormaliseDataImageMasked();
        }
    }
    else {
        if(!TreeToSegment->IsDataImageNormalised() && ! segment_param->flag_inDC) {
        cout<< "We have to normalise Image"<<endl;
        TreeToSegment->NormaliseDataImage();
        }
    }

    cout<< "Image is already normalised" << endl;
    std::cout<<TreeToSegment->GetDataImage()->nx<<endl;
    
    TreeToSegment->CreateAllocateAndInitializeParameters(0);
    
    if(!segment_param->flag_manual_priors){// case there is no priors for the general classes
        if(segment_param->OutliersMod==1){ // case we adopt the trilayered outlier model
            // first create children corresponding to outliers class
            float InlierWeight=1-segment_param->OutliersWeight;
            TreeToSegment->CreateAndAddChildWeight(segment_param,InlierWeight,0);
            nifti_image * ConstantOutlierPriors=TreeToSegment->BuildConstantPriors(1-InlierWeight);
            TreeToSegment->CreateAndAddChildPriors(segment_param,ConstantOutlierPriors,2);
            for (int c=0; c<segment_param->numb_classes; c++) {
                std::cout<<"Child number "<< c << "added"<<endl;
                nifti_image * ConstantPriors=TreeToSegment->BuildConstantPriors((float)(1.0/segment_param->numb_classes)*InlierWeight);
                TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,ConstantPriors,1);
            }
        }
        else{
            if(!segment_param->flag_Outliers){ // no priors and no outlier model
                for (int c=0; c<segment_param->numb_classes; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    nifti_image * ConstantPriors=TreeToSegment->BuildConstantPriors((float)(1.0/segment_param->numb_classes));
                    TreeToSegment->CreateAndAddChildPriors(segment_param,ConstantPriors,1);
                }
            }
            else{ // no priors but outlier model at same level
                float WeightToUse=1-segment_param->OutliersWeight;
                for (int c=0; c<segment_param->numb_classes; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    nifti_image * ConstantPriors=TreeToSegment->BuildConstantPriors((float)(1.0/segment_param->numb_classes)*WeightToUse);
                    TreeToSegment->CreateAndAddChildPriors(segment_param,ConstantPriors,1);
                }
                nifti_image * ConstantOutlierPriors=TreeToSegment->BuildConstantPriors(1-WeightToUse);
                TreeToSegment->CreateAndAddChildPriors(segment_param,ConstantOutlierPriors,2); // Uniform distribution for outlier class at initialisation
            }
        }
    }
    else{ // Normal initialisation
        int numbPriors=VectorPriors.size();
        //NormalisationPriorsList(TreeToSegment->GetData(), PriorsVector);
        //cout<<"Address of Mask"<<TreeToSegment->GetMask()<<endl;
        //depends if we have an outlier class to consider or not
        switch (segment_param->OutliersMod) {
            case 0:{ // Most classical case with priors but no outliers model
                
                for (int c=0; c<numbPriors; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    TreeToSegment->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                }
                TreeToSegment->NormalisePriors();
                TreeToSegment->NormalisePriorsAdapted();
            }
                break;
            case 1:{
                // first create children corresponding to outliers class
                float InlierWeight=1-segment_param->OutliersWeight;
                TreeToSegment->CreateAndAddChildWeight(segment_param,InlierWeight,0);
                nifti_image * ConstantOutlierPriors=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                TreeToSegment->CreateAndAddChildPriors(segment_param,ConstantOutlierPriors,2);
                for (int c=0; c<numbPriors; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                }
                if (segment_param->flag_OutBrain) {
                    nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                }
                TreeToSegment->GetChild(0)->NormalisePriors();
                TreeToSegment->GetChild(0)->NormalisePriorsAdapted();
            }
                
                break;
            case 2:{ // Outliers model 2
                for (int c=0; c<numbPriors; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    TreeToSegment->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                }
                if (segment_param->flag_OutBrain) {
                    nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                    TreeToSegment->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                }
                TreeToSegment->NormalisePriors();
                TreeToSegment->NormalisePriorsAdapted();
                if(segment_param->flag_Outliers){
                    // In case there is outlier model at same level, need to set constant priors for outlier class
                    nifti_image * PriorsConstant=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                    //                nifti_set_filenames(PriorsConstant, (g_scratch_dir + "/ConstantPriors.nii.gz").c_str(), 0, 0);
                    //                nifti_image_write(PriorsConstant);
                    
                    TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstant,2);
                }
            }
                break;
            case 3:{ // Outliers model 3
                nifti_image * PriorsConstantOut=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                nifti_image * PriorsConstantIn=TreeToSegment->BuildConstantPriors(1-segment_param->OutliersWeight);
//                                TreeToSegment->SaveTmpResult(static_cast<float*>( VectorPriors[3]->data), g_scratch_dir + "/PrePriorsTest.nii.gz");
                if (segment_param->flag_IOAverage) {
                    nifti_image_free(PriorsConstantOut);
                    nifti_image_free(PriorsConstantIn);
                    nifti_image * PriorsConstantInInit=ReadFromFilename(segment_param->filename_IOAverage[0]);
                    nifti_image * PriorsConstantOutInit=ReadFromFilename(segment_param->filename_IOAverage[1]);
                    PriorsConstantIn=TreeToSegment->AddNii(PriorsConstantInInit, -1*segment_param->OutliersWeight);
                    PriorsConstantOut=TreeToSegment->AddNii(PriorsConstantOutInit, segment_param->OutliersWeight);
                    nifti_set_filenames(PriorsConstantOut, (g_scratch_dir + "/PriorsOutInit").c_str(), 0, 0);
                    nifti_set_filenames(PriorsConstantIn, (g_scratch_dir + "/PriorsInInit").c_str(), 0, 0);
                    nifti_image_write(PriorsConstantIn);
                    nifti_image_write(PriorsConstantOut);
                    
                }
                TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantIn,0);
                TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantOut,0);
                vector<nifti_image *> DoubledVectorPriors;
                for (int c=0; c<numbPriors; c++) {
                    nifti_image * CopiedPriors=nifti_copy_nim_info(VectorPriors[c]);
                    CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                    float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                    float * PriorsToCopy_PTR=static_cast<float *>(VectorPriors[c]->data);
                    int numbvox=CopiedPriors->nvox;
                    for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                        *CopiedPriors_PTR=*PriorsToCopy_PTR;
                    }
                    DoubledVectorPriors.push_back(CopiedPriors);
                    
                }
//                TreeToSegment->SaveTmpResult(static_cast<float*>( VectorPriors[3]->data), g_scratch_dir + "/PrePriorsTest.nii.gz");
                for (int c=0; c<numbPriors; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,DoubledVectorPriors[c],2);
                }
//                nifti_image * TestPre=VectorPriors[3];
//                nifti_set_filenames(TestPre, (g_scratch_dir + "/PrePrePriorsTest.nii.gz").c_str(), 0, 0);
//                nifti_image_write(TestPre);
//                TreeToSegment->SaveTmpResult(static_cast<float*>( DoubledVectorPriors[3]->data), g_scratch_dir + "/PrePriorsTest.nii.gz");
                if (segment_param->flag_OutBrain) {
                    nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                        nifti_image * CopiedPriors=nifti_copy_nim_info(OutBrainPriors);
                        CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                        float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                        float * PriorsToCopy_PTR=static_cast<float *>(OutBrainPriors->data);
                        int numbvox=CopiedPriors->nvox;
                        for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                            *CopiedPriors_PTR=*PriorsToCopy_PTR;
                        }
                    
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,CopiedPriors,2);
                }
//                nifti_image * Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
//                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
//                nifti_image_write(Test);
//                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                TreeToSegment->NormalisePriors();
                TreeToSegment->NormalisePriorsAdapted();
//                Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
//                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
//                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
//                nifti_image_write(Test);
            }
                break;
                
            case 5:{ // Outliers model 5 (outliers priors from combination of inlier priors)
                nifti_image * PriorsConstantOut=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                nifti_image * PriorsConstantIn=TreeToSegment->BuildConstantPriors(1-segment_param->OutliersWeight);
                TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantIn,0);
                TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantOut,0);
                vector<nifti_image *> DoubledVectorPriors;

                for (int c=0; c<numbPriors; c++) {
                    nifti_image * CopiedPriors=nifti_copy_nim_info(VectorPriors[c]);
                    CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                    float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                    float * PriorsToCopy_PTR=static_cast<float *>(VectorPriors[c]->data);
                    int numbvox=CopiedPriors->nvox;
                    for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                        *CopiedPriors_PTR=*PriorsToCopy_PTR;
                    }
                    DoubledVectorPriors.push_back(CopiedPriors);
                    
                }
                for (int c=0; c<numbPriors; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
//                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,DoubledVectorPriors[c],2);
                }
                int numbOutliersPriors=segment_param->OutliersCombined.size(); // checked before that it is valid
                for(int o=0;o<numbOutliersPriors;o++){
                    nifti_image * OutlierPrior = CreateNewPrior(DoubledVectorPriors,segment_param->OutliersCombined[o]);
                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param, OutlierPrior, 2);
                }
                if (segment_param->flag_OutBrain) {
                    nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                    nifti_image * CopiedPriors=nifti_copy_nim_info(OutBrainPriors);
                    CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                    float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                    float * PriorsToCopy_PTR=static_cast<float *>(OutBrainPriors->data);
                    int numbvox=CopiedPriors->nvox;
                    for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                        *CopiedPriors_PTR=*PriorsToCopy_PTR;
                    }
                    
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,CopiedPriors,2);
                }
                                nifti_image * Test=TreeToSegment->GetChild(1)->GetChild(0)->GetPriorsDirect();
                                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                               nifti_image_write(Test);
                //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                TreeToSegment->NormalisePriors();
                TreeToSegment->NormalisePriorsAdapted();
                //                Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
                //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                //                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                //                nifti_image_write(Test);
            }
                break;
            case 6:{ // Outliers model 6 (outliers priors in priors_out)

                nifti_image * PriorsConstantOut=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                nifti_image * PriorsConstantIn=TreeToSegment->BuildConstantPriors(1-segment_param->OutliersWeight);
                TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantIn,0);
                TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantOut,0);

                for (int c=0; c<numbPriors; c++) { // creation of children inliers with priors
                    std::cout<<"Child number "<< c << "added"<<endl;
                    
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                    //                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,DoubledVectorPriors[c],2);
                }
                int numbOutliersPriors=VectorPriorsOut.size(); // checked before that it is valid
                for(int o=0;o<numbOutliersPriors;o++){
                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param, VectorPriorsOut[o], 2);
                }
                if (segment_param->flag_OutBrain) {
                    nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                    nifti_image * CopiedPriors=nifti_copy_nim_info(OutBrainPriors);
                    CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                    float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                    float * PriorsToCopy_PTR=static_cast<float *>(OutBrainPriors->data);
                    int numbvox=CopiedPriors->nvox;
                    for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                        *CopiedPriors_PTR=*PriorsToCopy_PTR;
                    }
                    
                    TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                    TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,CopiedPriors,2);
                }
                nifti_image * Test=TreeToSegment->GetChild(1)->GetChild(0)->GetPriorsDirect();
                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                nifti_image_write(Test);
                //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                TreeToSegment->NormalisePriors();
                TreeToSegment->NormalisePriorsAdapted();
                //                Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
                //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                //                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                //                nifti_image_write(Test);
            }
                break;
                
                case 7:{ // Outliers model 7 (Priors for 3rd level separating brain and non brain part)
                    
                    nifti_image * PriorsConstantOut=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                    nifti_image * PriorsConstantIn=TreeToSegment->BuildConstantPriors(1-segment_param->OutliersWeight);
                    TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantIn,0);
                    TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantOut,0);
                    vector<nifti_image *> VectorBNBPriors=ReadFromFilenamesVector(segment_param->filename_priors_bnb);
                    for (int l=0; l<2; l++) {
                        for (int bnb=0; bnb<2; bnb++) {
                            if (l==0) {
                                TreeToSegment->GetChild(l)->CreateAndAddChildPriors(segment_param, VectorBNBPriors[bnb], 0);
                            }
                            else{
                                nifti_image * CopyBNB=TreeToSegment->GetChild(0)->GetChild(bnb)->CopyPriors();
                                TreeToSegment->GetChild(l)->CreateAndAddChildPriors(segment_param, CopyBNB, 0);
                            }
                            for (int c=0; c<numbPriors; c++) { // creation of children inliers with priors
                                std::cout<<"Child number "<< c << "added"<<endl;
                                int Dist=l==0?1:2; // Giving the appropriate distribution type for the children
                                if (l==0 && bnb==0) {
                                    TreeToSegment->GetChild(l)->GetChild(bnb)->CreateAndAddChildPriors(segment_param,VectorPriors[c],Dist);
                                }
                                else{
                                    nifti_image * CopyPriors=TreeToSegment->GetChild(0)->GetChild(0)->GetChild(c)->CopyPriors();
                                    TreeToSegment->GetChild(l)->GetChild(bnb)->CreateAndAddChildPriors(segment_param, CopyPriors, Dist);
                                }
                            }
                        }
                    }
                    if (segment_param->flag_OutBrain) {
                        nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                        nifti_image * CopiedPriors=nifti_copy_nim_info(OutBrainPriors);
                        CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                        float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                        float * PriorsToCopy_PTR=static_cast<float *>(OutBrainPriors->data);
                        int numbvox=CopiedPriors->nvox;
                        for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                            *CopiedPriors_PTR=*PriorsToCopy_PTR;
                        }
                        for (int bnb=0; bnb<2; bnb++) {
                            TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                            TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,CopiedPriors,2);
                        }

                    }
//                    nifti_image * Test=TreeToSegment->GetChild(1)->GetChild(0)->GetPriorsDirect();
//                    nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
//                    nifti_image_write(Test);
                    //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                    TreeToSegment->NormalisePriors();
                    TreeToSegment->NormalisePriorsAdapted();
                    //                Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
                    //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                    //                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                    //                nifti_image_write(Test);
                }
                break;
                case 8:{ // Outliers model 8:
                    nifti_image * PriorsConstantOut=TreeToSegment->BuildConstantPriors(segment_param->OutliersWeight);
                    nifti_image * PriorsConstantIn=TreeToSegment->BuildConstantPriors(1-segment_param->OutliersWeight);
//                    TreeToSegment->SaveTmpResult(static_cast<float*>( VectorPriors[3]->data), g_scratch_dir + "/PrePriorsTest.nii.gz");
                    TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantIn,0);
                    TreeToSegment->CreateAndAddChildPriors(segment_param,PriorsConstantOut,0);
                    if (segment_param->flag_bnbComb) {
                        int numbBrainPriors=segment_param->vecBP.size();
                        bool * BP=new bool[segment_param->numb_classes]; // boolean array telling which prior should be considered under brain mask and therefore classified as it...
                        for (int c=0; c<segment_param->numb_classes; c++) {
                            BP[c]=0;
                        }
                        for (int nbp=0; nbp<numbBrainPriors; nbp++) {
                            if (segment_param->vecBP[nbp]>segment_param->numb_classes) {
                                cout << "Impossible priors stated in bnb Comb"<<endl;
                                return NULL;
                            }
                            BP[segment_param->vecBP[nbp]]=1;

                        }
                        vector<nifti_image*>VectorBNBPriors=ReadFromFilenamesVector(segment_param->filename_priors_bnb);
                        for (int c=0; c<segment_param->numb_classes; c++) {
                            TreeToSegment->MultiplyNii(VectorPriors[c], VectorBNBPriors[!BP[c]]);
                            nifti_image * CopiedPriors=nifti_copy_nim_info(VectorPriors[c]);
                            CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                            float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                            float * PriorsToCopy_PTR=static_cast<float *>(VectorPriors[c]->data);
                            int numbvox=CopiedPriors->nvox;
                            for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                                *CopiedPriors_PTR=*PriorsToCopy_PTR;
                            }
                            TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                            TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,CopiedPriors,2);
                        }
                    }
                    else{
                    vector<nifti_image *> DoubledVectorPriors;
                    for (int c=0; c<numbPriors; c++) {
                        nifti_image * CopiedPriors=nifti_copy_nim_info(VectorPriors[c]);
                        CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                        float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                        float * PriorsToCopy_PTR=static_cast<float *>(VectorPriors[c]->data);
                        int numbvox=CopiedPriors->nvox;
                        for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                            *CopiedPriors_PTR=*PriorsToCopy_PTR;
                        }
                        DoubledVectorPriors.push_back(CopiedPriors);
                        
                    }
//                    TreeToSegment->SaveTmpResult(static_cast<float*>( VectorPriors[3]->data), g_scratch_dir + "/PrePriorsTest.nii.gz");
                    for (int c=0; c<numbPriors; c++) {
                        std::cout<<"Child number "<< c << "added"<<endl;
                        
                        TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                        TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,DoubledVectorPriors[c],2);
                    }
//                    nifti_image * TestPre=VectorPriors[3];
//                    nifti_set_filenames(TestPre, (g_scratch_dir + "/PrePrePriorsTest.nii.gz").c_str(), 0, 0);
//                    nifti_image_write(TestPre);
//                    TreeToSegment->SaveTmpResult(static_cast<float*>( DoubledVectorPriors[3]->data), g_scratch_dir + "/PrePriorsTest.nii.gz");
                    }
                    if (segment_param->flag_OutBrain) {
                        nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                        nifti_image * CopiedPriors=nifti_copy_nim_info(OutBrainPriors);
                        CopiedPriors->data = (void *) calloc(CopiedPriors->nvox, sizeof(float));
                        float * CopiedPriors_PTR=static_cast<float *>(CopiedPriors->data);
                        float * PriorsToCopy_PTR=static_cast<float *>(OutBrainPriors->data);
                        int numbvox=CopiedPriors->nvox;
                        for (int i=0; i<numbvox; i++,CopiedPriors_PTR++,PriorsToCopy_PTR++) {
                            *CopiedPriors_PTR=*PriorsToCopy_PTR;
                        }
                        
                        TreeToSegment->GetChild(0)->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                        TreeToSegment->GetChild(1)->CreateAndAddChildPriors(segment_param,CopiedPriors,2);
                    }
                    //                nifti_image * Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
                    //                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                    //                nifti_image_write(Test);
                    //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                    TreeToSegment->SetFlagOutliers(3); // Afterwards behaves almost as with OM3. Only Change will occur with the adaptation of the priors since the mask will be applied again after smoothing
                    TreeToSegment->NormalisePriors();
                    TreeToSegment->NormalisePriorsAdapted();
                    //                Test=TreeToSegment->GetChild(0)->GetPriorsDirect();
                    //                TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
                    //                nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
                    //                nifti_image_write(Test);
                }
                break;

                
            default:{
                for (int c=0; c<numbPriors; c++) {
                    std::cout<<"Child number "<< c << "added"<<endl;
                    TreeToSegment->CreateAndAddChildPriors(segment_param,VectorPriors[c],1);
                }
                if (segment_param->flag_OutBrain) {
                    nifti_image * OutBrainPriors=TreeToSegment->CreateOutBrainPriors(segment_param);
                    TreeToSegment->CreateAndAddChildPriors(segment_param,OutBrainPriors,1);
                }

                TreeToSegment->NormalisePriors();
                TreeToSegment->NormalisePriorsAdapted();
            }
                break;
        }
    }
//    nifti_image * Test=TreeToSegment->GetChild(1)->GetChild(0)->GetPriorsDirect();
//    nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
//    nifti_image_write(Test);
    TreeToSegment->NormalisePriors();
//    TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
    TreeToSegment->NormalisePriorsAdapted();
//    Test=TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsDirect();
//    nifti_set_filenames(Test, (g_scratch_dir + "/Test.nii.gz").c_str(), 0, 0);
//    nifti_image_write(Test);
//    TreeToSegment->SaveTmpResult(TreeToSegment->GetChild(0)->GetChild(0)->GetPriorsAdaptedDirect(), g_scratch_dir + "/Test.nii.gz");
    std::cout<<"Priors normalised now "<< TreeToSegment->ArePriorsNormalised()<<endl;
    std::cout<<"Priors adapted also normalised "<<TreeToSegment->ArePriorsAdaptedNormalised()<<endl;
    if(segment_param->flag_DP==1){
        vector<int *> CountChildrenInput=ReadCountfromFiles(segment_param->filename_DP);
        float * DPChildrenInput=NULL;
        int numbchild=CountChildrenInput.size();
        if(numbchild==TreeToSegment->GetNumberChildren()){
            TreeToSegment->SetFlagDistClassInd(1);
            DPChildrenInput=new float[MaxSupport*numbchild];
            for(int i=0;i<MaxSupport*numbchild;i++){
                DPChildrenInput[i]=0;
            }
            for(int c=0;c<numbchild;c++){
                float * DP_tmp=PoissonRegressionKernelSmoothing(CountChildrenInput[c],segment_param->smoothing_order,segment_param);
                for(int i=0;i<MaxSupport;i++){
                    DPChildrenInput[c*MaxSupport+i]=DP_tmp[i];
                }
                if(DP_tmp!=NULL){
                    delete[] DP_tmp;
                    DP_tmp=NULL;
                }
            }
        }
        else{
            vector<int*> CountChildrenInput=ReadCountFromDistributionTextFile(segment_param->filename_DP[0]);
            if(!segment_param->flag_DistClassInd){
                TreeToSegment->SetFlagDistClassInd(0);
                vector<int*> DistributionChildrenInput=ReadDistributionFromTextFile(segment_param->filename_DP[0]);
                float * CountHistogram= MultivariateJointHistogram(CountChildrenInput,DistributionChildrenInput);
                vector<int> DimToUse;
                int numbClasses=CountChildrenInput.size();
                for(int c=0;c<numbClasses;c++){
                    DimToUse.push_back(MaxSupport);
                }
                //                float gauss_std=1;
                float * DPChildrenInput=GaussianBlurring(CountHistogram, segment_param->DPGauss_std, DimToUse);
                int SizeHistogram=pow(MaxSupport,numbClasses);
                // Normalisation of the gaussian filtered distribution
                float sumHistBlurred=0;
                for(int i=0;i<SizeHistogram;i++){
                    sumHistBlurred+=DPChildrenInput[i];
                }
                for(int i=0;i<SizeHistogram;i++){
                    DPChildrenInput[i]/=sumHistBlurred;
                }
                delete [] CountHistogram;
                for(int c=0;c<numbClasses;c++){
                    delete[] DistributionChildrenInput[c];
                    delete[] CountChildrenInput[c];
                }
                TreeToSegment->SetDPChildren(DPChildrenInput);
            }
            if(segment_param->flag_DistClassInd){
                TreeToSegment->SetFlagDistClassInd(1);
                int Csize=CountChildrenInput.size();
                float * DPChildrenInput=new float[Csize*MaxSupport];
                for(int c=0;c<Csize;c++){
                    float * DPChildrenInput_tmp=PoissonRegressionKernelSmoothing(CountChildrenInput[c],segment_param->smoothing_order,segment_param);
                    for(int i=0;i<MaxSupport;i++){
                        DPChildrenInput[c*MaxSupport+i]=DPChildrenInput_tmp[i];
                    }
                    delete [] DPChildrenInput_tmp;
                    DPChildrenInput_tmp=NULL;
                    delete [] CountChildrenInput[c];
                    CountChildrenInput[c]=NULL;
                }
                CountChildrenInput.clear();
                TreeToSegment->SetDPChildren(DPChildrenInput);
            }
        }
    }

    TreeToSegment->SetFlagCovPriors( segment_param->CovPriorsType);
    TreeToSegment->SetFlagMeanPriors(segment_param->MeanPriors);
    float * PriorsCov=NULL;
    switch (TreeToSegment->GetFlagCovPriors()) {
        case 2:
            PriorsCov=TreeToSegment->GetPriorsCovMatrixGeneral();
            break;
        case 3:
            PriorsCov=NoiseCovPriorsOtsu(TreeToSegment->GetDataImage());
            break;
        default:
            break;
    }
    
    TreeToSegment->SetCovPriorsMatrix(PriorsCov);
    
    TreeToSegment->InitialiseBasicTree(segment_param);
    
    return TreeToSegment;
}




// Returns the index of the first non NULL pointer in the ImagesVector considered
int FirstNotNULL(vector<nifti_image *> ImagesVector){
    if (ImagesVector.size()==0){
        std::cout<< "No image in vector"<<endl;
        return -1;
    }
    int finn = ImagesVector.size();
    for (int i=ImagesVector.size()-1; i>=0; i--) {
        if (ImagesVector[i]!=NULL){
            finn=i;
        }
    }
    return finn;
}

// Returns the index of the last non NULL pointer to a nifti_image in the ImageVector considered
int LastNotNULL(vector<nifti_image *> ImagesVector){
    if (ImagesVector.size()==0) {
        return -2;
    }
    int linn=0;
    int sizeVector=ImagesVector.size();
    for (int i=0; i<sizeVector; i++) {
        if (ImagesVector[i]!=NULL) {
            linn=i;
        }
    }
    return linn;
}

// Read the file where the count of the occurences of number of classes are kept and give it back under a vector of int arrays
vector<int *> ReadDistributionFromTextFile(string DistributionTextFile){
    ifstream DistributionText (DistributionTextFile.c_str());
    vector<int *> DistributionVector;
    if(!DistributionText){
        cout<<"Pb in opening file"<<endl;
        return DistributionVector;
    }
    else {
        string line;
        
        while(getline(DistributionText,line)){
            vector<int> LineContent;
            istringstream in(line);
            int num=0;
            while(in >> num){
                LineContent.push_back(num);
            }
            // convert the vector obtained to an int array
            int lenline=LineContent.size();
            int * DistributionArray=new int[lenline];
            for(int i=0;i<lenline;i++){
//                cout<< LineContent[i]<<" ";
                DistributionArray[i]=LineContent[i];
            }
//            cout<<endl;
            if(lenline>0){
                DistributionVector.push_back(DistributionArray);
            }
        }
        return DistributionVector;
    }
}

// Read the file where the count of the occurences of number of classes are kept and give it back under a vector of int arrays
vector<int *> ReadCountFromDistributionTextFile(string DistributionTextFile){
    ifstream DistributionText (DistributionTextFile.c_str());
    vector<int *> CountsVector;
    if(!DistributionText){
        cout<<"Pb in opening file"<<endl;
        return CountsVector;
    }
    else {
        string line;
        while(getline(DistributionText,line)){
            int * Counts=new int[MaxSupport];
            for(int i=0;i<MaxSupport;i++){
                Counts[i]=0;
            }
            vector<int> LineContent;
            istringstream in(line);
            int num=0;
            while(in >> num){
                Counts[num-1]++;
            }
            int lenline=line.length();
            if(lenline>0){
                //                for(int i=0;i<MaxSupport;i++){
                //                    cout<<Counts[i]<<" ";
                //                }
                //                cout<<endl;
                CountsVector.push_back(Counts);
            }
        }
        return CountsVector;
    }
    
}


TreeEM * ReadTreeFromTextFileWithAdapt(string TreeTextFile, SEG_PARAMETERS * segment_param, vector<nifti_image *> PriorsGCVector, vector<nifti_image*>PriorsIOVector, vector<string> AdaptTextFile){
//    First take care of possible adaptation with AdaptTextFile
    int typeAdaptTextFile=AdaptTextFile.size();
    
    float ConstAdapt[MaxNumbModal];
    float LinearAdapt[MaxNumbModal];
    float MinMaxStore1[MaxNumbModal*2];
    float MinMaxStore2[MaxNumbModal*2];
    for(int m=0;m<MaxNumbModal;m++){
        MinMaxStore1[m*2]=0;
        MinMaxStore1[m*2+1]=-1;
        MinMaxStore2[m*2]=0;
        MinMaxStore2[m*2+1]=-1;
    }
    bool flag_AdaptNeeded=0;

        switch (typeAdaptTextFile) {
                
            case 1:{ // Case where the coefficients are already written in the file
                ifstream textAdapt (AdaptTextFile[0].c_str());
                if(!textAdapt){
                    std::cout<<"could not open the text file properly for the params adaptation!"<<endl;
                }
                else{
                std::string lineAdapt;
                while(getline(textAdapt,lineAdapt)){
                    istringstream in(lineAdapt);
                    std:: string type;
                    in >> type;
                    if(type == "Linear"){
                        flag_AdaptNeeded=1;
                        float LinearWeight=0;
                        int m=0;
                        while (in >> LinearWeight) {
                            LinearAdapt[m]=LinearWeight;
                            m++;
                        }
                    }
                    else if(type== "Constant"){
                        flag_AdaptNeeded=1;
                        float ConstWeight=0;
                        int m=0;
                        while (in >> ConstWeight) {
                            ConstAdapt[m]=ConstWeight;
                            m++;
                        }
                    }
                }
                }
            }
                break;
            case 2:{ // Case where the minimum and the maximum along with the modality are given in the text file. The text file gives : Name of Modality, min and max for this modality on the inlier part,
                ifstream textAdapt1 (AdaptTextFile[0].c_str());
                ifstream textAdapt2 (AdaptTextFile[1].c_str());
                if(!textAdapt1 || !textAdapt2){
                    std::cout<<"could not open the text file properly for the params adaptation!"<<endl;
                }
                std::string lineAdapt;
                while(getline(textAdapt1,lineAdapt)){
                    istringstream in(lineAdapt);
                    std:: string Modality;
                    in >> Modality;
                    float min;
                    float max;
                    in >> min;
                    in >> max;
                    if (Modality=="FLAIR") {
                        MinMaxStore1[4]=min;
                        MinMaxStore1[5]=max;
                    }
                    if(Modality == "T2"){
                        MinMaxStore1[2]=min;
                        MinMaxStore1[3]=max;
                    }
                    else if(Modality== "T1"){
                        MinMaxStore1[0]=min;
                        MinMaxStore1[1]=max;
                    }
                }
                while(getline(textAdapt2,lineAdapt)){
                    istringstream in(lineAdapt);
                    std:: string Modality;
                    in >> Modality;
                    float min;
                    float max;
                    in >> min;
                    in >> max;
                    if (Modality=="FLAIR") {
                        MinMaxStore2[4]=min;
                        MinMaxStore2[5]=max;
                    }
                    if(Modality == "T2"){
                        MinMaxStore2[2]=min;
                        MinMaxStore2[3]=max;
                    }
                    else if(Modality== "T1"){
                        MinMaxStore2[0]=min;
                        MinMaxStore2[1]=max;
                    }
                }
//                Initialise Linear and Const with 1 or 0
                for (int m=0; m<MaxNumbModal; m++) {
                    LinearAdapt[m]=1;
                    ConstAdapt[m]=1;
                }
//                Once read, build the Linear and Const weight assuming that there is compatibility in the
                for (int m=0; m<MaxNumbModal; m++) {
                    if (MinMaxStore1[2*m+1]>0 && MinMaxStore2[2*m+1]>0 ) { // Meaning that there is the same modality m on which to make linear combination, transforming 2 into 1. 2 correspond then to the file of the mean while 1 correspond to file of the subject
                        LinearAdapt[m]=(MinMaxStore1[2*m+1]-MinMaxStore1[2*m])/(MinMaxStore2[2*m+1]-MinMaxStore2[2*m]);
                        cout<<"for modality "<<m<<" linear is "<<LinearAdapt[m]<<endl;
                        ConstAdapt[m]=MinMaxStore1[2*m]-MinMaxStore2[2*m]*LinearAdapt[m];
                        flag_AdaptNeeded=1;
                    }
                }
            }
                break;
            default:
                break;
        }
    

    
    if(segment_param==NULL){
        segment_param=new SEG_PARAMETERS();
        segment_param->quantMax=1;
        segment_param->quantMin=0;
    }
//    int numbmodal=0;
    //        int numbGchild=0;
    
    
    
    // Second, build Tree
    
    ifstream text (TreeTextFile.c_str());
    //    text.open(TreeTextFile,ios::in);
    if(!text){
        std::cout<<"could not open the text file properly !"<<endl;
        return NULL;
    }
    else{
        TreeEM * TreeResult=new TreeEM();
        nifti_image * MRFImage=NULL;
        nifti_image * ImageCorrectedToSegment=NULL;
        string DataCorrectedFilename=segment_param->filename_inDC;
        ImageCorrectedToSegment=ReadFromFilename(DataCorrectedFilename);
        TreeResult->SetData(ImageCorrectedToSegment);
//        TreeResult->NormaliseDataImageNonLog();
        string MaskFilename=segment_param->filename_mask;
        nifti_image* MaskToUse=ReadFromFilename(MaskFilename);
        segment_param->flag_mask=1;
        segment_param->filename_mask=MaskFilename;
        TreeResult->SetMask(MaskToUse);
        TreeResult->CreateAllocateAndInitializeParameters(0);
        if (ImageCorrectedToSegment!=NULL) {
            int * L2S_PTR=TreeResult->GetL2S();
            float * DataCorrected_PTR=static_cast<float *>(ImageCorrectedToSegment->data);
            int numelmasked=TreeResult->GetNumberMaskedElements();
            int numel=TreeResult->GetNumberElements();
            int numbmodal=TreeResult->GetNumberModalities();
            float * DataCorrectedToSet=new float[numbmodal*numelmasked];
            for (int m=0; m<numbmodal; m++) {
                int j=0;
                for (int i=0; i<numel; i++) {
                    if (L2S_PTR[i]>=0) {
                        DataCorrectedToSet[j+m*numelmasked]=DataCorrected_PTR[m*numel+i];
                        j++;
                    }
                }
            }
            TreeResult->SetDataBFCorrected(DataCorrectedToSet);
            delete [] DataCorrectedToSet;
            DataCorrectedToSet=NULL;
            TreeResult->UncorrectData(segment_param->bias_order);

//            nifti_image * TestUncorrected=TreeResult->GetDataImage();
//            nifti_set_filenames(TestUncorrected, (g_scratch_dir + "/Uncorrected.nii.gz").c_str(), 0, 0);
//            nifti_image_write(TestUncorrected);
            cout<<"Uncorrected saved"<<endl;                     cout<<"Testing"<<endl;
            TreeResult->SetFlagCovPriors( segment_param->CovPriorsType);
            TreeResult->SetFlagMeanPriors(segment_param->MeanPriors);
            float * PriorsCov=NULL;
            switch (TreeResult->GetFlagCovPriors()) {
                case 2:
                    PriorsCov=TreeResult->GetPriorsCovMatrixGeneral();
                    break;
                case 3:
                    PriorsCov=NoiseCovPriorsOtsu(TreeResult->GetDataImage());
                    break;
                default:
                    break;
            }
            
            TreeResult->SetCovPriorsMatrix(PriorsCov);
        }
        
        std::string line;
        if(segment_param==NULL){
            segment_param=new SEG_PARAMETERS();
            segment_param->quantMax=1;
            segment_param->quantMin=0;
        }
//        int numbmodal=0;
        //        int numbGchild=0;
        while(getline(text,line)){
            istringstream in(line);
            std:: string type;
            int level;
            int numbElements=0;
            in >> type;
            in >> level;
            in >> numbElements;
            if(level==1){
                segment_param->numb_classes=numbElements;
            }
//            bool pastRoot=0;
            int l=0;
            while(!line.empty()){
                //            while(type !="Level" || !pastRoot){
//                pastRoot=1;
                //            in.clear();
                getline(text,line);
                l++;
                //                cout<<line;
                istringstream in2(line);
                
                //                cout<< line.empty() <<endl;
                if(!line.empty()){
                    in2 >> type;
                    cout << type;
                }
                else{
                    type.clear();
                }
                if(level==0){ // Parsing of the root level
                    
                    if(type=="OutliersMode"){ // Parsing of the outliers mode
                        int OutliersMode=0;
                        in2 >> OutliersMode;
                        segment_param->OutliersMod=OutliersMode;
                        segment_param->flag_Outliers=(OutliersMode>0);
                        TreeResult->SetFlagOutliers(OutliersMode);
                    }
                    else if(type == "MRFImage"){
                        string MRFfilename;
                        in2>> MRFfilename;
                        MRFImage=ReadFromFilename(MRFfilename);
                    }
//                    else if(type =="GMatrixIn"){
//                        string GMatrixInFilename;
//                        in2>> GMatrixInFilename;
//                        segment_param->flag_GMatrixIn=1;
//                        segment_param->filename_GMatrix=GMatrixInFilename;
//                    }
//                    else if(type =="GMatrixPost"){
//                        string GMatrixPostFilename;
//                        in2 >> GMatrixPostFilename;
//                        segment_param->flag_GMatrixPost=1;
//                        segment_param->filename_GMatrixPost=GMatrixPostFilename;
//                    }
                    else if(type =="AtlasWeight"){
                        float AtlasWeight;
                        while (in2 >> AtlasWeight) {
                            segment_param->AtlasWeight.push_back(AtlasWeight);
                        }
                        //                        in2 >> AtlasWeight;
                        //                        segment_param->AtlasWeight=AtlasWeight;
                    }
                    else if(type =="AtlasSmoothing"){
                        float AtlasSmoothing;
                        while (in2 >> AtlasSmoothing){
                            segment_param->AtlasSmoothing.push_back(AtlasSmoothing);
                        }
                    }
                    
                }
                else if(level>=1){
                    for(int c=0;c<numbElements;c++){
                        cout<<c<<" child in creation"<<endl;
                        TreeEM * ChildToAdd=ReadChildLevel(text,TreeResult,segment_param);
                        if(ChildToAdd->GetParent()->GetParameters()==NULL){
                            ChildToAdd->GetParent()->SetParameters(new Parameters());
                        }
                        if(ChildToAdd->GetDistributionType()==1 &&segment_param->flag_CovPriors){
                            int numbmodal=TreeResult->GetNumberModalities();
                            int numbmodalSq=numbmodal*numbmodal;
                            float * CovPriorsMatrix=new float[numbmodalSq];
                            for(int i=0;i<numbmodalSq;i++){
                                CovPriorsMatrix[i]=0;
                            }
                            if (ChildToAdd->GetParametersValue() !=NULL && segment_param->flag_meanPriors) {
                                float * ParametersValue=ChildToAdd->GetParametersValue();
                                for (int i=0; i<numbmodalSq; i++) {
                                    CovPriorsMatrix[i]=ParametersValue[numbmodal+i];
                                }
                            }
                            ChildToAdd->SetFlagCovPriorsParameters(1);
                            ChildToAdd->SetCovPriorsMatrix(CovPriorsMatrix);
                            if(CovPriorsMatrix !=NULL){
                                delete [] CovPriorsMatrix;
                                CovPriorsMatrix=NULL;
                            }
                        }

                        ChildToAdd->GetParent()->AddChild(ChildToAdd);
                    }
                }
            }
        }
        
        
        if (flag_AdaptNeeded) {
            TreeResult->AdaptParametersToAffineTransform(LinearAdapt,ConstAdapt);
        }
        vector<TreeEM *> LeavesVector=TreeResult->GetAllLeaves();
        int CountDeleted=0;
        int numbLeaves=LeavesVector.size();
        for (int l=0; l<numbLeaves; l++) {
            bool flag_CompatibleParams=LeavesVector[l]->CheckCompatibilityParams();
            flag_CompatibleParams=1;
            if (!flag_CompatibleParams && TreeResult->GetFlagOutliers()==3) {
                vector<int> HierarchyLeaf=LeavesVector[l]->GetHierarchyVector();
                if (HierarchyLeaf.size()==2) {
                    cout<<"Pb in Params adaptation"<<endl;
                }
                else if(HierarchyLeaf[0]==0){
                    cout<<"Pb in params adaptation normal inliers"<<endl;
                }
                else{
                    delete LeavesVector[l];
                    LeavesVector[l]=NULL;
                    CountDeleted++;
                }
            }
        }
        TreeResult->DeleteAllNULLChildren();
        TreeResult->CollapseOnlyChildTot();
        cout<<"All children added"<<endl;
        
        vector<TreeEM *> VectorLeaves=TreeResult->GetAllLeaves();
        int numbleaves=VectorLeaves.size();
        for (int l=0; l<numbleaves; l++) {
            if(VectorLeaves[l]->GetDistributionType()==1 &&segment_param->MeanPriors>0.01){
                int numbmodal=TreeResult->GetNumberModalities();
                float * MeanPriorsMatrix2=new float[numbmodal];
                for(int m=0;m<numbmodal;m++){
                    MeanPriorsMatrix2[m]=0;
                }
                float * PossibleMean=VectorLeaves[l]->GetMean();
                if (PossibleMean!=NULL) {
                    cout<<"Possible to input mean for leave "<<l<< " !"<<endl;
                    for(int i=0;i<numbmodal;i++){
                        MeanPriorsMatrix2[i]=PossibleMean[i];
                    }
                    VectorLeaves[l]->SetMeanPriorsMatrix(MeanPriorsMatrix2);

                    cout<<"Mean added"<<endl;
                }
                if(MeanPriorsMatrix2!=NULL){
                    delete [] MeanPriorsMatrix2;
                    MeanPriorsMatrix2=NULL;
                }
                MeanPriorsMatrix2=NULL;
                
            }
        }

        

        //        cout<<"MRFimage ...";
        // If MRF modification available, setting at the leaves
        if(MRFImage!=NULL){
            // First check if size of MRFImage compatible with number of leaves
            vector<TreeEM*> LeavesVector=TreeResult->GetAllLeaves();
            int numbLeaves=LeavesVector.size();
            int numelmasked=TreeResult->GetNumberMaskedElements();
            int numel=TreeResult->GetNumberElements();
            float * MRF_PTR=static_cast<float*>(MRFImage->data);
            if(numbLeaves==MRFImage->dim[3]){
                for(int l=0;l<numbLeaves;l++){
                    int * S2L_PTR=TreeResult->GetS2L();
                    float * MRFToSet=new float[numelmasked];
                    float * MRF_PTRtmp=&MRF_PTR[l*numel];
                    for(int i=0;i<numelmasked;i++,S2L_PTR++){
                        MRFToSet[i]=MRF_PTRtmp[*S2L_PTR];
                    }
                    LeavesVector[l]->SetMRF(MRFToSet);
                    if(MRFToSet!=NULL){
                        delete [] MRFToSet;
                        MRFToSet=NULL;
                    }
                }
            }
            nifti_image_free(MRFImage);
        }
        cout<<"... taken care of"<<endl;
        cout<< "Image is already normalised" << endl;
        
        if(!segment_param->flag_MRF){
            vector<TreeEM*> LeavesVector=TreeResult->GetAllLeaves();
            int numbLeaves=LeavesVector.size();
            for(int l=0;l<numbLeaves;l++){
                if(LeavesVector[l]->GetMRF()!=NULL){
                    LeavesVector[l]->SetMRF(NULL);
                }
            }
        }
//        Setting the priors for the general tissue classes if the priors exist
        int sizePriorsGC=PriorsGCVector.size();
        if (sizePriorsGC==TreeResult->GetNumberGeneralClasses()) {
            vector<TreeEM*> GeneralClassesVector=TreeResult->GetGeneralClassesVector();
            int numbGC=GeneralClassesVector.size();
            for (int c=0; c<numbGC; c++) {
                GeneralClassesVector[c]->SetPriors(PriorsGCVector[c]);
                GeneralClassesVector[c]->SetPriorsAdapted(static_cast<float *>(PriorsGCVector[c]->data));
                if (TreeResult->GetFlagOutliers()==3) {
                    nifti_image * CopyPriors=CopyFloatNiiImage(PriorsGCVector[c]);
                    TreeResult->GetNodeOutlier()->GetChild(c)->SetPriors(CopyPriors);
                    TreeResult->GetNodeOutlier()->GetChild(c)->SetPriorsAdapted(static_cast<float *>(PriorsGCVector[c]->data));
                    
                }
            }
        }
        
        // Setting the priors for the inlier outlier separation
        nifti_set_filenames(PriorsIOVector[0], (g_scratch_dir + "/TestPriorsStrange.nii.gz").c_str(), 0, 0);
//        nifti_image_write(PriorsIOVector[0]);
        if (PriorsIOVector.size()==2 && TreeResult->GetFlagOutliers()==3) {
            TreeResult->GetNodeInlier()->SetPriors(PriorsIOVector[0]);
//            TreeEM * TestTreeInlier=TreeResult->GetNodeInlier();
            nifti_set_filenames(TreeResult->GetNodeInlier()->GetPriors(), (g_scratch_dir + "/TestPriorsStrange2.nii.gz").c_str(), 0, 0);
//            nifti_image_write(TreeResult->GetNodeInlier()->GetPriors());
            TreeResult->GetNodeOutlier()->SetPriors(PriorsIOVector[1]);
            TreeResult->GetNodeInlier()->SetPriorsAdapted(static_cast<float *>(PriorsIOVector[0]->data));
            TreeResult->GetNodeOutlier()->SetPriorsAdapted(static_cast<float*>(PriorsIOVector[1]->data));
        }
        TreeResult->SavePriorsAdaptedHierarchy(segment_param);
        
        
        return TreeResult;
    }
}

// Recreate a tree from text file
TreeEM * ReadTreeFromTextFile(string TreeTextFile, SEG_PARAMETERS * segment_param,char * ChangePath=NULL){
    
    cout<<"Name of text file is "<<TreeTextFile<<endl;
    ifstream text (TreeTextFile.c_str());
    //    text.open(TreeTextFile,ios::in);
    if(!text){
        std::cout<<"could not open the text file properly !"<<endl;
        return NULL;
    }
    else{
        TreeEM * TreeResult=new TreeEM();
        nifti_image * MRFImage=NULL;
        nifti_image * ImageCorrectedToSegment=NULL;
        std::string line;
        if(segment_param==NULL){
            segment_param=new SEG_PARAMETERS();
            segment_param->quantMax=1;
            segment_param->quantMin=0;
        }
//        else{
//        segment_param->filename_input.clear();
//        segment_param->filename_priors.clear();
//        segment_param->flag_manual_priors=0;
//        }
        int numbmodal=0;
//        int numbGchild=0;
        while(getline(text,line)){
            istringstream in(line);
            std:: string type;
            int level;
            int numbElements=0;
            in >> type;
            in >> level;
            in >> numbElements;
            if(level==1){
                segment_param->numb_classes=numbElements;
            }
//            bool pastRoot=0;
            int l=0;
            while(!line.empty()){
                //            while(type !="Level" || !pastRoot){
//                pastRoot=1;
                //            in.clear();
                getline(text,line);
                l++;
//                cout<<line;
                istringstream in2(line);
                
//                cout<< line.empty() <<endl;
                if(!line.empty()){
                    in2 >> type;
                    cout << type;
                }
                else{
                    type.clear();
                }
                if(level==0){ // Parsing of the root level
//                    cout<<type;
                    if(type=="Images"){ // Parsing of the images
                        
                        in2 >> numbmodal;
                        vector<string> ImageFiles;
                        //                segment_param->filename_input.clear();
                        string ImageFilename_tmp;
                        for(int i=0;i<numbmodal;i++){
                            getline(text,line);
                            istringstream in3(line);
                            in3>>ImageFilename_tmp;
                            
                            
                            if (ChangePath!=NULL) {
                                int IndexDC=ImageFilename_tmp.find_last_of('/');
                                string FilenameDC_e=ImageFilename_tmp.substr(IndexDC+1,ImageFilename_tmp.length());
                                ImageFilename_tmp=ChangePath+FilenameDC_e;
                            }
                            //                    cout<<ImageFilename_tmp<<endl;
                            // char * ImageFilename=const_cast<string>(ImageFilename_tmp.c_str());
//                            cout<<ImageFilename_tmp<<endl;
                            if (ImageFiles.size()!=0) {
                                if (ImageFilename_tmp!=ImageFiles[0]) {
                                    ImageFiles.push_back(ImageFilename_tmp);
                                }
                                else{
                                    cout<<"same name for file"<<endl;
                                }
                            }

                            else {
                                ImageFiles.push_back(ImageFilename_tmp);
                            }
                            
                            //                    cout<<ImageFiles[i]<<endl;
                        }
//                        segment_param->filename_input=ImageFiles;
//                        //                cout<<segment_param->filename_input[0]<<endl;
//                        segment_param->flag_input=1;
//                        segment_param->numbmodal=numbmodal;
//                        //                cout<<ImageFiles[0]<<endl;
//                        vector<nifti_image*> ImagesToSegment=ReadFromFilenamesVector(segment_param->filename_input);
                        vector<nifti_image*> ImagesToSegment=ReadFromFilenamesVector(ImageFiles);
                        nifti_image * DataImage=CreateDataImage(ImagesToSegment);
                        TreeResult->SetData(DataImage);
                        
                    }
                    else if(type == "BF"){
                        int BForder;
                        in2>>BForder;
                        int numbBF=(int)((BForder+1)*(BForder+2)*(BForder+3))/6;
                        float * BFcoeffs=new float[numbBF * numbmodal];
                        for(int m=0;m<numbmodal;m++){
                            getline(text,line);
                            istringstream in4(line);
                            for(int j=0;j<numbBF;j++){
                                in4 >> BFcoeffs[m*numbBF+j];
                            }
                            
                        }
                        TreeResult->SetBFCoeffs(BFcoeffs);
                        if (BFcoeffs!=NULL && numbBF*numbmodal>0) {
                            delete [] BFcoeffs;
                            BFcoeffs=NULL;
                        }
                    if (ImageCorrectedToSegment!=NULL) {
                        int * L2S_PTR=TreeResult->GetL2S();
                        float * DataCorrected_PTR=static_cast<float *>(ImageCorrectedToSegment->data);
                        int numelmasked=TreeResult->GetNumberMaskedElements();
                        int numel=TreeResult->GetNumberElements();
                        int numbmodal=TreeResult->GetNumberModalities();
                        float * DataCorrectedToSet=new float[numbmodal*numelmasked];
                        for (int m=0; m<numbmodal; m++) {
                            int j=0;
                            for (int i=0; i<numel; i++) {
                                if (L2S_PTR[i]>=0) {
                                    DataCorrectedToSet[j+m*numelmasked]=DataCorrected_PTR[m*numel+i];
                                    j++;
                                }
                            }
                        }
                        TreeResult->SetDataBFCorrected(DataCorrectedToSet);
                        delete [] DataCorrectedToSet;
                        DataCorrectedToSet=NULL;
//                        float * DataTest=static_cast<float *>(TreeResult->GetDataImage()->data);
//                        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+TreeResult->GetNumberElements()]<<endl;
                        TreeResult->UncorrectData(BForder);
//                        DataTest=static_cast<float *>(TreeResult->GetDataImage()->data);
//                        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+TreeResult->GetNumberElements()]<<endl;
//                        nifti_image * TestUncorrected=TreeResult->GetDataImage();
//                        nifti_set_filenames(TestUncorrected, (g_scratch_dir + "/Uncorrected.nii.gz").c_str(), 0, 0);
//                        nifti_image_write(TestUncorrected);
                        cout<<"Uncorrected saved"<<endl;
//                        DataTest=static_cast<float *>(TreeResult->GetDataImage()->data);
//                        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+TreeResult->GetNumberElements()]<<endl;
//                        cout<<"Testing"<<endl;
                    }
//                        TreeResult->MakeDataBFCorrected();
                    }
                    else if(type=="DataCorrected"){
                        string DataCorrectedFilename;
                        in2 >> DataCorrectedFilename;
                        if (ChangePath!=NULL) {
                            int IndexDC=DataCorrectedFilename.find_last_of('/');
                            string FilenameDC_e=DataCorrectedFilename.substr(IndexDC+1,DataCorrectedFilename.length());
                            DataCorrectedFilename=ChangePath+FilenameDC_e;
                        }
                        
                        ImageCorrectedToSegment=ReadFromFilename(DataCorrectedFilename);
                        TreeResult->SetData(ImageCorrectedToSegment);
//                        float * DataTest=static_cast<float *>(TreeResult->GetDataImage()->data);
//                        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+TreeResult->GetNumberElements()]<<endl;
                        // When resetting, need to clear BF coefficients and
                    }
                    else if(type=="Mask"){ // Parsing of the mask
                        string MaskFilename;
                        in2 >> MaskFilename;
                        if (ChangePath!=NULL) {
                            int IndexDC=MaskFilename.find_last_of('/');
                            string FilenameMask_e=MaskFilename.substr(IndexDC+1,MaskFilename.length());
                            MaskFilename=ChangePath+FilenameMask_e;
                        }
                        nifti_image* MaskToUse=ReadFromFilename(MaskFilename);
                        segment_param->flag_mask=1;
                        segment_param->filename_mask=MaskFilename;
                        TreeResult->SetMask(MaskToUse);
                        TreeResult->CreateAllocateAndInitializeParameters(0);
                    }
                    else if(type=="OutliersMode"){ // Parsing of the outliers mode
                        int OutliersMode=0;
                        in2 >> OutliersMode;
                        segment_param->OutliersMod=OutliersMode;
                        segment_param->flag_Outliers=(OutliersMode>0);
                        TreeResult->SetFlagOutliers(OutliersMode);
                    }
                    else if(type == "MRFImage"){
                        string MRFfilename;
                        in2>> MRFfilename;
                        // FIX: Apply ChangePath correction to MRFImage filename (bug fix 08/05/2026)
                        // This ensures the MRF file path is resolved relative to the scratch directory
                        if (ChangePath!=NULL) {
                            int IndexMRF=MRFfilename.find_last_of('/');
                            string FilenameMRF_e=MRFfilename.substr(IndexMRF+1,MRFfilename.length());
                            MRFfilename=ChangePath+FilenameMRF_e;
                        }
                        MRFImage=ReadFromFilename(MRFfilename);
                    }
                    else if(type =="GMatrixIn"){
                        string GMatrixInFilename;
                        in2>> GMatrixInFilename;
                        // FIX: Apply ChangePath correction to GMatrixIn filename (bug fix 08/05/2026)
                        // This ensures the GMatrix file path is resolved relative to the scratch directory
                        if (ChangePath!=NULL) {
                            int IndexGMatrixIn=GMatrixInFilename.find_last_of('/');
                            string FilenameGMatrixIn_e=GMatrixInFilename.substr(IndexGMatrixIn+1,GMatrixInFilename.length());
                            GMatrixInFilename=ChangePath+FilenameGMatrixIn_e;
                        }
                        segment_param->flag_GMatrixIn=1;
                        segment_param->filename_GMatrix=GMatrixInFilename;
                    }
                    else if(type =="GMatrixPost"){
                        string GMatrixPostFilename;
                        in2 >> GMatrixPostFilename;
                        // FIX: Apply ChangePath correction to GMatrixPost filename (bug fix 08/05/2026)
                        // This ensures the GMatrix file path is resolved relative to the scratch directory
                        if (ChangePath!=NULL) {
                            int IndexGMatrixPost=GMatrixPostFilename.find_last_of('/');
                            string FilenameGMatrixPost_e=GMatrixPostFilename.substr(IndexGMatrixPost+1,GMatrixPostFilename.length());
                            GMatrixPostFilename=ChangePath+FilenameGMatrixPost_e;
                        }
                        segment_param->flag_GMatrixPost=1;
                        segment_param->filename_GMatrixPost=GMatrixPostFilename;
                    }
                    else if(type =="AtlasWeight"){
                        float AtlasWeight;
                        while (in2 >> AtlasWeight) {
                            segment_param->AtlasWeight.push_back(AtlasWeight);
                        }
//                        in2 >> AtlasWeight;
//                        segment_param->AtlasWeight=AtlasWeight;
                    }
                    else if(type =="AtlasSmoothing"){
                        float AtlasSmoothing;
                        while (in2 >> AtlasSmoothing){
                        segment_param->AtlasSmoothing.push_back(AtlasSmoothing);
                        }
                    }

                }
                else if(level>=1){
//                    float * DataTest=static_cast<float *>(TreeResult->GetDataImage()->data);
//                    cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+TreeResult->GetNumberElements()]<<endl;
                    for(int c=0;c<numbElements;c++){
                        cout<<c<<" child in creation"<<endl;
                        TreeEM * ChildToAdd=ReadChildLevel(text,TreeResult,segment_param);
                        if(ChildToAdd->GetParent()->GetParameters()==NULL){
                            ChildToAdd->GetParent()->SetParameters(new Parameters());
                        }
                        
                        
                        ChildToAdd->GetParent()->AddChild(ChildToAdd);
                    }
                }
                
            }
        }

        
        cout<<"All children added"<<endl;
//        cout<<"MRFimage ...";
        // If MRF modification available, setting at the leaves
        if(MRFImage!=NULL){
            // First check if size of MRFImage compatible with number of leaves
            vector<TreeEM*> LeavesVector=TreeResult->GetAllLeaves();
            int numbLeaves=LeavesVector.size();
            int numelmasked=TreeResult->GetNumberMaskedElements();
            int numel=TreeResult->GetNumberElements();
            float * MRF_PTR=static_cast<float*>(MRFImage->data);
            if(numbLeaves==MRFImage->dim[3]){
                for(int l=0;l<numbLeaves;l++){
                    int * S2L_PTR=TreeResult->GetS2L();
                    float * MRFToSet=new float[numelmasked];
                    float * MRF_PTRtmp=&MRF_PTR[l*numel];
                    for(int i=0;i<numelmasked;i++,S2L_PTR++){
                        MRFToSet[i]=MRF_PTRtmp[*S2L_PTR];
                    }
                    LeavesVector[l]->SetMRF(MRFToSet);
                    if(MRFToSet!=NULL){
                        delete [] MRFToSet;
                        MRFToSet=NULL;
                    }
                }
            }
            nifti_image_free(MRFImage);
        }
        cout<<"... taken care of"<<endl;
//        cout<<"Quantilization..."<<endl;
//        if (segment_param!=NULL) {
//            cout<<segment_param->quantMin<<" "<<segment_param->quantMax;
//        }
//        TreeResult->QuantilizeDataImage(segment_param);
//        cout<<"...taken care of";
//        if (!TreeResult->IsDataImageNormalised()) {
//            cout<< "We have to normalise Image"<<endl;
//            TreeResult->NormaliseDataImage();
//        }
        cout<< "Image is already normalised" << endl;
        
        // WARNING !!! temporary only !!!
//        TreeResult->NormalisePriors();
//        TreeResult->UpdateNonNormResp(segment_param);
//        TreeResult->UpdateNormResp();
        
        // Clearing MRF if not to be used after initialisation
        if(!segment_param->flag_MRF){
            vector<TreeEM*> LeavesVector=TreeResult->GetAllLeaves();
            int numbLeaves=LeavesVector.size();
            for(int l=0;l<numbLeaves;l++){
                if(LeavesVector[l]->GetMRF()!=NULL){
                    LeavesVector[l]->SetMRF(NULL);
                }
            }
        }
        
        if(segment_param->flag_DP==1){
            
            int sizeDP=segment_param->filename_DP.size();
            vector<int *> CountChildrenInput;
            if(sizeDP==TreeResult->GetNumberChildren()){
                TreeResult->SetFlagDistClassInd(1);
                CountChildrenInput=ReadCountfromFiles(segment_param->filename_DP);
            }
            if(sizeDP==1){
                
                CountChildrenInput=ReadCountFromDistributionTextFile(segment_param->filename_DP[0]);
            }
            if(!segment_param->flag_DistClassInd){
                TreeResult->SetFlagDistClassInd(0);
                vector<int*> DistributionChildrenInput=ReadDistributionFromTextFile(segment_param->filename_DP[0]);
                float * CountHistogram= MultivariateJointHistogram(CountChildrenInput,DistributionChildrenInput);
                vector<int> DimToUse;
                int numbClasses=CountChildrenInput.size();
                for(int c=0;c<numbClasses;c++){
                    DimToUse.push_back(MaxSupport);
                }
                //                float gauss_std=1;
                float * DPChildrenInput=GaussianBlurring(CountHistogram, segment_param->DPGauss_std, DimToUse);
                // Normalisation of the gaussian filtered distribution
                int SizeHistogram=pow(MaxSupport,numbClasses);
                float sumHistBlurred=0;
                for(int i=0;i<SizeHistogram;i++){
                    sumHistBlurred+=DPChildrenInput[i];
                }
                for(int i=0;i<SizeHistogram;i++){
                    DPChildrenInput[i]/=sumHistBlurred;
                }
                delete [] CountHistogram;
                for(int c=0;c<numbClasses;c++){
                    delete[] DistributionChildrenInput[c];
                    delete[] CountChildrenInput[c];
                }
                TreeResult->SetDPChildren(DPChildrenInput);
            }
            else{
                TreeResult->SetFlagDistClassInd(1);
                int Csize=CountChildrenInput.size();
                float * DPChildrenInput=new float[Csize*MaxSupport];
                for(int c=0;c<Csize;c++){
                    float * DPChildrenInput_tmp=PoissonRegressionKernelSmoothing(CountChildrenInput[c],segment_param->smoothing_order,segment_param);
                    for(int i=0;i<MaxSupport;i++){
                        DPChildrenInput[c*MaxSupport+i]=DPChildrenInput_tmp[i];
                    }
                    delete [] DPChildrenInput_tmp;
                    DPChildrenInput_tmp=NULL;
                    //                DPChildrenInput.push_back(PoissonRegressionKernelSmoothing(CountChildrenInput[c],segment_param->smoothing_order,segment_param));
                    //                segment_param->BWfactor=2;
                    //                DPTest.push_back(PoissonRegressionKernelSmoothing(CountChildrenInput[c],segment_param->smoothing_order,segment_param));
                    delete [] CountChildrenInput[c];
                    CountChildrenInput[c]=NULL;
                }
                CountChildrenInput.clear();
                TreeResult->SetDPChildren(DPChildrenInput);
            }
        }
        if (segment_param!=NULL) {
            delete segment_param;
            segment_param=NULL;
        }
//        float * DataTest=static_cast<float *>(TreeResult->GetDataImage()->data);
//        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+TreeResult->GetNumberElements()]<<endl;
        return TreeResult;
    }
    
}

TreeEM * ReadChildLevel(istream & text,TreeEM * TreeResult, SEG_PARAMETERS * segment_param){
    TreeEM * ChildToAdd=new TreeEM();
    int numbmodal=TreeResult->GetNumberModalities();
    string line;
    getline(text,line);
//    cout<<line<<endl;
    istringstream in(line);
    string first;
    in >> first;
    if(first=="Class"){ // Find parent in the tree hierarchy
        vector<int> Hierarchy;
        int num;
        while(in >> num){
            Hierarchy.push_back(num);
        }
        int level=Hierarchy.size();
        TreeEM * ParentOfAddedChild=TreeResult;
        if(level>1){
            for(int l=0;l<level-1;l++){
                ParentOfAddedChild=ParentOfAddedChild->GetChild(Hierarchy[l]);
            }
        }
        ChildToAdd->SetParent(ParentOfAddedChild);
    }
    
//    nifti_image* PriorImage=NULL;
    int Distribution=0;
    float weight=0;
    string PriorFilename;
    
    getline(text,line);
    
    while(!line.empty()){
        istringstream in2(line);
//        cout<<line<<endl;
        string type;
        if(!line.empty()){
            in2 >>type;
        }
        else{
            type.clear();
        }
        if(type == "Prior"){ // WARNING Temporarily not set
            in2>> PriorFilename;
            //PriorImage=ReadFromFilename(const_cast<char*>(PriorFilename.c_str()));
//            PriorImage=ReadFromFilename(PriorFilename);
////            ChildToAdd->SetPriors(PriorImage);
//            //            segment_param->filename_priors.push_back(const_cast<char*>(PriorFilename.c_str()));
//            segment_param->filename_priors.push_back(PriorFilename);
//            //            segment_param->flag_manual_priors++;
//            segment_param->flag_manual_priors=1;
        }
//        else if(type == "AdaptedPriors"){ // WARNING Temporarily not set
//            string AdaptedPriorFilename;
//            in2>> AdaptedPriorFilename;
//            nifti_image * PriorAdapted=ReadFromFilename(AdaptedPriorFilename);
//            ChildToAdd->SetPriorsAdapted(static_cast<float*>(PriorAdapted->data));
////            nifti_image_free(PriorAdapted);
////            ChildToAdd->SetPriors(PriorAdapted);
////            ChildToAdd->SetPriorsAdapted(static_cast<float*>(PriorAdapted->data));
//            segment_param->filename_PriorsAdaptedOut.push_back(AdaptedPriorFilename);
//            cout<<segment_param->filename_PriorsAdaptedOut[0]<<endl;
//        }
        else if(type == "Weight"){
            in2 >> weight;
            ChildToAdd->SetNormWeight(weight);
        }
        else if (type == "UniformDist"){
            Parameters * ParametersFinal=new Parameters();
            ParametersFinal->DistributionType=2;
            ChildToAdd->SetParameters(ParametersFinal);
            delete ParametersFinal;
            ParametersFinal=NULL;
        }
        else if(type == "Leaf"){
            in2>>Distribution;
            Parameters * ParametersFinal=new Parameters();
            switch(Distribution){
                default : {
                    
                    ParametersFinal->DistributionType=Distribution;
                    
                    ParametersFinal->SizeParameters=TreeResult->CalculateSizeParameters(Distribution);
                    int sizeParams=ParametersFinal->SizeParameters;
                    ParametersFinal->ValueParameters=new float[numbmodal+numbmodal*numbmodal];
                    float * ParametersValue = new float[numbmodal+numbmodal*numbmodal];
                    for (int m=0; m<sizeParams; m++) {
                        ParametersValue[m]=0;
                        ParametersFinal->ValueParameters[m]=0;
                    }
                    
                    for(int i=0;i<2;i++){
                        getline(text,line);
                        istringstream in3(line);
                        in3>>type;
                        if(type=="Mean"){
                            for(int m=0;m<numbmodal;m++){
                                in3>>ParametersValue[m];
                            }
                        }
                        else if(type == "Variance"){
                            for(int m=0;m<numbmodal*numbmodal;m++){
                                in3>>ParametersValue[numbmodal+m];
                            }
                            
                        }
                    }
                    
                    for (int m=0; m<sizeParams; m++) {
                        ParametersFinal->ValueParameters[m]=ParametersValue[m];
                    }
                    delete [] ParametersValue;
                    ParametersValue=NULL;
                }
            }
            ChildToAdd->SetParameters(ParametersFinal);
            delete ParametersFinal;
            ParametersFinal=NULL;
        }
        else if(type == "Mixture"){
            Parameters * ParametersFinal=new Parameters();
            ChildToAdd->SetParameters(ParametersFinal);
            delete ParametersFinal;
            ParametersFinal=NULL;
        }
        getline(text,line);
    }
    return ChildToAdd;
}

vector<int> ReadReclassifLesionFromFile(string DecisionTextFile){
    vector<int> DecisionVector;
    cout<<"Name of text file is "<<DecisionTextFile<<endl;
    ifstream text (DecisionTextFile.c_str());
    //    text.open(TreeTextFile,ios::in);
    if(!text){
        std::cout<<"could not open the text file properly !"<<endl;
        return DecisionVector;
    }
    else{
        std::string line;
    while(getline(text,line)){
        istringstream in(line);
        std:: string initLabel;
        int label;
        int decision=0;
        in >> initLabel;
        in >> label;
        in >> decision;
        DecisionVector.push_back(decision);
    }
    }
    return DecisionVector;
}

// From Text file of reclassification decision about lesions return a vector of boolean asserting if the labels should be reclassified or not
map<int,int> ReadReclassifDecisionFromFile(string DecisionTextFile){

    map<int,int> DecisionVector;
    cout<<"Name of text file is "<<DecisionTextFile<<endl;
    ifstream text (DecisionTextFile.c_str());
    //    text.open(TreeTextFile,ios::in);
    if(!text){
        std::cout<<"could not open the text file properly !"<<endl;
        return DecisionVector;
    }
    else{
        std::string line;
    while(getline(text,line)){
        istringstream in(line);
        std:: string initLabel;
        int label;
        int decision=0;
//        in >> initLabel;
        in >> label;
//        DecisionVector.push_back(label);
        in >> decision;
        cout << label << " " << decision << endl;
//        DecisionVector.push_back(decision);
        DecisionVector.insert( pair<int,int>(label,decision) );
    }
    }
    return DecisionVector;
    
}

// From 2 4D images with presumably the same 4D dimension, smoothing of each image in 4D and averaging given the float weight applied to S1 (1-w applied to S2) No normalisation at this stage. We assume compatible dimensions
vector<nifti_image *> SmoothAndAverageImagesForPriors(nifti_image * Summarised1, nifti_image * Summarised2, float Weight, SEG_PARAMETERS * segment_param,int Level=0){
    vector<nifti_image *> VectorSAImages;
    vector<nifti_image*> VectorTempSImages1;
    vector<nifti_image*> VectorTempSImages2;
    if (Summarised1==NULL && Summarised2==NULL) {
        return VectorSAImages;
    }
    if (Summarised2!=NULL) {
        VectorTempSImages2=Demerge4DImage(Summarised2);
        cout<<"Image 2 demerged"<<endl;
        int numbimages=VectorTempSImages2.size();
        for (int n=0; n<numbimages; n++) {
            float * ToSmooth=static_cast<float *>(VectorTempSImages2[n]->data);
            vector<int> DimVector;
            int numel=VectorTempSImages2[n]->nvox;
            for (int d=0; d<3; d++) {
                DimVector.push_back(VectorTempSImages2[n]->dim[d+1]);
            }
            float Smoothing=segment_param->AtlasSmoothing[0];
            int sizeAtlasAveraging=segment_param->AtlasAveraging.size();
            if (segment_param->AtlasAveraging.size()>0){
                if (Level<sizeAtlasAveraging & Level>0){
                    Smoothing=segment_param->AtlasAveraging[Level];
                    cout << "Smoothing is "<<Smoothing << "for level "<<Level<<endl;
                }
                else{
                    Smoothing=segment_param->AtlasAveraging[sizeAtlasAveraging-1];
                }
            }
            float * SmoothedTemp=GaussianBlurring(ToSmooth, Smoothing, DimVector,0);
            float MaxValue=-1;
            float MaxPre=-1;
            for (int i=0; i<numel ; i++) {
                if (ToSmooth[i]>MaxPre) {
                    MaxPre=ToSmooth[i];
                }
                ToSmooth[i]=SmoothedTemp[i];
                if (SmoothedTemp[i]>MaxValue) {
                    MaxValue=SmoothedTemp[i];
                }
            }
//            cout<<"max for smoothed 2 is "<<MaxValue<<" and for init is "<<MaxPre<<endl;
            delete [] SmoothedTemp;
            SmoothedTemp=NULL;
        }
//        nifti_set_filenames(VectorTempSImages2[1], g_scratch_dir + "/SmoothingInit.nii.gz", 0, 0);
//        nifti_image_write(VectorTempSImages2[1]);
    }
    if(Summarised1!=NULL){
       VectorTempSImages1=Demerge4DImage(Summarised1);
        int numbimages=VectorTempSImages1.size();
        for (int n=0; n<numbimages; n++) {
            float * ToSmooth=static_cast<float *>(VectorTempSImages1[n]->data);
            vector<int> DimVector;
            int numel=VectorTempSImages1[n]->nvox;
            for (int d=0; d<3; d++) {
                DimVector.push_back(VectorTempSImages1[n]->dim[d+1]);
            }
            float * SmoothedTemp=GaussianBlurring(ToSmooth, segment_param->AtlasSmoothing[0], DimVector,0);
            float MaxValue=-1;
            float MaxPre=-1;
            for (int i=0; i<numel ; i++) {
                if (ToSmooth[i]>MaxPre) {
                    MaxPre=ToSmooth[i];
                }
                ToSmooth[i]=SmoothedTemp[i];
                if (SmoothedTemp[i]>MaxValue) {
                    MaxValue=SmoothedTemp[i];
                }
            }
//             cout<<"max for smoothed 1 is "<<MaxValue<<" and for init is "<<MaxPre<<endl;
            delete [] SmoothedTemp;
            SmoothedTemp=NULL;
        }
//        nifti_set_filenames(VectorTempSImages1[1], g_scratch_dir + "/SmoothingInit.nii.gz", 0, 0);
//        nifti_image_write(VectorTempSImages1[1]);
    }
//    Doing the linear combination while possible and taking only one if only one available
    int numb1=VectorTempSImages1.size();
    int numb2=VectorTempSImages2.size();
    for (int n=numb2; n<numb1; n++) {
        VectorTempSImages2.push_back(NULL);
    }
    for (int n=numb1; n<numb2; n++) {
        VectorTempSImages1.push_back(NULL);
    }
    
    int numbfin=VectorTempSImages1.size();
    vector<float> VectorCombWeight;
    VectorCombWeight.push_back(Weight);
    VectorCombWeight.push_back(1-Weight);
    for (int n=0; n<numbfin; n++) {
        vector<nifti_image *> VectorCombImages;
        VectorCombImages.push_back(VectorTempSImages1[n]);
        VectorCombImages.push_back(VectorTempSImages2[n]);
        nifti_image* LCImage=LinearCombWeighted(VectorCombImages,VectorCombWeight);
        VectorSAImages.push_back(LCImage);
    }
    
//    Free the temporary images in both VectorTemp
    for (int n=0; n<numbfin; n++) {
        if (VectorTempSImages2[n]!=NULL) {
            nifti_image_free(VectorTempSImages2[n]);
            VectorTempSImages2[n]=NULL;
        }
        if (VectorTempSImages1[n]!=NULL) {
            nifti_image_free(VectorTempSImages1[n]);
            VectorTempSImages1[n]=NULL;
        }
    }
    return VectorSAImages;
}

nifti_image * LinearCombWeighted(vector<nifti_image *> VectorImagesComb, vector<float> VectorWeightComb){
    int fnn=FirstNotNULL(VectorImagesComb);
    int sizeVIC=VectorImagesComb.size();
    if (fnn==sizeVIC) {
        return NULL;
    }
    int numbI=VectorImagesComb.size();
    int numbW=VectorWeightComb.size();
    int numbfin=numbI<=numbW?numbI:numbW;
    vector<bool> ImagesNULL;
    vector<bool> WeightZero;
    float sumWeights=0;
    for (int n=0; n<numbfin; n++) {
        bool WeightFlag=VectorWeightComb[n]!=0?1:0;
        sumWeights+=VectorWeightComb[n];
        bool ImageFlag=VectorImagesComb[n]!=NULL?1:0;
        ImagesNULL.push_back(ImageFlag);
        WeightZero.push_back(WeightFlag);
    }
    if (sumWeights==0) {
        return NULL;
    }
    else{
        for (int n=0; n<numbfin; n++) {
            VectorWeightComb[n]/=sumWeights;
        }
    }
    nifti_image * ResultComb=nifti_copy_nim_info(VectorImagesComb[fnn]);
    ResultComb->data=(void *)calloc(ResultComb->nvox, sizeof(float));
    int numel=ResultComb->nvox;
    float * ResultCombData=static_cast<float*>(ResultComb->data);
    for (int i=0; i<numel; i++) {
        ResultCombData[i]=0;
    }
    for (int n=0; n<numbfin; n++) {
        if (ImagesNULL[n]&&WeightZero[n]) {
            float * ImageToCombData=static_cast<float *>(VectorImagesComb[n]->data);
            for (int i=0; i<numel; i++) {
                ResultCombData[i]+=VectorWeightComb[n]*ImageToCombData[i];
            }
        }
    }
    return ResultComb;
}

// Copy in a vector of nii images the images demerged from the 4D stacked one
vector<nifti_image *> Demerge4DImage(nifti_image * Image4D){
    vector<nifti_image*> VectorImages;
    if (Image4D==NULL) {
        return VectorImages;
    }
    int numbimages=Image4D->nu*Image4D->nt;
    int numel=Image4D->nx*Image4D->ny*Image4D->nz;
    float * Image4DData=static_cast<float *>(Image4D->data);
    for (int v=0; v<numbimages; v++) {
//        Preparing image vector
        nifti_image * ImageVector=nifti_copy_nim_info(Image4D);
        ImageVector->dim[0]=3;
        ImageVector->dim[4]=1;
        ImageVector->dim[5]=1;
        nifti_update_dims_from_array(ImageVector);
        ImageVector->data=(void *)calloc(ImageVector->nvox, sizeof(float));
        float * ImageData=static_cast<float *>(ImageVector->data);
//        Filling image and pushing it to vector
        for (int i=0; i<numel; i++) {
            ImageData[i]=Image4DData[v*numel+i];
        }
        VectorImages.push_back(ImageVector);
        
    }
    return VectorImages;
}

nifti_image * CreateNewPrior(vector<nifti_image*> DoubledVectorPriors,vector<int> ListToCombine){
    //Checking that there are priors to combined
    int numbPriors=DoubledVectorPriors.size();
    if (numbPriors==0) {
        cout<<"No priors to combine"<<endl;
        return NULL;
    }
    //Checking if ListToCombine is non void
    int numbToCombine=ListToCombine.size();
    if (numbToCombine==0) {
        cout<<"Nothing in the combining list"<<endl;
        return NULL;
    }
    // First checking that the combination offered in ListToCombine is compatible with the list of nifti_images
    for(int i=0;i<numbToCombine;i++){
        if (ListToCombine[i]>=numbPriors){
            cout<<"Impossible combination of priors"<<endl;
            return NULL;
        }
        for(int j=0;j<numbToCombine;j++){
            if (ListToCombine[i]==ListToCombine[j] && i!=j) {
                cout<<"Impossible to combine twice same priors"<<endl;
                return NULL;
            }
        }
    }
    // We assume that the dimensions between the priors vectors are compatible
    nifti_image * CombinedPriors=nifti_copy_nim_info(DoubledVectorPriors[0]);
    CombinedPriors->data=(void *)calloc(CombinedPriors->nvox, sizeof(float));
    vector<float*> DataCombiningPriors; // Constructing vector containing all priors to combine
    for (int i=0; i<numbToCombine; i++) {
        float * Priors_PTR=static_cast<float*>(DoubledVectorPriors[ListToCombine[i]]->data);
        DataCombiningPriors.push_back(Priors_PTR);
    }
    // Doing the proper addition
    float * CombinedPriors_PTR=static_cast<float*>(CombinedPriors->data);
    int numel = CombinedPriors->nvox;
    for (int i=0; i<numel; i++) {
        CombinedPriors_PTR[i]=0;
    }
    for (int i=0; i<numel; i++, CombinedPriors_PTR++) {
        for (int c=0; c<numbToCombine; c++) {
            *CombinedPriors_PTR+=*DataCombiningPriors[c];
            DataCombiningPriors[c]++;
        }
    }
    return CombinedPriors;
}


vector<string> ReadListFilenamesFromTextFile(string TextFile){
    vector<string> FilenameVector;
    ifstream text (TextFile.c_str());
    if(!text){
        std::cout<<"could not open the text file properly for the params adaptation!"<<endl;
    }
    else{
        std::string line;
        while(getline(text,line)){
            istringstream in(line);
            std:: string filename;
            in >> filename;
            if(strlen(filename.c_str())>0)
            FilenameVector.push_back(filename);
        }
    }
    return FilenameVector;

}

/*************************************** FUNCTIONS NEEDED WHEN USING PROGRESSIVE MODEL COMPLEXIFICATION WITH ADDED MODALITIES ****************/


nifti_image * CreateUncorrectedBFImage(nifti_image * DataCorrected, float * BFCoeffs, int BForder,nifti_image * Mask){
    nifti_image * CorrectedFirst=NULL;
    if (DataCorrected!=NULL) {
        CorrectedFirst=CopyFloatNiiImage(DataCorrected);
    }
    return CorrectedFirst;
}

nifti_image * CopyFloatNiiImage(nifti_image * FloatNiiToCopy){
    if (FloatNiiToCopy==NULL) {
        return NULL;
    }
    nifti_image * CopiedImage=nifti_copy_nim_info(FloatNiiToCopy);
    CopiedImage->data=(void *)calloc(CopiedImage->nvox, sizeof(float));
    float * CopiedImageData=static_cast<float *>(CopiedImage->data);
    float * FloatNiiToCopyData=static_cast<float *>(FloatNiiToCopy->data);
    int numel=CopiedImage->nvox;
    for (int i=0; i<numel; i++) {
        CopiedImageData[i]=FloatNiiToCopyData[i];
    }
    return CopiedImage;
    
}

#endif


