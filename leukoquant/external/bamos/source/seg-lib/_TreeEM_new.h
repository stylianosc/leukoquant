#pragma once

//#ifndef _Seg_InitAndRun_h
//#define _Seg_InitAndRun_h
//#endif

//
//  TreeEM.h
//  TreeSeg
//
//  Created by Carole Sudre on 08/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 07/05/2026.
//

//#ifndef TreeSeg_TreeEM_h
//#define TreeSeg_TreeEM_h

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <limits>

#include "nifti1.h"
#include "nifti1_io.h"
#include "_seg_matrix.h"
#include "_seg_tools.h"
//#include "_seg_EM.h"
#include "_SVDCov.h"
#include "_EigenvaluesCov.h"



using namespace std;

class TreeEM;

#define PrecisionTYPE double
#define MaxNumbLeaves 30
#define MaxSupport 15
#define MaxNumbModal 5
#define MaxBForder 5
#define WeightThreshold 10E-6
#define MaxNumbBF ((MaxBForder+1)*(MaxBForder+2)*(MaxBForder+3))/6


typedef struct SEG_PARAMETERS{
    int  bias_order;
    PrecisionTYPE Bias_threshold;
    PrecisionTYPE ConvThreshold;
    PrecisionTYPE DPGauss_std;
    vector <float> AtlasWeight;
    vector <float> AtlasSmoothing;
    vector <float> AtlasAveraging;
    PrecisionTYPE VLkappa; // Used only if VL implementation in place
    int KernelSize;
    float quantMin;
    float quantMax;
    float OutliersWeight;
    float UnifSplitWeight;
    float WeightMini;
    float AcceptanceThreshold;
    float DistInitUnif;
    float Mahal;
    float weightIOAverage;
    float weightGCAverage;
    float MeanPriors;
    int OutliersMod;
    int  maxIteration;
    int maxIterationBF;
    int minIteration;
    int maxNumbLeaves;
    int numbDP;
    int numbCount;
    int choiceInitSplit;
    int uniformTypeChange;
    int VarianceInitUnif;
    int class_keptpriorsmax;
//    int  verbose_level;
    int numb_classes;
    int numb_classes_out;
    int numbmodal;
    int smoothing_order;
    PrecisionTYPE BWfactor;
    bool flag_Outliers;
    bool flag_manual_priors;
    bool flag_manual_priors_out;
    bool flag_manual_priors_bnb;
    bool flag_OutBrain;
    bool flag_input;
    bool flag_mask;
    bool flag_NormMask;
    bool flag_out;
    bool flag_Bias;
    bool flag_bc_out;
    bool flag_bf_out;
    bool flag_data_out;
    bool flag_EMfirst_out;
    bool flag_BiASM;
    bool flag_DP;
    bool flag_DistClassInd;
    bool flag_Count;
    bool flag_Countmod;
    bool flag_intxt;
    bool flag_CEM;
    bool flag_TypicalityAtlas;
    bool flag_CommonChanges;
    bool flag_MRF;
    bool flag_MRFPost;
    bool flag_MRFIn;
    bool flag_GMatrix;
    bool flag_GMatrixIn;
    bool flag_GMatrixPost;
    bool flag_PriorsAdaptedOut;
    bool flag_MRFOut;
    bool flag_inDC; // Input of already data corrected
    bool flag_inDCFile;
    bool flag_DeleteUnderWeight;
    bool flag_CovPriors;
    int flag_OutlierAtlas;
    bool flag_optMRFOut;
    bool flag_bnbComb;
    bool flag_DGMPrior;
    bool flag_JuxtaCorrection;
    bool flag_RefinedSOrdering;
    bool flag_KMeansModif;
    bool flag_meanPriors;
    bool flag_progressiveBFC; // Stating if BF correction is progressive in term of number of polynomial used
    bool flag_unifTot; // if using kmeans initialisation for the splitting of uniforms assess if we will test also the part with bigger variance in later stages of split trials
    bool flag_savePriors; // Determines if priors will be saved as output.
    bool flag_IOAverage; // If we use two preliminary IO segmentation for atlases
    bool flag_GCAverage; // If we use two preliminary segmentations for atlases;
    bool flag_AdaptTransform; // If there is an histogram matching to apply on parameters when building again
    bool flag_BoostAtlas; // If we want to boost the outliers atlas after the initial segmentation
    bool flag_FurtherAddedModa;
    int CovTest;
    int CovPriorsType;
    int CovPriorsSplit;
    int CovPriorsMerge;
    int BICFP;
    int AcceptanceType;
    int SMOrder;
    int PriorsKept;
    int ProgressivePriors;
    int choiceInitSplitUnifKMeans;
    int SplitAccept;
    int MaxRunEM;
    int IndexCSF;
    int IndexWM;
    int IndexGM;
    int IndexOut;
    vector<int> vecBP;
    vector<int> Modalities;
    vector<int> NumberAddedModalities;
    vector<vector<int> > OutliersCombined;
    vector<string> filename_input;
    vector<string> filename_out;
    vector<string> filename_DP;
    vector<string> filename_Count;
    string filename_datacorrected;
    string filename_intxt;
    string filename_correction;
    string filename_mask;
    string filename_datatxt;
    string filename_EMfirst;
    string filename_GMatrixPost;
    string filename_GMatrix;
    string filename_DGMPrior;
    string filename_inDC;
    vector<string> filename_priors;
    vector<string> filename_priors_out;
    vector<string> filename_priors_bnb;
    vector<string> filename_PriorsAdaptedOut;
    vector<string> filename_IOAverage;
    vector<string> filename_GCAverage;
    string filename_MRFOut;
    vector<string> filename_AdaptTransform;

//    SEG_PARAMETERS(){
//        filename_datacorrected=NULL;
//        filename_correction=NULL;
//        filename_mask=NULL;
//    }

//    ~SEG_PARAMETERS(){
//        if(filename_datacorrected!=NULL){
//            delete[] filename_datacorrected;
//            filename_datacorrected=NULL;
//        }
//        if(filename_correction!=NULL){
//            delete[] filename_correction;
//            filename_correction=NULL;
//        }
//        if(filename_mask!=NULL){
//            delete[] filename_mask;
//            filename_mask=NULL;
//        }
//    }

  }SEG_PARAMETERS;



typedef struct Parameters{
    int DistributionType;
    int SizeParameters;
    bool PriorsCovFlag;
    float * PriorsCov;
    float * PriorsMean;
    float * ValueParameters;
    Parameters(){
        DistributionType=0;
        SizeParameters=0;
        ValueParameters=NULL;
        PriorsCovFlag=0;
        PriorsCov=NULL;
        PriorsMean=NULL;
    }
~Parameters(){
    if(ValueParameters!=NULL){
        delete [] ValueParameters;
        ValueParameters=NULL;
    }
    if(PriorsCov!=NULL){
        delete [] PriorsCov;
        PriorsCov=NULL;
    }
    if(PriorsMean!=NULL){
        delete [] PriorsMean;
        PriorsMean=NULL;
    }

}
} Parameters;

typedef struct GenKLD{
    float KLDTot;
//    float * KLD;
    GenKLD(){
        KLDTot=-1;
//        KLD=NULL;
    }
    ~GenKLD(){
        KLDTot=-1;
//        if(KLD!=NULL){
//            delete[] KLD;
//            KLD=NULL;
//        }
    }

    GenKLD* CopyGenKLD(){
        GenKLD * CopiedGenKLD=new GenKLD();
        CopiedGenKLD->KLDTot=this->KLDTot;
//        CopiedGenKLD->KLD=new float[MaxNumbModal];
//        if(this->KLD!=NULL){
//            for(int m=0;m<MaxNumbModal;m++){
//                CopiedGenKLD->KLD[m]=this->KLD[m];
//            }
//        }
//        else{
//            CopiedGenKLD->KLD=NULL;
//        }
        return CopiedGenKLD;
    }
}GenKLD;

typedef struct MergeKLD{
    TreeEM * Parent;
    int Child1;
    int Child2;
    float KLDTot;
//    float * KLD;

    MergeKLD(){
        Parent=NULL;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
//        KLD=NULL;
//        KLD=new float[MaxNumbModal];
//        for(int m=0;m<MaxNumbModal;m++){
//            KLD[m]=0;
//        }
    }

    ~MergeKLD(){
        Parent=NULL;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
//        if(KLD!=NULL){
//            delete [] KLD;
//            KLD=NULL;
//        }
    }

    MergeKLD * CopyMergeKLD(){
        MergeKLD * CopiedMergeKLD=new MergeKLD();
        CopiedMergeKLD->KLDTot=this->KLDTot;
        CopiedMergeKLD->Parent=this->Parent;
        CopiedMergeKLD->Child1=this->Child1;
        CopiedMergeKLD->Child2=this->Child2;
//        int numbmodal=CopiedSplitKLD->Parent->GetNumberModalities();
//        CopiedMergeKLD->KLD=new float[MaxNumbModal];
//        if(this->KLD!=NULL){
//        for(int m=0;m<MaxNumbModal;m++){
//            CopiedMergeKLD->KLD[m]=this->KLD[m];
//        }
//        }
//        else{
//            CopiedMergeKLD->KLD=NULL;
//        }
        return CopiedMergeKLD;
    }
}MergeKLD;

typedef struct MergeKLD_b{
    vector<int> HierarchyParent;
    int Child1;
    int Child2;
    float KLDTot;
    //    float * KLD;

    MergeKLD_b(){
        vector<int> HierarchyParent;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
        //        KLD=NULL;
        //        KLD=new float[MaxNumbModal];
        //        for(int m=0;m<MaxNumbModal;m++){
        //            KLD[m]=0;
        //        }
    }

    ~MergeKLD_b(){
        HierarchyParent.clear();
//        HierarchyParent.shrink_to_fit();
//        Child1=-1;
//        Child2=-1;
//        KLDTot=-1;
        //        if(KLD!=NULL){
        //            delete [] KLD;
        //            KLD=NULL;
        //        }
    }

    MergeKLD_b * CopyMergeKLD_b(){
        MergeKLD_b * CopiedMergeKLD=new MergeKLD_b();
        CopiedMergeKLD->KLDTot=this->KLDTot;
        int sizeHierarch=this->HierarchyParent.size();
        for(int i=0;i<sizeHierarch;i++){
            CopiedMergeKLD->HierarchyParent.push_back(this->HierarchyParent[i]);
        }
        CopiedMergeKLD->Child1=this->Child1;
        CopiedMergeKLD->Child2=this->Child2;
        //        int numbmodal=CopiedSplitKLD->Parent->GetNumberModalities();
        //        CopiedMergeKLD->KLD=new float[MaxNumbModal];
        //        if(this->KLD!=NULL){
        //        for(int m=0;m<MaxNumbModal;m++){
        //            CopiedMergeKLD->KLD[m]=this->KLD[m];
        //        }
        //        }
        //        else{
        //            CopiedMergeKLD->KLD=NULL;
        //        }
        return CopiedMergeKLD;
    }
}MergeKLD_b;

typedef struct SplitKLD{
    TreeEM * Parent;
    int ChildToSplit;
    int TypeChange; // 0 for change to Gaussian, 1 for classical split, 2 for change from 1 uniform to 2 gaussians
    int InitSplitUnif;
    float KLDTot;
    float * KLD;

    SplitKLD(){
        Parent=NULL;
        ChildToSplit=-1;
        KLDTot=-1;
        TypeChange=1;
        InitSplitUnif=0;
//        KLD=NULL;
//        KLD=new float[MaxNumbModal];
//        for(int m=0;m<MaxNumbModal;m++){
//            KLD[m]=0;
//        }
    }

    ~SplitKLD(){
        Parent=NULL;
//        if(KLD!=NULL){
//            delete [] KLD;
//            KLD=NULL;
//        }
    }

    SplitKLD * CopySplitKLD(){
        SplitKLD * CopiedSplitKLD=new SplitKLD();
        CopiedSplitKLD->KLDTot=this->KLDTot;
        CopiedSplitKLD->Parent=this->Parent;
        CopiedSplitKLD->ChildToSplit=this->ChildToSplit;
        CopiedSplitKLD->TypeChange=this->TypeChange;
        CopiedSplitKLD->InitSplitUnif=this->InitSplitUnif;
//        int numbmodal=CopiedSplitKLD->Parent->GetNumberModalities();
//        CopiedSplitKLD->KLD=new float[MaxNumbModal];
//        if(this->KLD!=NULL){
//        for(int m=0;m<MaxNumbModal;m++){
//            CopiedSplitKLD->KLD[m]=this->KLD[m];
//        }
//        }
//        else{
//            CopiedSplitKLD->KLD=NULL;
//        }
        return CopiedSplitKLD;
    }

} SplitKLD;

typedef struct SplitKLD_b{
    vector<int> HierarchyParent;
    int ChildToSplit;
    int TypeChange; // 0 for change to Gaussian, 1 for classical split, 2 for change from 1 uniform to 2 gaussians
    int InitSplitUnif; // Type of initialisation strategy when splitting uniform 0 for normal with kmeans and minimal variance, 1 for kmeans with maximal variance, ....
    float KLDTot;

    SplitKLD_b(){
        vector<int> HierarchyParent;
        ChildToSplit=-1;
        KLDTot=-1;
        TypeChange=1;
        InitSplitUnif=0;
        //        KLD=NULL;
        //        KLD=new float[MaxNumbModal];
        //        for(int m=0;m<MaxNumbModal;m++){
        //            KLD[m]=0;
        //        }
    }

    ~SplitKLD_b(){
        HierarchyParent.clear();
        HierarchyParent.swap(HierarchyParent);

//        HierarchyParent.shrink_to_fit();
        //        if(KLD!=NULL){
        //            delete [] KLD;
        //            KLD=NULL;
        //        }
    }

    SplitKLD_b * CopySplitKLD_b(){
        SplitKLD_b * CopiedSplitKLD=new SplitKLD_b();
        CopiedSplitKLD->KLDTot=this->KLDTot;
        int sizeHierarch=this->HierarchyParent.size();
        for (int i=0;i<sizeHierarch;i++){
            CopiedSplitKLD->HierarchyParent.push_back(this->HierarchyParent[i]);
        }
        CopiedSplitKLD->ChildToSplit=this->ChildToSplit;
        CopiedSplitKLD->TypeChange=this->TypeChange;
        CopiedSplitKLD->InitSplitUnif=this->InitSplitUnif;
        //        int numbmodal=CopiedSplitKLD->Parent->GetNumberModalities();
        //        CopiedSplitKLD->KLD=new float[MaxNumbModal];
        //        if(this->KLD!=NULL){
        //        for(int m=0;m<MaxNumbModal;m++){
        //            CopiedSplitKLD->KLD[m]=this->KLD[m];
        //        }
        //        }
        //        else{
        //            CopiedSplitKLD->KLD=NULL;
        //        }
        return CopiedSplitKLD;
    }

} SplitKLD_b;

class TreeEM{
    TreeEM * Parent; // NULL if Root of the Tree
    vector<TreeEM *> Children; // Children of the tree : NULL if leaf
    float * DPChildren;// DP for number of subclasses per child only stored at root
    bool FlagDistClassInd;
    int FlagUTC;
    int FlagOutliers; // Only need to be stored at root.
    int FlagCovPriors; // Only need to be stored at root.
    float FlagMeanPriors; // Only need to be stored at root
    nifti_image * DataImage;
    nifti_image * Mask;
    nifti_image * Priors;
    float * PriorsAdapted;
    float * PartPriorsAdapted;
    int * S2L;
    int * L2S;
//    PrecisionTYPE * NonNormResp; // Non normalized responsability
    float * NormResp; // Normalised responsability
//    PrecisionTYPE * Distribution;// Probabilistic distribution
    float NonNormWeight; // Non normalised weight
    float NormWeight; // Normalised weight
    Parameters * ParametersDistribution; // Parameters corresponding to the distribution ( non null pointer for the value only if leaf of the tree and therefore non mixture distribution)
    bool * SplitCheck;
    bool * MergeCheck;
    float * BFCoeffs;
    //PrecisionTYPE * BFBasisFunctions;
//    PrecisionTYPE * BFCorrection;
    float * DataBFCorrected;
    float IndFactor;
    int NumberMaskedElements;
    int * HardSeg;
    float * MRF;
    float * GMatrix;

    public :

    float * LogGaussBlur(float * GaussianBlur,int TotalSize);
    inline float pow_int(const float x,
                                    int exp);
    vector<nifti_image *> ReadFromFilenamesVectorProgressive(vector<string> Filenames, int numbElementsToRead, int begin=0);
    int FirstNotNULL(vector<nifti_image *> ImagesVector);
    int LastNotNULL(vector<nifti_image *> ImagesVector);
    nifti_image * CreateDataImage(vector<nifti_image*> ImagesToSegment);
    void invertMatrix(float * MatrixToInvert,int sizeMatrix);
    float * ProductMatrix(float *MatrixA, float M, int *Size);
    float * DiffMatrix(float * A, float * B, int* Size);
    float * AddMatrix(float * A, float * B, int* Size);
    float NormInfMatrix(float * A, int* Size);
    float * MatrixSquareRoot(float * Matrix, int Size);
    float * ExpMatrix(float * mat, int maxit, int Size);
    float * LogMatrix(float * Matrix, int Size);
    float * WeightedLogMean(vector<float *> MatrixVector, vector<float> WeightVector, int Size);
    float * InvertMatrix(float * MatrixToInvert, int SizeMatrix);
    float determinant(float * data,int size);
    float * ProductMatrix(float * MatrixL, float * MatrixR, int * Size);
    float * TransposeMatrix(float * MatrixToTranspose, int SizeM,int SizeN);
    float Determinant_lib(float **a,int n);
    // Constructors, copy and destructor
    TreeEM();
    TreeEM(nifti_image * Data);

    TreeEM * CopyTree(TreeEM * ParentNew);
    nifti_image * CopyImage();
    nifti_image * CopyMask();
    nifti_image * CopyPriors();
    float * CopyPriorsAdapted();
    float * CopyPartPriorsAdapted();
    int * CopyS2L();
    int * CopyL2S();
//    PrecisionTYPE * CopyDistribution();
//    PrecisionTYPE * CopyNonNormResp();
    float * CopyNormResp();
    Parameters * CopyParameters();
    bool * CopyMergeCheck();
    bool * CopySplitCheck();
    float * CopyBFCorrection();
    float * CopyDataBFCorrected();
    float * CopyBFCoeffs();
    float CopyIndFactor();
    float * CopyDPChildren();
    int * CopyHardSeg();
    ~TreeEM();
    void MakeEmpty();

    // Get methods
    TreeEM * GetParent();
    vector<TreeEM * >GetChildren();
    int GetNumberChildren();
    int GetNumberGeneralClasses();
    TreeEM * GetChild(int c);
    nifti_image * GetDataDirect();
    nifti_image * GetDataImage();
    float GetMaxDataModal(int m);
    float GetMaxDataMaskedModal(int m);
    float GetMinDataModal(int m);
    float GetMinDataMaskedModal(int m);
    float * GetExtremaDataSeg(int ExtremumType, float threshold=0.5);
    float * GetDataModalMasked(int m);
    float * GetDataShifted(int dim,int dir);
    float * GetDataShiftedCorrected(int dim, int dir);
    float * GetDataCorrelation(int dim);
    float GetDataCorrelationMod(int dim, int* modalChoice);
    float MakeIndFactorMod(int dim);
    float MakeIndFactorTotMod();
    float MakeIndFactor(int dim);
    float MakeIndFactorTot();
    int GetNumberElements();
    int GetNumberModalities();
    int GetNumberMaskedElements();
    int MakeNumberMaskedElements();
//    PrecisionTYPE * GetDistribution();
    int GetNumbbins();
nifti_image * GetMask();
    nifti_image * GetMaskDirect();
    nifti_image * GetPriors();
    nifti_image * GetPriorsDirect();
    float * GetPriorsAdaptedDirect();
    float * GetPriorsAdapted();
    float * GetPartPriorsAdapted_bis();
    float * GetPartPriorsAdaptedDirect();
    float * GetPartPriorsAdaptedMasked();
    float * GetPartPriorsAdapted();
int * GetS2L();
    int * GetS2LDirect();
int * GetL2S();
    int * GetL2SDirect();
    int * GetHardSeg();
    int * GetHardSegDirect();
    void MakeHardSeg();
    void MakeL2S();
    void MakeS2L();
//PrecisionTYPE * GetNonNormResp();
float * GetNormResp();
float  GetNonNormWeight();
float  GetNormWeight();
    float GetPartNormWeight();
    float GetPartWeightLevel2();
    float GetPartWeightLevel3();
    float GetPartWeightLevel(int Level);
   float GetPartNormWeightAbovePriors();
    float GetPartNormWeightAboveTopPriors();
    float GetCompleteWeight();
Parameters * GetParameters();
int GetDistributionType();
    bool GetPriorsCovFlag();
    float * GetPriorsCovMatrix();
    float * GetPriorsCovMatrixGeneral();
    float * GetMean();
    float * GetVariance();
    float * GetMeanVarianceGeneralClasses();
    float * GetMeanDirect();
    void GetMeanDirect_bis(float * MeanResult);
    float ** GetMeanDirectLabel(int * Label, int numbLabels);
    float * GetDiagVarianceDirect_corr();
    float * GetDiagVarianceDirect();
    float * GetVarianceDirect();
    float * GetVarianceDirect(float * MeanDirect);
    float ** GetVarianceDirectLabel(int * Label, int numbLabels);
    float * GetParametersValue();
    int GetSizeParameters();
    vector<nifti_image*> GetPriorsVector();
    vector<TreeEM*> GetPriorsNodeVector();
    vector<TreeEM *> GetAllPriorsNodeVectorParent();
    vector<TreeEM *> GetDeepestPriorNodeVector();
    void GetAllPriorsNodeVector(vector<TreeEM *> & PriorNodeVector);
    vector<TreeEM*> GetGeneralClassesVector();
    vector<TreeEM*> GetMainNodesVector();
    TreeEM * GetNodeOutlier();
    TreeEM * GetNodeInlier();
    vector<TreeEM *> GetOutliersMainNodesVector();
    vector<TreeEM *> GetOutliersMainNodesVectorChangeableAtlases(SEG_PARAMETERS * segment_param);
    vector<float *> GetPriorsAdaptedVector();
    vector<float *> GetPartPriorsAdaptedVector();
    vector<float *> GetPartPriorsAdaptedVectorParent();
    vector<float *> GetAllPartPriorsAdaptedVector();
    vector<float *> GetPriorsAdaptedVectorParent();
    vector<float *> GetAllPriorsAdaptedVector();
    void GetLeaves(vector<TreeEM*> & LeavesVector);
    vector<TreeEM *> GetAllLeavesAnatomicalClass(int anatClass);
    vector<TreeEM *> GetAllGaussianLeavesAnatomicalClass(int anatClass);
    vector<TreeEM *> GetAllLeaves();
    vector<TreeEM *> GetAllNodes();
    vector<TreeEM *> GetAllGaussianLeaves();
    vector<TreeEM *> GetAllDirectLeaves();
    int GetNumberDirectLeaves();
    int GetAnatClass();
    TreeEM * GetCorrespondingUniform();
    bool * GetSplitCheck();
    bool * GetMergeCheck();
    void GetNumberLeaves(int & numbleaves);
    int GetNumberAllLeaves();
    int GetNumberFreeParameters();
    float GetLogLikelihood();
    float GetLogLikelihoodCEM(int c);
    float * GetBFCoeffs();
    float * GetBFCoeffsDirect();
    //PrecisionTYPE * GetBFBasisFunctionsDirect();
    //PrecisionTYPE * GetBFBasisFunctions();
    float * GetBFCorrection();
    float * GetBFCorrectionDirect();
    float * GetDataBFCorrected();
    float * GetDataBFCorrectedDirect();
float GetIndFactor();
float * GetDPChildren();
float * GetDPChildrenDirect();
bool GetFlagDistClassInd();
int GetFlagOutliers();
    int GetFlagUTC();
    int GetFlagCovPriors();
    float GetFlagMeanPriors();
    float * GetPriorsCovMatrixTotal();
    float * GetPriorsMeanTotal();
    float * GetPriorsMeanMatrix();
    float * MakeScalingVariance();
    TreeEM * GetMiniVariance();
    enum TreeComponent {
        ROOT = 0,
        BRANCH=1,
        LEAF = 2,
        INITIALNODE=3
    };
    enum PartialResultType {
        DISTRIBUTION=0,
        NONNORMRESP=1,
        NORMRESP=2
    };
    enum OperationType {
        SPLIT=2,
        MERGE=0
    };




    float * GetPartialResult(PartialResultType ResultType,SEG_PARAMETERS * segment_param);

    // Set methods
    void SetParent(TreeEM * Parent);
    void ReinitialiseTree();
    void SetChildren(vector<TreeEM *> Children);
    void AddChild(TreeEM * ChildAdded);
    void ChangeData(nifti_image * DataInput);
    void SetData(nifti_image * data);
    void NormaliseDataImage();
    void NormaliseDataImageNonLog();
    void NormaliseDataImageMasked();
    void QuantilizeDataImage(SEG_PARAMETERS * segment_param);
    void ConvertDataImageToFloat();
//    void SetDistribution(float * Distribution);
    void SetMask(nifti_image * MaskImage);
    void SetNumberMaskedElements(int NumberMaskedElementsInput);
    void MakeMaskBinary();
    void SetPriors(nifti_image * Priors);
    void MakePriorsFloat();
    void MakePriorsProbabilityType();
    void SetPriorsAdapted(float * PriorsAdapted);
    void SetPartPriorsAdapted(float * PartPriorsAdapted);
    void MakePriorsAdaptedProbabilityType();
    void SetS2L(int * S2L);
    void SetL2S(int * L2S);
//    void SetNonNormResp(PrecisionTYPE * NonNormResp);
    void SetNormResp(float * NormResp);
    void InitialiseNormRespWithPriors();
    void InitialiseNormResp();
    void SetNonNormWeight(float NonNormWeight);
    void SetNormWeight(float  NormWeight);
    void SetParameters(Parameters * ParametersToSet);
    void AdaptParametersToAffineTransform(float * LinearWeight, float * ConstWeight);
    void SetDistributionType(int DistributionType);
    void SetBFCoeffsSeparated(float * BFCoeffsInput, vector<int> BFOrderPerModality);
    void SetBFCoeffs(float * BFCoeffs);
    void SetBFCorrection(float * BFCorrection);
    void SetDataBFCorrected(float * DataBFCorrectedInput);
    void MakeParametersMixture();
    void MakeTreeBasic();
    void SetIndFactor(float IndFactorInput);
    void SetDPChildren(float* DPChildrenInput);
    void SetFlagDistClassInd(bool flag_DistClassInd);
    void SetFlagOutliers(int flag_Outliers);
    void SetFlagUTC(int flag_UTC);
    void SetFlagCovPriorsParameters(bool flag_CovPriors);
    void SetFlagCovPriors(int flag_CovPriors);
    void SetFlagMeanPriors(float flag_MeanPriors);
    void SetHardSeg(int * HardSegInput);
    void ModifyCovPriors(SEG_PARAMETERS * segment_param);
    void UpdateCovPriors(SEG_PARAMETERS * segment_param);
    float * CreateCovPriorsFromExistingCovMatrices();
    float * CreateCovPriorsForInverseWishart();
    void SetCovPriorsMatrix(float * NewCovPriorsMatrix);
    void SetMeanPriorsMatrix(float * NewMeanPriorsMatrix);

// EM related methods
    void UpdateParameters();
    void UpdateParametersCEM(int c);
    void UpdateHardSeg();
    void UpdateGaussianParameters();
    void UpdateGaussianParametersCEM(int c);
    void UpdateBFCoeffs();
    void UpdateBFCoeffsSeparated(vector<int> BFOrderperModality, int IndexModalitiesBFTogether);
    void UpdateBFCorrection();
    void UpdateDataBFCorrected();
    void UpdateDataBFCorrectedSeparated(vector<int> BFOrderperModality);
    void UpdateGaussianMean();
    void UpdateGaussianMeanCEM(int c);
    void UpdateGaussianVariance();
    void UpdateGaussianVariancePriors();
    void UpdateGaussianMeanPriors();
    void UpdateGaussianVarianceCEM(int c);
    void UpdatePriorsAdapted(SEG_PARAMETERS * segment_param);
    void UpdatePartPriorsAdapted();
//    void UpdateDistribution();

//    void SumToDistributionWeightedDistributionOfChild(int c);

//    void ReputDistributionToZero();
    void UpdateGaussianDistribution();

    void UpdateNonNormResp(SEG_PARAMETERS * segment_param);
    void UpdateNonNormRespCEM();
    void UpdateNormRespRoot();

//    void ReputNonNormRespToZero();
//    void CalculatingNonNormResp();

//    void SumToNonNormRespWeightedNonNormRespOfChild(int c);

//    void MultiplyCurrentNonNormRespByNormWeight();

    void UpdateNormResp();
    void UpdateTypicality(SEG_PARAMETERS * segment_param);
    void MakeNonNormWeightedSum(float * NonNormSum,SEG_PARAMETERS * segment_param);
    void MakeNonNormWeightedSumCEM(float * NonNormSum);
    void MakeWeightedDist(float * WeigtedDist);
    void MakeWeightedDistPriors(float * WeightedDist);
    float * MakeUniformDistribution();
    float * MakeGaussianDistribution();
//    void NormaliseNonNormResp();

    void DivideNonNormRespByRootNonNormResp();

    void UpdateNonNormWeights();
    void UpdateNonNormWeightsCEM(int c);
    void UpdateNormWeights();
    float * NormalisedPriorsForCCL();
    float CompleteLogLikelihoodFinBis();
    float CompleteLogLikelihoodFin(SEG_PARAMETERS * segment_param);
    float CompleteLogLikelihood();
    void EMCompleteIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, SEG_PARAMETERS * segment_param);
    vector<Parameters *> KeepingCopyOfParameters();
    vector<float> KeepingCopyNormWeight();
    vector<float *> KeepingCopyNormResp();
    vector<float *> KeepingCopyMRF();
    float * KeepingCopyBFCoeffs(vector<int> BForderPerModality);
    void ComeBackInTime(vector<Parameters *> BackupParams, vector<float> BackupNormWeight,SEG_PARAMETERS * segment_param);
    void ComeBackInTime(vector<Parameters *> BackupParams, vector<float> BackupNormWeight, vector<float *> BackupNormResp, vector<float *> BackupMRF,SEG_PARAMETERS * segment_param);
    void CEMPartialIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, int c);
    void CEMCompleteIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param);
    void RunFullEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, SEG_PARAMETERS * segment_param);
    void EMCompleteIterationBFSeparated(vector<int> BFOrderperModality, int IndexModalitiesBFTogether,float & LogLikelihood, float & OldLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param);
    void RunFullEMSeparatedBF(vector<int> BFOrderperModality, int IndexModalityBFTogether, float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param);
    void RunPartialCEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, int c);
    void RunFullCEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param);
    void ReplugParameters(TreeEM * TreeToReplug);
    TreeEM * RunFullBiASM(SEG_PARAMETERS * segment_param);
    TreeEM * RunFullBiASM_bis(SEG_PARAMETERS * segment_param);
    TreeEM * RunFullBiASM_ter(SEG_PARAMETERS * segment_param);
    TreeEM * RunCompleteBaMoS(SEG_PARAMETERS * segment_param);
    TreeEM * RunBaMoSLoop(SEG_PARAMETERS * segment_param);
    //Tree related methods
    void CreateAndAddChild();
    void CreateAndAddChildPriors(SEG_PARAMETERS * segment_param,nifti_image * Priors,int DistributionType);
    void CreateAndAddChildWeight(SEG_PARAMETERS *segment_param, float Weight, int DistributionType);
//int GetLevel();
    void CollapseOnlyChild();
    void ReplaceUnifDist(SEG_PARAMETERS * segment_param);
    void CollapseOnlyChildTot();
    void CollapseChild(TreeEM * ChildToCollapse);
    void AllLeavesToChildrenLevel();
    int Index(TreeEM * ElementToFind);
    void FindAllAtCertainLevel(vector<TreeEM*> & ResultVector, int Level);
    void FindDepth( int & ResDepth);
    TreeEM * FindRoot();
    void DividePartNonNormResp(PrecisionTYPE * SumNonNormResp);
    void InitialiseBasicTree(SEG_PARAMETERS * segment_param);
    void ClearSMChecks();
    void PutAllLeavesToChildrenLevel();
    void ModifyNormWeightsForChildrenLevel();
    void ModifyNormWeightsForMainNodesChildrenLevel();
    void PutAllLeavesToMainNodesChildrenLevel(SEG_PARAMETERS * segment_param);
    void RecPutAllLeavesToChildrenLevel();
    void DeleteUnderWeight(SEG_PARAMETERS * segment_param);
    void DeleteAllNULLChildren();

// Tests on tree
bool IsEqual(TreeEM * TreeEqualityCheck);
bool IsLeaf();
    bool IsBranch();
    bool IsRoot();
void DoesBelongToCurrentTree(TreeEM * Element, bool & test);
int WhatTypeOfTreeComponent();
    bool IsTreeValid(int LevelCheck=0);
    bool IsTreeStructureValid();
    bool IsInitialisationValid();
    bool AreWeightsValid();
    bool IsNormRespValid();
    bool AreNormRespValid();
    bool CheckAcceptanceWeights(vector<SplitKLD *>SplitTry,SEG_PARAMETERS * segment_param);
    bool CheckAcceptanceWeights(vector<SplitKLD_b *>SplitTry,SEG_PARAMETERS * segment_param);
    bool CheckCompatibilityParams();
    bool CheckUnifDistPosition(SEG_PARAMETERS * segment_param);
    bool IsNormRespValidLeaves();
    bool IsNormRespValidGeneral();
    bool IsBasicTree();
    bool IsNumbChildOKWithDistType();
    bool IsDataImageNormalised();
    bool IsDataFloat();
    bool IsDataIsotropic();
    float DistanceFactor(int dim);
    bool IsPriorsFloat();
    bool AreBFCoeffsDirectUseable();
    bool IsStructureSimilar(TreeEM * TreeTest);
    bool IsNewModelAccepted(TreeEM * CopiedTree,OperationType Operation,SEG_PARAMETERS * segment_param);
    // Parameters structure

    int CalculateSizeParameters();
    int CalculateSizeParameters(int DistributionType);
    void CreateAllocateAndInitializeParameters(int DistributionType);
    bool CheckForValidityOfParametersStructure();
    bool CheckForSizeParametersValidity();
    bool CheckForTreeValidity();

    bool ArePriorsNormalised();
    bool ArePriorsAdaptedNormalised();
    bool IsPriorsVectorFilledWithNULL();
    bool IsPriorsAdaptedVectorFilledWithNULL();
    bool IsOneOfThePriorsNULL();
    bool IsOneOfThePriorsAdaptedNULL();
    void NormalisePriors();
    void NormalisePriorsAdapted();
    void NormalisePriorsOverFlown();
    void NormalisePriorsAdaptedOverFlown();
    bool IsPriorsProbabilityType();
    bool IsPriorsAdaptedProbabilityType();
    bool CheckForValidityOfParametersStructure(Parameters * ParametersToCheck);
    bool CheckForSizeParametersValidity( Parameters * ParametersToSet);
    bool IsMRFZero();

    /**************** METHODS FOR SM OPERATIONS *********************/
    // Functions needed for calculation of KLD

    vector<float *> GetDistHistogram();
    float * GetDistHistogramTotal();
    vector<float *> GetUniformDistributionHist();
    vector<float *> GetGaussianDistributionHist();
    float * GetUniformDistributionHistTotal();
    float * GetGaussianDistributionHistTotal();
    vector<float *> GetDataHistogram();
    float * GetDataHistogramTotal();
    float * GetDataHistogramTotalAversion();
    float GetDistCompKLDTotal();
    float * GetDistCompKLD();

    float CriterionCalculation();
    float CriterionCalculationDPInd();
    float CriterionCalculationDPNonInd(int ChildModified);
    float CriterionCalculationSplitDP();
    float CriterionCalculationMergeDP();

    // Functions needed for merging operations
    GenKLD * GetGenKLD();
    MergeKLD * GetKLDMerge(int Child1,int Child2);
    MergeKLD_b * GetKLDMerge_b(int Child1, int Child2);
    vector<MergeKLD *> GetVectorKLDMergeChildren();
    vector<MergeKLD_b *> GetVectorKLDMergeChildren_b();
    bool IsMergeChecked();
    vector<MergeKLD *> GetVectorKLDMergeLeaves();
    vector<MergeKLD_b *> GetVectorKLDMergeLeaves_b();
    MergeKLD * GetToMerge();
    void MergeOperation(int Child1,int Child2,SEG_PARAMETERS * segment_param);
    Parameters * ParametersForMerging(int Child1,int Child2);
    vector<MergeKLD*> GetMergeMoreVertical(int numbclasses);
    vector<MergeKLD_b *> GetMergeMoreVertical_b(int numbclasses);
    vector<int> OrderingMergingChildren();
    vector<MergeKLD*> OrderingMergingLeaves();
    vector<MergeKLD_b *> OrderingMergingLeaves_b();
    vector<MergeKLD*> GetMergeOrderVertical();
    vector<MergeKLD_b *> GetMergeOrderVertical_b();
    TreeEM * RunMergeOperation(MergeKLD * MergeTry, bool & AcceptanceDecision,SEG_PARAMETERS * segment_param);
    TreeEM * RunMergeOperation(vector<MergeKLD_b *> MergeTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param);
    TreeEM * RunMergeOperation(vector<MergeKLD *> MergeTry, bool & AcceptanceDecision,SEG_PARAMETERS * segment_param);
    float * GetKLDModel();
    float * GetSSDModel();
    // Functions needed for splitting operations
    SplitKLD * GetKLDSplit(int ChildToSplit);
    SplitKLD_b * GetKLDSplit_b(int ChildInput);
    vector<SplitKLD *> GetVectorKLDSplitChildren(SEG_PARAMETERS * segment_param);
    vector<SplitKLD_b *> GetVectorKLDSplitChildren_b(SEG_PARAMETERS * segment_param);
    vector<SplitKLD *> GetVectorKLDSplitLeaves(SEG_PARAMETERS * segment_param);
    vector<SplitKLD_b *> GetVectorKLDSplitLeaves_b(SEG_PARAMETERS * segment_param);
    vector<SplitKLD_b *> GetVectorKLDSplitLeavesUniform_b(SEG_PARAMETERS * segment_param);
    bool IsSplitChecked();
    SplitKLD * GetToSplit(SEG_PARAMETERS * segment_param);

    SplitKLD_b * GetToSplit_b(SEG_PARAMETERS * segment_param);
    vector<int> OrderingSplittingChildren();
    vector<SplitKLD*> OrderingSplittingLeaves(SEG_PARAMETERS* segment_param);
    vector<SplitKLD_b *> OrderingSplittingLeaves_b(SEG_PARAMETERS * segment_param);
    vector<SplitKLD_b *> OrderingSplittingLeavesUniform_b(SEG_PARAMETERS * segment_param, int InitSplitUnif);
    vector<SplitKLD*> GetSplitOrder();
    vector<SplitKLD_b*> GetSplitOrder_b();
    vector<SplitKLD*> GetSplitMoreVertical(int numbclasses,SEG_PARAMETERS* segment_param);
    vector<SplitKLD_b*> GetSplitMoreVertical_b(int numbclasses, SEG_PARAMETERS * segment_param);
    int * Combination(int n,int k);
    int * CombinationBis(int * TabSize,int numbclasses);
    vector<SplitKLD*> GetToSplitMore(vector<int> ClassesToSplit,SEG_PARAMETERS * segment_param);
    vector<SplitKLD_b *> GetToSplitMore_b(vector<int> ClassesToSplit, SEG_PARAMETERS * segment_param);
    vector<SplitKLD*> GetSplitOrderVertical(SEG_PARAMETERS * segment_param);
    vector<SplitKLD_b *> GetSplitOrderVertical_b(SEG_PARAMETERS * segment_param);
    inline int Factorial(int x);
    int NumbComb(int k, int n);
    void SplitOperation(SplitKLD * SplitTry, int choiceInit, SEG_PARAMETERS* segment_param);
    void SplitOperation(SplitKLD_b * SplitTry, int choiceInit,SEG_PARAMETERS * segment_param);
    void SplitOperation(int ChildToSplit,int choiceInit,SEG_PARAMETERS* segment_param);
    Parameters ** ParametersDoubleForSplittingUniform(int choiceInit);
    Parameters ** ParametersForSplitting(int choiceInit);
    Parameters * ParametersForSplittingUniform3(SEG_PARAMETERS * segment_param,int InitChoiceSplitUnifKMeans);
    int * KMeansInitialisationForUniform(float ** MeanK, int numbLabels,SEG_PARAMETERS * segment_param);
    void MultiplyNii(nifti_image * ImageToMultiply, float MultFactor);
    void MultiplyNii(nifti_image * ImageToMultiply, nifti_image * MultiplicativeImage);
    nifti_image * AddNii(nifti_image * ImageToChange, float FloatToAdd);
    nifti_image * CreateNormaliseOppositeImage(nifti_image * PriorsOutliers);
    float * MeanForSplitInitialisation(float MaxEigen, float * MaxVector);
    float * MeanForSplitInitialisationUniform(float MaxEigen, float* MaxVector);
    void SortingSVDEigen(float * SVDEigen, int * Index);
    float * VarianceForSplitInitialisation(float MaxEigen, float * MaxVector);
    float * VarianceForSplitInitialisationUniform(float MaxEigen, float * MaxVector);
    TreeEM * RunSplitOperation(SplitKLD * SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param);
    bool AreSplitPossible(vector<SplitKLD *> SplitTry, SEG_PARAMETERS * segment_param);
    TreeEM * RunSplitOperation(vector<SplitKLD *> SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param);
    TreeEM* RunSplitOperation_b(vector<SplitKLD_b *> SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param);
    bool AreSplitPossible(vector<SplitKLD_b *> SplitTry, SEG_PARAMETERS * segment_param);
    float * AdaptPriors(SEG_PARAMETERS * segment_param);
    void AdaptPriorsChildren(SEG_PARAMETERS * segment_param);
    void AdaptPriorsAllLevels(SEG_PARAMETERS * segment_param);
    void AdaptPriorsOM3PK8(SEG_PARAMETERS * segment_param);
    void AdaptPriorsAllLevelsOM3PK5(SEG_PARAMETERS * segment_param);
    void AdaptPriorsOM7PK5(SEG_PARAMETERS * segment_param);
    void FullAdaptPriors(SEG_PARAMETERS * segment_param);
    nifti_image * BuildConstantPriors(float WeightInput);
    nifti_image * BuildOutliersPriors(SEG_PARAMETERS * segment_param);
    float * BuildAtypicalityMap(SEG_PARAMETERS * segment_param);
    nifti_image * BuildOutliersPriors_bis(float MahalDistance);
    nifti_image * CreatePriorsFromAdaptedPriors();
    nifti_image * CreateOutBrainPriors(SEG_PARAMETERS * segment_param);
    vector<nifti_image *> CreatePriorsFromSplittingUniform(SEG_PARAMETERS * segment_param,int InitSplitUnif);
    nifti_image * TransformNormRespIntoPriors(SEG_PARAMETERS * segment_param);
    nifti_image * TransformArrayIntoPriors(float * NewSeg, SEG_PARAMETERS * segment_param);
    // Methods for bias field correction
    float * MakeBasisFunctions();
    float * MakeWMatrixChildren();
    float * MakeWMatrixLeaves();
    float * MakeRMatrixChildren();
    float * MakeAtWAMatrixChildren();
    float * MakeAtWRVectorChildren();
    float * MakeInvAtWAMatrixChildren();
    float * MakeFinalBFCoeffsChildren();
    float * MakeBFCoeffsDirect();
    float * MakeBFCoeffsSeparated(vector<int> BFOrderperModality, int IndexModalityTogether);
    float * MakeBFCorrection(int BForder);
    float * GenerateBFCorrection(int BForder, int numbmodal, float * BFCoeffs_newB, nifti_image* DataRef=NULL);
    void UncorrectData(int BForder);
    float * MakeBFCorrectionSeparated(vector<int> BFOrderperModality);
    float * MakeDataBFCorrected();
    float * MakeDataBFCorrectedSeparated(vector<int> BFOrderperModality);

    // Saving Results
    void PrintGMatrix();
    void SaveAllClasses(string filenameOut, SEG_PARAMETERS * segment_param);
    void SaveGeneralClasses(string filenameOut, SEG_PARAMETERS * segment_param);
    void SaveBFBasisFunctions(string filenameBF);
    void SaveBFBasisFunctions(int BFnumber, string filenameBF);
    void SaveBFCorrection(string filenameBF);
    void SaveBFCorrectedData(string filenameBF);
    void SavePriorsAdapted(SEG_PARAMETERS* segment_param, string filename);
    void SavePriorsAdaptedHierarchy(SEG_PARAMETERS * segment_param);
    void SaveAllLesionClasses(SEG_PARAMETERS * segment_param);
    void SaveLesionClass(SEG_PARAMETERS * segment_param);


    bool * GetCSFOutIndicVL(SEG_PARAMETERS * segment_param);
    bool * GetMSLesIndicVL(SEG_PARAMETERS * segment_param);
    float * GetMSOutlierBeliefs(SEG_PARAMETERS * segment_param);
    float * GetCSFOutlierBeliefs(SEG_PARAMETERS * segment_param);
    float * GetOutlierBeliefs(SEG_PARAMETERS * segment_param);
    nifti_image * CreateTypicalityAtlas(SEG_PARAMETERS * segment_param);
    nifti_image * Sum4DNii(nifti_image * Image4D);
    float * MakeTypicalityFloatClass(int IndexClass, SEG_PARAMETERS * segment_param);
    nifti_image * CreateTypicalityAtlas4D(SEG_PARAMETERS * segment_param);
    float * MakeTypicalityWeight(SEG_PARAMETERS * segment_param);
    vector<TreeEM *> GetWMLesions(SEG_PARAMETERS * segment_param);
    vector<float *> GetMeanGeneralClassesVector();
    vector<float *> GetVarianceGeneralClassesVector();
    int * GetModalities(SEG_PARAMETERS * segment_param);
    int GetLevel();
    int GetNumberLevels();
    vector<TreeEM*> GetAllTreesFromLevel(int l);
    int FindIndex(TreeEM * TreeToLook);
    TreeEM * FindUnifDist(SEG_PARAMETERS * segment_param);
    vector<TreeEM *> GetUniformLeavesVector();
    vector<vector<int> > FindAllUniformDistHierarchies(SEG_PARAMETERS * segment_param);
    int FindGeneralClass();
    TreeEM* FindGeneralClassPriors();
    int * GetHierarchy();
    vector<int> GetHierarchyVector();
    TreeEM * FindFromHierarchy(vector<int> HierarchyVector);
    TreeEM * FindMainNode();

    void SaveTreeInTextFile(string filenameOut,SEG_PARAMETERS * segment_param);
    void ModifyCountFiles(vector<string> CountFiles);

    float * GaussianBlurring(float * CountHistogram, float gauss_std, vector<int> dim);

    // MRF functions
    void ClearMRFModel();
    float * CopyGMatrix();
    float * CopyMRF();
    float * GetMRF();
    float * GetGMatrix();
    float * GetGMatrixDirect();
    void SetGMatrix(float * GMatrixToSet);
    void SetMRF(float * MRFToSet);
    void UpdateMRF(SEG_PARAMETERS * segment_param);
    bool PrepareGInfoFromFile(string GMatrixFilename);
    float * CreateGMatrixFromInfo(bool optMRFOut,string GMatrixFilename);
    float * PrepareGMatrixFromFile(string GMatrixFilename, bool optMRFOut);
    float * GMulSumNeighbNormRespExp(SEG_PARAMETERS * segment_param);
    float * GetNormRespShiftedSumN();
    float * GetSumNeighboursNormResp_bis();
    float * GetSumNeighboursNormResp(SEG_PARAMETERS * segment_param);
    int * CorrespondingCoord(int Index, int * Shift);
    float * GetNormRespShifted(int dim, int dir,SEG_PARAMETERS * segment_param);
    int * MakeHardSegCombined();
    int * MakeHardSegLeaves();
    bool * GetSegOutliers();
    int GetNumberInliers();
    int * HardSegLeavesShifted(int dim,int dir);
    int * SumHardSegNeighboursLeaves();
    float * MakeHistogramHardSegLeaves();
    float * HistogramSoftSegMRF(SEG_PARAMETERS * segment_param);
    float * MRFOptMakeLogRatio(SEG_PARAMETERS * segment_param);
    float * MRFOptSolveLS(SEG_PARAMETERS * segment_param);
    int * MRFOptAMatrix();

    void SaveTmpResultMasked(float * ResultToSave,string filename);
    void SaveTmpResult(float * ResultToSave,string filename);
    void SaveMRFImage(SEG_PARAMETERS * segment_param);
    void SavePriorsAdapted(SEG_PARAMETERS* segment_param);

    nifti_image * ReadFromFilename(string Filename);
    vector<nifti_image *> ReadFromFilenamesVector(vector<string> FilenamesVector);


    void ReinitialisationForModalityChange();
    bool HaveSimilarStructure(TreeEM * ExistingTree);
    void ResuscitateUniform(SEG_PARAMETERS * segment_param);
    void RepopulateNormResp(TreeEM * ExistingTree);
    void ProgressivePriorsAdaptation(SEG_PARAMETERS * segment_param);
    TreeEM * BuildTreeWithAddedModalityFromExistingModel(TreeEM * ExistingTree, SEG_PARAMETERS * segment_param);

    float GetMahalDist(float * ValueArray, float * MeanArray, float * InvertVarianceArray, int numbmodal);

    template <class NewTYPE>
    int seg_changeDatatype(nifti_image *image)
    {
        switch(image->datatype){
        case DT_BINARY:
            seg_changeDatatype1<NewTYPE,bool>(image);
            break;
        case NIFTI_TYPE_UINT8:
            seg_changeDatatype1<NewTYPE,unsigned char>(image);
            break;
        case NIFTI_TYPE_INT8:
            seg_changeDatatype1<NewTYPE,char>(image);
            break;
        case NIFTI_TYPE_UINT16:
            seg_changeDatatype1<NewTYPE,unsigned short>(image);
            break;
        case NIFTI_TYPE_INT16:
            seg_changeDatatype1<NewTYPE,short>(image);
            break;
        case NIFTI_TYPE_UINT32:
            seg_changeDatatype1<NewTYPE,unsigned int>(image);
            break;
        case NIFTI_TYPE_INT32:
            seg_changeDatatype1<NewTYPE,int>(image);
            break;
        case NIFTI_TYPE_FLOAT32:
            seg_changeDatatype1<NewTYPE,float>(image);
            break;
        case NIFTI_TYPE_FLOAT64:
            seg_changeDatatype1<NewTYPE,double>(image);
            break;
        default:
            fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tThe initial image data type (%d) is not supported\n",image->datatype);
            exit(1);
        }
        return 1;
    }

    template <class T> T* MakeLong(T *ArrayToLong, int * L2S, int numel){
        T* ArrayLonged=new T[numel];
        for(int i=0;i<numel;i++){
            if(L2S[i]<0)
            ArrayLonged[i]=0;
            else{
                ArrayLonged[i]=ArrayToLong[L2S[i]];
            }
        }
        return ArrayLonged;

    }

    template <class T> T* MakeLongPad(T *ArrayToLong, T ValuePadding, int * L2S, int numel){
        T* ArrayLonged=new T[numel];
        for(int i=0;i<numel;i++){
            if(L2S[i]<0)
                ArrayLonged[i]=ValuePadding;
            else{
                ArrayLonged[i]=ArrayToLong[L2S[i]];
            }
        }
        return ArrayLonged;

    }

    template< class T> T* MakeSmall(T * ArrayToSmall, int * S2L, int numelmasked){
        if (ArrayToSmall==NULL){
            return NULL;
        }
        T* ArraySmalled=new T[numelmasked];
        for(int i=0;i<numelmasked;i++){
            ArraySmalled[i]=ArrayToSmall[S2L[i]];
        }
        return ArraySmalled;
    }

    template <class T> bool * OpposeBoolArray(T * ArrayToOppose,int numel){
        if(ArrayToOppose==NULL){
            return NULL;
        }
        bool * OpposedArray=new bool[numel];
        for(int i=0;i<numel;i++){
            OpposedArray[i]=!ArrayToOppose[i];
        }
        return OpposedArray;
    }

    template <class T1, class T2> T2* TranscribeArray(T1 * ArrayToTranscribe, int numel){
        if(ArrayToTranscribe==NULL){
            return NULL;
        }
        T2 * TranscribedArray=new T2[numel];
        for(int i=0;i<numel;i++){
            TranscribedArray[i]=(T2)ArrayToTranscribe[i];
        }
        return TranscribedArray;
    }

    template <class T> T GetMaxArray(T* ArrayToMax, int numel){
        if (ArrayToMax==NULL){
            return -std::numeric_limits<T>::max();
        }

        T MaxResult=-std::numeric_limits<T>::max();
        for(int i=0;i<numel;i++){
            if (MaxResult<ArrayToMax[i]){
            MaxResult=ArrayToMax[i];
            }
        }
        return MaxResult;
    }

    template <class T> int GetLabelSize(T ValueToCount, T* ArrayToCheck, int numel){
        if(ArrayToCheck==NULL){
            return -1;
        }
        int CountValue=0;
        for(int i=0;i<numel;i++){
            if(ArrayToCheck[i]==ValueToCount){
                CountValue++;
            }
        }
        return CountValue;
    }


    template<class T> T* Repmat(int * Pattern, int size, int L, int C){
        int sizeTot=size*L*C;
        int sizeMin=size*L;
        T * Result=new T[sizeTot];
        for (int s=0; s<size; s++) {
            for (int l=0; l<L; l++) {
                Result[s*L+l]=Pattern[s];
            }
        }
        for(int c=0;c<C;c++){
            for(int s=0;s<sizeMin;s++){
                Result[c*sizeMin+s]=Result[s];
            }
        }
        return Result;
    }

    template<class T> vector<int> GetIndicesValue(T * ArrayToCheck, T ValueToFind, int numel){
        vector<int> IndicesValue;
        if(ArrayToCheck==NULL){
            return IndicesValue;
        }
        for(int i=0;i<numel;i++){
            if(ArrayToCheck[i]==ValueToFind){
                IndicesValue.push_back(i);
            }
        }
        return IndicesValue;
    }



    template <class NewTYPE, class DTYPE>
    int seg_changeDatatype1(nifti_image *image)
    {
        // the initial array is saved and freeed
//        DTYPE *initialValue = (DTYPE *)malloc(image->nvox*sizeof(DTYPE));
        DTYPE * initialValue=new DTYPE[image->nvox];
        DTYPE *ImagePTR=static_cast<DTYPE *>(image->data);
//        memcpy(initialValue, image->data, image->nvox*sizeof(DTYPE));
        int numel=image->nvox;
//        int numel=image->nx*image->ny*image->nz;
        for(int i=0;i<numel;i++){
            initialValue[i]=ImagePTR[i];
        }

        // the new array is allocated and then filled
        if(sizeof(NewTYPE)==sizeof(unsigned char)) image->datatype = NIFTI_TYPE_UINT8;
        else if(sizeof(NewTYPE)==sizeof(float)) image->datatype = NIFTI_TYPE_FLOAT32;
        else if(sizeof(NewTYPE)==sizeof(double)) image->datatype = NIFTI_TYPE_FLOAT64;
        else{
            fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tOnly change to unsigned char, float or double are supported\n");
            exit(1);
        }
        free(image->data);
        image->nbyper = sizeof(NewTYPE);
        image->data = (void *)calloc(image->nvox,sizeof(NewTYPE));
        NewTYPE *dataPtr = static_cast<NewTYPE *>(image->data);
        if(sizeof(NewTYPE)==sizeof(unsigned char)){
            for(unsigned int i=0; i<image->nvox; i++) dataPtr[i] = (unsigned char)(round(initialValue[i]));
        }
        else{
            for(unsigned int i=0; i<image->nvox; i++) dataPtr[i] = (NewTYPE)(initialValue[i]);
        }
//        nifti_set_filenames(image,"/Users/Carole/Documents/PhD/TestImage.nii.gz",0,0);
//        nifti_image_write(image);

//        free(initialValue);
        delete [] initialValue;
        return 1;
    }
//    METHODS TO CORRECT FOR WRONG HOLES IN WM

//    ComponentLabel methods
    void    RecurseCompLabel(float * Data, int * PositionTmp, int Label, int * ComponentsLabel,int Neigh,int * Dim,int * Shift, float thresh,int& Times);
    int * ComponentLabeling(float * ImageToLabel, int Neigh, int * Dim, int * Shift, float thresh);
    void RelabelComponents(int * ComponentsLabel, int * Dim);
    void RepassComponentsLabel(int * ComponentsLabel, int Neigh,int * Dim,int * Shift);
    int * GetListNeighbours(int CurrentIndex, int * Dim, int * Shift,int NeighbourhingType);
    bool * CreateLesionBool(int * LesionLabeling, int Label,int numel );
    vector<int > GetIndicesBorderLesion (int * OrderedLabels, int Label, int * Dim, int * Shift);
    float * ProportionNeighboursLesion(int * LesionLabeling, int Label, int * HardSegLong,SEG_PARAMETERS * segment_param);
    void  GetListNeighbours_bis(int * ListNeighbours, int CurrentIndex, int * Dim, int * Shift,int NeighbouringType);
    float * GetLesionGMCSFSeg(int * HardSeg, int IndexClass, SEG_PARAMETERS * segment_param);

nifti_image *  CorrectInitialEMResultForWrongWMHoles(SEG_PARAMETERS * segment_param);
    void CorrectNormResp(float * Correction, int InitialClass, int FinalClass, int CorrectionType);
   void CorrectNormRespToUniform(float * Correction, int GeneralInit, int UniformFin, int CorrectionType);

 void   CorrectWMOforCSFInclusion(SEG_PARAMETERS * segment_param);
    int * GetCorrespondingModality(SEG_PARAMETERS * segment_param);

    bool * GetWMDGMSegBool(int * HardSegLeaves,SEG_PARAMETERS* segment_param);


    void RebuildNormRespFromLeaves();
    void NullAllNonLeavesNormResp();
    float * SumNormRespChildren();

    int * GetCorrespondanceOrderedVolume(int * VolumeLabels, int maxLabel);
    int * OrderedVolumeLabel(int * ComponentLabels, int MiniSize, int numel);
    int * GetVolumeLabels(int * ComponentLabel, int numel);
    bool IsOutsideMask(vector<int> ListIndices);
};



