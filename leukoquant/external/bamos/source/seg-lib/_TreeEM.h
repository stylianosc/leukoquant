//
//  TreeEM.h
//  TreeSeg
//
//  Created by Carole Sudre on 08/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 07/05/2026.
//

#ifndef TreeSeg_TreeEM_h
#define TreeSeg_TreeEM_h

#include <iostream>
#include <vector>
#include "nifti1.h"
#include "nifti1_io.h"

using namespace std;

class TreeEM;

typedef struct Parameters{
    int DistributionType;
    int SizeParameters;
    float * ValueParameters;
    Parameters(){
        DistributionType=0;
        SizeParameters=0;
        ValueParameters=NULL;
    }
~Parameters(){
    if(ValueParameters!=NULL){
        delete [] ValueParameters;
        ValueParameters=NULL;
    }
}
} Parameters;

typedef struct MergeKLD{
    TreeEM * Parent;
    int Child1;
    int Child2;
    float KLDTot;
    float * KLD;
    
    MergeKLD(){
        Parent=NULL;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
        KLD=NULL;
    }
    
    ~MergeKLD(){
        Parent=NULL;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
        if(KLD!=NULL){
            delete [] KLD;
            KLD=NULL;
        }
    }
}MergeKLD;

typedef struct SplitKLD{
    TreeEM * Parent;
    int ChildToSplit;
    float KLDTot;
    float * KLD;
    
    SplitKLD(){
        Parent=NULL;
        ChildToSplit=-1;
        KLDTot=-1;
        KLD=NULL;
    }
    
    ~SplitKLD(){
        if(KLD!=NULL){
            delete [] KLD;
            KLD=NULL;
        }
    }
    
} SplitKLD;

class TreeEM{
    TreeEM * Parent; // NULL if Root of the Tree
    vector<TreeEM *> Children; // Children of the tree : NULL if leaf
    nifti_image * DataImage;
    nifti_image * Mask;
    nifti_image * Priors;
    int * S2L;
    int * L2S;
    float * NonNormResp; // Non normalized responsability
    float * NormResp; // Normalised responsability
    float * Distribution;// Probabilistic distribution
    float NonNormWeight; // Non normalised weight
    float NormWeight; // Normalised weight
    Parameters * ParametersDistribution; // Parameters corresponding to the distribution ( non null pointer for the value only if leaf of the tree and therefore non mixture distribution)
    bool * SplitCheck;
    bool * MergeCheck;
    float * BFCoeffs;
    //float * BFBasisFunctions;
    float * BFCorrection;
    public :
    
    // Constructors, copy and destructor
    TreeEM();
    TreeEM(nifti_image * Data);
    
    TreeEM * CopyTree(TreeEM * ParentNew);
    nifti_image * CopyImage();
    nifti_image * CopyMask();
    nifti_image * CopyPriors();
    int * CopyS2L();
    int * CopyL2S();
    float * CopyDistribution();
    float * CopyNonNormResp();
    float * CopyNormResp();
    Parameters * CopyParameters();
    bool * CopyMergeCheck();
    bool * CopySplitCheck();
    float * CopyBFCorrection();
    float * CopyBFCoeffs();
    ~TreeEM();
    void MakeEmpty();
    
    // Get methods
    TreeEM * GetParent();
    vector<TreeEM * >GetChildren();
    int GetNumberChildren();
    TreeEM * GetChild(int c);
    nifti_image * GetDataDirect();
    nifti_image * GetDataImage();
    float GetMaxDataModal(int m);
    float GetMinDataModal(int m);
    float * GetDataModalMasked(int m);
    float * GetDataShifted(int dim,int dir);
    float * GetDataCorrelation(int dim);
    float GetIndFactor(int dim);
    float GetIndFactorTot();
    int GetNumberElements();
    int GetNumberModalities();
    int GetNumberMaskedElements();
    float * GetDistribution();
    int GetNumbbins();
nifti_image * GetMask();
    nifti_image * GetMaskDirect();
    nifti_image * GetPriors();
    nifti_image * GetPriorsDirect();
int * GetS2L();
    int * GetS2LDirect();
int * GetL2S();
    int * GetL2SDirect();
    void MakeL2S();
    void MakeS2L();
float * GetNonNormResp();
float * GetNormResp();
float  GetNonNormWeight();
float  GetNormWeight();
    float GetPartNormWeight();
Parameters * GetParameters();
int GetDistributionType();
    float * GetMean();
    float * GetVariance();
    float * GetMeanDirect();
    float * GetDiagVarianceDirect();
    float * GetParametersValue();
    int GetSizeParameters();
    vector<nifti_image*> GetPriorsVector();
    void GetLeaves(vector<TreeEM*> & LeavesVector);
    vector<TreeEM *> GetAllLeaves();
    vector<TreeEM *> GetAllDirectLeaves();
    int GetNumberDirectLeaves();
    bool * GetSplitCheck();
    bool * GetMergeCheck();
    void GetNumberLeaves(int & numbleaves);
    int GetNumberAllLeaves();
    int GetNumberFreeParameters();
    float GetLogLikelihood();
    float * GetBFCoeffs();
    float * GetBFCoeffsDirect();
    //float * GetBFBasisFunctionsDirect();
    //float * GetBFBasisFunctions();
    float * GetBFCorrection();
    float * GetBFCorrectionDirect();

    enum TreeComponent {
        ROOT = 0,
        BRANCH=1,
        LEAF = 2,
        INITIALNODE=3
    };
    enum PartialResultType {
        DISTRIBUTION=0,
        NONNORMRESP=1,
        NORMRESP=2,
    };
    
    

    
    float * GetPartialResult(PartialResultType ResultType);
    
    // Set methods
    void SetParent(TreeEM * Parent);
    void ReinitialiseTree();
    void SetChildren(vector<TreeEM *> Children);
    void AddChild(TreeEM * ChildAdded);
    void SetData(nifti_image * data);
    void NormaliseDataImage();
    void ConvertDataImageToFloat();
    void SetDistribution(float * Distribution);
    void SetMask(nifti_image * MaskImage);
    void MakeMaskBinary();
    void SetPriors(nifti_image * Priors);
    void MakePriorsFloat();
    void MakePriorsProbabilityType();
    void SetS2L(int * S2L);
    void SetL2S(int * L2S);
    void SetNonNormResp(float * NonNormResp);
    void SetNormResp(float * NormResp);
    void InitialiseNormRespWithPriors();
    void SetNonNormWeight(float NonNormWeight);
    void SetNormWeight(float  NormWeight);
    void SetParameters(Parameters * ParametersToSet);
    void SetDistributionType(int DistributionType);
    void SetBFCoeffs(float * BFCoeffs);
    void SetBFCorrection(float * BFCorrection);
    void MakeParametersMixture();
    

// EM related methods
    void UpdateParameters();
    void UpdateGaussianParameters();
    void UpdateBFCoeffs();
    void UpdateBFCorrection();
    void UpdateGaussianMean();
    void UpdateGaussianVariance();
    void UpdateDistribution();

    void SumToDistributionWeightedDistributionOfChild(int c);

    void ReputDistributionToZero();
    void UpdateGaussianDistribution();

    void UpdateNonNormResp();

    void ReputNonNormRespToZero();
    void CalculatingNonNormResp();

    void SumToNonNormRespWeightedNonNormRespOfChild(int c);

    void MultiplyCurrentNonNormRespByNormWeight();

    void UpdateNormResp();
    void NormaliseNonNormResp();

    void DivideNonNormRespByRootNonNormResp();

    void UpdateNonNormWeights();
    void UpdateNormWeights();
    float CompleteLogLikelihood();
    void EMCompleteIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration);
    void RunFullEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration);
    TreeEM * RunFullBiASM();
 
    //Tree related methods
    void CreateAndAddChild();
    void CreateAndAddChild(nifti_image * Priors,int DistributionType);
int GetLevel();
    void CollapseOnlyChild();
    void CollapseOnlyChildTot();
    void CollapseChild(TreeEM * ChildToCollapse);
    void AllLeavesToChildrenLevel();
    int Index(TreeEM * ElementToFind);
    void FindAllAtCertainLevel(vector<TreeEM*> & ResultVector, int Level);
    void FindDepth( int & ResDepth);
    TreeEM * FindRoot();
    void DividePartNonNormResp(float * SumNonNormResp);
    void InitialiseBasicTree();
    void ClearSMChecks();
    void PutAllLeavesToChildrenLevel();
    void ModifyNormWeightsForChildrenLevel();
    void RecPutAllLeavesToChildrenLevel();

// Tests on tree
bool IsEqual(TreeEM * TreeEqualityCheck);
bool IsLeaf();
    bool IsBranch();
    bool IsRoot();
void DoesBelongToCurrentTree(TreeEM * Element, bool & test);
int WhatTypeOfTreeComponent();
    bool IsTreeValid();
    bool IsBasicTree();
    bool IsNumbChildOKWithDistType();
    bool IsDataImageNormalised();
    bool IsDataFloat();
    bool IsPriorsFloat();
    bool AreBFCoeffsDirectUseable();
    // Parameters structure
    
    int CalculateSizeParameters();
    int CalculateSizeParameters(int DistributionType);
    void CreateAllocateAndInitializeParameters(int DistributionType);
    bool CheckForValidityOfParametersStructure();
    bool CheckForSizeParametersValidity();
    
    bool ArePriorsNormalised();
    bool IsPriorsVectorFilledWithNULL();
    bool IsOneOfThePriorsNULL();
    void NormalisePriors();
    bool IsPriorsProbabilityType();
    bool CheckForValidityOfParametersStructure(Parameters * ParametersToCheck);
    bool CheckForSizeParametersValidity( Parameters * ParametersToSet);
    
    /**************** METHODS FOR SM OPERATIONS *********************/
    // Functions needed for calculation of KLD

    vector<float *> GetDistHistogram();
    float * GetDistHistogramTotal();
    vector<float *> GetGaussianDistributionHist();
    float * GetGaussianDistributionHistTotal();
    vector<float *> GetDataHistogram();
    float * GetDataHistogramTotal();
    
    float GetDistCompKLDTotal();
    float * GetDistCompKLD();
    
    float CriterionCalculation();
    
    // Functions needed for merging operations 
    
    MergeKLD * GetKLDMerge(int Child1,int Child2);
    vector<MergeKLD *> GetVectorKLDMergeChildren();
    bool IsMergeChecked();
    vector<MergeKLD *> GetVectorKLDMergeLeaves();
    MergeKLD * GetToMerge();
    void MergeOperation(int Child1,int Child2);
    Parameters * ParametersForMerging(int Child1,int Child2);
    TreeEM * RunMergeOperation();
    
    // Functions needed for splitting operations
    SplitKLD * GetKLDSplit(int ChildToSplit);
    vector<SplitKLD *> GetVectorKLDSplitChildren();
    vector<SplitKLD *> GetVectorKLDSplitLeaves();
    bool IsSplitChecked();
    SplitKLD * GetToSplit();
    void SplitOperation(int ChildToSplit);
    Parameters ** ParametersForSplitting();
    float * MeanForSplitInitialisation();
    float * VarianceForSplitInitialisation();
    TreeEM * RunSplitOperation();
    
    // Methods for bias field correction
    float * MakeBasisFunctions();
    float * MakeWMatrixChildren();
    float * MakeWMatrixLeaves();
    float * MakeRMatrixChildren();
    float * MakeAtWAMatrixChildren();
    float * MakeAtWRVectorChildren();
    float * MakeInvAtWAMatrixChildren();
    float * MakeFinalBFCoeffsChildren();
    float * MakeBFCorrection();
    float * MakeDataBFCorrected();
    
    // Saving Results
    void SaveAllClasses(char * filenameOut);
    void SaveGeneralClasses(char * filenameOut);
    void SaveBFBasisFunctions(char * filenameBF);
    void SaveBFBasisFunctions(int BFnumber, char * filenameBF);
    void SaveBFCorrection(char * filenameBF);
    void SaveBFCorrectedData(char * filenameBF);
    void SaveTreeInTextFile(char * filenameOut);

};





#endif
