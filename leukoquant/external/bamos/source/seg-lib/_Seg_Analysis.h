#ifndef _SEG_ANALYSIS_H
#define _SEG_ANALYSIS_H


#include "nifti1_io.h"
#include "_TreeEM_new.h"
#include "_Seg_InitAndRun.h"
#include <iostream>
#include<sstream>
#include<string>
#include<set>
#include<list>
#include<deque>
#include "math.h"
#include<numeric>
#include<stdlib.h>
#include<algorithm>
//#include <boost/format.hpp>
using namespace std;






enum log_level_t {
    LOG_NOTHING=0,
    LOG_CRITICAL=1,
    LOG_ERROR=2,
    LOG_WARNING=3,
    LOG_INFO=4,
    LOG_DEBUG=5
};

log_level_t GetLevelMap(int level_int);

static log_level_t GLOBAL_LEVEL=LOG_NOTHING;
//namespace log_impl {
//class formatted_log_t {
//public:
//    formatted_log_t( log_level_t level, const wchar_t* msg ) : fmt(msg), level(level) {}
//    ~formatted_log_t() {
//        // GLOBAL_LEVEL is a global variable and could be changed at runtime
//        // Any customization could be here
//        if ( level <= GLOBAL_LEVEL ) wcout << level << L" " << fmt << endl;
//    }
//    template <typename T>
//    formatted_log_t& operator %(T value) {
//        fmt % value;
//        return *this;
//    }
//protected:
//    log_level_t     level;
//    boost::wformat      fmt;
//};
//}//namespace log_impl
//// Helper function. Class formatted_log_t will not be used directly.
//template <log_level_t level>
//log_impl::formatted_log_t log(const wchar_t* msg) {
//    return log_impl::formatted_log_t( level, msg );
//}


struct SEG_ANALYSIS{
    bool flag_VesselSeed;
    bool flag_CorrectConnect; // For consistency of correction in longitudinal analysis
    bool flag_infarcts;
    bool flag_PV;
    bool flag_layers;
    bool flag_layercreation;
    bool flag_simple;
    bool flag_WMCheck;
    bool flag_ITCheck;
    bool flag_CSFCheck;
    bool flag_test;
    bool flag_HS;
    bool flag_DCh;
    bool flag_DCs;
    bool flag_TP;
    bool flag_TN;
    bool flag_FN;
    bool flag_FP;
    bool flag_Acc;
    bool flag_Sens;
    bool flag_Spec;
    bool flag_PSI;
    bool flag_Corr;
    bool flag_FPR;
    bool flag_TPR;
    bool flag_VD;
    bool flag_DE;
    bool flag_OE;
    bool flag_TLLSoft;
    bool flag_TLLHard;
    bool flag_fileOut;
    bool flag_outTxt;
    bool flag_numbActive;
    bool flag_in;
    bool flag_data;
    bool flag_BF;
    bool flag_removespurious;
    bool flag_outliers;
    bool flag_KLD;
    bool flag_SSD;
    bool flag_TextFile;
    bool flag_SegTot;
    bool flag_Les;
    bool flag_LaplaceSolImage;
    bool flag_LapAnalysis;
    bool flag_refLes;
    bool flag_compJoint;
    bool flag_compMult;
    bool flag_mask;
    bool flag_maskMatch;
    bool flag_WMCard;
    bool flag_Distance;
    bool flag_NW;
    bool flag_changePath;
    bool flag_FLAIRonly;
    bool flag_correct;
    bool flag_RuleTextFile;
    bool flag_connect;
    bool flag_CorrectIV;
    bool flag_AvDist;
    bool flag_Parc;
    int flag_VentrSeg;
    bool flag_Hemi;
    bool flag_Analysis;
    bool flag_AnalysisConnect;
    bool flag_GivenModa;
    bool flag_GivenMahalRule;
    bool flag_AcceptedGM;
    bool flag_Quadrant;
    bool flag_MNITransform;
    bool flag_MNITemplate;
    bool flag_inLes;
    bool flag_outWMI;
    bool flag_getCH;
    bool flag_BFGen;
    bool flag_AffGen;
    bool flag_CodedOutliers; // If we perform the coding of the ouliers
    bool flag_OtherSeg;
    bool flag_VesselRuleTextFile;
    bool flag_vecIOLap; // If the In and Out segmentations to use for the laplace analysis are provided
    bool flag_nameLap;
    bool flag_outConnect;
    bool flag_TO; // If we want to consider all outliers for the building of the lesion instead of selecting special classes
    bool flag_IO; // Inliers/Outliers separation from Tree
    bool flag_inSum; // On if there is file for summarised segmentation.
    bool flag_inOptionText;
    bool flag_refLesConnect; // On if there is a name file for the connected components of the GT
    bool flag_GT; // ON if there is analysis to do on the reference image to use
    bool flag_GTSave; // if ON means that we do and save the GT analysis..
    bool flag_inConnect; // if already connected components of lesion are provided.
    bool flag_PrintConnect; // if we have to print in 3 directions all modalities the selected area.
    bool flag_ConnectRefAnalysis;
    bool flag_Evaluation;
    bool flag_Saving; // To save everything...
    bool flag_mean; // for mean of vector of images
    bool flag_inFPTPReclassif;
    bool flag_inWMSeg;
    bool flag_acceptJCMS;
    bool flag_correctCerr;
    float flag_Secondary;
    bool flag_inDataComp;
    bool flag_inLesCorr;
    bool flag_inVentricleSeg;
    bool flag_inPriorsICSF;
    bool flag_inPriorsECSF;
    bool flag_inPriorsOut;
    bool flag_inPriorsDGM;
    bool flag_inPriorsWM;
    bool flag_inPriorsCGM;
    bool flag_checkInliers;
    bool flag_juxtaCorrection;
    bool flag_RG; // Flag for region growing;
    bool flag_SPCorrection;
    bool flag_segGIF;
    bool flag_inArtefact;
    bool flag_WMMaps;
    bool flag_inMahal;
    bool flag_Vesselness;
    bool flag_outputDir;
    bool flag_LS;
    bool flag_sc;
    bool flag_Euc;
    bool flag_GetConnect;
    bool flag_LesionReduction;
    bool flag_GMPre; // If there is a file gathering islands of GM already selected as potential problematic regions
    bool flag_CSFPre; // Idem but for CSF
    bool flag_LCRuleTextFile;
    bool flag_inOutText;
    bool flag_inExFiles;
    bool flag_inROIFiles;
    bool flag_getCoM;
    bool flag_oldLesion;
    bool flag_inQuadrant;
    bool flag_Gaussian;
    bool flag_BFText;
    bool flag_inLobes;
    bool flag_inToAnalyse;
    bool flag_RuleCorr;
    bool flag_GIFPriors;
    bool flag_Dot;
    bool flag_Grad;
    bool flag_getEntropy;
    float flag_minSize;
    float dist_ring;
    bool flag_MapLap;
    bool flag_inAuthorised;
    vector<int> vecModalities;
    vector<int> vecRuleMahal;
    vector<int> vecCheckModalities;
    vector<float> Gaussian_Zscore;
    int Opt_JointLC; // Option for taking care of lesion origins separately or not 1 for joint 0 for separated
    int index_Gaussian; // Index of the class to consider when building Gaussian noise
    int numbMult;
    int flag_segWeighted;
    int ExploRadius;
    int LesionRuleType;
    int LesionUniformRuleType;
    int flag_correctionLevel;
    log_level_t verbose_level;
    bool flag_CST;
    int IndexWM;
    int IndexGM;
    int IndexCSF;
    int IndexOut;
    int IndexCGM;
    int IndexDGM;
    int IndexICSF;
    int IndexECSF;
    int IndexBrainstem;
    int Neigh;
    int MiniSize;
    int numbCoded; // Number of modalities to code for the outliers
    vector<int> numbLaminae;
    vector<float> BFParam;
    vector<float> AffParam;
    vector<float> vecLesionReduction;
    vector<float> vecLesionReductionGeneral; // Contains the options for the threshold, the step size and the blurring applied
    int flag_segType;
    float val_ext;
    int thresh_cortex;
    bool flag_CorExt;
    bool flag_LapIn;

    float flag_LesWMI; // Mahalanobis distance that will give a weight of 1 for the WMI part
    bool flag_match;
    bool flag_prion;
    bool flag_getBB;
    int orderMatch;
    float ThresholdSP;
    float thresh_quant;
    float weightThreshold;
    float thresh_RG; // Float to have for region growing
    int weightCompClass;
    int GTAnalysisType;
    float MahalThresholdCheck;
    int MaxValueEntr;
    int RadEntropy;
    float EdginessThreshold;
    vector<char*>filename_In;
    char* filename_compJoint;
    vector<char*> filename_compMult;
    char * filename_inAuthorised;
    char * filename_Hemi;
    char * filename_dot1;
    char * filename_dot2;
    char * filename_Out;
    char * filename_OutTxt;
    char * filename_mask;
    char * filename_maskMatch;
    char * filename_Dist;
    char * filename_Ref;
    char * filename_InLes;
    char * filename_TextFile;
    char * filename_SegTot;
    char * filename_RuleTextFile;
    char * filename_RuleCorr;
    char * filename_VesselRuleTextFile;
    char * filename_WMCard;
    char * filename_test;
    char * filename_Parc;
    char * filename_Artefact;
    char * filename_LaplaceSolImage;
    char * filename_RefConnect;
    char * filename_changePath;
    char * filename_inConnect;
    char * filename_inLesCorr;
    char * filename_inLes;
    char * filename_inDataComp;
    char * filename_inWMSeg;
    char * filename_inPriorsICSF;
    char * filename_inPriorsECSF;
    char * filename_inPriorsOut;
    char * filename_inPriorsDGM;
    char * filename_inPriorsCGM;
    char * filename_inPriorsWM;
    char * filename_inVentricleSeg;
    char * filename_inSum;
    char * filename_MNITransform;
    char * filename_MNITemplate;
    char * inOptionText;
    char * filename_saveMatch;
    char * filename_saveMean;
    char * filename_inFPTPReclassif;
    char * nameLap;
    char * filename_GMPre;
    char * filename_CSFPre;
    char * filename_LCRuleTextFile;
    char * filename_inOutText;
    char * filename_inQuadrant;
    char * filename_inLacMask;
    char * filename_inToAnalyse;
    char * filename_inMahal;
    char * filename_inCST;
    char * name_outputDir;
    char * filename_GIFPriors;
    vector<float> vecThresh;
    vector<char *> filename_inExFiles;
    vector<char *> filename_inROIFiles;
    vector<char*> filename_Data;
    vector<char*> filename_BF;
    vector<char*> filename_Outliers;
    vector<int> vecElementIn;
    vector<int> vecElementOut;
    vector<int> vecLeavesToAdd;
    vector<char*> filename_RefImages;
    vector<char *> filename_ImagesMean;
    vector<char*> filename_FloatImages;
    vector<char *> filenames_IOLap;
    vector<char *> filename_ImagesGaussian;
    vector<char *> filename_vectorRef;
    vector<char *> filename_vectorLes;
    vector<char *> filename_vectorNames;
    vector<char *> filename_vectorMask;
    vector<char *> filename_vectorOut;
    vector<char *> filename_vectorVentr;
    vector<char *> filename_inLobes;
    char * filename_BFText;
};

typedef struct Point {
  float x, y, z;
  Point(){
      float x;
      float y;
      float z;
  }
  Point(float X, float Y, float Z){
      x=X;
      y=Y;
      z=Z;
  }

  Point operator-(Point p) const {
    return Point(x - p.x, y - p.y, z - p.z);
  }

  Point operator+(Point p) const {
    return Point(x + p.x, y + p.y, z + p.z);
  }
  
  Point operator*(float f){
      return Point(x*f, y*f,z*f);
  }

  Point transcribe(float* p){
      return Point(p[0], p[1],p[2]);
  }

  Point cross(Point p) const {
    return Point(
      y * p.z - p.y * z,
      z * p.x - p.z * x,
      x * p.y - p.x * y
    );
  }

  float dot(Point p) const {
    return x * p.x + y * p.y + z * p.z;
  }

  float norm() const {
    return sqrt(x*x + y*y + z*z);
  }

  Point normalise(){
     return Point(x/norm(),y/norm(),z/norm());
  }
}Point;


typedef struct Face {
  vector< Point> v;
  Face(){
      vector< Point> v;
  }

  Face(vector < Point> V){
      v=V;
  }


  Face( Point V1,  Point V2,  Point V3){
      v.push_back(V1);
      v.push_back(V2);
      v.push_back(V3);
  }

  bool isDegenerate(){
      Point dir1 = v[1]-v[0];
      Point dir2 = v[2]-v[0];
      Point dir3 = v[2]-v[1];
      if(dir1.norm()==0 || dir2.norm()==0 || dir3.norm()==0){
          return true;
      }
      return false;
  }

  float area(){
      Point dir1 = v[1]-v[0];
      Point dir2 = v[2]-v[0];
      return 0.5*dir1.cross(dir2).norm();
  }
  Point greatestside(){
//      float G=0;
      Point S1=v[1]-v[0];
      Point S2=v[2]-v[0];
      Point S3=v[2]-v[1];
      if (S1.norm()<S2.norm()){
          if(S2.norm()<S3.norm()){
              return S3;
          }
          return S2;
      }
      else {
          if(S1.norm()<S3.norm()){
              return S3;
          }
          return S1;
      }
  }

  int greatestside_ind(){
//      float G=0;
      Point S1=v[1]-v[0];
      Point S2=v[2]-v[0];
      Point S3=v[2]-v[1];
      if (S1.norm()<S2.norm()){
          if(S2.norm()<S3.norm()){
              return 0;
          }
          return 1;
      }
      else {
          if(S1.norm()<S3.norm()){
              return 0;
          }
          return 2;
      }
  }

   Point normal() const {
     Point dir1 = v[1] - v[0];
     Point dir2 = v[2] - v[0];
     Point n  = dir1.cross(dir2);
    float d = n.norm();
     Point ToReturn;
    ToReturn.x=n.x/d;
    ToReturn.y=n.y/d;
    ToReturn.z=n.z/d;
    return  ToReturn;
  }
}Face;

typedef struct RuleCorr{
    vector<int> ClassesComparison;
    vector<int> SignComparison;
    vector<float> ZScoreComparisonBuffer;
    vector<float> ZScoreComparisonMax;
    RuleCorr(){
        vector<int > ClassesComparison;
        vector<int > SignComparison;
        vector<float > ZScoreComparisonBuffer;
        vector<float> ZScoreComparisonMax;
    }
    ~RuleCorr();
}RuleCorr;

typedef struct SimpleRule{
    vector<int> Modalities;
    vector<int> Origins;
    vector<int *> ClassesComparison;
    vector<int *> ComparisonType;
    vector<float *> ZScoreComparison;
    SimpleRule(){
        vector<int> Modalities;
        vector<int> Origins;
        vector<int *> ClassesComparison;
        vector<int *> ComparisonType;
        vector<float *> ZScoreComparison;
    }
    ~SimpleRule(){
        int CCsize=ClassesComparison.size();
        int CTsize=ComparisonType.size();
        int ZCsize=ZScoreComparison.size();
        for(int s=0;s<CCsize;s++){
            if(ClassesComparison[s]!=NULL){
            delete [] ClassesComparison[s];
            ClassesComparison[s]=NULL;
            }
        }
        for(int s=0;s<CTsize;s++){
            if(ComparisonType[s]!=NULL){
            delete [] ComparisonType[s];
            ComparisonType[s]=NULL;
            }
        }
        for(int s=0;s<ZCsize;s++){
            if(ZScoreComparison[s]!=NULL){
            delete [] ZScoreComparison[s];
            ZScoreComparison[s]=NULL;
            }
        }
    }
}SimpleRule;

typedef struct EvaluationReport{
    float DSC;
    float OER;
    float DE;
    float AvDist;
    float VD;
    float TPR;
    float TPR_card;
    float FPR;
    float FPR_card;
    float FNR_card;
    float VolumeRef;
    float VolumeSeg;
    float LabelRef;
    float LabelSeg;
    int FP;
    int FN;
    int TP;
    int DEFP;
    int DEFN;
    int OEFN;
    int OEFP;
    EvaluationReport(){
        VolumeRef=-1;
        VolumeSeg=-1;
        LabelRef=-1;
        LabelSeg=-1;
        DSC=-1;
        OER=-1;
        DE=-1;
        AvDist=-1;
        FNR_card=-1;
        VD=-1;
        TPR=-1;
        TPR_card=-1;
        FPR=-1;
        FPR_card=-1;
        FP=-1;
        FN=-1;
        TP=-1;
        DEFP=-1;
        DEFN=-1;
    }
    
};

typedef struct Rule{
    int * Modalities;
    bool * Acceptance;
    vector<int> Questionable;
    vector<int> QuestionableUniform;
    vector<int> RefinedCheck;
    vector<int> RefinedCheck2;
    vector<int> RefinedCheckUniform;
    vector<int> RefinedCheckUniform2;
    vector< vector <int> > CheckSuspiciousInliers;
    int * SignComparison;
    int * CorrespondingGClassComparison;
    bool * AcceptanceUniform;
    int * SignComparisonUniform;
    int * CorrespondingGClassComparisonUniform;
    vector <int> InliersRule;
    Rule(){
        Modalities=NULL;
        Acceptance=NULL;
        SignComparison=NULL;
        CorrespondingGClassComparison=NULL;
        AcceptanceUniform=NULL;
        SignComparisonUniform=NULL;
        CorrespondingGClassComparisonUniform=NULL;

    }
    ~Rule(){
        if(Modalities!=NULL){
            delete [] Modalities;
            Modalities=NULL;
        }
        if(Acceptance!=NULL){
            delete [] Acceptance;
            Acceptance=NULL;
        }
        if(SignComparison!=NULL){
            delete[] SignComparison;
            SignComparison =NULL;
        }
        if(CorrespondingGClassComparison!=NULL){
            delete [] CorrespondingGClassComparison;
            CorrespondingGClassComparison = NULL;
        }
        if(AcceptanceUniform!=NULL){
            delete [] AcceptanceUniform;
            AcceptanceUniform=NULL;
        }
        if(SignComparisonUniform!=NULL){
            delete[] SignComparisonUniform;
            SignComparisonUniform =NULL;
        }
        if(CorrespondingGClassComparisonUniform!=NULL){
            delete [] CorrespondingGClassComparisonUniform;
            CorrespondingGClassComparisonUniform = NULL;
        }
    }
}Rule;

enum LesionClass {
    NL=0,
    PV=10,
    PVNGM=11,
    PVD=12,
    SC1=13,
    SC2=14,
    SP=-10,
    SCS=6,
    SCS2=-6,
    CS=-2,
    IV=-11,
    TV=-12,
    FV=-13,
    Out=-14,
    RecFP=-15,
    RecTP=15
};

enum OutlierType{
    Iron=-1,
    ECSFRemain=-2,
    NoiseFlat=-3,
    VentricleBorder=-4,
    CorticalSheet=-5,
    ICSF=-6,
    SuspECSF=-7,
    OutLac=-8,
    CerebellumSuspect=-18,
    
    DGMLacune=1,
    WMDGMLacune=2,
    WMLacune=3,
    SuspIron=-9,
    SuspLacune=5,
    Lacune=6,
    Infarct=7,
    
//    Then treating options of WMH
    PeriV=10,
    PeriVNGM=11,
    PeriVD=12,
    SubCort1=13,
    SubCort2=14,
    NonClassified=15,
    JuxtaCorrMS=17,
    SuspCS=16,
    SuspCS2=-16,
    FalseWMH=-17,
    SeptPell=-10,
    IntV=-11,
    ThirdV=-12,
    FourthV=-13,
    OutWMH=-14,
    ReclassifFP=-15,
    ReclassifTP=15,
    ReclassifNonContrast=-20
};

enum LesionSubcategories{
    WMH=30,
    PVWMH_PureWM=31,
    PVWMH_BG=32,
    DWMH_PureWM=33,
    DWMH_Mixed=34,
    
    Lacunes_SCR=1,
    Lacunes_SCE=2,
    Lacunes_BGR=3,
    Lacunes_BGE=4,
    Lacunes_PVR=5,
    Lacunes_PVE=6,
    PVS_BGR=7,
    PVS_BGE=8,
    PVS_SCR=9,
    PVS_SCE=10,
    PVS_E=11,
    PVS_R=12,
    InfarctLesion=40,
    
    PVSLacSusp=13,
    
    IronDeposit=20,
    UndeterminedLac=21,
    UndeterminedWMH=22,
    NonLesion=0
};


typedef struct TextureDescriptors{
    float Contrast;
    float Homogeneity;
    float Entropy;
    float Energy;
    float Correlation;
    TextureDescriptors(){
        Contrast=0;
        Homogeneity=0;
        Entropy=0;
        Energy=0;
        Correlation=0;
    }
};

typedef struct TextureDescriptors1{
    float Mean;
    float Variance;
    float Skewness;
    float Kurtosis;
    TextureDescriptors1(){
        Mean=0;
        Variance=0;
        Skewness=0;
        Kurtosis=0;
    }
};


typedef struct Outlier{
    float Volume;
    float Surface;
    float SAV;
    float MeanOC;
    float MeanOCU;
    vector<float *> MahalVec;
    vector<float *> MahalVecQuant;
    float * Mean;
    float * Variance;
    int * BoundingBox;
    float * SVDEigen;
    float * RatioBoundingBox;
    float * ExtentVolumeBB;
    float * ExtentBB;
    float * ProportionOrigin;
    float * ProportionNeighbour;
    float PropBorderExtent;
    float Compactness;
    float DistanceGrav;
    bool DGMBelonging;
    float SPPotProp;
    float PropArtefact;
    float PropInterface;
    float PropCer;
    float PropCerEro;
    float PropICSF;
    float PropECSF;
    float PropCGM;
    float PropDGM;
    float PropWMHNeigh;
    float ProportionNeighbourArtefact;
    vector<float> VectorDiffGrav;
    LesionSubcategories LesionType;
    int LesionCode;
    vector<int> CentreGravity;
    vector<float> DistanceGMI;
    vector<float> DistanceWMI;
    vector<float> DistanceCSFI;
    vector<float> DistanceOutI;
    vector<float> DistanceSP;
    vector<float> DistanceVentricle;
    vector<TextureDescriptors *> Texture;
    vector<TextureDescriptors1 *> Texture1;
    vector <float *> Edginess;
    OutlierType OutlierClass;
    Outlier(){
        Volume=-1;
        Surface=-1;
        SAV=-1;
        Mean=NULL;
        Variance=NULL;
        BoundingBox=NULL;
        SVDEigen=NULL;
        ExtentBB=NULL;
        RatioBoundingBox=NULL;
        ProportionNeighbour=NULL;
        ExtentVolumeBB=NULL;
        ProportionOrigin=NULL;
        DGMBelonging=0;
        LesionType=NonLesion;
        LesionCode=-2;
        SPPotProp=0.51;
        PropICSF=0;
        PropECSF=0;
        PropDGM=0;
        PropCGM=0;
        PropWMHNeigh=0;
        PropCer=0;
        PropCerEro=0;
        PropArtefact=-1;
        ProportionNeighbourArtefact=-1;
        PropBorderExtent=0;
    }
    ~Outlier(){
        int sizeMahalVec=MahalVec.size();
        if(sizeMahalVec>0){
            for(int s=0;s<sizeMahalVec;s++){
                delete [] MahalVec[s];
                MahalVec[s]=NULL;
            }
        }
        int sizeEdginess=Edginess.size();
        if(sizeEdginess>0){
            for(int s=0;s<sizeEdginess;s++){
                delete [] Edginess[s];
                Edginess[s]=NULL;
            }
        }
        int sizeMahalVecQuant=MahalVecQuant.size();
        if(sizeMahalVecQuant>0){
            for(int s=0;s<sizeMahalVecQuant;s++){
                delete [] MahalVecQuant[s];
                MahalVecQuant[s]=NULL;
            }
        }
        if(Mean!=NULL){
            delete [] Mean;
            Mean=NULL;
        }
        if(Variance!=NULL){
            delete [] Variance;
            Variance=NULL;
        }
        if(BoundingBox!=NULL){
            delete [] BoundingBox;
            BoundingBox=NULL;
        }
        if(SVDEigen!=NULL){
            delete [] SVDEigen;
            SVDEigen=NULL;
        }
        if(ExtentVolumeBB!=NULL){
            delete [] ExtentVolumeBB;
            ExtentVolumeBB=NULL;
        }
        if(ExtentBB!=NULL){
            delete [] ExtentBB;
            ExtentBB=NULL;
        }
        if(RatioBoundingBox!=NULL){
            delete [] RatioBoundingBox;
            RatioBoundingBox=NULL;
        }
        if(ProportionOrigin!=NULL){
            delete [] ProportionOrigin;
            ProportionOrigin=NULL;
        }
        if(ProportionNeighbour!=NULL){
            delete [] ProportionNeighbour;
            ProportionNeighbour=NULL;
        }
        int sizeText=Texture.size();
        if(sizeText>0){
            for(int s=0;s<sizeText;s++){
                if(Texture[s]!=NULL){
                    delete Texture[s];
                    Texture[s]=NULL;
                }
            }
        }
        int sizeText1=Texture1.size();
        if(sizeText>0){
            for(int s=0;s<sizeText1;s++){
                if(Texture1[s]!=NULL){
                    delete Texture1[s];
                    Texture1[s]=NULL;
                }
            }
        }
    }
}Outlier;



typedef struct Lesion{
    float Volume;
    float Surface;
    float SAV;
    float * RatioToNormal;
    float * Mean;
    float * Variance;
    int * BoundingBox;
    float * RatioBoundingBox;
    float * ProportionOrigin;
    float DistanceCSF;
    float DistanceGM;
    float DistanceWM;
    float DistanceOut;
    float ProportionCSF;
    float ProportionGM;
    float ProportionWM;
    float ProportionOut;
    float ProportionGravCSF;
    float ProportionGravGM;
    float ProportionGravWM;
    float ProportionGravOut;
    float Compactness;
    float DistanceGrav;
    float DistanceVentricle;
    bool NeighbourWMI;
    bool NeighbourGMC;
    bool NeighbourCSF;
    bool NeighbourOut;
    bool NeighbourVentricle;
    bool DGMBelonging;
    float PropDGM;
    float SPPotProp;
    float PropICSF;
    vector<int> VectorDiffGrav;
    float DistanceAcceptedGM;
    LesionClass LesionType;
    int LesionCode;
    vector<int> CentreGravity;
    vector<float> DistanceGMI;
    vector<float> DistanceWMI;
    vector<float> DistanceCSFI;
    vector<float> DistanceOutI;
    vector<float> DistanceSP;
    vector<TextureDescriptors *> Texture;
    vector<TextureDescriptors1 *> Texture1;
    Lesion(){
        Volume=-1;
        Surface=-1;
        SAV=-1;
        RatioToNormal=NULL;
        Mean=NULL;
        Variance=NULL;
        BoundingBox=NULL;
        RatioBoundingBox=NULL;
        ProportionOrigin=NULL;
        DistanceCSF=-1;
        DistanceGM=-1;
        DistanceWM=-1;
        DistanceOut=-1;
        DistanceVentricle=-1;
        DistanceAcceptedGM=-1;
        DGMBelonging=0;
        NeighbourWMI=0;
        NeighbourGMC=0;
        NeighbourCSF=0;
        NeighbourVentricle=0;
        LesionType=Out;
        LesionCode=-2;
        SPPotProp=0.51;
        PropICSF=0;
    }
    ~Lesion(){
        if(RatioToNormal!=NULL){
            delete [] RatioToNormal;
            RatioToNormal=NULL;
        }
        if(Mean!=NULL){
            delete [] Mean;
            Mean=NULL;
        }
        if(Variance!=NULL){
            delete [] Variance;
            Variance=NULL;
        }
        if(BoundingBox!=NULL){
            delete [] BoundingBox;
            BoundingBox=NULL;
        }
        if(RatioBoundingBox!=NULL){
            delete [] RatioBoundingBox;
            RatioBoundingBox=NULL;
        }
        if(ProportionOrigin!=NULL){
            delete [] ProportionOrigin;
            ProportionOrigin=NULL;
        }
        int sizeText=Texture.size();
        if(sizeText>0){
            for(int s=0;s<sizeText;s++){
                if(Texture[s]!=NULL){
                    delete Texture[s];
                    Texture[s]=NULL;
                }
            }
        }
        int sizeText1=Texture1.size();
        if(sizeText>0){
            for(int s=0;s<sizeText1;s++){
                if(Texture1[s]!=NULL){
                    delete Texture1[s];
                    Texture1[s]=NULL;
                }
            }
        }
    }
}Lesion;

enum ResultType {
    TP = 0,
    TN=1,
    FP = 2,
    FN=3,
    OE=4,
    OEFN=5,
    OEFP=6,
    DEFP=7,
    DEFN=8
};



typedef struct LesionSimple{
    float Volume;
    float Surface;
    float SAV;
}LesionSimple;

typedef struct LesionIntensity{
    float Volume;
    vector<float> Param;
    vector<float> MeanRatio;
    vector<float> MahalDistance;
    vector<float> ExtremaMahal;
    vector<float> Quantilisation;
    vector<float> ProportionMahal;
}LesionIntensity;

typedef struct LesionRefIntensity{
    float Volume;
    float TrueWorkingVolume;
    vector<int> AssociatedLabel;
    float * TPMahal;
    float * FNMahal;
    float * FPMahal;
    float * SegMahal;
    float * RefMahal;
    float * ProportionCloseVentricle;

    float * PropMahalFP;
    float * PropMahalFN;
    float * PropMahalTP;
    float * PropMahalSeg;
    float * PropMahalRef;
    int TPl;
    int FPl;
    int FNl;
    int Segl;
    int Refl;
    float DSC;
    int FPNotAgree;
    int FNNotAgree;
    int TPNotAgree;
    int SegNotAgree;
    int RefNotAgree;
    ResultType FlagTypeConnectedComponent;
    LesionRefIntensity(){
        FlagTypeConnectedComponent=TN;
        Volume=-1;
        TrueWorkingVolume=-1;
        TPl=0;
        FNl=0;
        FPl=0;
        Refl=0;
        Segl=0;
        FPNotAgree=0;
        TPNotAgree=0;
        FNNotAgree=0;
        SegNotAgree=0;
        RefNotAgree=0;
        DSC=-1;
        ProportionCloseVentricle=new float[5];
        FPMahal=new float[6];
        FNMahal=new float[6];
        TPMahal=new float[6];
        SegMahal=new float [6];
        RefMahal=new float[6];
        PropMahalFP=new float[10];
        PropMahalFN=new float[10];
        PropMahalTP=new float[10];
        PropMahalSeg=new float[10];
        PropMahalRef=new float[10];
        for (int i=0; i<3; i++) {
            FPMahal[i]=pow(-1, i+1)*1000*i;
            FNMahal[i]=pow(-1, i+1)*1000*i;
            TPMahal[i]=pow(-1, i+1)*1000*i;
            SegMahal[i]=pow(-1, i+1)*1000*i;
            RefMahal[i]=pow(-1, i+1)*1000*i;
        }
        for (int i=3;i<6;i++){
            FPMahal[i]=0;
            FNMahal[i]=0;
            TPMahal[i]=0;
            SegMahal[i]=0;
            RefMahal[i]=0;
        }
        for(int i=0;i<5;i++){
            ProportionCloseVentricle[i]=0;
        }
    }
    ~LesionRefIntensity(){
        if(TPMahal!=NULL){
            delete [] TPMahal;
            TPMahal=NULL;
        }
        if(FNMahal!=NULL){
            delete [] FNMahal;
            FNMahal=NULL;
        }
        if(FPMahal!=NULL){
            delete [] FPMahal;
            FPMahal=NULL;
        }
        if(SegMahal!=NULL){
            delete [] SegMahal;
            SegMahal=NULL;
        }
        if(RefMahal!=NULL){
            delete [] RefMahal;
            RefMahal=NULL;
        }
        if(ProportionCloseVentricle!=NULL){
            delete [] ProportionCloseVentricle;
            ProportionCloseVentricle=NULL;
        }
        if(PropMahalFN!=NULL){
            delete [] PropMahalFN;
            PropMahalFN=NULL;
        }
        if(PropMahalFP!=NULL){
            delete [] PropMahalFP;
            PropMahalFP=NULL;
        }
        if(PropMahalTP!=NULL){
            delete [] PropMahalTP;
            PropMahalTP=NULL;
        }
        if(PropMahalSeg!=NULL){
            delete [] PropMahalSeg;
            PropMahalSeg=NULL;
        }
        if(PropMahalRef!=NULL){
            delete [] PropMahalRef;
            PropMahalRef=NULL;
        }
    }

}LesionRefIntensity;


typedef struct LeafID{
    vector<int> HierarchyVector;
    float Threshold;
    float RatioInside;
    float RatioGM;
    float RatioWM;
    float RatioOut;
    float RatioCSF;
    float NeighInside;
    float NeighCSF;
    float NeighGM;
    float NeighWM;
    float NeighOut;
    float PropInsideContribution;
    LeafID(){
        Threshold=0;
        RatioInside=0;
        RatioGM=0;
        RatioWM=0;
        RatioOut=0;
        RatioCSF=0;
        NeighInside=0;
        NeighCSF=0;
        NeighGM=0;
        NeighWM=0;
        NeighOut=0;
        PropInsideContribution=0;
    }
    ~LeafID(){
    }
}LeafID;

float * CreateVentriclesFromParc(float * ParcData, int numel);
float * CreateHemisphereFromParc(float * ParcData, int Side, int numel);
nifti_image * CombineVectorNifti(vector<nifti_image *> VectorSegCombine, float threshold);
nifti_image * HardSegmentationTemp(TreeEM * TreeToAnalyse, bool * BoolOppSeg, int Option, SEG_ANALYSIS * segment_analysis);
nifti_image * HardSegmentationCompare(nifti_image * SegToAnalyse, nifti_image * SummarisedSeg, int TypeHard, nifti_image * Mask);
nifti_image * HardSegmentationThresholdFromNormResp(float * NormResp, TreeEM* TreeToAnalyse, float Threshold);
nifti_image * HardSegmentationThreshold(nifti_image * SegToAnalyse,float Threshold,int Option=1);
nifti_image * WMIHardSegmentation(nifti_image * SummarisedSeg1, nifti_image * SegLesBasis,SEG_ANALYSIS * segment_analysis,int * L2SToUse);
nifti_image * HardSegmentationIndex(nifti_image * SoftSegmentation, int IndexHard, int * L2S);
nifti_image * HardSegmentationLesionReconstruct_bis(nifti_image * SegToAnalyse, nifti_image * SummarisedSeg, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
nifti_image * SoftSegmentationLesionReconstruct(nifti_image * SegToAnalyse, TreeEM * TreeToAnalyse);
nifti_image * HardSegmentationLesionReconstruct(nifti_image * SegToAnalyse, TreeEM * TreeToAnalyse);
nifti_image * HardSegmentation_bis(nifti_image * SoftSegmentation, nifti_image * Mask);
nifti_image* HardSegmentation(nifti_image * SoftSegmentation);
nifti_image * NeighbourhingWeight(nifti_image * SoftSeg,int ExploRadius);
vector<nifti_image*> DistanceFromClassesHS(nifti_image * SoftSeg1, nifti_image * SoftSeg2,int ExploRadius);
vector<int> IndicesToExplore(nifti_image* HS1,int CentralIndex, int ExploRadius);
vector<int> ValidIndicesFromExplo(nifti_image * HS1, int ClassToCheck, vector<int> IndicesToExplore);
float MinimumDistance(nifti_image * HS1,int CentralIndex,vector<int>ValidIndices);

float * HardDiceScores(nifti_image * SoftSegmentation1, nifti_image* SoftSegmentation2);
float * SoftDiceScores(nifti_image * SoftSegmentation1, nifti_image* SoftSegmentation2);

bool * SegBoolArray(nifti_image * BoolSeg, nifti_image * Mask);
bool * SegBoolArrayClass(nifti_image * BoolSeg, int c, nifti_image * Mask);
float * SegFloatArray(nifti_image * FloatSeg,  nifti_image * Mask);
float * SegFloatArrayClass(nifti_image * FloatSeg, int c, nifti_image * Mask);
float Correlation(nifti_image * ImageToCorrelate1, nifti_image * ImageToCorrelate2,nifti_image * Mask);

float ProbabilitySimilarityIndex(float * SoftSeg, bool * BinSegRef, int numelmasked);
float Sensitivity(bool * HardClass1, bool * HardClass2, int numelmasked);
float Specificity(bool * HardClass1, bool * HardClass2, int numelmasked);
float Accuracy(bool * HardClass1, bool * HardClass2, int numelmasked);
float FPR(bool * HardClass1, bool * HardClass2, int numelmasked);
float TPR(bool * HardClass1, bool * HardClass2, int numelmasked);
float DSC(bool * HardClass1, bool * HardClass2, int numelmasked);
vector<float> RateVector_card(int * LabelsSeg, int * LabelsRef, int numel);
EvaluationReport * CreateEvaluationReport(nifti_image * LesNii, nifti_image * RefNii, nifti_image * MaskNii);
float VD(bool * HardClass1, bool * HardClass2, int numelmasked);
float AverageDistanceMetric(bool * BoolSeg, bool * BoolRef, nifti_image * ImageRef);
float TLL(float * SoftLes, int numelmasked);
float TLL(bool * HardLes, int numelmasked);

int TruePositives_card(int * LabelsSeg, int * LabelsRef, int numel);
int TruePositives_bis(bool * HardClass1, bool * HardClass2, int numelmasked );
float * TruePositives(nifti_image * SoftSegmentation1, nifti_image* SoftSegmentation2);
ResultType * ClassifLabel(int * LabelSeg, int * LabelRef, int numel);
vector<LesionRefIntensity * > AnalysisPerRefConnectedIntensity_bis(int * LabelSeg, int * LabelRef, nifti_image * MahalDistMaps, SEG_ANALYSIS * segment_analysis,nifti_image * DistanceVentricleMap);
vector<LesionRefIntensity * > AnalysisPerRefConnectedIntensity(int * LabelSeg, int * LabelRef, nifti_image * MahalDistMaps, SEG_ANALYSIS * segment_analysis);
int * AnalysisPerRefConnected(int * LabelSeg, int * LabelRef, int numel);
vector<int>GetListCorrespondingLabels(int Label,int * LabelSeg,int * LabelRef,int numel);
int DetectionError(int * LabelSeg, int * LabelRef, int numel);
int DetectionErrorFP(int * LabelSeg, int * LabelRef, int numel);
int DetectionErrorFN(int * LabelSeg, int * LabelRef, int numel);
float GetThresholdForVolume(float * BlurredLesion,bool * MaskData, int  NumberNeeded,int numel);
nifti_image * ReductionLesion(nifti_image * LesNii, nifti_image * MaskNii, float Threshold, int NumberOut);
nifti_image * ReductionLesionVector(nifti_image * LesNii, nifti_image * MaskNii, vector<float> vecLesionReduction);
float * VolumeRange(nifti_image * RefNii,vector<float> ThreshVec,nifti_image *MaskNii);
int OutlineError(int * LabelSeg, int * LabelRef, int numel);
int FalsePositives_card(int * LabelsSeg, int * LabelsRef, int numel);
int FalsePositives_bis(bool * HardClass1, bool * HardClass2, int numelmasked );
float * FalsePositives(nifti_image * SoftSegmentation1, nifti_image* SoftSegmentation2);
int TrueNegatives_bis(bool * HardClass1, bool * HardClass2, int numelmasked );
float * TrueNegatives(nifti_image * SoftSegmentation1, nifti_image* SoftSegmentation2);
int FalseNegatives_bis(bool * HardClass1, bool * HardClass2, int numelmasked );
float * FalseNegatives(nifti_image * SoftSegmentation1, nifti_image* SoftSegmentation2);

int NumbActiveVoxels(nifti_image * Image);
int NumbActiveVoxels_bis(nifti_image * Mask);

nifti_image * PrepareCompareSoftSegmentations(vector<nifti_image*> SegmentationsToJoin, nifti_image * Mask);
nifti_image * PrepareCompareSoftSegmentationOutlierMult(vector<nifti_image*> PartialSeg,nifti_image *Mask);
void PrepareCompareSoftSegmentationOutlierJoint(nifti_image* SegToCompare,nifti_image *Mask);

nifti_image * Binarisation(nifti_image * Mask);
int CalculateNumberMaskedElements(nifti_image * Mask);
bool CheckCompatibleDimensions(nifti_image * Image1, nifti_image * Image2);
bool CheckCompatibilityVectorNii(vector<nifti_image*> VectorNii);
nifti_image * Floatisation(nifti_image * Image);
void MaskImage(nifti_image * ImageToMask, nifti_image * Mask);
void MakeSoftSegmentation(nifti_image* SoftSegmentation);

// Quantitative Dice Scores

// Tests on images
bool IsResultSegmentation(nifti_image * SoftSegmentation);
nifti_image * MakeResultSegmentation(nifti_image * SoftSegmentation);


//KLD, BF and related
nifti_image * NormaliseImage(nifti_image * DataImage, nifti_image * Mask);
vector<float *> MakeDataHistogram(nifti_image* DataImage, nifti_image * Mask, nifti_image * Outliers);
vector<float *> MakeDistHistogram(int numbmodal,vector<float *> ParametersVector);
float * MakeDistHistogramTotal(int numbmodal,vector<float*> ParametersVector);
float * MakeGaussianDistributionHist(int numbmodal, float * Parameters,int modal);
float * MakeGaussianDistributionHistTotal(int numbmodal, float *Parameters);
float * MakeDataHistogramTotal(nifti_image * DataImage, nifti_image* Mask, nifti_image * Outliers);
float MakeKLDTot(nifti_image * DataImage, nifti_image* Mask,nifti_image* Seg, nifti_image*Outliers);
vector<float> MakeKLDVector(nifti_image * DataImage, nifti_image* Mask,nifti_image* Seg, nifti_image*Outliers);
vector<nifti_image *> ReadFromFilenamesVector(vector<string> FilenamesVector);
nifti_image * CreateBFCorrectedImage(nifti_image * DataImage, nifti_image* BFImage);
nifti_image * CreateDataImage(vector<nifti_image*> ImagesToSegment);

// Reconstruction from text file and total seg

nifti_image * ResultImageFromSeg(nifti_image * BinarySeg, nifti_image * RefSeg, nifti_image * Mask, ResultType RT);
vector<TreeEM *> InlierSuspiciousClasses(Rule * LesionRule, TreeEM * TreeToAnalyse, vector<int> Modalities, nifti_image * LesionImageTot );
nifti_image * SuspiciousImage(Rule * LesionRule, TreeEM * TreeToAnalyse, vector<int> Modalities, nifti_image * LesionImageTot);
vector<TreeEM *> FindLesionClasses(TreeEM * TreeToAnalyse,Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
vector<TreeEM*> FindLesionClassesInliersFLAIRonly(TreeEM * TreeToAnalyse, vector<int> Modalities, SEG_ANALYSIS * segment_analysis);
vector<TreeEM *> FindLesionClassesInliers(TreeEM * TreeToAnalyse,Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
vector <TreeEM *> FindClassesToInlier(TreeEM * TreeToAnalyse, vector<int> Modalities, vector<int> Rule);
vector<TreeEM *> FindCSFWMOutliers(TreeEM * TreeToAnalyse, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
nifti_image * FindCSFWMOutliersUniform(TreeEM * TreeToAnalyse, vector<int> Modalities, SEG_ANALYSIS * segment_analysis,nifti_image * CorrectionJuxta=NULL);
float * GetMeanLesion(nifti_image * LesionImageTot, TreeEM * TreeToAnalyse);
vector<int> GetModalitiesFromTextFile(SEG_ANALYSIS * segment_analysis);
int * GetCorrespondingModality( vector<int> Modalities, vector<int> CheckOfModality);
float * GetMeanDirectVector(vector<TreeEM *> LeavesTrulyLesion);
bool ToCheckLesionLeave(TreeEM * LesionClassToCheck, Rule * LesionRule,vector<float*> MeanGeneralClasses,vector<int> Modalities);
bool RefinedChecking(TreeEM * LeavesToCheck,Rule * LesionRule,float * MeanTrueLesions,TreeEM* TreeToAnalyse,vector<TreeEM *>LeavesTrulyLesions,vector<int> Modalities,SEG_ANALYSIS *segment_analysis);
void RebuildNormRespFromLeaves(TreeEM * TreeToRebuild);
float * SumNormRespChildren(TreeEM * TreeToRebuild);
void NullAllNonLeavesNormResp(TreeEM * TreeToRebuild);
float * CompareMeanLesion(vector<TreeEM *> LesionClasses, vector<int> CheckVector, TreeEM * TreeToAnalyse, vector<int> Modalities);
float * GetMeanToAssess(vector<TreeEM *> LeavesToConsiderVector, int modality);
float ExtremaFromArray(float * ArrayToAssess, int size, int extremaType);
bool CheckLeavesLesionRule(TreeEM * Leaf,Rule * LesionRule,vector<float*>MeanGeneralClasses);
nifti_image * PVSExtraction(RuleCorr * Rule, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis, vector<float *> MahalVec);
float * ProbaCorrWeight(TreeEM * TreeToAnalyse, RuleCorr* RuleCorrection, nifti_image* SegToAnalyse, vector<float*> MahalVec, SEG_ANALYSIS* segment_analysis);
int * ImproveCorrWeight(TreeEM * TreeToAnalyse, RuleCorr* RuleCorrection, nifti_image * SegToAnalyse, vector<float *> MahalVec, SEG_ANALYSIS * segment_analysis);
int * CreateCorrWeight(TreeEM * TreeToAnalyse, RuleCorr * RuleCorrection, nifti_image * SegToAnalyse, vector<float*> MahalVec, SEG_ANALYSIS * segment_analysis);
RuleCorr * BuildRuleCorrFromTextFile(SEG_ANALYSIS * segment_analysis);
SimpleRule * BuildLCRuleFromTextFile(SEG_ANALYSIS * segment_analysis);
Rule * BuildRuleFromTextFile(TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis,bool Type=1);
Rule * BuildingRule(TreeEM * TreeToAnalyse, int RuleType, int UniformRuleType, vector<int> Modalities);
nifti_image * ReconstructLesionImage(TreeEM * TreeToAnalyse,  Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
nifti_image * ReconstructLesionImageTot(TreeEM * TreeToAnalyse,  Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
nifti_image * ReconstructLesionImageTotWeighted(TreeEM *TreeToAnalyse, vector<TreeEM*> LesionClasses, nifti_image * LesionPartsFromUniform,vector<int>Modalities, SEG_ANALYSIS * segment_analysis,bool secondary=0,nifti_image* CorrectionJuxta=NULL);
float * GetWeightLesionClasses(vector<TreeEM *>LesionClasses, TreeEM * TreeToAnalyse,SEG_ANALYSIS *segment_analysis);
float GetExtremaMahalDistLabel(int * OrderedLabel, int Label, float * MeanSeg, float * InvertedCovariance, nifti_image* DataComp,int OptionExtrema);
float GetMahalDist(float * ValueArray, float * MeanArray, float * InvertVarianceArray, int numbmodal);
nifti_image * LesionFromUniform(TreeEM * TreeToAnalyse, Rule * LesionRule, vector<int> Modalities,vector<TreeEM *> LesionClasses,SEG_ANALYSIS * segment_analysis);
nifti_image * LesionFromTotOutliers(TreeEM * TreeToAnalyse, Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
nifti_image * ExtractParenchyma(nifti_image * SummarisedSeg, nifti_image * VentricleSeg, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
//nifti_image * ExtractWMDGM(nifti_image * SummarisedSeg,nifti_image * VentricleSegmentation, SEG_ANALYSIS * segment_analysis, TreeEM * TreeToAnalyse);
LeafID * LeafCharacterisation(TreeEM * TreeToAnalyse, TreeEM * LeafStudy, SEG_ANALYSIS * segment_analysis, nifti_image * SummarisedSeg,float Threshold);
nifti_image * DGMSegmentationPriors(nifti_image * SummarisedSeg, SEG_ANALYSIS * segment_analysis, nifti_image * DGMPriors,int * L2S);
bool * VentricleSegmentationFromParc(TreeEM* TreeToAnalyse,SEG_ANALYSIS * segment_analysis);
void CreateInterfaceGMCSF(nifti_image * GIFParcellation, bool * Interface);
bool * VentricleSegmentationDirect(TreeEM * TreeToAnalyse,  SEG_ANALYSIS * segment_analysis);
nifti_image * VentricleSegmentationPriors(nifti_image * SummarisedSeg, SEG_ANALYSIS * segment_analysis, TreeEM * TreeToAnalyse, nifti_image * ICSFPriors);
nifti_image * VentricleSegmentation(nifti_image * SummarisedSeg,SEG_ANALYSIS * segment_analysis, TreeEM * TreeToAnalyse);
//nifti_image * VentricleSegmentation_bis(nifti_image * SummarisedSeg, TreeEM * TreeToAnalyse, SEG_ANALYSIS *segment_analysis);
nifti_image * SummarizedSegmentation(TreeEM * TreeToAnalyse, vector<TreeEM *> LesionClasses, nifti_image * LesionUniform, vector<TreeEM*> CSFWMClasses, nifti_image* CSFWMUniform,SEG_ANALYSIS * segment_analysis,float * ICSF=NULL,nifti_image * TrueLesionSeg=NULL);
void WriteWMLIdentityCard(TreeEM * TreeToAnalyse, Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
float * ProportionIndicesSegmentation( vector<int> IndicesPath, nifti_image * HardSegmentation);
float * ProportionFromBoolVec(bool * BoolToAssess, vector<bool *> VectorBool,int numel);
float * ProportionNeighboursLesion(int * LesionLabeling, int Label, nifti_image *SummarisedDecoupled,SEG_ANALYSIS * segment_analysis);
nifti_image * EuclideanDistanceMap(nifti_image * ReferenceSeg, nifti_image * MaskSeg);
nifti_image * MahalDistMaps(nifti_image * TissueComp, nifti_image * SegToAnalyse, nifti_image * DataToCompare);
nifti_image * SeverityMap(TreeEM * TreeToAnalyse,nifti_image * SegToAnalyse, vector<int> Modalities, SEG_ANALYSIS * segment_analysis);
int * ComponentLabelingNotRefined(nifti_image * SegToAnalyse, int Neigh, float thresh=0.5);
int * ComponentLabeling(nifti_image * SegToAnalyse, int Neigh, float thresh=0.5);
int * ComponentLabeling(float * ImageToLabel, int Neigh, int * Dim, int * Shift, float thresh=0.5);
void RepassComponentsLabel(int * ComponentsLabel,int Neigh,int * Dim,int *Shift);
void GetListNeighbours_bis(int * ListNeighbours, int CurrentIndex, int * Dim, int * Shift,int NeighbouringType);
int * GetListNeighboursToLook(int CurrentIndex, int * Dim,int * Shift);
int * GetListNeighbours(int CurrentIndex, int * Dim, int * Shift, int NeighbourhingType);
void RecurseCompLabel(float * Data, int * PositionTmp, int Label, int * ComponentsLabel,int Neigh,int * Dim,int * Shift,float thresh,int& Times);
void RelabelComponents(int * ComponentsLabel, int * Dim);
nifti_image * ImageComponentLabeling(nifti_image * SegToAnalyse, SEG_ANALYSIS * segment_analysis);
void CorrectionCST(TreeEM* TreeToAnalyse, nifti_image * SegToAnalyse, SEG_ANALYSIS* segment_analysis);
int GetLabelSize(int Label, int * ComponentLabel, int numel);
int * RefinedComponentLabel(int * ComponentLabel, int MiniSize, int numel, int * CorrespondanceTrueLabel, int maxLabel);
int * OrderedVolumeLabel(int * ComponentLabel, float MiniSize, int numel,float VolumeVox=1,float MaxSize=1000000);
nifti_image * ImageRefinedLabel(nifti_image * SegToAnalyse, SEG_ANALYSIS * segment_analysis);
bool IsLesionAcceptable(int Label, int * ComponentLabel, int MiniSize, int numel);
int * GetVolumeLabels(int * ComponentLabel, int numel);
int GetLabelSize( int Label, int * ComponentLabel, int numel);
bool IsLabelPresent(int Label, int * ComponentLabel,int numel);
int * GetCorrespondanceOrderedVolume(int * VolumeLabels, int maxLabel);
int GetMaxLabel(int * ComponentLabel, int numel);
nifti_image * CopyFloatNii(nifti_image * FloatNiiToCopy);
nifti_image * MakeMeanFromVectorImages(vector<nifti_image *> VectorImagesToMean,bool * Mask);
nifti_image * CorrectionLesionFromClassif(nifti_image * ImageLesionClassif,nifti_image * SegToAnalyse);
nifti_image * ImageOrderedVolumeLabel(nifti_image * SegToAnalyse, SEG_ANALYSIS * segment_analysis);
nifti_image * LesionVoxelwise(TreeEM * TreeToAnalyse, Rule * LesionRule, vector<int> Modalities,SEG_ANALYSIS * segment_analysis);
nifti_image * LesionCorrectedTrueSeg_bis(nifti_image * LesionInitialSeg, nifti_image * SummarisedSeg, TreeEM * TreeToAnalyse, bool TypeSeg, SEG_ANALYSIS * segment_analysis,nifti_image *VentricleSeg );
nifti_image * LesionReconstructTrueSegOutlier_ter(int * OrderedLabelsLesion, nifti_image * LesionInitialSeg, nifti_image * SummarisedSeg, vector<Outlier *> VectorOutlier, TreeEM * TreeToAnalyse, int TypeSeg,SEG_ANALYSIS * segment_analysis);
nifti_image * LesionReconstructTrueSeg_ter(int * OrderedLabelsLesion, nifti_image * LesionInitialSeg, nifti_image * SummarisedSeg, vector<Lesion *> VectorLesion, TreeEM * TreeToAnalyse, int TypeSeg,SEG_ANALYSIS * segment_analysis);
nifti_image * LesionCorrectedTrueSeg(nifti_image * LesionInitialSeg, TreeEM * TreeToAnalyse, bool TypeSeg, SEG_ANALYSIS * segment_analysis );
bool * CreateObjectsFromParcelLesSeg(nifti_image * ParcellationImage, nifti_image * LesionTot, vector<int> LabelVector, int LesionHandling);
vector<int> GetIndicesBorderLesionExt (int * OrderedLabels, int Label, int * Dim, int * Shift);
vector<int> GetIndicesToExploreFloatDist(int CentralIndex, float MaxDist, int * Dim, int * Shift, float * PixDim, bool flag_plane);
vector<int> GetIndicesToExplore(int BorderLesionIndex, int * OrderedLabels, int Radius, int * Dim);
bool * CreateLesionBool(int * LesionLabeling, int Label,int numel );

/**************** FUNCTIONS TO POPULATE THE LESION STRUCTURE ***********************/
float * QuantiliseMask(bool * BinToQuantilise, nifti_image * MahalDistImage);
float * GetParamFromImages(nifti_image * Seg,nifti_image * Data,int Option);
float * GetParamLabel(int * OrderedLabels, nifti_image * DataCompImage, int Label);
float * GetParamLabel(int * OrderedLabels, TreeEM * TreeToAnalyse, int Label);
vector<int> GetIndicesMinPath(int IndexBegin, int IndexEnd, int * Dim,int * Shift ,float * PixDim);
vector<int> GetIndicesBorderLesion (int * OrderedLabels, int Label, int * Dim, int * Shift);
LesionIntensity * BuildLesionIntensityFromConnectedLabel(int * OrderedLabels, int Label, float * MeanWMSeg, float * CovarianceWMSeg, nifti_image * DataCompImage,nifti_image * MahalDistImage);
LesionSimple * BuildLesionSimpleFromConnectedLabels(int * OrderedLabels, nifti_image * InitialSeg, int Label,int BorderType);
Lesion * BuildLesionFromConnectedLabels(int * OrderedLabels, TreeEM * TreeToAnalyse, int Label, SEG_ANALYSIS * segment_analysis);
Lesion * BuildLesionFromConnectedLabels_bis(int * OrderedLabels, TreeEM * TreeToAnalyse, bool * HardSegSummarised, int Label,SEG_ANALYSIS * segment_analysis);
int * CreateInsideFromBorder(bool * Ring, bool* Centre, float * DistCentre, int* Dim, int* Shift, float* PixDim, bool flag_plane);
float GetDistanceSegToPoint(int Index, bool* Seg, int* Dim, int*Shift, float* PixDim);
float GetMinDistanceFromIndicesVectors(vector<int> IndicesSource, vector<int> IndicesReached,TreeEM * TreeToAnalyse,int * Dim,int * Shift);
float * ProportionCompareGravity(bool * SegStudied, int IndexCentreGrav, int * Dim);
void CorrectedCentreGrav(bool* Mask, int * CentreGravCoord, int * Dim, int * Shift, float * PixDim);
int GetCenterGravity(bool * Mask, int * Dim);
float GetDistanceLabelToPoint(int * OrderedLabels, int Label, TreeEM* TreeToAnalyse, int IndexPoint);
float GetDistanceBetweenSeg(bool * Seg1, bool * Seg2, int * Dim, int * Shift, float * PixDim);
float GetDistanceToSeg(int * OrderedLabels, int Label, TreeEM * TreeToAnalyse, bool * HardSeg);
float GetDistanceToNormalTissue(int * OrderedLabels, int Label, TreeEM * TreeToAnalyse, int TissueType, SEG_ANALYSIS * segment_analysis);
float GetSurface(int * OrderedLabels,TreeEM * TreeToAnalyse,int Label,int BorderType);
float GetSurfaceImage(int * OrderedLabel, nifti_image * ImageRef, int Label, int BorderType);
int * Get6Neighbours(int CurrentIndex, int * Dim, int * Shift);
float * PrintData(nifti_image * DataImage,int Modal, int * NBBToPrint,int DirectionAvoid);
void PrintThreshRange(float * Results, vector<float> VecThresh, ostream & TxtFile);
void PrintEvalReport(EvaluationReport * Eval, ostream & TxtFile);
void PrintOutlierCharacteristics(Outlier * PrintedOutlier,ostream& TxtFileLesion,int numbmodal, int * GravImage, float * PixDim, SEG_ANALYSIS * segment_analysis);
void PrintAnalysisPerRefIntensityConnectedMahalProp(vector<LesionRefIntensity *>VectorLesionIntensity, int Index, ostream & TxtFileLesionIntensity);
void PrintAnalysisPerRefIntensityConnected(vector<LesionRefIntensity *>VectorLesionIntensity, ostream & TxtFileLesionIntensity);
void PrintLesionForFile(Lesion * PrintedLesion, ostream& TxtFileLesion,int numbmodal, int * GravImage, float * PixDim);
void PrintLesionRefIntensityMahalProp(LesionRefIntensity * LesionToPrint, int Index, ostream & TxtFileLesion);
void PrintOutlierCharacteristicsForFile(Outlier * PrintedOutlier,ostream& TxtFileLesion,int numbmodal, int * GravImage, float * PixDim, SEG_ANALYSIS * segment_analysis);
void PrintLesionRefIntensity(LesionRefIntensity * LesionToPrint, ostream & TxtFileLesion);
void PrintAnalysisPerRefConnected(int * ResultAnalysisRefConnected, int MaxLabel,ostream & TxtResRefConnected);
void PrintCountVolperType(int * CountperType,float * VolumeperType,ostream& TxtFileLesion);
void PrintNormDistResult(float * ProportionLesVolDist, int * ExtentLesion, float * LayerIntensity, int MaxLabel, int numbmodal, ostream & TxtNormDistFile,SEG_ANALYSIS* segment_analysis);
void PrintLesionIntensity(LesionIntensity * LesionToPrint, ostream & TxtFileLesion);
void PrintLesion(Lesion * PrintedLesion, ostream& TxtFileLesion,int numbmodal, int * GravImage, float * PixDim);
void PrintLesionSimple(LesionSimple * PrintedLesion, ostream& TxtFileLesion);
void PrintLeafID(LeafID * LeafToPrint, ostream& TxtFileLeaf,int numbmodal);
float RadiusLesionCoalescence(int * OrderedLabels, int Label,bool * AuthorisedRegion, nifti_image * BasisImage, float & Radius, int & MetLabel,int maxRadius);
int * CountLesionNormDistance(nifti_image * LesionLabelCorr, float * LengthNormSol, int NumberLaminae, float * PixDim, int * L2S, int * S2L, int numelmasked);
float * LayersLesIntensityRelation(nifti_image * LesionSeg,float * LengthNormSol, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);

nifti_image * SummarisedSegDecoupled(nifti_image * SummarisedSeg, nifti_image * SegToAnalyse, int CorrespondingDecouplingIndex);
vector<LesionIntensity *> GetVectorLesionIntensity(int * OrderedLabels, float* ParamWM,float * InvertCovarianceWM,nifti_image * DataCompImage, nifti_image * MahalDistImage);
vector<LesionSimple *> GetVectorLesionSimple(int * OrderedLabels,nifti_image * InitialSeg);
vector<Lesion *> GetVectorLesion_bis(int * OrderedLabels,nifti_image * SummarisedSeg, TreeEM * TreeToAnalyse,SEG_ANALYSIS* segment_analysis,nifti_image * VentricleSeg);
vector<Lesion *> GetVectorLesion_ter(int * OrderedLabels,nifti_image * SummarisedSeg, TreeEM * TreeToAnalyse,SEG_ANALYSIS* segment_analysis,bool * VentricleBool, bool * DGMBool, bool * AcceptedGMBool,float * ICSF, float * DGM);
vector<Lesion *> GetVectorLesion_quat(int * OrderedLabels,nifti_image * SummarisedSeg, TreeEM * TreeToAnalyse,SEG_ANALYSIS* segment_analysis,bool * VentricleBool, bool * DGMBool, bool * AcceptedGMBool, bool * SPRegion, float * ICSFPriors, float * DGMPriors,int * OrderedLabelsNIV, nifti_image * SegToAnalyse, nifti_image * SegToAnalyseNIV);
void CorrectOrderedLabels(int LabelBase, int LabelChange, int * OrderedLabelBase, int * OrderedLabelChange, nifti_image * SegBase, nifti_image * SegChange);

vector<vector<Outlier*> > GetVectorOutliersMultiType(float* MultiDataValue, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis, float mult_modify=1);
vector<Outlier *> GetVectorOutliers_NIVCorr(int * OrderedLabels, int * OrderedLabelsNIV, nifti_image * SegToAnalyse, nifti_image * SegToAnalyseNIV,TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis, float mult_modify=1);
nifti_image * OutlierSelected(SimpleRule * Rule, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis );
vector<Outlier *> GetVectorOutliers(int * OrderedLabels,TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis, float mult_modify=1);
vector<Lesion *> GetVectorLesion(int * OrderedLabels,TreeEM * TreeToAnalyse,SEG_ANALYSIS* segment_analysis);
bool IsCoherentEnough(vector<int> IndicesLabel,int * Dim, int * Shift);
vector<int> GetListIndicesLabel(int * ComponentLabel, int * Dim, int Label);
void RefineNeighbourhood(int * ComponentLabel, int * Dim, int * Shift);
bool * BorderExtraction_bis(TreeEM * TreeToAnalyse, nifti_image * SummarisedSeg, SEG_ANALYSIS * segment_analysis, int TissueType);
bool * BorderExtractionThreshold(TreeEM * TreeToAnalyse, nifti_image * SummarisedSeg,vector<int> OutIndexVector,float Threshold);
bool * BorderExtraction(TreeEM * TreeToAnalyse,SEG_ANALYSIS * segment_analysis, int TissueType);
int * LesionClassification(int * ComponentLabel, vector<Lesion *> LesionToClassify,int * Dim);
int * OutlierClassification(int * ComponentLabel, vector<Outlier *> OutlierVectorToClassify,int * Dim);
OutlierType OutlierTypePreliminary(Outlier * OutlierToClassify,SEG_ANALYSIS * segment_analysis, float mult_modify=1);
OutlierType ClassifyOutlierWMH(Outlier * OutlierToClassify, SEG_ANALYSIS * segment_analysis, float mult_modify=1);
vector<int > SymmetryLabelsFP(vector<Outlier*> VectorOutlier, int Label, int * OrderedLabels, int* Dim, int Direction);
void CorrectionForSymmetry(vector<Outlier*> VectorOutliers, int Direction, float ThresholdDifference, float Range, int * OrderedLabels, int * Dim);
OutlierType ClassifyOutlierLacune(Outlier * OutlierToClassify,SEG_ANALYSIS * segment_analysis);
OutlierType ClassifyOutlierSC(Outlier * OutlierToClassify, SEG_ANALYSIS * segment_analysis);
void LesionCode(Lesion * LesionToClassify);
void LesionCodeImprovedT1T2(Lesion * LesionToCode);
void LesionCodeImproved(Lesion * LesionToCode);
int * CountLesionType(vector<Lesion*> LesionVector);
float * VolumeLesionType(vector<Lesion*> LesionVector);
float TotVolumeLesion(vector<Lesion*> LesionVector);


/***************  FUNCTIONS FOR LAPLACE THICKNESS SOLVING *******************************/
float * CreateRidge(float * DistData, int numel, int *Dim, int*Shift,float* PixDim, int omitted_dir);
float * ExtractionSkeletonBis(bool * LesionBool, int numel, int* Dim, int* Shift, float* PixDim,nifti_image * BasisImage);
float * ExtractionSkeleton(bool * LesionBool, int numel, int* Dim, int* Shift, float* PixDim,nifti_image * BasisImage);
float * NormDistWithinLesionRidge(bool * Lesion, bool * Ridge,int numel, vector<int> DimVector, int* Dim, int* Shift, float* PixDim);
float * NormDistWithinLesion(bool * Lesion,int numel,vector<int> DimVector,int* Dim,int* Shift,float* PixDim,int dilationNumber,bool * maxMask,nifti_image * BasisImage);
float * NormDistInside(int CoG,bool * Lesion, int numel, bool*maxMask, vector<int> DimVector, int* Dim, int* Shift, float* PixDim);
float * DistanceFromLobe(float * LaplaceSolution, bool * ObjectIn, bool* ObjectOut, int* L2S, int* S2L,int numelmasked, float * PixDim, int * Dim, int * Shift,nifti_image * TestSave);
float * NormalisedLaplaceLength(float * LaplaceSolution, bool * Border0, bool * Border1, int * L2S, int * S2L, int numelmasked, float * PixDim, int * Dim, int * Shift,nifti_image * TestSave=NULL);
int * CreateLayersFromNormDistSol(float * NormDistSol, int numbLaminae, int numelmasked, float epsilon=0);
bool * CreateBorderFromBool(bool * SegToConsider, int * Dim, int * Shift);
bool * CreateBorderExtFromBool(bool * SegToConsider, int * Dim,int * Shift,int*L2S);
bool * CreateBorderFromBool_L2S(bool * SegToConsider, int * Dim, int * Shift, int * L2S, int * S2L, int numelmasked);
void InitialiseDistanceGivenBorder_L2S(float * Distance, bool * Border, int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked,bool * MaskIn);
float * SolvingLaplaceEquation(int * LabelRegions, int * Dim, int * Shift, float * PixDim, int DistMax, int * L2S, int * S2L, int numelmasked,nifti_image * TestToSave=NULL);
int * CreateLaplaceLabelling(bool * SegMini, bool * SegMaxi, int * Dim);
void LengthSolvingWithLaplace(float * Distance, bool * Border, float * LaplaceSolution,int * Dim, int * Shift,float * PixDim,int * L2S, int * S2L, int numelmasked,nifti_image * TestSave=NULL, bool* Mask=NULL);
float * ShiftedVersion(float * ToShift, int * Trans, int * Dim, int * Shift, int * L2S, int * S2L, int numelmasked);
void UpdateFluxIteration(float * NewUpdated, float * OldFromWhichToCalculate, int * Dim, int * Shift, float * PixDim, int Dist, int * L2S, int * S2L, int numelmasked);
void RefinementUpdatedFlux(float * FluxToRefine, float * OldFluxValue, int * LabelRegions,int * Dim, int* Shift, float * PixDim,int Dist, int * L2S, int * S2L, int numelmasked);
void MappingSolvingWithLaplace(float * Mapping, bool * Border, float * LaplaceSolution,int * Dim, int * Shift,float * PixDim,int * L2S, int * S2L, int numelmasked,nifti_image * TestImage);
float * LaplaceMappingExt(float * LaplaceSolution, bool * ObjectOut, int * L2S, int * S2L, int numelmasked, float * PixDim, int * Dim, int * Shift,nifti_image * TestSave);
void InitialiseMappingGivenBorder_L2S(float * Mapping, bool * Border, int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked, bool * MaskIn, int d);
void MappingSolvingWithLaplace3D(float * Mapping, bool * Border, float * LaplaceSolution,int * Dim, int * Shift,float * PixDim,int * L2S, int * S2L, int numelmasked,nifti_image * TestImage);
float * PerformMapping(float * Mapping, nifti_image * MapImage, int* L2S);
float * CorrectMapping(float* Mapping, bool * Border, int * L2S, int* S2L, int* Dim, int* Shift, float* PixDim, int numelmasked);

void RecursiveConsideration(float * ShiftedVersion, bool * NotConsidered, int *Dim, int *Shift,int * L2S, int * S2L, int numelmasked);
vector<int> IndicesCrossedDirection( int SourceIndex, int TargetIndex, int ExploRadius, int * Dim, int * Shift);
float * GradientImage(float * ToGrad, int dim, int Type, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
float * Tangent(float * ToTangent,int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked );
float * NormedGradientSquared(float * ToGrad, int dim, int sign, int * Dim, int * Shift,float * PixDim,int * L2S, int * S2L, int numelmasked);
float * SecondCrossedDer(float * ToGrad, int dim1, int dim2, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
float * Hessian(float * ToGrad, int dim1, int dim2, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
float * Vesselness(float * ToGrad,int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked,int polarity);
float * NormNormedGradient(float * ToGrad, int Type, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
float * EdginessBorderLesion(bool * LesionBool, float * Data, int* Dim,int*Shift, float*PixDim,int numel,int numbmodal );
float * CreateLineFromEigenCoMBB(float * Eigen, int  CoM, int* BB,int* Dim,int* Shift);
float * CreateLineFromCoMBB(bool* Seg, int CoM, int* BB,float * MainVector,int* Dim,int* Shift);
nifti_image * CreateEntropy(nifti_image * SegToAnalyse, SEG_ANALYSIS * segment_analysis);

void RegionGrowing(float * Map, int* Dim, int* Shift, int Neigh, float * MapFilling, float Thresh, nifti_image* ImageSave=NULL);
float * SecondDer(float * ToGrad, int dim, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
int * InitialiseStatusFromDistance(float * PhiIn, int * Dim, bool * Mask=NULL);
int * InitialisationSignFunction(bool * SegInit, int * Dim, bool * Mask=NULL);
float * InitialiseDistanceFromSignFunction(int * SignFunction, int * Dim, int * Shift);
float * CurvatureCalculation_2(float * ImageToCurve, int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked);
float TotalEnergyCalculation(float * NewFlux, int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked);
vector<int> IndicesCrossedDirection( int SourceIndex, int TargetIndex, int ExploRadius, int * Dim, int * Shift);
vector<int> IndicesToExploreVector(int CentralIndex, int ExploRadius, int * Dim, int * Shift);
vector<int> RefineVectorList(vector<int> ListToConsider, vector<int> ToDeconsider);
int * CorrespondingCoordinates(int Index, int * Dim, int * Shift);
int * CorrespondingCoordinates_ter(int Index, int * Dim, int * Shift,int sizeDim);
void CorrespondingCoordinates_bis(int * ResultCoordinates, int Index, int * Dim, int * Shift);
nifti_image * CorrectionNii(nifti_image * NiiImage, bool * Correction);
void OROperationBool(bool * Array1, bool * Array2, bool * Result, int numel);
void ANDOperationBool(bool * Array1, bool * Array2, bool * Result, int numel);
void XOROperationBool(bool * Array1, bool * Array2, bool * Result, int numel);
void AddTo(float * ImageToModify, float * ImageToAdd, int numel,bool * Mask=NULL);
bool * OpposeBoolArray(bool * ArrayToOppose, int numelmasked);
void MultiplyFloatArrayBy(float * FloatArrayMultiplied, float MultiplicativeFloat, int numel,bool * Mask=NULL);
float * MultiplyElementwise(float * Array1, float * Array2, int numel,bool * Mask=NULL);
float * AddElementwise(float * Array1, float * Array2, int numel,bool * Mask=NULL);
int * AddElementwiseInt(int * Array1, int * Array2, int numel);
float * AddElementwiseMasked(float * Array1, float * Array2, bool * MaskArray, int numel);
float * AddElementwise(vector<float *> VectorArray, int numel,bool * Mask=NULL);
float * MinElementwise(vector<float *> VectorArray,int numel,bool * Mask=NULL);
float * MaxElementwise(vector<float *> VectorArray,int numel,bool * Mask=NULL);
float * MultiplyElementwise(vector<float *> VectorArray, int numel,bool * Mask=NULL);
float * PowerElementwise(float * ArrayToPower, float power, int numel,bool * Mask=NULL);
float * DivideElementwise(vector<float *> VectorArray, int numel,bool * Mask=NULL);
float * InverseElementwise(float * ArrayToInvert, int numel,bool * Mask=NULL);
bool ValidityCoordinates(int * CoordinatesToCheck, int * Dim);
void RecursiveConsideration(float * ShiftedVersion, bool * NotConsidered, int *Dim, int *Shift);


/**************************** FUNCTIONS FOR LEVEL SET AND DISTANCE FAST MARCHING SOLVING ***********************************/
typedef struct PARAM_LS {
    float EpsCurv;
    float WidthNB;
    float PropNBSusp;
    int numbLaminae;
    bool flag_EquiVolume;// Define if EquiVolume strategy must be chosen or not
    float MaxDist;//Distance maximum used for solving paraboloid fitting
    int MaxNC; //Number of times no change in level set 0 curve before accepting final result
}PARAM_LS;


vector<float *> CompleteEquiVolumeLSSolving(bool * SegMini, bool * SegMaxi, int * Dim, float * PixDim, PARAM_LS* Param_LS,int * L2S, int * S2L, int & numelmasked,float * PhiIn=NULL, float * PhiOut=NULL, nifti_image * TestToSave=NULL);

int * CreateSignFunction(bool * BinarySeg, int * Dim);
bool * CreateCriterionFromSign(int * Sign, int TypeCrit, int numel, bool * Mask);
bool CheckAndCorrectCrossing(bool * SegMini, bool * SegMaxi, int * Dim);
nifti_image * HolesJuxta(TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
void CorrectWMI(TreeEM* TreeToAnalyse,nifti_image* ToCorrect,nifti_image* BaseCorrection, SEG_ANALYSIS* segment_analysis);
bool * CheckAndFillHoles(bool * SegToCheck, int * Dim, int * Shift);
int * FindBoundingBox(bool * SegToCheck, int * Dim , int * Shift);
int * ListCoordinatesBoolSeg(bool *LesionBool, int *Dim, int *Shift);
float * ListCoordinatesTransformed(int *ListCoordinates, int LabelSize,float * PixDim);
float * DemeanListCoordinates(float * ListCoordinates, int numbActive, int numbDir, log_level_t LOG_VERB=LOG_NOTHING);
bool IsBoolOutsideBoundingBox(bool * LesionBool, int * BoundingBox, int * Dim);
bool * CreateBoolSegFromSumPriorsIdx(nifti_image * SummarisedSeg,vector<float *> PriorsVector ,vector<int> IndexVector,int * L2S,float PriorsThreshold );
bool * CreateBoolBoundingBox(int * BoundingBox, int * Dim);
bool CheckForHoles(bool * SegToCheck, int * Dim, int * Shift);
bool * CreateMask(bool * SegMini, bool * SegMaxi, int * Dim);

template <typename T1, typename T2> T1 * CorrectByMask(T1* ArrayToChange, T2* ArrayToUseForMasking, int numel, T2 Thresh=0){
    T1* CorrectedArray = new T1[numel];
    for (int i=0;i<numel;i++){
        CorrectedArray[i]=ArrayToUseForMasking[i]>Thresh?ArrayToChange[i]:0;
    }
    return CorrectedArray;
}

template <typename T> T * CopyArray(T * ArrayToCopy, int numel){
    T * ArrayCopied=new T [numel];
    for (int i=0; i<numel; i++) {
        ArrayCopied[i]=ArrayToCopy[i];
    }
    return ArrayCopied;
}

template <typename T1,typename T2> void AddElementwiseInPlace(T1* ArrayToAddTo, T2* ArrayToAdd,int numel){
    for(int i=0;i<numel;i++){
        ArrayToAddTo[i]+=(T1)ArrayToAdd[i];
    }
    return;
}

template <typename T1, typename T2> T2 * TranscribeArray(T1 * ArrayToTranscribe, int numel){
    T2 * ArrayTranscribed=new T2[numel];
    for (int i=0; i<numel; i++) {
        ArrayTranscribed[i]=(T2)ArrayToTranscribe[i];
    }
    return ArrayTranscribed;
}

template <typename T1, typename T2, typename T3> T3* MaskArray(T1* Array, T2* Mask, int numel){
    T3 * MaskedArray = new T3[numel];
    if (Mask==NULL){

        for(int i=0;i<numel;i++){
            MaskedArray[i]=(T3)Array[i];
        }
    }
    else{

    for(int i=0;i<numel;i++){
        MaskedArray[i]=Mask[i]>0?(T3)Array[i]:0;
    }
    }
    return MaskedArray;
}

template <typename T> T* AddArray(vector<T *> ArrayVector, int numel,bool * Mask=NULL){
    T * ArraySum=new T[numel];
    for (int i=0; i<numel; i++) {
        ArraySum[i]=0;
    }
    int numbToSum=ArrayVector.size();
    if (Mask!=NULL) {
        for (int i=0; i<numel; i++) {
            if (Mask[i]) {
                for(int n=0;n<numbToSum;n++){
                    ArraySum[i]+=ArrayVector[n][i];
                }
            }
        }
    }
    else{
    for (int i=0; i<numel; i++) {
        for(int n=0;n<numbToSum;n++){
            ArraySum[i]+=ArrayVector[n][i];
        }
    }
    }
    return ArraySum;
}


template <typename T> nifti_image* AddNii(vector<nifti_image *> NiiVector){
    int initInd=FirstNotNULL(NiiVector);
    int numel=NiiVector[initInd]->nvox;
    T * ArraySum=new T[numel];
    for (int i=0; i<numel; i++) {
        ArraySum[i]=0;
    }
    int numbToSum=NiiVector.size();
            for(int n=initInd;n<numbToSum;n++){
                if(NiiVector[n]!=NULL){
                    T* ArrayData=static_cast<T*>(NiiVector[n]->data);
    for (int i=0; i<numel; i++) {

            ArraySum[i]+=ArrayData[i];
        }
    }
    }
            nifti_image * AddedNii=CreateNiiFromArray(ArraySum,NiiVector[initInd],numel);
            delete[] ArraySum;
    return AddedNii;
}

template<typename T> int GetIndexNValue(T* ArrayToSort, int N, bool order, int numel){
    T* ArrayToSortB=CopyArray(ArrayToSort, numel);
    if (N>numel){
        return -1;
    }
    quickSort(ArrayToSortB, numel);
    T Value;
    if (order) {
        Value=ArrayToSortB[N];
    }
    else{
        Value=ArrayToSortB[numel-1-N];
    }
    for (int i=0; i<numel; i++) {
        if (ArrayToSort[i]==Value) {
            return i;
        }
    }
    return 0;
}

//template <typename T> vector<size_t> sort_indices(vector<T> &v){
//    vector<size_t> idx(v.size());
//    for(size_t i=0; i!=idx.size();++i) idx[i]=i;
//    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2){return v[i1]>v[i2];});
//    return idx;
//}

template <typename T> T* SubtractArray(T* Array1, T* Array2, int numel){
    T * ArraySub=new T[numel];
    for (int i=0; i<numel; i++) {
        ArraySub[i]=0;
    }
        for (int i=0; i<numel; i++) {
                ArraySub[i]=Array1[i]-Array2[i];
        }
    return ArraySub;
}

template <typename T1, typename T2> void MultiplyElementwiseChange(T1* Array1, T2* Array2, int numel){
    for (int i=0; i<numel; i++) {
        Array1[i]*=Array2[i];
    }
}

template <typename T1, typename T2> vector <T2> GetListValues(T1 *Array,int numel){
    set<T2> SetOfValues;
    for (int i=0; i<numel; i++) {
        if (Array[i]>0) {
            SetOfValues.insert((T2)Array[i]);
        }
    }

    vector<T2> Result(SetOfValues.begin(), SetOfValues.end());
    return Result;
}

template <typename T1, typename T2, typename T3> T3* MultiplyElementwiseChoose(T1 * Array1, T2* Array2, int numel){
    T3 * ResultArray=new T3[numel];
    for (int i=0; i<numel; i++) {
        ResultArray[i]=(T3)(Array1[i]*Array2[i]);
    }
    return ResultArray;
}

template <typename T> T* CreateShort(T * ArrayToShort,int * S2L,int numelmasked){
    T* ShortArray=new T[numelmasked];
    for (int i=0; i<numelmasked; i++) {
        ShortArray[i]=ArrayToShort[S2L[i]];
    }
    return ShortArray;
}

template<typename T> T* CreateLong(T * ArrayToLong, int * L2S, int numel){
    T* LongArray=new T[numel];
    for (int i=0; i<numel; i++) {
        LongArray[i]=0; // we pan the other values with 0
        if (L2S[i]>=0) {
            LongArray[i]=ArrayToLong[L2S[i]];
        }
    }
    return LongArray;
}


template<typename T> T* CreateLongPadding(T * ArrayToLong, T ValuePad, int * L2S, int numel){
    T* LongArray=new T[numel];
    for (int i=0; i<numel; i++) {
        LongArray[i]=ValuePad; // we pad the other values with the ValuePad
        if (L2S[i]>=0) {
            LongArray[i]=ArrayToLong[L2S[i]];
        }
    }
    return LongArray;
}

template<typename T> int CountNonZero(T * ArrayToCount, int numel){
    int CountNonZeroRes=0;
    for (int i=0; i<numel; i++) {
        if (fabs((float)ArrayToCount[i])>10E-6) {
            CountNonZeroRes++;
        }
    }
    return CountNonZeroRes;
}

template<typename T1, typename T2> T2 * ThresholdArray(T1* ArrayToThreshold, T1 Threshold, int numel){
    T2* ThresholdedArray=new T2[numel];
    for (int i=0; i<numel; i++) {
        if (ArrayToThreshold[i]>Threshold) {
            ThresholdedArray[i]=(T2)ArrayToThreshold[i];
        }
        else{
            ThresholdedArray[i]=0;
        }
    }
    return ThresholdedArray;
}


template<typename T> T* CreateLongPaddingMulti(T* ArrayToLong, T ValuePad, int * L2S, int numel, int numbmodal){
    T* LongArray=new T[numel*numbmodal];
    int * L2SThresh=ThresholdArray<int,int>(L2S, -1, numel);
    int numelmasked=CountNonZero(L2SThresh, numel)+1;
    cout << "Number numelmasked is "<< numelmasked<<endl;
    delete [] L2SThresh;
    for (int i=0; i<numel; i++) {
        for (int m=0;m<numbmodal;m++){
            LongArray[i+m*numel]=ValuePad; // we pad the other values with the ValuePad
            if (L2S[i]>=0) {
                LongArray[i+m*numel]=ArrayToLong[L2S[i]+m*numelmasked];
            }
        }
    }
    cout << "Finished Padding multi"<<endl;
    return LongArray;
}


template<typename T1, typename T2> T2 * EqualArray(T1* ArrayToEqual, T1 EqualValue, int numel){
    T2* ThresholdedArray=new T2[numel];
    for (int i=0; i<numel; i++) {
        if (ArrayToEqual[i]<=(T1)(EqualValue+0.001) && ArrayToEqual[i]>=(T1)(EqualValue-0.001)) {
            ThresholdedArray[i]=(T2)ArrayToEqual[i];
        }
        else{
            ThresholdedArray[i]=0;
        }
    }
    return ThresholdedArray;
}



template <typename T> T * ExpArray(T* ArrayToExp, int numel){
    T * ExponentiatedArray=new T[numel];
    
    cout << numel <<"for exponentiation "<<endl;
    cout<<ArrayToExp[0]<<endl;
    cout<<ArrayToExp[numel-1]<<endl;
    for (int i=0; i<numel; i++) {
            ExponentiatedArray[i]=exp(ArrayToExp[i]);
    }
    return ExponentiatedArray;
}

template<typename T1, typename T2> T2 * UpperThresholdArray(T1* ArrayToThreshold, T1 Threshold, int numel){
    T2* ThresholdedArray=new T2[numel];
    for (int i=0; i<numel; i++) {
        if (ArrayToThreshold[i]<Threshold) {
            ThresholdedArray[i]=(T2)ArrayToThreshold[i];
        }
        else{
            ThresholdedArray[i]=0;
        }
    }
    return ThresholdedArray;
}

template <typename T1, typename T2> float SumOverMask(T1 * ArrayToSum, T2 * MaskToUse, T2 Threshold, int numel){
    float Result=0;
    if(MaskToUse==NULL){
        for (int i=0; i<numel; i++) {
            Result+=(float)ArrayToSum[i];
        }
        return Result;
    }
    else{
        for (int i=0; i<numel; i++) {
            if (MaskToUse[i]>Threshold) {
                Result+=(float)ArrayToSum[i];
            }
        }
        return Result;
    }
}

template <typename T> bool IsCoordValid(T * Point, int * Dim){
    for (int d=0; d<3; d++) {
        if (Point[d]<0 || Point[d]> Dim[d]-1) {
            return 0;
        }
    }
    return 1;
}

template <typename  T1, typename T2> float GetDistanceBetweenPoints(T1 * Point1 , T2 * Point2, int * Dim, float * PixDim) {
    if (!IsCoordValid<T1>(Point1, Dim) || !IsCoordValid<T2>(Point2, Dim)){
        return -1;
    }
    float Dist=0;
    for (int d=0; d<3; d++) {
        Dist+=(Point1[d]-Point2[d])*(Point1[d]-Point2[d])*PixDim[d]*PixDim[d];
    }
    Dist=sqrtf(Dist);
    return Dist;
}

float GetDistanceBetweenPoints(int Point1 , int Point2, int * Dim, float * PixDim);


// Save an array converted to float in a corresponding nifti image.
template <typename T> nifti_image * CreateNiiFromArray(T * ArrayForNii, nifti_image * BasicFloatImage,int numel){
    // We assume that there are numel elements in ArrayForNii and that there is no access problem when creating the nii image.
    // Moreover, we assume that the BasicFloatImage is also a nii image using floats.
    
    int numbmodal=numel/(BasicFloatImage->nx*BasicFloatImage->ny*BasicFloatImage->nz);
    if (numbmodal<=0) {
        return NULL;
    }
    nifti_image * BuiltNii=nifti_copy_nim_info(BasicFloatImage);
    BuiltNii->dim[0]=numbmodal>1?4:3;
    BuiltNii->dim[4]=numbmodal;
    BuiltNii->dim[5]=1;
    nifti_update_dims_from_array(BuiltNii);
//    cout<<"Number of modalities is "<<numbmodal<<endl;
//    cout<<"Number of elements is "<<numel<<endl;
    BuiltNii->data=(void *)calloc(BuiltNii->nvox, sizeof(float));
    float * BuiltData=static_cast<float *>(BuiltNii->data);
    cout<<"Memory allocated"<<endl;
    for (int i=0; i<numel; i++) {
        BuiltData[i]=(float)ArrayForNii[i];
    }
    cout << "Built"<<endl;
    return BuiltNii;
}





typedef struct DLS{
    float Distance;
    int Index;
    int Status;
    DLS(){
        Distance=-1;
        Index=-1;
        Status=-1;
    }
    ~DLS(){}
    DLS * CopyDLS_ptr(){
        DLS * DLSCopied=new DLS();
        DLSCopied->Status=this->Status;
        DLSCopied->Distance=this->Distance;
        DLSCopied->Index=this->Index;
        return DLSCopied;
    }
    DLS CopyDLS(){
        DLS DLSCopied;
        DLSCopied.Distance=this->Distance;
        DLSCopied.Index=this->Index;
        DLSCopied.Status=this->Status;
        return DLSCopied;
    }
}DLS;

bool CheckAndCorrectCrossing(bool * SegMini, bool * SegMaxi, int * Dim);
bool * CreateMask(bool * SegMini, bool * SegMaxi, int * Dim);


void HeapSort(DLS * DLSArray,int n);
void BuildMaxHeap(DLS * DLSArray,int n);
void MaxHeapify(DLS * DLSArray,int i,int n);
void HeapSort(vector<DLS> & DLSVector,int n);
void BuildMaxHeap(vector<DLS>  & DLSVector,int n);
void MaxHeapify(vector<DLS> & DLSVector,int i,int n);
void HeapSort(vector<DLS *> & DLSVector,int n);
void BuildMaxHeap(vector<DLS * >  & DLSVector,int n);
void MaxHeapify(vector<DLS *> & DLSVector,int i,int n);
int FindRapidLowerBound(vector<DLS> DLSVector, DLS DLSToFind, int IndexLastInserted=0);
int FindRapidUpperBound(vector<DLS> DLSVector, DLS DLSToFind, int IndexLastInserted=0);
int FindRapidLowerBound(vector<DLS *> DLSVector, DLS * DLSToFind, int IndexLastInserted=0);
int FindRapidUpperBound(vector<DLS *> DLSVector, DLS * DLSToFind, int IndexLastInserted=0);
void RecodeDistance(float & ToRecode, float * PixDim,int * Dim, int * ShiftRecode);
bool NearBorderNB(bool * NarrowBand, bool * NarrowBandSuspicious, int numel);
DLS * InitialiseDLSArray(int * InitialStatus, float * InitialDistance, int * Dim);
vector<DLS> InitialiseDLSVector(int * InitialStatus, float * InitialDistance, int * Dim);
vector<DLS *> InitialiseDLSVectorPtr(int * InitialStatus, float * InitialDistance, int * Dim);
float * InitialiseObjectiveLS(float Alpha, float * AMap, float * PhiIn, float * PhiOut, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
float * InitialiseObjectiveLS(float * RhoMapToUse, float * PhiIn, float * PhiOut, int numelmasked);
void InitialiseState(int * StateToInitialise, float * DistanceInitialise, int * Dim,int *Shift, float * PixDim,bool * Mask=NULL);
float * InitialiseForceField(float * PhiAlpha, float * PhiObjective, int numelmasked, bool * NarrowBand);
bool * NarrowBandDefine(float * CurrentPhi, float WidthNB,int numelmasked);
void UpdateNarrowBand(bool * NarrowBand, float * CurrentPhi,float WidthNB, int numelmasked);
float * UpdateMatrix(float * ForceField, bool * NarrowBand, float * CurrentPhi, int * Dim,int * Shift,float * PixDim, float EpsCurv, int * L2S, int * S2L, int numelmasked);
void UpdateForceField(float * CurrentPhi, float *ForceField, float * InitialObjective, bool * NarrowBand, int numelmasked);
bool IsNearBorderNarrowBand(bool * NarrowBandSuspicious, int * Sign, int numelmasked);
bool TotalUpdateLS(float * CurrentPhi, int * Sign, bool * NarrowBand, float * ForceField, int * Dim, int * Shift, float * PixDim, float EpsCurv,int * L2S, int * S2L, int numelmasked, nifti_image * TestToSave=NULL);
int TopologicalNumber(int Index, bool * Criterion, int * Dim, int * Shift, int TypeNeighborhood,int * L2S, int * S2L, int numelmasked);
vector<int > GeodesicNeighborhood(int Index, int TypeNeighborhood, bool* Criterion, int LengthNeighborhood,int *Dim,int * Shift,int * L2S, int * S2L, int numelmasked);
bool IsIn(int Index, int * ListIndices, int LengthList);
bool IsIn(int Index, vector<int> ListIndices);
float SolveQuadraticForDistance(int Index, int * Status, float * DistanceFinal, int * Dim, int * Shift, float * PixDim);
float SolveQuadraticForDistance_short(int Index, int * Status, float * DistanceFinal, int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked);
nifti_image * EuclideanDistanceImage(nifti_image * BasisImage, bool * BasisObject, nifti_image * Mask);
void IterFastMarchingMap(int * Status, float * DistanceFinal, int * Dim, int * Shift, float * PixDim, int & Iteration, bool * Mask=NULL);
void IterFastMarchingPtr(int * Status, float * DistanceFinal, int * Dim, int * Shift, float * PixDim, int & Iteration, bool * Mask=NULL);
void IterFastMarching(int * Status, vector<DLS> & DLSVector, float * DistanceFinal, int * Dim, int * Shift, float * PixDim,int & Iteration, bool * Mask=NULL);
void RecurseFastMarching(int * Status, vector<DLS> & DLSVector, float * DistanceFinal, int * Dim, int * Shift, float * PixDim,int & Iteration, bool * Mask=NULL);
int FindProjection(int IndexProjected, float * ProjectionSurface, int * Dim, int * Shift, float * PixDim,int * L2S, int * S2L, int numelmasked);
int * CreateSignFromSignedDistance(float * PhiIn, float * PhiOut, int * Dim);
float * CurvatureParaboloidCalculation(int Index,float * Phi,int * Dim,int *Shift,float MaxDist,float *PixDim,int * L2S, int * S2L, int numelmasked,float Std=1);
int * ProjectingMapIndex(float * ProjectionSurface, int * Dim, int * Shift, float * PixDim, int * L2S, int * S2L, int numelmasked);
float * AMapCalculation(float * PhiIn, float * PhiOut, int * Dim,int * Shift, float * PixDim, float MaxDist, int * L2S, int * S2L, int numelmasked);
float * RhoMap(float * AMap, int * ProjIn, int * ProjOut, float * PhiIn, float * PhiOut, int * Dim, int * Shift, float * PixDim,float Alpha, int * L2S, int * S2L, int numelmasked);
bool * CreateBoolObjectFromParcellation(nifti_image * ParcellationImage, vector<int> ElementsToAssociate);

/*********************** FUNCTIONS FOR LESION FEATURES IN TERMS OF NORMALISED DISTANCE BETWEEN VENTRICLES AND CORTEX ********************/
void GetWeightedBin(float * WeightIndexBin,int index,int NumberLaminae,float *LengthNormSol_s, int numelmasked, int * L2S, int * S2L, nifti_image * BasisImage);
float * ProportionLesionVolumeDistanceWeighted(nifti_image * LesionSeg, float * LengthNormSol, int NumberLaminae,float * PixDim, int * L2S, int * S2L, int numelmasked);
float * ProportionLesionVolumeDistance(nifti_image * LesionSeg, float * LengthNormSol, int NumberLaminae,float * PixDim, int * L2S, int * S2L, int numelmasked);
nifti_image * OutliersLac(SimpleRule * Rule, TreeEM* TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
float * AnisotropicMorphologicalChange(float * CountHistogram1, int KernelSize, vector<int> dim, bool Value, float * PixDim);
float * Erosion_bis(float * CountHistogram1, int KernelSize, vector<int> dim, bool Value);

int * ExtentLesionNormDist(int * LesionLabel, float * SolNormDist, int NumberLaminae, int * L2S, int * S2L, int numelmasked );


bool * Dilation(bool * ToDilate, int Dilation, int * Dim, int * Shift);
int * MakeL2S(bool * Mask, int * Dim);
int * MakeS2L(bool * Mask, int * Dim,int & numelmasked);
void  SaveTmpResult(float * ResultToSave, string FilenameToSave, nifti_image * ImageToUse);
void  SaveTmpResult(vector<float *> ResultToSave, string FilenameToSave, nifti_image * ImageToUse);

template<typename T> void SaveTmpResult_short(T * ResultToSave, string FilenameToSave, nifti_image * ImageToUse,int * L2S){
    // Create Image to fill with ResultToSave
    nifti_image * ImageToFill=nifti_copy_nim_info(ImageToUse);
    ImageToFill->data=(void *)calloc(ImageToFill->nvox, sizeof(float));
    float * Filling_PTR=static_cast<float *>(ImageToFill->data);
    int numel=ImageToFill->nvox;
    int j=0;
    int * L2S_PTR=L2S;
    for (int i=0; i<numel; i++,L2S_PTR++) {
        if (*L2S_PTR>=0) {
            Filling_PTR[i]=(float)ResultToSave[j];
            j++;
        }
        else{
            Filling_PTR[i]=0;
        }
    }
    nifti_set_filenames(ImageToFill, FilenameToSave.c_str(), 0, 0);
    nifti_image_write(ImageToFill);
    nifti_image_free(ImageToFill);
    ImageToFill=NULL;
    return;
}
// Save TemporaryResult based on existing nifti_image float
template<typename T >void  SaveTmpResult(T * ResultToSave, string FilenameToSave, nifti_image * ImageToUse){
    // Create Image to fill with ResultToSave
    nifti_image * ImageToFill=nifti_copy_nim_info(ImageToUse);
    ImageToFill->data=(void *)calloc(ImageToFill->nvox, sizeof(float));
    float * Filling_PTR=static_cast<float *>(ImageToFill->data);
    int numel=ImageToFill->nvox;
    for (int i=0; i<numel; i++) {
        Filling_PTR[i]=(float)ResultToSave[i];
    }
    nifti_set_filenames(ImageToFill, FilenameToSave.c_str(), 0, 0);
    nifti_image_write(ImageToFill);
    nifti_image_free(ImageToFill);
    ImageToFill=NULL;
    return;
}

template<typename T> int GetIndexMax(T * ArrayToMax, int numel){
    T MaxFound=-10000000;
    int IndexMax=0;
    for (int i=0; i<numel; i++) {
        if (MaxFound<ArrayToMax[i]) {
            MaxFound=ArrayToMax[i];
            IndexMax=i;
        }
    }
    return IndexMax;
    
}
template<typename T> int GetIndexSecondMax(T * ArrayToMax, int numel, int ind_max){
    T MaxFound=-10000000;
    int IndexMax=0;
    for (int i=0; i<numel; i++) {
        if (MaxFound<ArrayToMax[i] && i!=ind_max) {
            MaxFound=ArrayToMax[i];
            IndexMax=i;
        }
    }
    return IndexMax;

}

template<typename T> int GetIndexMaxMask(T * ArrayToMax,bool * Mask, int numel){
    T MaxFound=-10000000;
    int IndexMax=0;
    for (int i=0; i<numel; i++) {
        if (Mask[i]){
        if (MaxFound<ArrayToMax[i]) {
            MaxFound=ArrayToMax[i];
            IndexMax=i;
        }
        }
    }
    return IndexMax;

}

template<typename T> int GetIndexMin(T * ArrayToMin, int numel){
    T MinFound=10000000;
    int IndexMin=0;
    for (int i=0; i<numel; i++) {
        if (MinFound>ArrayToMin[i]) {
            MinFound=ArrayToMin[i];
            IndexMin=i;
        }
    }
    return IndexMin;
    
}

template<typename T> int GetIndexMinMasked(T * ArrayToMin, bool * Mask, int numel){
    T MinFound=10000000;
    int IndexMin=-1;
    for (int i=0; i<numel; i++) {
        if (MinFound>ArrayToMin[i] && Mask[i]) {
            MinFound=ArrayToMin[i];
            IndexMin=i;
        }
    }
    return IndexMin;
    
}

template<typename T> void AbsoluteInPlace(T* Array, int numel){
    for(int i=0;i<numel;i++){
        Array[i]=(T)sqrt(Array[i]*Array[i]);
    }
    return;
}

template <typename T> void AddConstantInPlace(T* Array, T ToAdd, int numel){
    for(int i=0;i<numel;i++){
        Array[i]+=ToAdd;
    }
    return;
}


template <typename T> void MultiplyElementwiseInPlace(T* Array, T* ToMultiply, int numel){
    for(int i=0;i<numel;i++){
        Array[i]*=ToMultiply[i];
    }
    return;
}

template<typename T> T* Absolute(T* Array,int numel){
    T* Result=new T[numel];
    for(int i=0;i<numel;i++){
        Result[i]=(T)sqrt(Array[i]*Array[i]);
    }
    return Result;
}

template<typename T> T * Sign(T* Array,int numel){
    T* Result=new T[numel];
    for(int i=0;i<numel;i++){
        Result[i]=0;
        if(Array[i]>0){
            Result[i]=1;
        }
        else if(Array[i]<0){
            Result[i]=-1;
        }
    }
    return Result;
}

float GetSum(float * ArrayToSum, int numel);
float GetMax(float * ArrayToMax,int numel);
template<typename T> T GetMax(T * ArrayToMax, int numel){
    T MaxFound=-1000000;
    if (ArrayToMax==NULL){
        return MaxFound;
    }
    for (int i=0; i<numel; i++) {
        if (MaxFound<ArrayToMax[i]) {
            MaxFound=ArrayToMax[i];
        }
    }
    return MaxFound;
}

template<typename T1,typename T2> T1 GetMaxArrayMasked(T1 * ArrayToMin, T2* Mask, int numel){
    T1 MaxFound=-1000000;
    if (ArrayToMin==NULL){
        return MaxFound;
    }
    for (int i=0; i<numel; i++) {
        if(Mask == NULL){
            if (MaxFound<ArrayToMin[i]) {
                MaxFound=ArrayToMin[i];
            }
        }
        else{
        if (Mask[i]>0) {
            if (MaxFound<ArrayToMin[i]) {
                MaxFound=ArrayToMin[i];
            }
        }
        }
    }
    return MaxFound;
}

template<typename T1,typename T2> T1 * GetMaxDataMulti(T1 * ArrayToMin, T2* Mask, int numel,int numbmodal){
    T1 * MaxFound=new T1[numbmodal];
    for (int m=0; m<numbmodal; m++) {
        MaxFound[m]=-100000;
    }
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            for (int m=0; m<numbmodal; m++) {
                if (MaxFound[m]<ArrayToMin[m*numel+i]) {
                    MaxFound[m]=ArrayToMin[i+m*numel];
                }
            }

        }
    }
    return MaxFound;
}

template<typename T1,typename T2> T1 * GetMinElementwise(vector<T1*> ArrayComp,T2*Mask,int numel){
    T1 * MinOver=new T1[numel];
    bool flag_suppress=0;
    if (Mask==NULL){
        flag_suppress=1;
        Mask=new T2[numel];
        for(int i=0;i<numel;i++){
            Mask[i]=1;
        }
    }
    int Vsize=ArrayComp.size();
    if(Vsize==0){
        return NULL;
    }
    for(int i=0;i<numel;i++){
        T1 Min=numel;
        MinOver[i]=0;
        if(Mask[i]){
            for(int s=0;s<Vsize;s++){
                if (ArrayComp[s][i]<Min){
                    Min=ArrayComp[s][i];
                }
            }
            MinOver[i]=Min;
        }
    }
    if(flag_suppress){
        delete [] Mask;
        Mask=NULL;
    }
    return MinOver;
}

template<typename T1,typename T2> T1 * GetMinDataMulti(T1 * ArrayToMin, T2* Mask, int numel,int numbmodal){
    T1 * MaxFound=new T1[numbmodal];
    for (int m=0; m<numbmodal; m++) {
        MaxFound[m]=100000;
    }
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            for (int m=0; m<numbmodal; m++) {
                if (MaxFound[m]>ArrayToMin[m*numel+i]) {
                    MaxFound[m]=ArrayToMin[i+m*numel];
                }
            }

        }
    }
    return MaxFound;
}

template<typename T> T GetMin(T * ArrayToMin, int numel){
    T MinFound=1000000;
    for (int i=0; i<numel; i++) {
        if (MinFound>ArrayToMin[i]) {
            MinFound=ArrayToMin[i];
        }
    }
    return MinFound;
}

template<typename T1,typename T2> T1 GetMin(T1 * ArrayToMin, T2* Mask, int numel){
    T1 MinFound=1000000;
    if (ArrayToMin==NULL){
        return MinFound;
    }

    for (int i=0; i<numel; i++) {
        if (Mask==NULL){
            if (MinFound>ArrayToMin[i]) {
                MinFound=ArrayToMin[i];
            }
        }
        else{
        if (Mask[i]>0) {
            if (MinFound>ArrayToMin[i]) {
                MinFound=ArrayToMin[i];
            }
        }
        }
    }
    return MinFound;
}

//template<typename T1,typename T2> T1 * GetMinDataMulti(T1 * ArrayToMin, T2* Mask, int numel, int numbmodal){
//    T1 * MinFound=new T1[numbmodal];
//    for (int m=0; m<numbmodal; m++) {
//        MinFound[m]=100000;
//    }
//    for (int i=0; i<numel; i++) {
//        if (Mask[i]>0) {
//            for (int m=0; m<numbmodal; m++) {
//                if (MinFound[m]>ArrayToMin[i+numel*m]) {
//                    MinFound[m]=ArrayToMin[i+numel*m];
//                }
//            }

//        }
//    }
//    return MinFound;
//}

template <typename T1, typename T2> int GetCountUnder(T1* ArrayToStudy, T2* Mask, T1 Threshold, int numel){
    int CountTrue=0;
    int CountMask=0;
    for (int i=0; i<numel; i++) {
        if(Mask[i]!=0){
            if (ArrayToStudy[i]<Threshold) {
                CountTrue++;
            }
            CountMask++;
        }
    }
    //    cout <<"CountMask is "<<CountMask;
    return CountTrue;
}


template <typename T1, typename T2> int GetCountAbove(T1* ArrayToStudy, T2* Mask, T1 Threshold, int numel){
    int CountTrue=0;
    int CountMask=0;
    for (int i=0; i<numel; i++) {
        if(Mask[i]!=0){
            if (ArrayToStudy[i]>Threshold) {
                CountTrue++;
            }
            CountMask++;
        }
    }
//    cout <<"CountMask is "<<CountMask;
    return CountTrue;
}


template<typename T1,typename T2> float GetProportionAbove(T1* ArrayToStudy, T2* Mask, T1 Threshold,int numel){
    int CountTrue=0;
    int CountTotal=0;
    if (ArrayToStudy==NULL){
        return 0;
    }
    for (int i=0; i<numel; i++) {
        if (Mask==NULL){
            CountTotal++;
            if (ArrayToStudy[i]>Threshold) {
                CountTrue++;
            }
        }
        else{
        if(Mask[i]!=0){
            CountTotal++;
            if (ArrayToStudy[i]>Threshold) {
                CountTrue++;
            }
        }
        }
    }
    if (CountTotal>0) {
        return (float)CountTrue/(float)(1.0*CountTotal);
    }
    else{
        return 0;
    }
}

template<typename T1,typename T2> float GetProportionUnder(T1* ArrayToStudy, T2* Mask, T1 Threshold,int numel){
    int CountTrue=0;
    int CountTotal=0;
    for (int i=0; i<numel; i++) {
        if(Mask[i]!=0){
            CountTotal++;
            if (ArrayToStudy[i]<Threshold) {
                CountTrue++;
            }
        }
    }
    if (CountTotal>0) {
        return (float)CountTrue/(float)(1.0*CountTotal);
    }
    else{
        return 0;
    }
}

template<typename T1,typename T2> T1 GetMax(T1 * ArrayToMin, T2* Mask, int numel){
    T1 MaxFound=-1000000;
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            if (MaxFound<ArrayToMin[i]) {
                MaxFound=ArrayToMin[i];
            }
        }
    }
    return MaxFound;
}



template<typename T1,typename T2> T1 GetMeanData(T1 * ArrayToMean, T2* Mask , int numel){
    T1 MeanFound=0;
    T1 Normalisation=0;
    if (ArrayToMean==NULL) {
        return MeanFound;
    }
    for (int i=0; i<numel; i++) {
        if (Mask==NULL){
            MeanFound+=Mask[i]*ArrayToMean[i];
            Normalisation+=Mask[i];
        }
        else{
        if (Mask[i]>0) {
            MeanFound+=Mask[i]*ArrayToMean[i];
            Normalisation+=Mask[i];
        }
        }
    }
    if (Normalisation>0) {
        return MeanFound/Normalisation;
    }
    else{
        return 0;
    }
}

template<typename T1,typename T2> T1 * GetMeanDataMulti(T1 * ArrayToMean, T2* Mask, int numel, int numbmodal){
    T1 * MeanResults=new T1[numbmodal];
    for (int m=0; m<numbmodal; m++) {
        MeanResults[m]=0;
    }
//    cout<<GetMax(Mask,numel)<<endl;
    float Normalisation=0;
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            for (int m=0 ; m<numbmodal; m++) {
                MeanResults[m]+=ArrayToMean[numel*m+i]*Mask[i];
            }
            Normalisation+=Mask[i];
        }
    }
    if (Normalisation>0) {
        for (int m=0; m<numbmodal; m++) {
//            cout<<"MeanRes "<<MeanResults[m]<<" "<<Normalisation<<endl;
            MeanResults[m]/=(T1)Normalisation;
        }
        return MeanResults;
    }
    else{
        for (int m=0; m<numbmodal; m++) {
            MeanResults[m]/=Normalisation;
        }
        return MeanResults;
    }
}

template<typename T1,typename T2> vector<T1 *> GetQuantilesMulti(T1 * ArrayToAnalyse, T2 * Mask, int numel, int numbmodal, vector<float> QuantilesToCalculate){
    vector<T1*> QuantilesVector;
    int sizeQuantiles=QuantilesToCalculate.size();
    int CountNonZeroMask=0;
    if (Mask!=NULL){
     CountNonZeroMask=CountNonZero(Mask, numel);
    }

    else{
        CountNonZeroMask=numel;
    }
//        cout<<sizeQuantiles<<" "<<CountNonZeroMask<<endl;
    int * IndicesToChoose=new int[sizeQuantiles];
    for (int s=0; s<sizeQuantiles; s++) {
        int IndQuantMin=(int)round(QuantilesToCalculate[s]*CountNonZeroMask);
        IndQuantMin=IndQuantMin<=0?0:IndQuantMin;
        IndQuantMin=IndQuantMin>=CountNonZeroMask-1?(CountNonZeroMask-1):IndQuantMin;
        IndicesToChoose[s]=IndQuantMin;
//        cout << IndQuantMin << endl;
    }
    for (int m=0; m<numbmodal; m++) {
        T1 * QuantilesResults=new T1[sizeQuantiles];
        T1 * ArrayFinToSort=new T1[CountNonZeroMask];
        int j=0;
        for (int i=0; i<numel; i++) {
            if (Mask==NULL || Mask[i]>10E-6) {
                ArrayFinToSort[j]=ArrayToAnalyse[m*numel+i];
                j++;
            }
        }
        float * TranscribedArrayFin=TranscribeArray<T1, float>(ArrayFinToSort, CountNonZeroMask);
        HeapSort(TranscribedArrayFin, CountNonZeroMask-1);
        for (int s=0; s<sizeQuantiles; s++) {
            QuantilesResults[s]=(T1) TranscribedArrayFin[IndicesToChoose[s]];
        }
        
        delete [] ArrayFinToSort;
        delete [] TranscribedArrayFin;
        ArrayFinToSort=NULL;
        TranscribedArrayFin=NULL;
        QuantilesVector.push_back(QuantilesResults);
    }
    delete [] IndicesToChoose;
    IndicesToChoose=NULL;
    return QuantilesVector;
}

template<typename T1,typename T2> T1 GetVarianceData(T1 * ArrayToMean, T2* Mask , int numel){
    T1 MeanData=GetMeanData(ArrayToMean, Mask, numel);
    T1 Normalisation=0;
    T1 Variance=0;
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            Variance+=Mask[i]*pow(ArrayToMean[i]-MeanData,2);
            Normalisation+=Mask[i];
        }
    }
    if (Normalisation>0) {
        return Variance/(Normalisation);
    }
    else{
        return 0;
    }
    
}

template<typename T1, typename T2> T1 *GetVarianceDataMulti(T1 * ArrayToVariance, T2 * Mask,int numel, int numbmodal){
    T1 * MeanData=GetMeanDataMulti(ArrayToVariance, Mask, numel, numbmodal);
    T1 Normalisation=0;
    T1 * Variance=new T1 [numbmodal];
    for (int m=0; m<numbmodal; m++) {
        Variance[m]=0;
    }
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            for (int m=0; m<numbmodal; m++) {
            
            Variance[m]+=Mask[i]*pow(ArrayToVariance[m*numel+i]-MeanData[m],2);
            }
            Normalisation+=Mask[i];
        }
    }
    if (Normalisation>0) {
        for (int m=0; m<numbmodal; m++) {
            Variance[m]/=Normalisation;
        }
    }
    else{
        for (int m=0; m<numbmodal; m++) {
            Variance[m]=0;
        }
    }
    delete [] MeanData;
    MeanData=NULL;
    return Variance;
}

template<typename T1,typename T2> T1 NCCCalculation(T1 * Array1, T1 * Array2, T2* Mask, int numel){
    float NCC=0;
    float Mean1=GetMeanData(Array1, Mask, numel);
    float Mean2=GetMeanData(Array2, Mask, numel);
    float Variance1=GetVarianceData(Array1, Mask, numel);
    float Variance2=GetVarianceData(Array2, Mask, numel);
    int CountMask=0;
    for (int i=0; i<numel; i++) {
        if (Mask[i]) {
            NCC+=(Array1[i]-Mean1)*(Array2[i]-Mean2);
            CountMask++;
        }
    }
    NCC/=(sqrtf(Variance1*Variance2))*CountMask;
    return NCC;
}

template<typename T1,typename T2> T1 GetSkewnessData(T1 * ArrayToMean, T2* Mask , int numel){
    T1 MeanData=GetMeanData(ArrayToMean, Mask, numel);
    T1 VarianceData=GetVarianceData(ArrayToMean, Mask, numel);
    T1 Normalisation=0;
    T1 Skewness=0;
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            Skewness+=Mask[i]*pow((ArrayToMean[i]-MeanData)/sqrt(VarianceData),3);
            Normalisation+=Mask[i];
        }
    }
    if (Normalisation>0) {
        return Skewness/Normalisation;
    }
    else{
        return 0;
    }
    
}

template<typename T1,typename T2> T1 GetKurtosisData(T1 * ArrayToMean, T2* Mask , int numel){
    T1 MeanData=GetMeanData(ArrayToMean, Mask, numel);
    T1 VarianceData=GetVarianceData(ArrayToMean, Mask, numel);
    T1 Normalisation=0;
    T1 Kurtosis=0;
    for (int i=0; i<numel; i++) {
        if (Mask[i]>0) {
            Kurtosis+=Mask[i]*pow((ArrayToMean[i]-MeanData)/sqrt(VarianceData),4);
            Normalisation+=Mask[i];
        }
    }
    if (Normalisation>0) {
        return Kurtosis/Normalisation-3;
    }
    else{
        return 0;
    }
}

template<typename T1,typename T2> T1 * GetStatisticData(T1 * ArrayToMean, T2* Mask , int numel){
    T1 * ArrayStats=new T1[9];
    ArrayStats[0]=GetMeanData(ArrayToMean, Mask, numel);
    ArrayStats[1]=GetVarianceData(ArrayToMean, Mask, numel);
    ArrayStats[2]=GetMin(ArrayToMean, Mask, numel);
    ArrayStats[3]=GetMax(ArrayToMean, Mask, numel);
    vector<float> Quantiles;
    Quantiles.push_back(0.5);
    Quantiles.push_back(0.25);
    Quantiles.push_back(0.75);
    vector<T1 *> QuantilesResults=GetQuantilesMulti(ArrayToMean, Mask, numel, 1, Quantiles);
    ArrayStats[4]=QuantilesResults[0][0];
    ArrayStats[5]=QuantilesResults[0][1];
    ArrayStats[6]=QuantilesResults[0][2];
    delete [] QuantilesResults[0];
    QuantilesResults[0]=NULL;
    
    ArrayStats[7]=GetKurtosisData(ArrayToMean, Mask, numel);
    ArrayStats[8]=GetSkewnessData(ArrayToMean, Mask, numel);
    return ArrayStats;
    
}

float GetMin(float * ArrayToMin,int numel);

template<typename T> void SaveTmpResult_short(vector<T *> ResultToSave, string FilenameToSave, nifti_image * ImageToUse, int * L2S){
    // Create Image to fill with ResultToSave
    nifti_image * ImageToFill=nifti_copy_nim_info(ImageToUse);
    ImageToFill->dim[0]=ResultToSave.size()>1?4:3;
    ImageToFill->dim[4]=ResultToSave.size();
    nifti_update_dims_from_array(ImageToFill);
    int numbClasses=ResultToSave.size();
    ImageToFill->data=(void *)calloc(ImageToFill->nvox, sizeof(float));
    float * Filling_PTR=static_cast<float *>(ImageToFill->data);
    int numel=ImageToFill->nx*ImageToFill->ny*ImageToFill->nz;
    int j=0;
    for (int i=0; i<numel; i++) {
        if (L2S[i]>=0) {
        for (int n=0; n<numbClasses; n++) {
                Filling_PTR[i+n*numel]=(float)ResultToSave[n][j];
            }
            j++;
        }
    }
    nifti_set_filenames(ImageToFill, FilenameToSave.c_str(), 0, 0);
    nifti_image_write(ImageToFill);
    nifti_image_free(ImageToFill);
    ImageToFill=NULL;
    return;
}
template<typename T> void SaveTmpResult(vector<T *> ResultToSave, string FilenameToSave, nifti_image * ImageToUse){
    // Create Image to fill with ResultToSave
    nifti_image * ImageToFill=nifti_copy_nim_info(ImageToUse);
    ImageToFill->dim[0]=ResultToSave.size()>1?4:3;
    ImageToFill->dim[4]=ResultToSave.size();
    nifti_update_dims_from_array(ImageToFill);
    int numbClasses=ResultToSave.size();
    ImageToFill->data=(void *)calloc(ImageToFill->nvox, sizeof(float));
    float * Filling_PTR=static_cast<float *>(ImageToFill->data);
    int numel=ImageToFill->nx*ImageToFill->ny*ImageToFill->nz;
    for (int i=0; i<numel; i++) {
        for (int n=0; n<numbClasses; n++) {
            Filling_PTR[i+n*numel]=(float)ResultToSave[n][i];
        }
    }
    nifti_set_filenames(ImageToFill, FilenameToSave.c_str(), 0, 0);
    nifti_image_write(ImageToFill);
    nifti_image_free(ImageToFill);
    ImageToFill=NULL;
    return;
}

template<typename T> bool CheckContainsNaN(T * ArrayToCheck, int numel){
    for (int i=0; i<numel; i++) {
        if (ArrayToCheck[i]!=ArrayToCheck[i]) {
            return 1;
        }
    }
    return 0;
}

template<typename T> bool CheckContainsInf(T * ArrayToCheck, int numel){
    for (int i=0; i<numel; i++) {
        if (fabs(ArrayToCheck[i])>1E18) {
            return 1;
        }
    }
    return 0;
}

template<typename T> T* ExtractTPFromNii(nifti_image * NiiToExtract, int tp){
    if (NiiToExtract==NULL) {
        return NULL;
    }
    int nD=NiiToExtract->nt*NiiToExtract->nu;
    if (tp>=nD) {
        return NULL;
    }
    int numel=NiiToExtract->nx*NiiToExtract->ny*NiiToExtract->nz;
    T* ResultArray=new T[numel];
    float * DataToExtract=static_cast<float *>(NiiToExtract->data);
    float * DataToCopy=&DataToExtract[numel*tp];
    for (int i=0; i<numel; i++) {
        ResultArray[i]=(T)DataToCopy[i];
    }
    return ResultArray;
}

template<typename T> T * ErosionTemplate(T * CountHistogram1, int KernelSize, vector<int> dim, bool Value){
    
    int numel=1;
    int numbDim=dim.size();
    for (int d=0; d<numbDim; d++) {
        numel*=dim[d];
    }
    // To do filtering, kernel must be at least of size 2.
    KernelSize=KernelSize<3?3:KernelSize;
    
    // Construction of the kernel
    float * Kernel=new float[KernelSize];
//    int kernelradius=KernelSize/2;
    for(int i =0;i<KernelSize;i++){
        Kernel[i]=1;
    }
    
    
    // Construction of the Shift array
    int *  Shift=new int[numbDim];
    for(int c=0;c<numbDim;c++){
        if(c>0){
            Shift[c]=Shift[c-1]*dim[c-1];
        }
        else{
            Shift[c]=1;
        }
    }
    
    
    // Copying the initial array to erode
    int SizeHistogram=Shift[numbDim-1]*dim[numbDim-1];
    bool * ErodedResult=new bool[SizeHistogram];
    for (int i=0; i<SizeHistogram; i++) {
        ErodedResult[i]=(bool)CountHistogram1[i];
    }
    
    for(int c=0;c<numbDim;c++){ // To do the blurring in each direction
        for(int j=0;j<SizeHistogram/(dim[c]);j++){
            int i=0;
            int * RemainDim=new int[numbDim-1];
            int * IndRemainDim=new int[numbDim-1];
            int * RemainShift=new int[numbDim-1];
            int * FinalShift=new int[numbDim-1];
            for(int i=0;i<numbDim-1;i++){
                RemainShift[i]=1;
                FinalShift[i]=1;
            }
            for(int d=0;d<numbDim;d++){
                if(d!=c){
                    RemainDim[i]=dim[d];
                    IndRemainDim[i]=d;
                    i++;
                }
            }
            for(int i=0;i<numbDim-1;i++){
                if(i>0){
                    RemainShift[i]=RemainShift[i-1]*RemainDim[i-1];
                }
            }
            int t=j;
            int tmp=0;
            int * Index=new int[numbDim-1];
            for(int i=numbDim-2;i>=0;i--){
                tmp=t/RemainShift[i];
                t=t-tmp*RemainShift[i];
                Index[i]=tmp;
            }
            for(int i=0;i<numbDim-1;i++){
                if(c<IndRemainDim[i]){
                    FinalShift[i]=dim[c]*RemainShift[i];
                }
                else{
                    FinalShift[i]=RemainShift[i];
                }
            }
            int FinalIndex=0;
            for(int i=0;i<numbDim-1;i++){
                FinalIndex+=Index[i]*FinalShift[i];
            }
            for(int i=0;i<dim[c];i++){
                int FinalIndexB=FinalIndex+i*Shift[c];
                if (ErodedResult[FinalIndexB]==Value) {
                    for(int k=0;k<KernelSize;k++){
                        if(i-KernelSize/2+k>=0 && i-KernelSize/2+k<dim[c]){ // Only when remains in image considered for kernel. Main pb is that the kernel is then not normalised anymore
                            if ((bool)CountHistogram1[FinalIndexB+(-KernelSize/2+k)*Shift[c]]==!Value) {
                                ErodedResult[FinalIndexB]=!Value;
                            }
                        }
                    }
                }
            }
            delete [] Index;
            Index=NULL;
            if(RemainShift!=NULL){
                delete [] RemainShift;
                RemainShift=NULL;
            }
            if(RemainDim!=NULL){
                delete [] RemainDim;
                RemainDim=NULL;
            }
            if(IndRemainDim!=NULL){
                delete [] IndRemainDim;
                IndRemainDim=NULL;
            }
            if(FinalShift!=NULL){
                delete [] FinalShift;
                FinalShift=NULL;
            }
        }
        
    }
    delete [] Kernel;
    Kernel=NULL;
    delete [] Shift;
    Shift=NULL;
    T * FinalEroded=TranscribeArray<bool, T>(ErodedResult, numel);
    delete [] ErodedResult;
    ErodedResult=NULL;
    return FinalEroded;
}

/**************************METHODS FOR FACES POINTS VECTOR CONVEXHULL ************************************/
vector <Face> ConvexHullFrom3D(int * ListCoord, int numbBorder, int* Dim, float *PixDim );
float * NormalisedVector(float * Vector);
bool InsidePolyhedron(deque<Face> Faces,  Point Element);
void ReorientBB(vector<Point>ToReorient);
vector<Point> CreatePointForBB(float *BB);
float DistanceFromFace(Face FaceDistance, Point Element);
void ReorientFace(Face FaceToReorient,Point PointRef);
int DifferentInAdjacent3(Face Face1, Face Face2);
bool IsAdjacent(Face Face1, Face Face2);
deque <vector<int> > AssociatePointsToFace(deque<Face> FacesInList, vector<int> PointsIndex, vector<float*> CoordinatesPoints);
float DistanceFromFace(Face FaceDistance, Point Element);
void InsidePolyhedronList(deque<Face>Faces,vector<float *> ListBorderFloat,vector<bool> &OnlineBool,int numbBorder,int& numberOut);
float * ProjectionToPlane(float *Point1, float* VectNorm, float* Point3);
float * ProjectionToLine(float * Point1, float* Point2, float*Point3);
float DistanceFromLineFloat(float * Point1, float* Point2, float* Point3);
float DistanceFromLine(int * Point1, int* Point2, int* Point3, float* PixDim);
float * CrossProduct(float * Vect1, float* Vect2);
float * CreatePVBoxFromBBPoint(vector<Point> ListBB, int* Dim, int* Shift, float* PixDim);
vector <Point> MMBBFromConvexHull(vector<Face> FaceVector);
Point ProjectedElement(Face ProjectedFace, Point Element);
vector <int> PossibleIndices(Point Element, int* Shift, float * PixDim);
vector <float> PossiblePartial(Point Element, int* Shift, float* PixDim);
float * ProbaInFromFaces(deque<Face> FaceList, vector<Point> ListToCheck, int * Dim, int* Shift, float* PixDim);
Face FindClosestFace(deque<Face> FaceList,Point Element);
vector<Point> GetListPointsFromBoolSeg(bool * BBSet, int* Dim, int* Shift, float* PixDim);
deque<Face> CreateListFacesBB(vector<Point> ListBB);
vector<Point> RectangleFromFace(Face FaceToDraw);
bool IsInRectangle(vector<Point> ListP, Point Element);
vector<Point> AssociatedFurthest(vector<Face> FaceVector);
vector<float> FurthestDistance(vector<Face>FaceVector, vector<Point> OppPoints);
float *GetExtremaPoints(vector<Face> FaceVector);
vector<Point> ListPointsFromFaces(vector<Face> FaceVector);
float * CreateSegPoints(vector<Point>ListP,int* Dim,int *Shift,float* PixDim);
float *CreateCharMM(bool *LesionBool,vector<Point> ListBB, int* Dim,int *Shift, float* PixDim);

/******************************* METHODS TO TAKE CARE OF ANGLE DEFINITION ********************************/
float * ReadAffineFromTextFile(string TreeTextFile);
mat44 ReadAffineFromTextFilemat44(string TreeTextFile);
void TranscribeMat44(mat44 M, float * TranscribedM);
mat44 TranscribeFloatToMat44(float * MToTranscribe);
void TransformedIndexToTemplateXYZ(mat44 * Transformation, int Index, int * Dim, float * xyz);
float ThetaAngleFromXYZ(float * XYZ);
float PhiAngleFromXYZ(float * XYZ);
float * ThetaAngleBorder(bool * MaskToTheta, mat44 * Transfo, int * Dim);
float * PhiAngleBorder(bool * MaskToPhi, mat44 * Transfo, int * Dim);
nifti_image * CreateQuadrantResult(nifti_image * SegLesBasis, SEG_ANALYSIS * segment_analysis, TreeEM * TreeToAnalyse);
nifti_image * QuadrantTransformation(float * TransfoTot,float * CoordinatesCenterMNI, nifti_image * SegLesBasis, int * L2S);
int * CodePlaneFromQuadrant(nifti_image * QuadrantResult);
bool * PlaneSelectionFromCode(int * PlaneCode, nifti_image * QuadrantCode, int DirFixed);
float CodeQuadrant(float * CoordinatesResult, float * CoordinatesCenter);
int * SidesFromQuadrant(nifti_image * QuadrantImage, int Dir);
float * AddingNeighbourLesions(float * WeightedLesionImage, float * MahalToUse,TreeEM *TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
void CorrectingSegmentation(nifti_image * SegToAnalyse, bool * CorrectionToWithdraw);
bool * VentricleGMBorderCorrection(nifti_image * SegToAnalyse, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
bool * ExternalGMBorderWM(nifti_image * SegToAnalyse, TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
bool * ExternalGMBorderCorrection(nifti_image * SegToAnalyse,TreeEM * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
bool * PotentialSPRegion(nifti_image * QuadrantImage, bool * VentricleSeg, bool * CGMBool, float DistanceThreshold, nifti_image * Mask );
void CorrectDataAccordingToModa(TreeEM * TreeToAnalyse, vector<int> Modalities);


/************ METHODS TO TAKE CARE OF HISTOGRAM MATCHING **********/
float * GetMaskHistogram(float * Data, bool * Mask, float Min, float Max, float sizeBin,int numbmodal,int numel);
vector<float> CompleteHistogramMatching(SEG_ANALYSIS * segment_analysis ,bool * Mask);

float * HistogramFitting(float * RefData, float * FloatData, int numel, int FittingOrder);

void PrintingCoeffsHistMatch(vector<float> CoeffsHist, SEG_ANALYSIS *segment_analysis, ostream& TxtFile);
float * ApplyingPolyfitPiecewise(vector<vector<float> >AffineCoeffs,float * BlurredLesion,vector<float>TFTV,int numel);
float * ApplyingPolyfitWT(vector<float> CoeffMatching, float * ImageToModify, float Threshold,int numel);
nifti_image * ApplyingPolyfitMatching(vector<float> CoeffMatching, nifti_image * ImageToModify,bool * MaskToApply);
TextureDescriptors1 * CreateTextureDescriptors1(bool * SegToUse, float * DataInit, TreeEM * TreeToAnalyse);
TextureDescriptors * CreateTextureDescriptors(bool * SegToUse, float * DataInit, vector<float *> ShiftedData, float* PixDimShift, TreeEM * TreeToAnalyse);
float * MakeHistogram(float * DataImage, int sizeSeg, int numbmodal);
float * GetGLMC(bool * SegToCheck, float * DataCorr, float * DataShifted, TreeEM * TreeToAnalyse);
float * GetParamFromHistogram(float * Histogram, vector<int> DimHist);





/***********METHODS to get LACUNES and OTHER OUTLIERS **************/

float * AddNormResp(vector<TreeEM*> VectorToAdd);
nifti_image * CodedOutliers(TreeEM  * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
nifti_image * CodedOutliersMask(TreeEM  * TreeToAnalyse, bool * OutlierMaskShort, SEG_ANALYSIS * segment_analysis);
nifti_image * OutliersOrigin(TreeEM  * TreeToAnalyse, SEG_ANALYSIS * segment_analysis);
float * OutlierSelection(TreeEM* TreeToAnalyse, float * OutlierCoding, float * OutlierOrigin, vector<float> CodeSelection, vector<float> OriginSelection, bool * RejectedZone);
bool CheckForFPDGM(float PropDGM, vector<float *> MahalVec, int numbmodal);
LesionSubcategories OutlierRefinement(Outlier * OutlierInit, SEG_ANALYSIS * segment_analysis);
vector<TreeEM *> SelectionPVClasses(TreeEM* TreeToAnalyse, int Index1, int Index2);
int * TreatingPreLabel(int *ComponentPV, bool *RejectedZonePV,int * Dim, int * Shift);
int * TreatingLabelFromSetDiscard(int * ComponentPV, vector<int> LabelsDiscarded, int numel);
vector<int> LabelToDiscard(int * ComponentPV, bool * RejectedZonePV, int * Dim, int * Shift);

float * CreateSphereAroundVox(int Index, int * Dim , int * Shift, float* PixDim, float Distance);
float * CreateCircleAroundVox(int Index, int * Dim, int * Shift, float * PixDim, float Distance, int DimFixed);
void ExtractWriteCoMComponent(nifti_image * ImageToLabel, int Neigh, string Filename);
nifti_image * FinalROIPlacement(vector<nifti_image *> VectorNiiEx, vector<nifti_image *> VectorNiiROI, nifti_image * ImageRef, nifti_image * ROIRef,nifti_image * Mask, string TxtFileCoMFinalString);
vector<int> ReadClassifFromTextFile(SEG_ANALYSIS * segment_analysis);
int ReadClassif(string line);
vector<Outlier *> ReadOutliersFromTextFile(SEG_ANALYSIS * segment_analysis, TreeEM * TreeToAnalyse);
Outlier * ReadOutlier(istream & text, TreeEM * TreeResult,SEG_ANALYSIS * segment_analysis);
vector<float *> LocalSummaryStats(nifti_image * WMNii, nifti_image * Lobes, nifti_image * Layers , nifti_image * MahalNii);
float * LocalSummaryLesion(nifti_image * Lesion, nifti_image * Quadrant, nifti_image * Layers);
float * LocalSummaryMahal(nifti_image * Lesion, nifti_image * Quadrant, nifti_image * Layers, nifti_image * MahalImage);
int * LocalSummaryLesionNumber(nifti_image * Lesion, nifti_image * Quadrant, nifti_image * Layers);
float * LocalSummaryLesion_bis(nifti_image * Lesion, nifti_image * Lobes, nifti_image * Layers);
float * LocalSummaryMahal_bis(nifti_image * Lesion, nifti_image * Lobes, nifti_image * Layers, nifti_image * MahalImage);
int * LocalSummaryLesionNumber_bis(nifti_image * Lesion, nifti_image * Lobes, nifti_image * Layers);
void PrintLocalSummary(float * LocalSummaryToPrint,ostream& TxtFileLS);
void PrintLocalSummaryStats(vector<float *> LocalSummaryToPrint,ostream& TxtFileLS);
void PrintLocalSummaryNumber(int * LocalSummaryToPrint, ostream& TxtFileLS);
void PrintLocalSummaryMahal(float * LocalSummaryToPrint, int numbMahal, ostream& TxtFileLS);
void PrintLocalSummaryTot(float * LocalSummaryVol, int * LocalSummaryNumber, float * LocalSummaryMahal, int numbMahal, ostream& TxtFileLS);
void PrintLocalSummaryTot_bis(float * LocalSummaryVol, int * LocalSummaryNumber, float * LocalSummaryMahal, int numbMahal, int numbTot, ostream& TxtFileLS);


/************* Vessels Heat Map ************/
bool * CreateSeedsForVessel(nifti_image * SeedHeatMap, nifti_image* MahalNii,float Thresh);
nifti_image * VesselFromSeeds(bool* Seeds,nifti_image* MahalNii,nifti_image* Mask,float * HeatData,float * VesselLong);
int *CreateBoolVesselFromSeed(int Index,int * Seeds,nifti_image * MahalNii,bool * Mask,float * HeatData,float * VesselLong);
int CountNeighbourBinary(int * Vessel,int index ,int * Dim,int * Shift,int Neigh);
bool * CreateSeedsForVessel_bis(nifti_image * SeedHeatMap,bool * SeedPrevious, nifti_image* MahalNii,float Thresh);
nifti_image * CreateCodePVSFromMahal(nifti_image * MahalNii);
/************** Random Noise Simulator **************/
nifti_image * RandomGaussianNoiseNii(vector<float> Mean, vector<float> Std, nifti_image * BasicImage, nifti_image * Mask);
float randn (float mu, float sigma);
float randu(float min, float max);
float * randuSeries(float min, float max, int numb);



//class mystreambuf: public std::streambuf
//{
//};
//mystreambuf nostreambuf;
//std::ostream nocout(&nostreambuf);
inline void  log_verb(log_level_t x, string msg, log_level_t GL=GLOBAL_LEVEL){
    if (x >= GL){
        cout << msg << endl;
        return;
    }

}

#endif // _SEG_ANALYSIS_H
