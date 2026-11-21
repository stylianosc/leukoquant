#pragma once

// [leukoquant patch] Must come before any NiftyReg/NiftySeg headers below:
// glibc's <dirent.h> declares its own DT_UNKNOWN (directory entry type),
// which collides with NiftyReg's nifti1.h #define DT_UNKNOWN (NIfTI
// datatype code) -- whichever one the preprocessor sees first "wins" and
// corrupts the other's declaration as a raw macro substitution. Getting
// dirent.h fully parsed first, before either nifti header, avoids this
// regardless of the NiftyReg-vs-NiftySeg include order below (which is
// controlled separately, for a different collision).
#include <dirent.h>

#include "rapidxml.hpp"
#include "rapidxml_print.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <cmath>
#include <cstdlib>

// [leukoquant patch] NiftyReg's nifti1_io.h and NiftySeg's own vendored
// copy of the same header share a guard name (added to NiftyReg's copy to
// stop them redefining the same structs when both are linked together --
// see niftilib/nifti1_io.h in the vendored NiftyReg source for details).
// Whichever one is #included FIRST is the one that actually runs; the
// second is just skipped. NiftyReg's copy also does real RNifti setup
// (RN_DEBUG, Rf_warning fallback macros, etc.) that NiftySeg's plain copy
// knows nothing about, so NiftyReg's needs to go first here -- if
// NiftySeg's _seg_tools.h (below) pulled its copy in first instead, that
// RNifti setup would be silently skipped and later fail to compile.
#ifdef _WIN32
#pragma push_macro("_WIN32")
#undef _WIN32
#endif
#pragma pop_macro("_WIN32")
#include "_reg_f3d2.h"
#include "_reg_aladin.h"
#include "_reg_aladin_sym.h"
#include "_reg_tools.h"
#include "AladinContent.h"

#include "_seg_tools.h"
#include "_seg_EM.h"



enum unsigned_int {LABEL,PROB,LABELSEG,PROBSEG,INTEN,INTENSORT,INTENMAP,LIKE,VEL};

using namespace std;
using namespace rapidxml;

typedef std::map <unsigned char, float> DataPointMap;
typedef std::pair <unsigned char, float> DataPointPair;

char * xmlReadFile(string filename);
xml_document<> * xmlReadFromCharPtr(char * filename);



//namespace patch
//{
//    template < typename T > std::string to_string( const T& n )
//    {
//        std::ostringstream stm ;
//        stm << n;
//        return stm.str() ;
//    }
//}


// Reference Counter
class RefCounter
{
private:
    int count;
public:
    void AddRef()
    {
        count++;
    }
    int Release()
    {
        return --count;
    }
    int getCount()
    {
        return this->count;
    }
};


// SMART Pointer

template <class T>
class smartPointer
{
private:
    T*    pData; // Generic pointer to be stored
    RefCounter* reference;
public:
    smartPointer();

    smartPointer(T* pValue);

    smartPointer(const smartPointer<T>& sp);

    T* cPTR();

    ~smartPointer();

    T& operator* ();

    T* operator-> ();

    void printCounter();

    smartPointer<T>& operator = (const smartPointer<T>& sp);
};

// Nifti_dataset class

class niftiObject
{
protected:
    string filename;

    nifti_image * file;
public:
    niftiObject();

    niftiObject(nifti_image * nifti_img);

    virtual ~niftiObject();

    virtual nifti_image * getPtr();

    string getFilename();

    void setFilename(string _filename);
};

class niftiDynamicObject : public niftiObject
{
public:

    niftiDynamicObject();

    niftiDynamicObject(nifti_image * nifti_img, bool _own);

    niftiDynamicObject(string filename,
                       bool read_data);


    nifti_image * getPtr();

    nifti_image * getPtrNoData();

    bool getLoadStatus();

    void load_data();

    void unload_data();

};

class geoObject : public niftiObject
{

public:

    geoObject();

    geoObject(nifti_image * nifti_img);

    void addIfSmaller(nifti_image * _input);

};


// Data structure to store the results dynamically.

class mapImage
{
protected:
    vector<DataPointMap *> Data;
public:
    mapImage(int size);

    mapImage(int _size,
             int * maskImage);

    ~mapImage();

    void addToPair(int index,
                   unsigned char label,
                   float value);

    void assignToPair(int index,
                      unsigned char label,
                      float value);

    void assignSortToPair(int index,
                          unsigned char maxItemsInSort,
                          float intensity,
                          float distance);

    float getValue(int index,
                   unsigned char label);

    float getValueOrNaN(int index,
                        unsigned char label);

    unsigned char getKeyWithMaxVal(int index);

    void mask(int * maskImage);
};


// Base database class

class Database
{
protected:
    vector< smartPointer<niftiDynamicObject> > datasets;
    string path;
    string description;
    bool initialised;

public:
    Database();
    virtual ~Database();

    virtual void init(string _folderpath,
                      string targetfilename);

    virtual void init(string _folderpath,
                      string _description,
                      string targetfilename);

    void setDescription(string _discription);

    string getDescription();

    string getPath();

    size_t getNumberOfDatasets();

    void unloadAllData(bool _verbose);

    void unloadDataFromDataset(unsigned int numb);

    nifti_image* getDatasetPtr(unsigned int numb);

    nifti_image *getDatasetNoDataPtr(unsigned int numb);

    string getDatasetName(unsigned int numb);

    void setDataset(unsigned int numb,
                    smartPointer<niftiDynamicObject> newDataset);

};


// Derived database class with extra properties

class infoDatabase : public Database
{
protected:
    string fname;
    string nt;
    vector<string> extrafiles;
    unsigned int type;
    smartPointer<niftiDynamicObject> target_nifti;
    int * target_nifti_int_mask;

public:

    using Database::init;
    void init(string _folderpath,
              smartPointer< Database>rootdb,
              string target_filename);

    void init(smartPointer<Database> rootdb);

    void setType(string _type);

    unsigned int getType();

    void setFname(string _fname);

    string getFname();

    void setExtraFiles(vector<string> _extrafiles);

    vector<string> getExtraFiles();

    int getNt();

    void unloadTarget();

    void setTarget(string _target_filename);

    nifti_image * getTargetPtr();

    void setMaskPtr(int * _mask_PTR);

    int * getMaskPtr();

    ~infoDatabase();
};


// Main GIF class

class GIF
{
protected:
    // Vars
    char * xmlFileCharPtr;
    smartPointer<Database> rootdb;
    smartPointer<infoDatabase> cppdb;
    vector< smartPointer<mapImage> > resultsInfoDB;
    vector< smartPointer<geoObject> > resultsGeoDB;
    vector< smartPointer<infoDatabase> > infoDB;
    vector< smartPointer<infoDatabase> > geoDB;
    vector< smartPointer<mat44> > affDB;

    smartPointer<niftiDynamicObject> target;
    smartPointer<niftiDynamicObject> targetLaplacian;
    smartPointer<niftiDynamicObject> groupMask;

    smartPointer<mat44> DatabaseToTargetAffine;
    smartPointer< xml_document<> > xml;
    string basePath;
    string resultsPath;
    string inputAffine;

    int * targetMask;
    int sFormAligned;
    bool saveTransformations;
    bool updateLab;
    bool saveGeo;
    bool regStopAfterAffine;
    int verbose;
    int ompj;
    float kernelstdDef;
    float kernelstdSsd;
    float kernelstdLncc;
    float kernelstdT1DTI;
    float weightLSSD;
    float weightDEF;
    float weightLNCC;
    float weightT1DTI;
    float temperature;
    unsigned char sortMaxItems;
    float sortBeta;

    int nrr_metric;
    bool regRunSecondLevel;
    float regBendingEnergyWeight;
    float regJacobianLogWeight;
    int regSpeed;
    int regPlatf;
    float probthresh;

    float segRFsmo;
    float segRFalpha;
    float segMRFbeta;
    float segMinNumIter;

    // General Functions
    void RegisterSubjectAffine(unsigned int subnumb);

    void RegisterSubjectAffineImprove(unsigned int subnumb);

    void RegisterSubjectNRR(unsigned int subnumb);

    void RegisterAllSubjectsAndGetAffine();

    void RegisterAllSubjectsAndGetCPP();

    void AverageAffine();

    void PropagateSubjectInformation(unsigned int subnumb);

    void PropagateAllSubjectsInformation();

    // Flow Functions
    void informationFlow(unsigned int databasenumb,
                         unsigned int subnumb,
                         nifti_image * DeformationField,
                         nifti_image * DistanceImage);

    void informationFlowLabel(unsigned int databasenumb,
                              unsigned int subnumb,
                              nifti_image * DeformationField,
                              nifti_image * DistanceImage);

    void informationFlowProb(unsigned int databasenumb,
                             unsigned int subnumb,
                             nifti_image * DeformationField,
                             nifti_image * DistanceImage);

    void informationFlowInten(unsigned int databasenumb,
                              unsigned int subnumb,
                              nifti_image * DeformationField,
                              nifti_image * DistanceImage);

    void informationFlowIntenSort(unsigned int databasenumb,
                                  unsigned int subnumb,
                                  nifti_image * DeformationField,
                                  nifti_image * DistanceImage);

    void informationFlowIntenMAP(unsigned int databasenumb,
                                   unsigned int subnumb,
                                   nifti_image *DeformationField);

    void informationFlowLikelihood(unsigned int databasenumb,
                                   unsigned int subnumb,
                                   nifti_image *DeformationField);

    void informationFlowVelocityField(unsigned int databasenumb,
                                      unsigned int subnumb,
                                      nifti_image *DeformationField,
                                      nifti_image * DistanceImage);

    // Output functions
    void processAndSaveResults();

    void processAndSaveResultingFlow(unsigned int databasenumb);

    void processAndSaveFlowAsLabel(unsigned int databasenumb);

    void processAndSaveFlowAsProb(unsigned int databasenumb);

    void processAndSaveFlowAsLabelSeg(unsigned int databasenumb);

    void processAndSaveFlowAsProbSeg(unsigned int databasenumb);

    void processAndSaveFlowAsInten(unsigned int databasenumb);

    void processAndSaveFlowAsIntenSort(unsigned int databasenumb);

    void processAndSaveFlowAsIntenMAP(unsigned int databasenumb);

    void processAndSaveFlowAsVel(unsigned int databasenumb);

    void postprocessNeuromorphometricsFix(float * ResultingLabelsPtr, nifti_image * ResultingLabels);

    bool postprocessXML(string saveFilename,
                             string baseXMLFilename,
                             vector< vector<int> > mapping,
                             nifti_image * BiasFieldCorrected,
                             nifti_image * TissueSegmentation, nifti_image *categoricalSegmentation,
                             mapImage * currentInfoDB);

    double postprocessEstimateLabelCNR(int label1,int label2,
                                       vector< vector<int> > mapping,
                                       nifti_image * BiasFieldCorrected,
                                       nifti_image * TissueSegmentation,
                                       mapImage * resultsInfoDB);

    double postprocessEstimateTissueCNR(int tissue1,int tissue2,
                                        nifti_image * BiasFieldCorrected,
                                        nifti_image * TissueSegmentation);

    double postprocessEstimateLabelSNR(int label, vector<vector<int> > mapping,
                                       nifti_image * BiasFieldCorrected,
                                       nifti_image * TissueSegmentation,
                                       mapImage * resultsInfoDB);

    double postprocessEstimateTissueSNR(int tissue,
                                        nifti_image * BiasFieldCorrected,
                                        nifti_image * TissueSegmentation);

    double postprocessEstimateLabelProbabilisticVolume(int label,
                                          vector<vector<int> > * mapping,
                                          nifti_image * TissueSegmentation,
                                          mapImage * resultsInfoDB);

    double postprocessEstimateLabelCategoricalVolume(int label,
                                          nifti_image * categoricalSegmentation);

    double postprocessEstimateTissueProbabilisticVolume(int tissue,
                                           nifti_image * TissueSegmentation);

    double postprocessEstimateTissueCategoricalVolume(int tissue,
                                           nifti_image * TissueSegmentation);

    // Tools -> Fetch data ptr (owned by the object)
    char getSformAlignedDatabase();

    bool getAndSaveTransformations();

    char getVerbose();

    nifti_image * getDatabaseTarget(unsigned int databaseNumb);

    // Tools -> Create data (owned by the caller)
    nifti_image * createResampledImageUsingDeformationField(nifti_image * ImageToResample,
                                                            nifti_image * TargetImage,
                                                            nifti_image * DeformationField,
                                                            unsigned int typeOfInterpolation);

    nifti_image * createTransportVectorUsingDeformationField(nifti_image *ImageToResample,
                                                             nifti_image *FloatingImage,
                                                             nifti_image *TargetImage,
                                                             nifti_image * DeformationField,
                                                             unsigned int typeOfInterpolation);

    nifti_image * createResampledImageUsingSForm(nifti_image * ImageToResample,
                                                 nifti_image * TargetImage,
                                                 unsigned int typeOfInterpolation);

    nifti_image * createResampledImageUsingAffine(nifti_image * ImageToResample,
                                                  nifti_image * TargetImage,
                                                  mat44 * InputAffine,
                                                  unsigned int typeOfInterpolation);

    nifti_image * createDeformationFieldFromGlobalTargetToDatabaseTarget(nifti_image * DeformationField,
                                                                         unsigned int databasenumb);

    nifti_image * createDeformationFieldFromSubjectToGeneralTarget(unsigned int subnumb);

    nifti_image * createDistanceFromSubject(unsigned int subnumb,
                                            nifti_image * Deformation);

    nifti_image * createDistanceF(unsigned int subnumb,
                                  nifti_image * Deformation);

    nifti_image * createDistanceL(unsigned int subnumb,
                                  nifti_image * Deformation);

    nifti_image * createDistanceLNCC(unsigned int subnumb,
                                     nifti_image * Deformation);

    nifti_image * createDistanceT1DTI(unsigned int subnumb,
                                      nifti_image *Deformation);

    nifti_image * createTargetMask();


    vector< vector<int> > readLabMappingXml(string _filename,
                                            unsigned char _maxsize);

    void allocateInfoGeoResStuctures();

    void freeMemory();


    void normaliseImage(nifti_image *img_to_normalise, nifti_image *reference_image, int * mask);

    void processAndSaveFlowAsLike(unsigned int databasenumb);

    void WriteAffineFile(mat44 *mat, const char *fileName);

    nifti_image *createLaplacian();

    void setRoiMask(string _filename);

    void updateTargetMaskWithGroupMask();

public:
    GIF(string _xmlFilename,
        string _target_filename);

    ~GIF();

    // Main control/execute function

    void Run();

    // Object setup
    void setTransformationsFolder(string _target_folder);

    void setRootTargetMask(string _target_filename);

    void setDatabaseTarget(string _database_fname,
                           string _target_filename);

    void setDatabaseTargetMask(string _database_fname,
                               string _target_filename);

    void setSformAlignedDatabase(char _sform_aligned);

    void setInputAffineFromDatabase(string _input_affine);

    void setResultsPath(string _path);

    void setVerbose(char _verbose);

    void set_kernelstd_def(float _value);

    void set_kernelstd_ssd(float _value);

    void set_kernelstd_lncc(float _value);

    void set_kernelstd_t1dti(float _value);

    void set_weight_def(float _value);

    void set_weight_ssd(float _value);

    void set_weight_lncc(float _value);

    void set_weight_t1dti(float _value);

    void set_temperature(float _value);

    void set_nrr_metric(char _value);

    void set_sort_beta(float _value);

    void set_sort_number(char _value);

    void set_regRunSecondLevel(bool _value);

    void set_regBendingEnergyWeight(float _value);

    void set_regJacobianLogWeight(float _value);

    void set_regSpeed(int _value);

    void set_regPlatf(int _value);

    void set_probthresh(float _value);

    void set_segRFsmo(float _value);

    void set_segRFalpha(float _value);

    void set_segMRFbeta(float _value);

    void set_segMinNumIter(float _value);

    void set_updateLab(bool _value);

    void set_saveGeo(bool _value);

    void set_regAffOnly(bool _value);

    void set_ompj(int ompj);

};

