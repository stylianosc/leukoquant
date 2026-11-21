#include <iostream>
#include <time.h>
#include <sstream>
#include <string>
#include <iostream>
#include "_seg_common.h"
#include "_seg_tools.h"

#include "_reg_f3d.h"
#include "_reg_globalTrans.h"
#include "_reg_localTrans.h"
#include "_reg_aladin_sym.h"
#include "_reg_tools.h"
#include "_reg_resampling.h"



using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s <FMRI_time_Series> <ventricle mask> <tmp_folder> <output.csv>\n\n",exec);
    printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

namespace patch
{
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}
}


float * subtractGroupwiseAndGetDistances(vector<mat44 *> affDB)
{

    // Percentage outliers
    double percent=0.2;


    cout << "\nFinding Group Average"<<endl;

    mat44 *inputMatrices=(mat44 *)malloc(affDB.size() * sizeof(mat44));
    mat44 averageMatrix;

    // All the input matrices are log-ed
    for(size_t m=0; m<affDB.size(); ++m)
    {
        if(affDB[m]!=NULL)
        {
            inputMatrices[m] = reg_mat44_logm(affDB[m]);
        }
    }
    float *distances=(float *)malloc(affDB.size()*3 * sizeof(float));
    float *weight=(float *)malloc(affDB.size() * sizeof(float));
    float *weight2=(float *)malloc(affDB.size() * sizeof(float));

    for(size_t m=0; m<affDB.size(); ++m)
    {
        weight[m]=1;
        weight2[m]=1;
    }

    // All the exponentiated matrices are summed up into one matrix
    //temporary double are used to avoid error accumulation

    for(int iter=0; iter<10; iter++)
    {
        double tempValue[16]= {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,0
                              };
        float sumdistance=0;
        for(size_t m=0; m<affDB.size(); ++m)
        {
            if(affDB[m]!=NULL)
            {
                tempValue[ 0]+=weight[m]*(double)inputMatrices[m].m[0][0];
                tempValue[ 1]+=weight[m]*(double)inputMatrices[m].m[0][1];
                tempValue[ 2]+=weight[m]*(double)inputMatrices[m].m[0][2];
                tempValue[ 3]+=weight[m]*(double)inputMatrices[m].m[0][3];
                tempValue[ 4]+=weight[m]*(double)inputMatrices[m].m[1][0];
                tempValue[ 5]+=weight[m]*(double)inputMatrices[m].m[1][1];
                tempValue[ 6]+=weight[m]*(double)inputMatrices[m].m[1][2];
                tempValue[ 7]+=weight[m]*(double)inputMatrices[m].m[1][3];
                tempValue[ 8]+=weight[m]*(double)inputMatrices[m].m[2][0];
                tempValue[ 9]+=weight[m]*(double)inputMatrices[m].m[2][1];
                tempValue[10]+=weight[m]*(double)inputMatrices[m].m[2][2];
                tempValue[11]+=weight[m]*(double)inputMatrices[m].m[2][3];
                tempValue[12]+=weight[m]*(double)inputMatrices[m].m[3][0];
                tempValue[13]+=weight[m]*(double)inputMatrices[m].m[3][1];
                tempValue[14]+=weight[m]*(double)inputMatrices[m].m[3][2];
                tempValue[15]+=weight[m]*(double)inputMatrices[m].m[3][3];
                sumdistance+=weight[m];
            }
        }
        // Average matrix is computed
        tempValue[ 0] /= (double)sumdistance;
        tempValue[ 1] /= (double)sumdistance;
        tempValue[ 2] /= (double)sumdistance;
        tempValue[ 3] /= (double)sumdistance;
        tempValue[ 4] /= (double)sumdistance;
        tempValue[ 5] /= (double)sumdistance;
        tempValue[ 6] /= (double)sumdistance;
        tempValue[ 7] /= (double)sumdistance;
        tempValue[ 8] /= (double)sumdistance;
        tempValue[ 9] /= (double)sumdistance;
        tempValue[10] /= (double)sumdistance;
        tempValue[11] /= (double)sumdistance;
        tempValue[12] /= (double)sumdistance;
        tempValue[13] /= (double)sumdistance;
        tempValue[14] /= (double)sumdistance;
        tempValue[15] /= (double)sumdistance;

        // The final matrix is exponentiated
        averageMatrix.m[0][0]=(float)tempValue[ 0];
        averageMatrix.m[0][1]=(float)tempValue[ 1];
        averageMatrix.m[0][2]=(float)tempValue[ 2];
        averageMatrix.m[0][3]=(float)tempValue[ 3];
        averageMatrix.m[1][0]=(float)tempValue[ 4];
        averageMatrix.m[1][1]=(float)tempValue[ 5];
        averageMatrix.m[1][2]=(float)tempValue[ 6];
        averageMatrix.m[1][3]=(float)tempValue[ 7];
        averageMatrix.m[2][0]=(float)tempValue[ 8];
        averageMatrix.m[2][1]=(float)tempValue[ 9];
        averageMatrix.m[2][2]=(float)tempValue[10];
        averageMatrix.m[2][3]=(float)tempValue[11];
        averageMatrix.m[3][0]=(float)tempValue[12];
        averageMatrix.m[3][1]=(float)tempValue[13];
        averageMatrix.m[3][2]=(float)tempValue[14];
        averageMatrix.m[3][3]=(float)tempValue[15];
        // Free the array containing the input matrices

        for(size_t tp=0; tp<affDB.size(); ++tp)
        {
            weight[tp]=10000;
            if(affDB[tp]!=NULL)
            {

                mat44 MinusGroupwise=reg_mat44_minus(&(inputMatrices[tp]),&averageMatrix);

                mat44 MinusPairwise;
                mat44 MinusTP1;
                if(tp>0){
                    if(iter==0){
                    MinusTP1=reg_mat44_minus(&(inputMatrices[tp]),&(inputMatrices[0]));
                    }
                    MinusPairwise=reg_mat44_minus(&(inputMatrices[tp]),&(inputMatrices[tp-1]));
                }
                else{
                    for(int i=0; i<4; ++i)
                    {
                        for(int j=0; j<4; ++j)
                        {
                            if(iter==0){MinusTP1.m[i][j]=0;}
                            MinusPairwise.m[i][j]=0;
                        }
                    }
                }


                mat44 MinusGroupwise_transpose;
                mat44 MinusTP1_transpose;
                mat44 MinusPairwise_transpose;
                for(int i=0; i<4; ++i)
                {
                    for(int j=0; j<4; ++j)
                    {
                        MinusGroupwise_transpose.m[i][j]=MinusGroupwise.m[j][i];
                        if(iter==0){
                            MinusTP1_transpose.m[i][j]=MinusTP1.m[j][i];
                        }
                        MinusPairwise_transpose.m[i][j]=MinusPairwise.m[j][i];
                    }
                }
                mat44 MTMGroupwise=reg_mat44_mul(&MinusGroupwise_transpose,&MinusGroupwise);
                mat44 MTMTP1;
                if(iter==0){
                    MTMTP1=reg_mat44_mul(&MinusTP1_transpose,&MinusTP1);
                }
                mat44 MTMPairwise=reg_mat44_mul(&MinusPairwise_transpose,&MinusPairwise);

                double traceGroupwise=0;
                double traceTP1=0;
                double tracePairwise=0;
                for(size_t i=0; i<4; ++i)
                {
                    traceGroupwise+=MTMGroupwise.m[i][i];
                    if(iter==0){
                    traceTP1+=MTMTP1.m[i][i];
                    }
                    tracePairwise+=MTMPairwise.m[i][i];
                }
                distances[tp]=(traceGroupwise);
                if(iter==0){
                distances[tp+affDB.size()]=(traceTP1);
                }
                distances[tp+2*affDB.size()]=(tracePairwise);
                weight[tp]=1/(sqrt(traceGroupwise));
                weight2[tp]=1/(sqrt(traceGroupwise));
            }
        }


        quickSort(weight2,affDB.size());
        for(size_t m=0; m<affDB.size(); ++m)
        {
            weight[m]=weight[m]>=weight2[(int)ceil((double)(affDB.size())*percent)];
        }


        // Demean all the input affine matrices
        for(size_t i=0; i<affDB.size(); ++i){
            inputMatrices[i] = inputMatrices[i] - averageMatrix;
        }
    }



    // Put the demeaned matrices back by exponentiating and copying
    for(size_t m=0; m<affDB.size(); ++m){
        inputMatrices[m] = reg_mat44_expm(&inputMatrices[m]);
    }
    for(size_t m=0; m<affDB.size(); ++m)
    {
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                affDB[m]->m[i][j]=inputMatrices[m].m[i][j];
            }
        }
    }

    free(inputMatrices);

    return distances;
}

int main(int argc, char **argv)
{


    if (argc < 4)
    {
        Usage(argv[0]);
        return 0;
    }
    if(strcmp(argv[1], "-help")==0 || strcmp(argv[1], "-Help")==0 ||
            strcmp(argv[1], "-HELP")==0 || strcmp(argv[1], "-h")==0 ||
            strcmp(argv[1], "--h")==0 || strcmp(argv[1], "--help")==0)
    {
        Usage(argv[0]);
        return 0;
    }


    nifti_image * InputImages=nifti_image_read(argv[1],true);
    if(InputImages!=NULL)
    {
        if(InputImages->datatype!=NIFTI_TYPE_FLOAT32)
            seg_changeDatatype<float>(InputImages);
        nifti_update_dims_from_array(InputImages);
    }
    else{
        std::cout<<"Could not load FMRI data"<<std::endl;
    }
    float * InputImagesDataPtr = static_cast<float *>(InputImages->data);



    nifti_image * Ventricles=nifti_image_read(argv[2],true);
    if(Ventricles!=NULL)
    {
        if(Ventricles->datatype!=NIFTI_TYPE_FLOAT32)
            seg_changeDatatype<float>(Ventricles);
        nifti_update_dims_from_array(Ventricles);
    }
    else{
        std::cout<<"Could not load Ventricles mask"<<std::endl;
    }
    float * VentriclesDataPtr = static_cast<float *>(Ventricles->data);


    string folder=argv[3];

    string resultFile=argv[4];

    std::vector< std::vector<float> > CSV_values;

    CSV_values.resize(InputImages->nt);
    for (int i = 0; i < InputImages->nt; ++i)
      CSV_values[i].resize(6);


    // Read in the current groupwise
    nifti_image * CurrentGroupwise;
    string groupwiseFilename=folder+"/groupwiseMean.nii.gz";
    FILE *groupwiseFile;
    if( (bool)(groupwiseFile=fopen(groupwiseFilename.c_str(), "r")))
    {
        fclose(groupwiseFile);
        CurrentGroupwise=nifti_image_read(groupwiseFilename.c_str(),true);
        if(CurrentGroupwise!=NULL)
        {
            if(CurrentGroupwise->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<float>(CurrentGroupwise);
            nifti_update_dims_from_array(CurrentGroupwise);
        }
        else{
            std::cout<<"Could not load CurrentGroupwise data"<<std::endl;
        }
    }
    // If it does not exist, assign it to the first time point of the FMRI time series
    else{
        CurrentGroupwise=nifti_copy_nim_info(InputImages);
        CurrentGroupwise->dim[0]=3;
        CurrentGroupwise->dim[4]=CurrentGroupwise->nt=0;
        nifti_update_dims_from_array(CurrentGroupwise);
        CurrentGroupwise->datatype = NIFTI_TYPE_FLOAT32;
        CurrentGroupwise->nbyper = sizeof(float);
        CurrentGroupwise->data=(void *)malloc(CurrentGroupwise->nvox * CurrentGroupwise->nbyper);
        float *CurrentGroupwiseDataPtr = static_cast<float *>(CurrentGroupwise->data);
        memcpy(CurrentGroupwiseDataPtr, &InputImagesDataPtr[0], CurrentGroupwise->nvox * CurrentGroupwise->nbyper);
        nifti_update_dims_from_array(CurrentGroupwise);
        nifti_set_filenames(CurrentGroupwise,groupwiseFilename.c_str(),0,0);
    }



    // Allocate space for the next Groupwise result
    nifti_image * NextGroupwiseMean=nifti_copy_nim_info(InputImages);
    NextGroupwiseMean->dim[0]=3;
    NextGroupwiseMean->dim[4]=NextGroupwiseMean->nt=0;
    nifti_update_dims_from_array(NextGroupwiseMean);
    NextGroupwiseMean->datatype = NIFTI_TYPE_FLOAT32;
    NextGroupwiseMean->nbyper = sizeof(float);
    NextGroupwiseMean->data=(void *)malloc(NextGroupwiseMean->nvox * NextGroupwiseMean->nbyper);
    float *NextGroupwiseMeanDataPtr = static_cast<float *>(NextGroupwiseMean->data);
    nifti_update_dims_from_array(NextGroupwiseMean);
    nifti_set_filenames(NextGroupwiseMean,groupwiseFilename.c_str(),0,0);
    for(int i=0; i<NextGroupwiseMean->nvox;i++){NextGroupwiseMeanDataPtr[i]=0;}

    // Allocate space for the next Groupwise result
    string groupwiseVarFilename=folder+"/groupwiseVariance.nii.gz";
    nifti_image * NextGroupwiseVariance=nifti_copy_nim_info(InputImages);
    NextGroupwiseVariance->dim[0]=3;
    NextGroupwiseVariance->dim[4]=NextGroupwiseVariance->nt=0;
    nifti_update_dims_from_array(NextGroupwiseVariance);
    NextGroupwiseVariance->datatype = NIFTI_TYPE_FLOAT32;
    NextGroupwiseVariance->nbyper = sizeof(float);
    NextGroupwiseVariance->data=(void *)malloc(NextGroupwiseVariance->nvox * NextGroupwiseVariance->nbyper);
    float *NextGroupwiseVarianceDataPtr = static_cast<float *>(NextGroupwiseVariance->data);
    nifti_update_dims_from_array(NextGroupwiseVariance);
    nifti_set_filenames(NextGroupwiseVariance,groupwiseVarFilename.c_str(),0,0);
    for(int i=0; i<NextGroupwiseVariance->nvox;i++){NextGroupwiseVarianceDataPtr[i]=0;}


    // Allocate space for the resampled data after the groupwise
    string ressampledFilename=folder+"/resampledInput.nii.gz";
    nifti_image * RessampledInputImages=nifti_copy_nim_info(InputImages);
    nifti_update_dims_from_array(RessampledInputImages);
    RessampledInputImages->datatype = NIFTI_TYPE_FLOAT32;
    RessampledInputImages->nbyper = sizeof(float);
    nifti_update_dims_from_array(RessampledInputImages);
    RessampledInputImages->data=(void *)malloc(RessampledInputImages->nvox * RessampledInputImages->nbyper);
    float *RessampledInputImagesDataPtr = static_cast<float *>(RessampledInputImages->data);
    nifti_set_filenames(RessampledInputImages,ressampledFilename.c_str(),0,0);



    // Allocate space for the resampled data after the groupwise
    string ressampledErrorFilename=folder+"/errorInInput.nii.gz";
    nifti_image * ErrorInInputImages=nifti_copy_nim_info(InputImages);
    nifti_update_dims_from_array(ErrorInInputImages);
    ErrorInInputImages->datatype = NIFTI_TYPE_FLOAT32;
    ErrorInInputImages->nbyper = sizeof(float);
    nifti_update_dims_from_array(ErrorInInputImages);
    ErrorInInputImages->data=(void *)malloc(ErrorInInputImages->nvox * ErrorInInputImages->nbyper);
    float *ErrorInInputImagesDataPtr = static_cast<float *>(ErrorInInputImages->data);
    nifti_set_filenames(ErrorInInputImages,ressampledErrorFilename.c_str(),0,0);




    // Allocate the Affine Database
    vector< mat44 * > affDB;


    // Find or load all the affines
    for(int tp=0; tp<InputImages->nt;tp++){

        // Get current time point
        nifti_image * CurrentTpImage=nifti_copy_nim_info(InputImages);
        CurrentTpImage->dim[0]=3;
        CurrentTpImage->dim[4]=CurrentTpImage->nt=0;
        nifti_update_dims_from_array(CurrentTpImage);
        CurrentTpImage->datatype = NIFTI_TYPE_FLOAT32;
        CurrentTpImage->nbyper = sizeof(float);
        CurrentTpImage->data=(void *)malloc(CurrentTpImage->nvox * CurrentTpImage->nbyper);
        float *CurrentTpImageDataPtr = static_cast<float *>(CurrentTpImage->data);
        memcpy(CurrentTpImageDataPtr, &InputImagesDataPtr[tp*CurrentTpImage->nvox], CurrentTpImage->nvox * CurrentTpImage->nbyper);
        nifti_update_dims_from_array(CurrentTpImage);
        string filenameTp=folder+"/"+patch::to_string(tp)+".nii.gz";
        nifti_set_filenames(CurrentTpImage,filenameTp.c_str(),0,0);
        //nifti_image_write(CurrentTpImage);

        // Allocate and set name of the affine file
        affDB.push_back(new mat44);
        string aff_name=folder+"/"+patch::to_string(tp)+".txt";


        // Check if affine exists
        FILE *aff;
        if( (bool)(aff=fopen(aff_name.c_str(), "r")))
        {
            fclose(aff);
            reg_tool_ReadAffineFile(affDB[tp],(char *)aff_name.c_str());
            cout<<"Aligning TP "<<tp<<" - Alredy registered"<<endl;
        }
        //Otherwise run affine
        else{

            std::cout<<"Aligning TP "<<tp;
            std::cout.flush();
            // Set up and run Affine
            reg_aladin_sym<float> *REG_AFFINE=new reg_aladin_sym<float>;
            REG_AFFINE->SetInputReference(CurrentGroupwise);
            REG_AFFINE->SetInputFloating(CurrentTpImage);
            //REG_AFFINE->SetInputMask();
            REG_AFFINE->SetMaxIterations(10);
            REG_AFFINE->SetNumberOfLevels(3);
            REG_AFFINE->SetLevelsToPerform(3);
            REG_AFFINE->SetReferenceSigma(0.0f);
            REG_AFFINE->SetFloatingSigma(0.0f);
            REG_AFFINE->SetAlignCentre(1);
            REG_AFFINE->SetPerformRigid(1);
            REG_AFFINE->SetPerformAffine(0);
            REG_AFFINE->SetBlockPercentage(80);
            REG_AFFINE->SetInlierLts(80);
            REG_AFFINE->SetInterpolation(1);
            REG_AFFINE->SetVerbose(0);
            REG_AFFINE->setGpuIdx(999);
            REG_AFFINE->Run();
            std::cout<<" - Done"<<std::endl;
            mat44 * affine_tmp=REG_AFFINE->GetTransformationMatrix();

            for(int i=0; i<4; i++)
            {
                for(int j=0; j<4; j++)
                {
                    affDB[tp]->m[i][j]=affine_tmp->m[i][j];
                }
            }
            reg_mat44_disp(affDB[tp],(char *)(string("Saving TP "+patch::to_string(tp)).c_str()));
            reg_tool_WriteAffineFile(affDB[tp],aff_name.c_str());
        }



        nifti_image_free(CurrentTpImage);

    }

    // Demean all the affines
    float * distances=subtractGroupwiseAndGetDistances(affDB);



    // Write distances to file
    for(int tp=0; tp<InputImages->nt;tp++){
        CSV_values[tp][0]=distances[tp];
        CSV_values[tp][1]=distances[tp+InputImages->nt];
        CSV_values[tp][2]=distances[tp+2*InputImages->nt];
    }



    // Resamples images and create mean




    for(int tp=0; tp<InputImages->nt;tp++){

        // Get current time point
        nifti_image * CurrentTpImage=nifti_copy_nim_info(InputImages);
        CurrentTpImage->dim[0]=3;
        CurrentTpImage->dim[4]=CurrentTpImage->nt=0;
        nifti_update_dims_from_array(CurrentTpImage);
        CurrentTpImage->datatype = NIFTI_TYPE_FLOAT32;
        CurrentTpImage->nbyper = sizeof(float);
        CurrentTpImage->data=(void *)malloc(CurrentTpImage->nvox * CurrentTpImage->nbyper);
        float *CurrentTpImageDataPtr = static_cast<float *>(CurrentTpImage->data);
        memcpy(CurrentTpImageDataPtr, &InputImagesDataPtr[tp*CurrentTpImage->nvox], CurrentTpImage->nvox * CurrentTpImage->nbyper);
        nifti_update_dims_from_array(CurrentTpImage);
        string filenameTp=folder+"/"+patch::to_string(tp)+".nii.gz";
        nifti_set_filenames(CurrentTpImage,filenameTp.c_str(),0,0);
        //nifti_image_write(CurrentTpImage);

        std::cout<<"Ressampling ["<<tp<<"]"<<std::endl;

        // Resample the image using the updated affine
        // First ger the deformation field
        nifti_image *defField = nifti_copy_nim_info(CurrentGroupwise);
        defField->dim[0]=defField->ndim=5;
        defField->dim[1]=defField->nx=CurrentGroupwise->nx;
        defField->dim[2]=defField->ny=CurrentGroupwise->ny;
        defField->dim[3]=defField->nz=CurrentGroupwise->nz;
        defField->dim[4]=defField->nt=1;
        defField->pixdim[4]=defField->dt=1.0;
        defField->dim[5]=defField->nu=CurrentGroupwise->nz>1?3:2;
        defField->dim[6]=defField->nv=1;
        defField->dim[7]=defField->nw=1;
        defField->nvox =(size_t)defField->nx*defField->ny*defField->nz*defField->nt*defField->nu;
        defField->scl_slope=1.f;
        defField->scl_inter=0.f;
        defField->datatype = NIFTI_TYPE_FLOAT32;
        defField->nbyper = sizeof(float);
        defField->data = (void *)calloc(defField->nvox, defField->nbyper);
        reg_affine_getDeformationField(affDB[tp],defField,false,NULL);


        // Then, allocate the resampled image
        nifti_image *warpedImage = nifti_copy_nim_info(CurrentGroupwise);
        warpedImage->dim[0]=warpedImage->ndim=CurrentTpImage->dim[0];
        warpedImage->dim[4]=warpedImage->nt=CurrentTpImage->dim[4];
        warpedImage->cal_min=CurrentTpImage->cal_min;
        warpedImage->cal_max=CurrentTpImage->cal_max;
        warpedImage->scl_slope=CurrentTpImage->scl_slope;
        warpedImage->scl_inter=CurrentTpImage->scl_inter;
        reg_tools_changeDatatype<float>(CurrentTpImage);
        warpedImage->nbyper = CurrentTpImage->nbyper;
        warpedImage->nvox = (size_t)warpedImage->dim[1] * (size_t)warpedImage->dim[2] *
                (size_t)warpedImage->dim[3] * (size_t)warpedImage->dim[4];
        warpedImage->data = (void *)calloc(warpedImage->nvox, warpedImage->nbyper);
        float *warpedImageDataPtr = static_cast<float *>(warpedImage->data);

        // actually resample the data
        reg_resampleImage(CurrentTpImage,warpedImage,defField,NULL,1,0);



        for(int i=0; i<NextGroupwiseMean->nvox;i++){
            NextGroupwiseMeanDataPtr[i]+=warpedImageDataPtr[i];
            NextGroupwiseVarianceDataPtr[i]+=warpedImageDataPtr[i]*warpedImageDataPtr[i];
            RessampledInputImagesDataPtr[i+tp*CurrentTpImage->nvox]+=warpedImageDataPtr[i];
        }

        nifti_image_free(CurrentTpImage);
        nifti_image_free(warpedImage);
        nifti_image_free(defField);
    }



    // Estimate the mean and variance
    for(int i=0; i<NextGroupwiseMean->nvox;i++){
        NextGroupwiseMeanDataPtr[i]/=InputImages->nt;
        NextGroupwiseVarianceDataPtr[i]=NextGroupwiseVarianceDataPtr[i]/InputImages->nt-powf(NextGroupwiseMeanDataPtr[i],2);

    }
    nifti_image_write(RessampledInputImages);
    nifti_image_write(NextGroupwiseMean);
    nifti_image_write(NextGroupwiseVariance);


    nifti_image *Mask=nifti_copy_nim_info(NextGroupwiseMean);
    Mask->dim[0]=3;
    Mask->dim[4]=Mask->nt=0;
    nifti_update_dims_from_array(Mask);
    Mask->datatype = NIFTI_TYPE_FLOAT32;
    Mask->nbyper = sizeof(float);
    Mask->data=(void *)malloc(Mask->nvox * Mask->nbyper);
    float *MaskDataPtr = static_cast<float *>(Mask->data);
    memcpy(MaskDataPtr, &NextGroupwiseMeanDataPtr[0], Mask->nvox * Mask->nbyper);
    nifti_update_dims_from_array(CurrentGroupwise);
    string filenameMask=folder+"/brainMask.nii.gz";
    nifti_set_filenames(Mask,filenameMask.c_str(),0,0);
    ImageSize * CurrSize = new ImageSize [1]();
    CurrSize->numel=(long)(Mask->nx*Mask->ny*Mask->nz);
    CurrSize->xsize=Mask->nx;
    CurrSize->ysize=Mask->ny;
    CurrSize->zsize=Mask->nz;
    CurrSize->usize=(Mask->nu>1)?Mask->nu:1;
    CurrSize->tsize=(Mask->nt>1)?Mask->nt:1;
    otsu(MaskDataPtr,NULL, CurrSize);

    nifti_image_write(Mask);


    // Resamples Mean and Variance Back and create mean
    for(int tp=0; tp<InputImages->nt;tp++){

        // Get current time point
        nifti_image * CurrentTpImage=nifti_copy_nim_info(InputImages);
        CurrentTpImage->dim[0]=3;
        CurrentTpImage->dim[4]=CurrentTpImage->nt=0;
        nifti_update_dims_from_array(CurrentTpImage);
        CurrentTpImage->datatype = NIFTI_TYPE_FLOAT32;
        CurrentTpImage->nbyper = sizeof(float);
        CurrentTpImage->data=(void *)malloc(CurrentTpImage->nvox * CurrentTpImage->nbyper);
        float *CurrentTpImageDataPtr = static_cast<float *>(CurrentTpImage->data);
        memcpy(CurrentTpImageDataPtr, &InputImagesDataPtr[tp*CurrentTpImage->nvox], CurrentTpImage->nvox * CurrentTpImage->nbyper);
        nifti_update_dims_from_array(CurrentTpImage);
        string filenameTp=folder+"/"+patch::to_string(tp)+".nii.gz";
        nifti_set_filenames(CurrentTpImage,filenameTp.c_str(),0,0);
        //nifti_image_write(CurrentTpImage);

        std::cout<<"Estimating Error ["<<tp<<"]"<<std::endl;

        // Resample the image using the updated affine
        // First ger the deformation field
        nifti_image *defField = nifti_copy_nim_info(CurrentGroupwise);
        defField->dim[0]=defField->ndim=5;
        defField->dim[1]=defField->nx=CurrentGroupwise->nx;
        defField->dim[2]=defField->ny=CurrentGroupwise->ny;
        defField->dim[3]=defField->nz=CurrentGroupwise->nz;
        defField->dim[4]=defField->nt=1;
        defField->pixdim[4]=defField->dt=1.0;
        defField->dim[5]=defField->nu=CurrentGroupwise->nz>1?3:2;
        defField->dim[6]=defField->nv=1;
        defField->dim[7]=defField->nw=1;
        defField->nvox =(size_t)defField->nx*defField->ny*defField->nz*defField->nt*defField->nu;
        defField->scl_slope=1.f;
        defField->scl_inter=0.f;
        defField->datatype = NIFTI_TYPE_FLOAT32;
        defField->nbyper = sizeof(float);
        defField->data = (void *)calloc(defField->nvox, defField->nbyper);
        mat44 inverseAFF = reg_mat44_inv(affDB[tp]);
        reg_affine_getDeformationField(&inverseAFF,defField,false,NULL);


        // Then, allocate the resampled mean image
        nifti_image *warpedMean = nifti_copy_nim_info(CurrentGroupwise);
        warpedMean->dim[0]=warpedMean->ndim=CurrentTpImage->dim[0];
        warpedMean->dim[4]=warpedMean->nt=CurrentTpImage->dim[4];
        warpedMean->cal_min=CurrentTpImage->cal_min;
        warpedMean->cal_max=CurrentTpImage->cal_max;
        warpedMean->scl_slope=CurrentTpImage->scl_slope;
        warpedMean->scl_inter=CurrentTpImage->scl_inter;
        reg_tools_changeDatatype<float>(CurrentTpImage);
        warpedMean->nbyper = CurrentTpImage->nbyper;
        warpedMean->nvox = (size_t)warpedMean->dim[1] * (size_t)warpedMean->dim[2] *
                (size_t)warpedMean->dim[3] * (size_t)warpedMean->dim[4];
        warpedMean->data = (void *)calloc(warpedMean->nvox, warpedMean->nbyper);
        float *warpedMeanDataPtr = static_cast<float *>(warpedMean->data);

        // actually resample the mean data
        reg_resampleImage(NextGroupwiseMean,warpedMean,defField,NULL,1,0);


        // Then, allocate the resampled var image
        nifti_image *warpedVariance = nifti_copy_nim_info(CurrentGroupwise);
        warpedVariance->dim[0]=warpedVariance->ndim=CurrentTpImage->dim[0];
        warpedVariance->dim[4]=warpedVariance->nt=CurrentTpImage->dim[4];
        warpedVariance->cal_min=CurrentTpImage->cal_min;
        warpedVariance->cal_max=CurrentTpImage->cal_max;
        warpedVariance->scl_slope=CurrentTpImage->scl_slope;
        warpedVariance->scl_inter=CurrentTpImage->scl_inter;
        reg_tools_changeDatatype<float>(CurrentTpImage);
        warpedVariance->nbyper = CurrentTpImage->nbyper;
        warpedVariance->nvox = (size_t)warpedVariance->dim[1] * (size_t)warpedVariance->dim[2] *
                (size_t)warpedVariance->dim[3] * (size_t)warpedVariance->dim[4];
        warpedVariance->data = (void *)calloc(warpedVariance->nvox, warpedVariance->nbyper);
        float *warpedVarianceDataPtr = static_cast<float *>(warpedVariance->data);

        // actually resample the var data
        reg_resampleImage(NextGroupwiseVariance,warpedVariance,defField,NULL,1,0);



        // Then, allocate the resampled ventricle mask image
        nifti_image *warpedVentricles = nifti_copy_nim_info(CurrentGroupwise);
        warpedVentricles->dim[0]=warpedVentricles->ndim=CurrentTpImage->dim[0];
        warpedVentricles->dim[4]=warpedVentricles->nt=CurrentTpImage->dim[4];
        warpedVentricles->cal_min=CurrentTpImage->cal_min;
        warpedVentricles->cal_max=CurrentTpImage->cal_max;
        warpedVentricles->scl_slope=CurrentTpImage->scl_slope;
        warpedVentricles->scl_inter=CurrentTpImage->scl_inter;
        reg_tools_changeDatatype<float>(CurrentTpImage);
        warpedVentricles->nbyper = CurrentTpImage->nbyper;
        warpedVentricles->nvox = (size_t)warpedVentricles->dim[1] * (size_t)warpedVentricles->dim[2] *
                (size_t)warpedVentricles->dim[3] * (size_t)warpedVentricles->dim[4];
        warpedVentricles->data = (void *)calloc(warpedVentricles->nvox, warpedVentricles->nbyper);
        float *warpedVentriclesDataPtr = static_cast<float *>(warpedVentricles->data);

        // actually resample the ventricle mask data
        reg_resampleImage(Ventricles,warpedVentricles,defField,NULL,0,0);


        for(int i=0; i<NextGroupwiseMean->nvox;i++){
            ErrorInInputImagesDataPtr[i+tp*CurrentTpImage->nvox]=MaskDataPtr[i]*(CurrentTpImageDataPtr[i]-warpedMeanDataPtr[i])/warpedVarianceDataPtr[i];
        }

        float SNR_signal=0;
        int SNR_sum=0;
        for(int i=0; i<NextGroupwiseMean->nvox;i++){
            SNR_signal+=warpedVentriclesDataPtr[i]*CurrentTpImageDataPtr[i];
            SNR_sum+=warpedVentriclesDataPtr[i];
        }

        SNR_signal=SNR_signal/SNR_sum;
        float SNR_noise=0;
        for(int i=0; i<NextGroupwiseMean->nvox;i++){
            SNR_noise+=warpedVentriclesDataPtr[i]*powf(CurrentTpImageDataPtr[i]-SNR_signal,2);
        }
        SNR_noise=SNR_noise/SNR_sum;
        CSV_values[tp][3]=SNR_signal/sqrt(SNR_noise);


        nifti_image_free(CurrentTpImage);
        nifti_image_free(warpedMean);
        nifti_image_free(warpedVariance);
        nifti_image_free(warpedVentricles);
        nifti_image_free(defField);
    }

    nifti_image_write(ErrorInInputImages);




    for(int tp=0; tp<InputImages->nt;tp++){
        float mean=0;
        for(int iz=0; iz<ErrorInInputImages->nz;iz+=1){
            for(int iy=0; iy<ErrorInInputImages->ny;iy+=1){
                for(int ix=0; ix<ErrorInInputImages->nx;ix+=1){
                    int index=ix+iy*ErrorInInputImages->nx+iz*ErrorInInputImages->nx*ErrorInInputImages->ny;
                    mean+=MaskDataPtr[index]*ErrorInInputImagesDataPtr[index+tp*NextGroupwiseMean->nvox];

                }
            }
        }
        mean=mean/((float)(NextGroupwiseMean->nvox));

       CSV_values[tp][4]=mean;
    }



    for(int tp=0; tp<InputImages->nt;tp++){
        float mean0=0;
        for(int iz=0; iz<ErrorInInputImages->nz;iz+=2){
            for(int iy=0; iy<ErrorInInputImages->ny;iy+=1){
                for(int ix=0; ix<ErrorInInputImages->nx;ix+=1){
                    int index=ix+iy*ErrorInInputImages->nx+iz*ErrorInInputImages->nx*ErrorInInputImages->ny;
                    mean0+=MaskDataPtr[index]*ErrorInInputImagesDataPtr[index+tp*NextGroupwiseMean->nvox];

                }
            }
        }
        mean0=mean0/((float)(NextGroupwiseMean->nvox)/2);

        float mean1=0;
        for(int iz=1; iz<ErrorInInputImages->nz;iz+=2){
            for(int iy=0; iy<ErrorInInputImages->ny;iy+=1){
                for(int ix=0; ix<ErrorInInputImages->nx;ix+=1){
                    int index=ix+iy*ErrorInInputImages->nx+iz*ErrorInInputImages->nx*ErrorInInputImages->ny;
                    mean1+=MaskDataPtr[index]*ErrorInInputImagesDataPtr[index+tp*NextGroupwiseMean->nvox];

                }
            }
        }
        mean1=mean1/((float)(NextGroupwiseMean->nvox)/2);

        CSV_values[tp][5]=fabs(mean1-mean0);
    }



    ofstream myfile;
    myfile.open(resultFile.c_str());
    myfile<<"Movement - Distance To groupwise";
    myfile<<",Movement - Distance from first timepoint";
    myfile<<",Movement - Distance between adjacent timepoints";
    myfile<<",SNR per time point";
    myfile<<",Intensity - Mean z-scored signal deviation";
    myfile<<",Intensity - Mean z-scored inter-slice error";
    myfile<<std::endl;

    for(int tp=0; tp<InputImages->nt;tp++){
        for(int item=0; item<6;item++){
            if(item==0){
            myfile<<CSV_values[tp][item];
            }
            if(item>0){
            myfile<<","<<CSV_values[tp][item];
            }
        }
        myfile<<std::endl;
    }




    myfile.close();
    nifti_image_free(Mask);
    nifti_image_free(RessampledInputImages);
    nifti_image_free(ErrorInInputImages);
    nifti_image_free(NextGroupwiseMean);
    nifti_image_free(CurrentGroupwise);
    nifti_image_free(InputImages);
    return 0;
}


