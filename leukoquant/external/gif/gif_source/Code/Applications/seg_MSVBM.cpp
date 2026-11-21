
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include "_seg_common.h"
#include "_reg_localTransformation.h"
#include "_reg_resampling.h"
#include "_reg_maths.h"
#include "_reg_tools.h"
#include "_seg_tools.h"
#include "_seg_Topo.h"

using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("\n  BrainMesh Segmentation :\n  Usage ->\t%s <mandatory_options>\n\n",exec);
    printf("  * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-in <filename> \t\t| Input anatomical image filename\n");
    printf("\t-mask <filename>\t| Mask over the input image\n");
    printf("\t-outTmap <filename>\t\t| Output prior segemntation filename\n");
    printf("\t-image_db <path>\t| Folder path to database with anatomical images\n");
    printf("\t-cpp_db <path>\t\t| Folder path to database with cpp files\n");
    printf("\t-seg_db <nc> <path>\t\t| Folder path to database with the current probabilistic segmentations\n");
    printf("\t-gstring <str> <str> \t\t| The groupsize per arm. (Roughly 1/6 of the tot_num_subj in the 2 groups)\n");
    printf("\n  * * * * * * * * * * * * * * * * * Options * * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-v <int>\t\t| Verbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
    printf("\t-dnm \t\t| Do not modulate\n");
    printf("\t-fwhm <float> \t\t| FWHM smoothing (default=6mm)\n");
    printf("\t-temp <float>\t\t| Heat kernel temperature [default=0.5]\n");
    printf("\t-weigh <float> <float>\t| Heat kernel weight for intensity and displacement [default=(0.7,0.3)]\n");
    printf("\t-kstd <int> <int>\t| Local intensity and displacement metric kernel [default=(5,2)]\n");
    printf("\n  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    return;
}


int main(int argc, char **argv)
{

#ifndef NDEBUG
    printf("[NiftySeg DEBUG] NiftySeg has been compiled in DEBUG mode\n");
#endif
    if (argc < 1)
    {
        Usage(argv[0]);
        return 0;
    }

    time_t start,end;
    time(&start);
    string target_subject_file;
    string OutputTmapFilename;
    string cpp_folder;
    string image_folder;
    string seg_folder;
    string target_subject_mask;
    string grp0_str;
    string grp1_str;
    int verbose=0;
    bool modulate=1;

    int  kernelstd_img=5;
    int kernelstd_dis=2;
    float weightIntensity=0.7f;
    float weightDisplacement=0.3f;
    float temperature=0.5f;
    int float_nt=0;
    float fwhm=2.0f;




    /* read the input parameter */

    for(int i=1; i<argc; i++)
    {
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc)
        {
            target_subject_file=argv[++i];
        }
        else if(strcmp(argv[i], "-mask") == 0 && (i+1)<argc)
        {
            target_subject_mask=argv[++i];
        }
        else if(strcmp(argv[i], "-outTmap") == 0 && (i+1)<argc)
        {
            OutputTmapFilename = argv[++i];
        }
        else if(strcmp(argv[i], "-cpp_db") == 0 && (i+1)<argc)
        {
            cpp_folder = argv[++i];
        }
        else if(strcmp(argv[i], "-seg_db") == 0 && (i+2)<argc)
        {
            float_nt=atoi(argv[++i]);
            seg_folder = argv[++i];
        }
        else if(strcmp(argv[i], "-image_db") == 0 && (i+1)<argc)
        {
            image_folder = argv[++i];
        }
        else if(strcmp(argv[i], "-temp") == 0 && (i+1)<argc)
        {
            temperature = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-fwhm") == 0 && (i+1)<argc)
        {
            fwhm = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-dnm") == 0 && (i)<argc)
        {
            modulate=0;
        }
        else if(strcmp(argv[i], "-gstring") == 0 && (i+2)<argc)
        {
            grp0_str = argv[++i];
            grp1_str = argv[++i];
        }
        else if(strcmp(argv[i], "-weigh") == 0 && (i+2)<argc)
        {
            weightIntensity = atof(argv[++i]);
            weightDisplacement = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-kstd") == 0 && (i+2)<argc)
        {
            kernelstd_img = atoi(argv[++i]);
            kernelstd_dis = atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-v") == 0 && (i+1)<argc)
        {
            verbose = atoi(argv[++i]);
        }
        else
        {
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }


    if(image_folder.size()==0 || cpp_folder.size()==0 || OutputTmapFilename.size()==0|| target_subject_file.size()==0||seg_folder.size()==0)
    {
        cout << "ERROR: parameter ";
        if( target_subject_file.size()==0)
            cout << "in, ";
        if(OutputTmapFilename.size()==0)
            cout << "outSeg, ";
        if(image_folder.size()==0)
            cout << "image_db, ";
        if(cpp_folder.size()==0)
            cout << "cpp, ";
        if(seg_folder.size()==0)
            cout << "seg_db, ";

        cout << " is/are missing";
        Usage(argv[0]);
        return 0;
    }

    if(float_nt<1)
    {
        cout << "The number of tissue classes is smaller than 1. Check -seg_db <nc> param.";
        Usage(argv[0]);
        return 0;
    }

    string subject=target_subject_file.substr((target_subject_file.find_last_of("/")+1),(target_subject_file.find(".nii"))-(target_subject_file.find_last_of("/")+1));
    if(verbose>0)
        cout <<"Propagating info to subject: "<<subject<<endl;

    // reading in the subject image
    nifti_image * SubjectImage=nifti_image_read(target_subject_file.c_str(),true);
    if(SubjectImage == NULL)
    {
        fprintf(stderr,"* Error when reading the segmentation image: %s\n",target_subject_file.c_str());
        return 1;
    }
    if(SubjectImage->datatype!=NIFTI_TYPE_FLOAT32)
        seg_changeDatatype<float>(SubjectImage);
    // reading in the subject image's mask
    nifti_image * SubjectImageMask=NULL;
    int * SubjectImageMaskPtr=NULL;


    if(target_subject_mask.size()>0)
    {
        SubjectImageMask=nifti_image_read(target_subject_mask.c_str(),true);
        if(SubjectImageMask == NULL)
        {
            fprintf(stderr,"* Error when reading the target mask image: %s\n",target_subject_file.c_str());
            return 1;
        }
        if(SubjectImageMask->datatype!=NIFTI_TYPE_UINT8)
            seg_changeDatatype<unsigned char>(SubjectImageMask);
        unsigned char * SubjectImageMaskOriPtr = static_cast<unsigned char *>(SubjectImageMask->data);

        unsigned char * SubjectImageMask_tmp_Ptr=new unsigned char [SubjectImageMask->nvox];
        SubjectImageMaskPtr=new int [SubjectImageMask->nvox];

        ImageSize * CurrSize = new ImageSize [1]();
        CurrSize->numel=(int)(SubjectImageMask->nx*SubjectImageMask->ny*SubjectImageMask->nz);
        CurrSize->xsize=SubjectImageMask->nx;
        CurrSize->ysize=SubjectImageMask->ny;
        CurrSize->zsize=SubjectImageMask->nz;
        CurrSize->usize=(SubjectImageMask->nu>1)?SubjectImageMask->nu:1;
        CurrSize->tsize=(SubjectImageMask->nt>1)?SubjectImageMask->nt:1;

        Close_Forground_ConnectComp<unsigned char,unsigned char>(static_cast<void*>(SubjectImageMaskOriPtr),SubjectImageMask_tmp_Ptr,CurrSize);
        for(unsigned int index=0; index<SubjectImageMask->nvox; index++)
        {
            SubjectImageMaskPtr[index]=SubjectImageMask_tmp_Ptr[index];
        }
        delete [] SubjectImageMask_tmp_Ptr;
        nifti_image_free(SubjectImageMask);
    }


    // Allocating space for the final tmap
    nifti_image * ResultTmap =NULL;
    ResultTmap = nifti_copy_nim_info(SubjectImage);
    ResultTmap->dim[0]=4;
    ResultTmap->dim[4]=ResultTmap->nt=float_nt;
    ResultTmap->datatype=NIFTI_TYPE_FLOAT32;
    ResultTmap->cal_max=1;
    nifti_set_filenames(ResultTmap,OutputTmapFilename.c_str(),0,0);
    nifti_update_dims_from_array(ResultTmap);
    nifti_datatype_sizes(ResultTmap->datatype,&ResultTmap->nbyper,&ResultTmap->swapsize);
    ResultTmap->data = (void *) calloc(ResultTmap->nvox, sizeof(float));
    if(ResultTmap->data == NULL)
    {
        fprintf(stderr,"* Error when allocating ResultTmap. Not enough memory");
        return 1;
    }

    // Allocating space for the resampled T1's
    nifti_image * ResampledAnatomical = nifti_copy_nim_info(SubjectImage);
    ResampledAnatomical->dim[0]=3;
    ResampledAnatomical->dim[4]=1;
    ResampledAnatomical->datatype=NIFTI_TYPE_FLOAT32;
    nifti_set_filenames(ResampledAnatomical,"Resampled_Anatomical_and_LSSD.nii.gz",0,0);
    nifti_update_dims_from_array(ResampledAnatomical);
    nifti_datatype_sizes(ResampledAnatomical->datatype,&ResampledAnatomical->nbyper,&ResampledAnatomical->swapsize);
    ResampledAnatomical->data = (void *) calloc(ResampledAnatomical->nvox, sizeof(float));
    if(ResampledAnatomical->data == NULL)
    {
        fprintf(stderr,"* Error when allocating ResampledAnatomical. Not enough memory");
        return 1;
    }

    // Allocating space for the Deformation fields
    nifti_image * Deformation = nifti_copy_nim_info(SubjectImage);
    Deformation->dim[0]=Deformation->ndim=5;
    Deformation->dim[1]=Deformation->nx=SubjectImage->nx;
    Deformation->dim[2]=Deformation->ny=SubjectImage->ny;
    Deformation->dim[3]=Deformation->nz=SubjectImage->nz;
    Deformation->dim[4]=Deformation->nt=1;
    Deformation->pixdim[4]=Deformation->dt=1.0;
    if(SubjectImage->nz>1) Deformation->dim[5]=Deformation->nu=3;
    else Deformation->dim[5]=Deformation->nu=2;
    Deformation->pixdim[5]=Deformation->du=1.0;
    Deformation->dim[6]=Deformation->nv=1;
    Deformation->pixdim[6]=Deformation->dv=1.0;
    Deformation->dim[7]=Deformation->nw=1;
    Deformation->pixdim[7]=Deformation->dw=1.0;
    Deformation->nvox=Deformation->nx*Deformation->ny*Deformation->nz*Deformation->nt*Deformation->nu;
    if(sizeof(SegPrecisionTYPE)==8) Deformation->datatype = NIFTI_TYPE_FLOAT64;
    else Deformation->datatype = NIFTI_TYPE_FLOAT32;
    Deformation->nbyper = sizeof(SegPrecisionTYPE);
    Deformation->data = (void *)calloc(Deformation->nvox, Deformation->nbyper);
    if(Deformation->data == NULL)
    {
        fprintf(stderr,"* Error when allocating Deformation. Not enough memory");
        return 1;
    }

    // Allocate space for the Jacobian
    nifti_image *jacobianImage = nifti_copy_nim_info(SubjectImage);
    jacobianImage->cal_min=0;
    jacobianImage->cal_max=0;
    jacobianImage->scl_slope = 1.0f;
    jacobianImage->scl_inter = 0.0f;
    jacobianImage->datatype = NIFTI_TYPE_FLOAT32;
    jacobianImage->nbyper = sizeof(float);
    jacobianImage->data = (void *)calloc(jacobianImage->nvox, jacobianImage->nbyper);


    // Allocating space for the resampled segmentations

    nifti_image * ResampProbability=NULL;
    ResampProbability = nifti_copy_nim_info(SubjectImage);
    ResampProbability->dim[0]=ResampProbability->ndim=4;
    ResampProbability->dim[1]=ResampProbability->nx=SubjectImage->nx;
    ResampProbability->dim[2]=ResampProbability->ny=SubjectImage->ny;
    ResampProbability->dim[3]=ResampProbability->nz=SubjectImage->nz;
    ResampProbability->dim[4]=ResampProbability->nt=float_nt;
    ResampProbability->pixdim[4]=ResampProbability->dt=1.0;
    ResampProbability->dim[5]=ResampProbability->nu=1;
    ResampProbability->pixdim[5]=ResampProbability->du=1.0;
    ResampProbability->dim[6]=ResampProbability->nv=1;
    ResampProbability->pixdim[6]=ResampProbability->dv=1.0;
    ResampProbability->dim[7]=ResampProbability->nw=1;
    ResampProbability->pixdim[7]=ResampProbability->dw=1.0;
    ResampProbability->nvox=ResampProbability->nx*ResampProbability->ny*ResampProbability->nz*ResampProbability->nt*ResampProbability->nu;
    nifti_update_dims_from_array(ResampProbability);
    nifti_set_filenames(ResampProbability,"ResampProbability.nii.gz",0,0);
    ResampProbability->datatype = NIFTI_TYPE_FLOAT32;
    ResampProbability->nbyper = sizeof(float);
    ResampProbability->data = (void *)calloc(ResampProbability->nvox, ResampProbability->nbyper);
    if(ResampProbability->data == NULL)
    {
        fprintf(stderr,"* Error when allocating ResampProbability. Not enough memory");
        return 1;
    }

    ImageSize * CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz);
    CurrSizes->xsize=SubjectImage->nx;
    CurrSizes->ysize=SubjectImage->ny;
    CurrSizes->zsize=SubjectImage->nz;
    CurrSizes->usize=1;
    CurrSizes->tsize=float_nt;
    CurrSizes->numclass=1;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;









    // ***********************************************************************
    // ***********************************************************************
    // *************************      GROUP 0       **************************
    // ***********************************************************************
    // ***********************************************************************


    // *********** alocate weights for grp 0 ************
    float * Weights_group0=NULL;
    Weights_group0=new float [(SubjectImage->nvox)];
    if(Weights_group0 == NULL)
    {
        fprintf(stderr,"* Error when allocating Weights_group0. Not enough memory");
        return 1;
    }
    for(unsigned int index=0; index<((SubjectImage->nvox)); index++)
    {
        Weights_group0[index]=0;
    }

    // *********** alocate probs for grp 0 ************
    float * Probs_group0=NULL;
    Probs_group0=new float [(SubjectImage->nvox)*(float_nt)];
    if(Probs_group0 == NULL)
    {
        fprintf(stderr,"* Error when allocating Probs_group0. Not enough memory");
        return 1;
    }
    for(unsigned int index=0; index<((SubjectImage->nvox)*(float_nt)); index++)
    {
        Probs_group0[index]=0;
    }


    // GET FILES FROM GROUP 0
    vector<string> files_cpp = vector<string>();
    string group_string = grp0_str;
    get_all_files_that_match_string(cpp_folder,files_cpp,group_string);
    if(verbose>0)
    {
        cout << "Found "<<files_cpp.size()<<" cpp files for group "<<grp0_str<<endl;
        if(verbose>1)
        {
            for(unsigned int filenumb=0; filenumb<files_cpp.size(); filenumb++)
                cout << "\t -> "<<files_cpp[filenumb]<<endl;
        }
        cout<<"\nExtracting metrics and sampling from the database:"<<endl;
    }

    //for each file it found do
    for (unsigned int i = 0; i < files_cpp.size(); i++)
    {
        //for (int i = 0;i<5;i++) {
        string current_cpp_file=files_cpp[i];
        string current_patient_name=current_cpp_file.substr(current_cpp_file.find_last_of("/")+1,current_cpp_file.size()-(current_cpp_file.find_last_of("/")+1));
        current_patient_name=current_patient_name.substr(4,current_patient_name.find(".nii",0)-4);
        string segmentation_file_name=seg_folder;
        segmentation_file_name+="/";
        segmentation_file_name+=current_patient_name;
        segmentation_file_name+=".nii.gz";
        string Image_file_name=image_folder;
        Image_file_name+="/";
        Image_file_name+=current_patient_name;
        Image_file_name+=".nii.gz";

        FILE * F1=fopen(current_cpp_file.c_str(),"rb");
        FILE * F2=fopen(Image_file_name.c_str(),"rb");
        FILE * F3=fopen(segmentation_file_name.c_str(),"rb");
        if(verbose>0)
            cout << "Subj #"<<i<<
                 ":\n\tPID: "<< current_patient_name<<
                 "\n\tCPP: "<< current_cpp_file<< (F1!=0?" - Found":" - Not Found")<<
                 "\n\tIMG: "<< Image_file_name << (F2!=0?" - Found":" - Not Found")<<
                 "\n\tPROP:"<< segmentation_file_name<< (F3!=0?" - Found":" - Not Found\n")<<
                 endl;

        // If all files necessary for the current subject were found, do:
        bool dotheyexist= ((F1!=0) && (F2!=0) && (F3!=0));
        if(F1!=0)
        {
            fclose(F1);
        }
        if(F2!=0)
        {
            fclose(F2);
        }
        if(F3!=0)
        {
            fclose(F3);
        }

        if (dotheyexist)
        {

            //  ********************  Load CPP image  *******************
            //
            nifti_image * InputCPP=nifti_image_read(current_cpp_file.c_str(),true);
            if(InputCPP == NULL)
            {
                fprintf(stderr,"* Error when reading the cpp image: %s\n",current_cpp_file.c_str());
                break;
            }
            if(InputCPP->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<SegPrecisionTYPE>(InputCPP);

            //  ********************  Load Anatomical image  *******************
            //
            nifti_image * InputAnatomical=nifti_image_read(Image_file_name.c_str(),true);
            if(InputAnatomical == NULL)
            {
                nifti_image_free(InputAnatomical);
                fprintf(stderr,"* Error when reading the Anatomical image: %s\n",Image_file_name.c_str());
            }
            if(InputAnatomical->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<float>(InputAnatomical);

            //  ********************  Load Seg image  *******************
            //
            nifti_image * InputProp=nifti_image_read(segmentation_file_name.c_str(),true);
            if(InputProp == NULL)
            {
                fprintf(stderr,"* Error when reading the segmentation image: %s\n",segmentation_file_name.c_str());
            }
            if(InputProp->datatype!=NIFTI_TYPE_FLOAT32)
            {
                seg_changeDatatype<float>(InputProp);
            }

            //  ********************  START ALGORITHM  *******************
            if(verbose>1)
                cout<< "\n\tGET DEFORMATION FIELD"<<endl;
            if( InputCPP->intent_code==NIFTI_INTENT_VECTOR &&
                    strcmp(InputCPP->intent_name,"NREG_VEL_STEP")==0)
            {
                reg_spline_getDefFieldFromVelocityGrid(InputCPP,
                        Deformation,true);
            }
            else
            {
                reg_spline_getDeformationField(InputCPP,
                                               Deformation,
                                               SubjectImageMaskPtr, // mask
                                               false, //composition
                                               true // bspline
                                              );
            }

            if(verbose>1)
                cout<< "\tGET DISPLACEMENT FIELD FROM DEFORMATION FIELD"<<endl;
            nifti_image *Displacement = nifti_copy_nim_info(Deformation);
            Displacement->data = (void *)malloc(Deformation->nvox * Deformation->nbyper);
            memcpy(Displacement->data, Deformation->data, Deformation->nvox*Deformation->nbyper);

            reg_getDisplacementFromDeformation(Displacement);
            Displacement->nt=1;
            Displacement->dim[0]=Displacement->ndim=3;
            Displacement->dim[4]=1;
            Displacement->nvox=Displacement->nx*Displacement->ny*Displacement->nz;
            nifti_update_dims_from_array(Displacement);
            //cout << Displacement->nvox<<","<<Displacement->nx<<","<<Displacement->ny<<","<<Displacement->nz<<","<<Displacement->nt<<endl;
            SegPrecisionTYPE * DisplacementPtr = static_cast<SegPrecisionTYPE *>(Displacement->data);
            for(unsigned int index=0; index<Displacement->nvox; index++)
            {
                DisplacementPtr[index]=sqrt(pow(DisplacementPtr[index],2)+pow(DisplacementPtr[index+Displacement->nvox],2)+pow(DisplacementPtr[index+Displacement->nvox],2));
            }


            if(verbose>1)
                cout<< "\tGET SMOTH DISPLACEMENT FIELD FROM ORIGINAL DISPLACEMENT FIELD"<<endl;

            nifti_image *Displacement_Smooth = nifti_copy_nim_info(Displacement);
            Displacement_Smooth->data = (void *)malloc(Displacement->nvox * Displacement->nbyper);
            memcpy(Displacement_Smooth->data, Displacement->data, Displacement->nvox*Displacement->nbyper);
            float radius=20.0f;
            GaussianSmoothing(Displacement_Smooth,NULL, radius);
            if(verbose>1)
                cout<< "\tSUBTRACT SMOTH DISPLACEMENT FIELD FROM ORIGINAL DISPLACEMENT FIELD IN ORDER TO GET ONLY YHE HIGH-FREQUENCY DISPLACEMENT FIELD"<<endl;

            SegPrecisionTYPE * Displacement_SmoothPtr = static_cast<SegPrecisionTYPE *>(Displacement_Smooth->data);


            for(unsigned int index=0; index<Displacement_Smooth->nvox; index++)
            {
                DisplacementPtr[index]=(DisplacementPtr[index] - Displacement_SmoothPtr[index]);
            }

            radius=20.0f;
            GaussianSmoothing(Displacement,NULL, radius);

            for(unsigned int index=0; index<Displacement_Smooth->nvox; index++)
            {
                DisplacementPtr[index]=pow(DisplacementPtr[index],2);
            }

            if(verbose>1)
                cout<< "\tRESAMPLE ANATOMICAL IMAGE TO GET IMAGE SIMILARITY"<<endl;

            reg_resampleImage(InputAnatomical,
                              ResampledAnatomical,
                              Deformation,
                              SubjectImageMaskPtr,
                              1,
                              0);

            if(verbose>1)
                cout<< "\tCALCULATING METRIC"<<endl;
            SegPrecisionTYPE * ResampledAnatomicalPtr = static_cast<SegPrecisionTYPE *>(ResampledAnatomical->data);
            SegPrecisionTYPE * SubjectImagePtr = static_cast<SegPrecisionTYPE *>(SubjectImage->data);
            float allmeanNew=0;
            float allmeanInput=0;
            float allstdNew=0;
            float allstdInput=0;
            int maskcount=0;
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    allmeanNew+=ResampledAnatomicalPtr[index];
                    allmeanInput+=SubjectImagePtr[index];
                    maskcount++;
                }
            }

            allmeanNew=allmeanNew/(float)(maskcount);
            allmeanInput=allmeanInput/(float)(maskcount);
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    allstdNew+=(ResampledAnatomicalPtr[index]-allmeanNew)*(ResampledAnatomicalPtr[index]-allmeanNew);
                    allstdInput+=(SubjectImagePtr[index]-allmeanInput)*(SubjectImagePtr[index]-allmeanInput);
                }
            }
            allstdNew=allstdNew/(float)(maskcount);
            allstdInput=allstdInput/(float)(maskcount);
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    ResampledAnatomicalPtr[index]=pow(((ResampledAnatomicalPtr[index]-allmeanNew)/sqrt(allstdNew))-((SubjectImagePtr[index]-allmeanInput)/sqrt(allstdInput)),2);
                }
                else
                {
                    ResampledAnatomicalPtr[index]=0.0f;
                }
            }

            radius=kernelstd_img;
            GaussianSmoothing(ResampledAnatomical,NULL, radius);

            // ResampledAnatomicalPtr will be  the current subject's weight.
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    ResampledAnatomicalPtr[index]=(weightIntensity*ResampledAnatomicalPtr[index]+weightDisplacement*DisplacementPtr[index]);
                    ResampledAnatomicalPtr[index]=exp(-(ResampledAnatomicalPtr[index]/temperature));
                }
                else
                {
                    ResampledAnatomicalPtr[index]=0.0f;
                }
            }

            // Resample the probabilities, calculate the jacobian and modulate probabilities by the jacobian

            reg_resampleImage(InputProp,ResampProbability,Deformation,SubjectImageMaskPtr,1,0);
            float * ResampProbabilityPtr = static_cast<float *>(ResampProbability->data);
            if(modulate)
            {
                reg_defField_getJacobianMap(Deformation,jacobianImage);

                float * ResampJacPtr = static_cast<float *>(jacobianImage->data);

                for(int nt=0; nt<ResampProbability->nt; nt++)
                {
                    for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
                    {
                        if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                        {
                            ResampProbabilityPtr[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)]=ResampJacPtr[index]*ResampProbabilityPtr[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)];
                        }
                    }
                }
            }

            Gaussian_Filter_4D(ResampProbabilityPtr,fwhm/2.354820f,CurrSizes);


            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                for(int nt=0; nt<ResampProbability->nt; nt++)
                {
                    Probs_group0[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)]+=
                        ResampledAnatomicalPtr[index]
                        *ResampProbabilityPtr[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)];
                }
                Weights_group0[index]+=ResampledAnatomicalPtr[index];
            }

            cout << "\tUPDATING WEIGHTS, PROBS and INDECES"<<endl;

            // Free images
            nifti_image_free(Displacement_Smooth);
            nifti_image_free(Displacement);
            nifti_image_free(InputCPP);
            nifti_image_free(InputProp);
            nifti_image_free(InputAnatomical);

        }
    }

    files_cpp.clear();



    // ***********************************************************************
    // ***********************************************************************
    // *************************      GROUP 1       **************************
    // ***********************************************************************
    // ***********************************************************************


    // *********** alocate weights for grp 1 ************
    float * Weights_group1=NULL;
    Weights_group1=new float [(SubjectImage->nvox)];
    if(Weights_group1 == NULL)
    {
        fprintf(stderr,"* Error when allocating Weighted_vote_matrix_group1. Not enough memory");
        return 1;
    }
    for(unsigned int index=0; index<((SubjectImage->nvox)); index++)
    {
        Weights_group1[index]=0;
    }

    // *********** alocate probs for grp 1 ************
    float * Probs_group1=NULL;
    Probs_group1=new float [(SubjectImage->nvox)*(float_nt)];
    if(Probs_group1 == NULL)
    {
        fprintf(stderr,"* Error when allocating Weighted_vote_matrix_group1. Not enough memory");
        return 1;
    }
    for(unsigned int index=0; index<((SubjectImage->nvox)*(float_nt)); index++)
    {
        Probs_group1[index]=0;
    }


    // GET FILES FROM GROUP
    files_cpp.clear();
    group_string.clear();
    group_string= grp1_str;
    get_all_files_that_match_string(cpp_folder,files_cpp,group_string);
    if(verbose>0)
    {
        cout << "Found "<<files_cpp.size()<<" cpp files for group "<<grp1_str<<endl;
        if(verbose>1)
        {
            for(unsigned int filenumb=0; filenumb<files_cpp.size(); filenumb++)
                cout << "\t -> "<<files_cpp[filenumb]<<endl;
        }
        cout<<"\nExtracting metrics and sampling from the database:"<<endl;
    }

    //for each file it found do
    for (unsigned int i = 0; i < files_cpp.size(); i++)
    {
        string current_cpp_file=files_cpp[i];
        string current_patient_name=current_cpp_file.substr(current_cpp_file.find_last_of("/")+1,current_cpp_file.size()-(current_cpp_file.find_last_of("/")+1));
        current_patient_name=current_patient_name.substr(4,current_patient_name.find(".nii",0)-4);
        string segmentation_file_name=seg_folder;
        segmentation_file_name+="/";
        segmentation_file_name+=current_patient_name;
        segmentation_file_name+=".nii.gz";
        string Image_file_name=image_folder;
        Image_file_name+="/";
        Image_file_name+=current_patient_name;
        Image_file_name+=".nii.gz";

        FILE * F1=fopen(current_cpp_file.c_str(),"rb");
        FILE * F2=fopen(Image_file_name.c_str(),"rb");
        FILE * F3=fopen(segmentation_file_name.c_str(),"rb");
        if(verbose>0)
            cout << "Subj #"<<i<<
                 ":\n\tPID: "<< current_patient_name<<
                 "\n\tCPP: "<< current_cpp_file<< (F1!=0?" - Found":" - Not Found")<<
                 "\n\tIMG: "<< Image_file_name << (F2!=0?" - Found":" - Not Found")<<
                 "\n\tPROP:"<< segmentation_file_name<< (F3!=0?" - Found":" - Not Found\n")<<
                 endl;

        // If all files necessary for the current subject were found, do:
        bool dotheyexist= ((F1!=0) && (F2!=0) && (F3!=0));
        if(F1!=0)
        {
            fclose(F1);
        }
        if(F2!=0)
        {
            fclose(F2);
        }
        if(F3!=0)
        {
            fclose(F3);
        }

        if (dotheyexist)
        {

            //  ********************  Load CPP image  *******************
            //
            nifti_image * InputCPP=nifti_image_read(current_cpp_file.c_str(),true);
            if(InputCPP == NULL)
            {
                fprintf(stderr,"* Error when reading the cpp image: %s\n",current_cpp_file.c_str());
                break;
            }
            if(InputCPP->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<SegPrecisionTYPE>(InputCPP);

            //  ********************  Load Anatomical image  *******************
            //
            nifti_image * InputAnatomical=nifti_image_read(Image_file_name.c_str(),true);
            if(InputAnatomical == NULL)
            {
                nifti_image_free(InputAnatomical);
                fprintf(stderr,"* Error when reading the Anatomical image: %s\n",Image_file_name.c_str());
            }
            if(InputAnatomical->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<float>(InputAnatomical);

            //  ********************  Load Seg image  *******************
            //
            nifti_image * InputProp=nifti_image_read(segmentation_file_name.c_str(),true);
            if(InputProp == NULL)
            {
                fprintf(stderr,"* Error when reading the segmentation image: %s\n",segmentation_file_name.c_str());
            }
            if(InputProp->datatype!=NIFTI_TYPE_FLOAT32)
            {
                seg_changeDatatype<float>(InputProp);
            }

            //  ********************  START ALGORITHM  *******************
            if(verbose>1)
                cout<< "\n\tGET DEFORMATION FIELD"<<endl;
            if( InputCPP->intent_code==NIFTI_INTENT_VECTOR &&
                    strcmp(InputCPP->intent_name,"NREG_VEL_STEP")==0)
            {
                reg_spline_getDefFieldFromVelocityGrid(InputCPP,
                        Deformation,true);
            }
            else
            {
                reg_spline_getDeformationField(InputCPP,
                                               Deformation,
                                               SubjectImageMaskPtr, // mask
                                               false, //composition
                                               true // bspline
                                              );
            }

            if(verbose>1)
                cout<< "\tGET DISPLACEMENT FIELD FROM DEFORMATION FIELD"<<endl;
            nifti_image *Displacement = nifti_copy_nim_info(Deformation);
            Displacement->data = (void *)malloc(Deformation->nvox * Deformation->nbyper);
            memcpy(Displacement->data, Deformation->data, Deformation->nvox*Deformation->nbyper);

            reg_getDisplacementFromDeformation(Displacement);
            Displacement->nt=1;
            Displacement->dim[0]=Displacement->ndim=3;
            Displacement->dim[4]=1;
            Displacement->nvox=Displacement->nx*Displacement->ny*Displacement->nz;
            nifti_update_dims_from_array(Displacement);
            //cout << Displacement->nvox<<","<<Displacement->nx<<","<<Displacement->ny<<","<<Displacement->nz<<","<<Displacement->nt<<endl;
            SegPrecisionTYPE * DisplacementPtr = static_cast<SegPrecisionTYPE *>(Displacement->data);
            for(unsigned int index=0; index<Displacement->nvox; index++)
            {
                DisplacementPtr[index]=sqrt(pow(DisplacementPtr[index],2)+pow(DisplacementPtr[index+Displacement->nvox],2)+pow(DisplacementPtr[index+Displacement->nvox],2));
            }


            if(verbose>1)
                cout<< "\tGET SMOTH DISPLACEMENT FIELD FROM ORIGINAL DISPLACEMENT FIELD"<<endl;

            nifti_image *Displacement_Smooth = nifti_copy_nim_info(Displacement);
            Displacement_Smooth->data = (void *)malloc(Displacement->nvox * Displacement->nbyper);
            memcpy(Displacement_Smooth->data, Displacement->data, Displacement->nvox*Displacement->nbyper);
            float radius=20.0f;
            GaussianSmoothing(Displacement_Smooth,NULL, radius);
            if(verbose>1)
                cout<< "\tSUBTRACT SMOTH DISPLACEMENT FIELD FROM ORIGINAL DISPLACEMENT FIELD IN ORDER TO GET ONLY YHE HIGH-FREQUENCY DISPLACEMENT FIELD"<<endl;

            SegPrecisionTYPE * Displacement_SmoothPtr = static_cast<SegPrecisionTYPE *>(Displacement_Smooth->data);


            for(unsigned int index=0; index<Displacement_Smooth->nvox; index++)
            {
                DisplacementPtr[index]=(DisplacementPtr[index] - Displacement_SmoothPtr[index]);
            }

            radius=kernelstd_dis;
            GaussianSmoothing(Displacement,NULL, radius);

            for(unsigned int index=0; index<Displacement_Smooth->nvox; index++)
            {
                DisplacementPtr[index]=pow(DisplacementPtr[index],2);
            }

            if(verbose>1)
                cout<< "\tRESAMPLE ANATOMICAL IMAGE TO GET IMAGE SIMILARITY"<<endl;

            reg_resampleImage(InputAnatomical,
                              ResampledAnatomical,
                              Deformation,
                              SubjectImageMaskPtr,
                              1,
                              0);

            if(verbose>1)
                cout<< "\tCALCULATING METRIC"<<endl;
            SegPrecisionTYPE * ResampledAnatomicalPtr = static_cast<SegPrecisionTYPE *>(ResampledAnatomical->data);
            SegPrecisionTYPE * SubjectImagePtr = static_cast<SegPrecisionTYPE *>(SubjectImage->data);
            float allmeanNew=0;
            float allmeanInput=0;
            float allstdNew=0;
            float allstdInput=0;
            int maskcount=0;
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    allmeanNew+=ResampledAnatomicalPtr[index];
                    allmeanInput+=SubjectImagePtr[index];
                    maskcount++;
                }
            }

            allmeanNew=allmeanNew/(float)(maskcount);
            allmeanInput=allmeanInput/(float)(maskcount);
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    allstdNew+=(ResampledAnatomicalPtr[index]-allmeanNew)*(ResampledAnatomicalPtr[index]-allmeanNew);
                    allstdInput+=(SubjectImagePtr[index]-allmeanInput)*(SubjectImagePtr[index]-allmeanInput);
                }
            }
            allstdNew=allstdNew/(float)(maskcount);
            allstdInput=allstdInput/(float)(maskcount);
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    ResampledAnatomicalPtr[index]=pow(((ResampledAnatomicalPtr[index]-allmeanNew)/sqrt(allstdNew))-((SubjectImagePtr[index]-allmeanInput)/sqrt(allstdInput)),2);
                }
                else
                {
                    ResampledAnatomicalPtr[index]=0.0f;
                }
            }

            radius=kernelstd_img;
            GaussianSmoothing(ResampledAnatomical,NULL, radius);

            // ResampledAnatomicalPtr will be  the current subject's weight.
            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                {
                    ResampledAnatomicalPtr[index]=(weightIntensity*ResampledAnatomicalPtr[index]+weightDisplacement*DisplacementPtr[index]);
                    ResampledAnatomicalPtr[index]=exp(-(ResampledAnatomicalPtr[index]/temperature));
                }
                else
                {
                    ResampledAnatomicalPtr[index]=0.0f;
                }
            }

            // Resample the probabilities, calculate the jacobian and modulate probabilities by the jacobian

            reg_resampleImage(InputProp,ResampProbability,Deformation,SubjectImageMaskPtr,1,0);
            float * ResampProbabilityPtr = static_cast<float *>(ResampProbability->data);
            if(modulate)
            {
                reg_defField_getJacobianMap(Deformation,jacobianImage);

                float * ResampJacPtr = static_cast<float *>(jacobianImage->data);

                for(int nt=0; nt<ResampProbability->nt; nt++)
                {
                    for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
                    {
                        if(SubjectImageMaskPtr==NULL || SubjectImageMaskPtr[index]>0)
                        {
                            ResampProbabilityPtr[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)]=ResampJacPtr[index]*ResampProbabilityPtr[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)];
                        }
                    }
                }
            }

            Gaussian_Filter_4D(ResampProbabilityPtr,fwhm/2.354820f,CurrSizes);

            for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
            {
                for(int nt=0; nt<ResampProbability->nt; nt++)
                {
                    Probs_group1[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)]+=
                        ResampledAnatomicalPtr[index]
                        *ResampProbabilityPtr[index+nt*(SubjectImage->nx*SubjectImage->ny*SubjectImage->nz)];
                }
                Weights_group1[index]+=ResampledAnatomicalPtr[index];
            }

            // Free images
            nifti_image_free(Displacement_Smooth);
            nifti_image_free(Displacement);
            nifti_image_free(InputCPP);
            nifti_image_free(InputProp);
            nifti_image_free(InputAnatomical);

        }
    }





    cout <<endl<<endl<<"\tFINISHING CALCULATION"<<endl;
    float * TmapPTR = static_cast<float *>(ResultTmap->data);
    for(int nt=0; nt<ResultTmap->nt; nt++)
    {
        for(int index=0; index<SubjectImage->nx*SubjectImage->ny*SubjectImage->nz; index++)
        {
            double mean_grp0=Probs_group0[index+(nt)*SubjectImage->nvox]/Weights_group0[index];
            double mean_grp1=Probs_group1[index+(nt)*SubjectImage->nvox]/Weights_group1[index];
            //TmapPTR[index+nt*SubjectImage->nvox]=(mean_grp0);
            TmapPTR[index+nt*SubjectImage->nvox]=(mean_grp0-mean_grp1);
        }
    }
    cout <<"\tSAVING"<<endl;
    delete [] Probs_group0;
    delete [] Probs_group1;
    delete [] Weights_group0;
    delete [] Weights_group1;

    nifti_image_write(ResultTmap);
    nifti_image_free(ResultTmap);
    nifti_image_free(ResampProbability);
    nifti_image_free(Deformation);
    nifti_image_free(SubjectImage);

    time(&end);
    int minutes = (int)floorf(float(end-start)/60.0f);
    int seconds = (int)(end-start - 60*minutes);
    if(verbose>0)
        cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;

    return 0;
}

