#include <time.h>
#include <string>
#include <iostream>
#include "nifti1_io.h"
#include "_reg_tools.h"
#ifdef _OPENMP
#include "omp.h"
#endif

using namespace std;
#define SegPrecisionTYPE float
#define MSVBM_MaxSampleSize 200
#define MSVBM_Permutations 10000

void Usage(char *exec)
{
    printf("\n  BrainMesh Segmentation :\n  Usage ->\t%s <mandatory_options>\n\n",exec);
    printf("  * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-inSVBM <filename> \t\t| 4D density differences of all subjects\n");
    printf("\t-inVBM <filename_grp0> <filename_grp1>\t| 4D densities of all subjects for group 0 and group 1\n");
    printf("\t-outPmap <filename>\t\t| Output prior segemntation filename\n");
    printf("\t-mask <filename>\t| Mask over the input image\n");
    printf("\t-v <int>\t\t| Verbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
    printf("\n  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    return;
}

float PermuteTest(float & SVBM_Vector,int SVBM_image_size,float & VBM_Vector,int VBM_Vector_grp0_size,int VBM_Vector_grp1_size, float * sVBM_PVAL_Pos, float * VBM_PVAL_Pos,float * sVBM_PVAL_Neg, float * VBM_PVAL_Neg);

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
    string SVBM_file;
    string VBM_file_grp0;
    string VBM_file_grp1;
    string OutputPmap_file;
    string Mask_file;
    int verbose=0;


    // read the input parameter

    for(int i=1; i<argc; i++)
    {
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-inSVBM") == 0 && (i+1)<argc)
        {
            SVBM_file=argv[++i];
        }
        else if(strcmp(argv[i], "-inVBM") == 0 && (i+1)<argc)
        {
            VBM_file_grp0=argv[++i];
            VBM_file_grp1=argv[++i];
        }
        else if(strcmp(argv[i], "-outPmap") == 0 && (i+1)<argc)
        {
            OutputPmap_file=argv[++i];
        }
        else if(strcmp(argv[i], "-mask") == 0 && (i+1)<argc)
        {
            Mask_file=argv[++i];
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


    // reading in the VBM and sVBM images
    nifti_image * SVBM_image=nifti_image_read(SVBM_file.c_str(),true);
    if(SVBM_image == NULL)
    {
        fprintf(stderr,"* Error when reading the image: %s\n",SVBM_file.c_str());
        return 1;
    }
    if(verbose>1)
        cout << "SVBM read ["<< SVBM_image->nvox<<"] -> Done"<<endl;
    float * SVBM_image_ptr = static_cast<float *>(SVBM_image->data);




    nifti_image * VBM_image_grp0=nifti_image_read(VBM_file_grp0.c_str(),true);
    if(VBM_image_grp0 == NULL)
    {
        fprintf(stderr,"* Error when reading the image: %s\n",VBM_file_grp0.c_str());
        return 1;
    }
    if(verbose>1)
        cout << "VBM_image_grp0 ["<< VBM_image_grp0->nvox<<"] read -> Done"<<endl;
    float * VBM_image_grp0_ptr = static_cast<float *>(VBM_image_grp0->data);




    nifti_image *VBM_image_grp1=nifti_image_read(VBM_file_grp1.c_str(),true);
    if(VBM_image_grp1 == NULL)
    {
        fprintf(stderr,"* Error when reading the image: %s\n",VBM_file_grp1.c_str());
        return 1;
    }
    if(verbose>1)
        cout << "VBM_image_grp1 ["<< VBM_image_grp1->nvox<<"] read -> Done"<<endl;
    float * VBM_image_grp1_ptr = static_cast<float *>(VBM_image_grp1->data);




    // reading in the subject image's mask
    nifti_image * Mask_image=NULL;
    if(Mask_file.size()>0)
    {
        Mask_image=nifti_image_read(Mask_file.c_str(),true);
        if(Mask_image == NULL)
        {
            fprintf(stderr,"* Error when reading the target mask image: %s\n",Mask_file.c_str());
            return 1;
        }
        if(Mask_image->datatype!=NIFTI_TYPE_UINT8)
            reg_tools_changeDatatype<unsigned char>(Mask_image);
    }
    unsigned char * Mask_image_ptr = static_cast<unsigned char *>(Mask_image->data);



    // Allocating space for the final Pmap
    nifti_image * ResultPmap=NULL;
    ResultPmap = nifti_copy_nim_info(SVBM_image);
    ResultPmap->dim[0]=4;
    ResultPmap->dim[4]=4;
    ResultPmap->datatype=NIFTI_TYPE_FLOAT32;
    ResultPmap->cal_max=1;
    nifti_set_filenames(ResultPmap,OutputPmap_file.c_str(),0,0);
    nifti_update_dims_from_array(ResultPmap);
    nifti_datatype_sizes(ResultPmap->datatype,&ResultPmap->nbyper,&ResultPmap->swapsize);
    ResultPmap->data = (void *) calloc(ResultPmap->nvox, sizeof(float));
    if(ResultPmap->data == NULL)
    {
        fprintf(stderr,"* Error when allocating ResultTmap. Not enough memory");
        return 1;
    }
    float * ResultPmap_ptr = static_cast<float *>(ResultPmap->data);
    if(verbose>1)
        cout << "ResultPmap alocated ["<< ResultPmap->nvox<<"] read -> Done"<<endl;


    for(int index=0; index<(SVBM_image->nx*SVBM_image->ny*SVBM_image->nz); index++)
    {
        float mean=0;
        float stdev=0;
        for(int n=0; n<SVBM_image->nt; n++)
            mean+=SVBM_image_ptr[index+n*(SVBM_image->nx*SVBM_image->ny*SVBM_image->nz)];
        mean/=(float)SVBM_image->nt;

        for(int n=0; n<SVBM_image->nt; n++)
            stdev+=pow(SVBM_image_ptr[index+n*(SVBM_image->nx*SVBM_image->ny*SVBM_image->nz)]-mean,2);
        stdev=sqrt(stdev/(float)SVBM_image->nt);

        ResultPmap_ptr[index]=pow(mean,3)/pow(stdev,3);

    }
    nifti_image_write(ResultPmap);
    nifti_image_free(ResultPmap);
    exit(1);


    int index,index2,n=0;
    int voxelNumber=ResultPmap->nx*ResultPmap->ny*ResultPmap->nz;
    float SVBM_Vector[MSVBM_MaxSampleSize];
    float VBM_Vector[MSVBM_MaxSampleSize];

    int activeVoxelNumber=0;
    for(index=0; index<voxelNumber; index++)
        if(Mask_image_ptr[index]>0)
            activeVoxelNumber++;

    printf("Number of active index: %i\n", activeVoxelNumber);
    int *activeVoxel=new int [activeVoxelNumber];
    for(index=0; index<voxelNumber; index++)
        if(Mask_image_ptr[index]>0)
            activeVoxel[n++]=index;




#ifdef _OPENMP
    int tnum=omp_get_max_threads();
    #pragma omp parallel for default(none) \
    private(index,index2, n, SVBM_Vector, VBM_Vector) \
    shared(ResultPmap,ResultPmap_ptr,activeVoxel,activeVoxelNumber, \
           SVBM_image,SVBM_image_ptr,VBM_image_grp0,VBM_image_grp1, \
           VBM_image_grp0_ptr,VBM_image_grp1_ptr,voxelNumber,tnum)
#endif // _OPENMP
    for(index=0; index<activeVoxelNumber; index++)
    {
        index2=activeVoxel[index];
        for(n=0; n<SVBM_image->nt; n++)
            SVBM_Vector[n]=SVBM_image_ptr[index2+n*voxelNumber];
        for(n=0; n<VBM_image_grp0->nt; n++)
            VBM_Vector[n]=VBM_image_grp0_ptr[index2+n*voxelNumber];
        for(n=0; n<VBM_image_grp1->nt; n++)
            VBM_Vector[n+VBM_image_grp0->nt]=VBM_image_grp1_ptr[index2+n*voxelNumber];

        PermuteTest(SVBM_Vector[0],SVBM_image->nt,VBM_Vector[0],VBM_image_grp0->nt,VBM_image_grp1->nt, &ResultPmap_ptr[index2], &ResultPmap_ptr[index2+voxelNumber], &ResultPmap_ptr[index2+2*voxelNumber], &ResultPmap_ptr[index2+3*voxelNumber]);

#ifdef _OPENMP
        if(omp_get_thread_num()==0)
            if(floor(((float)(index)/1000.0f))==((float)(index)/1000.0f))
                printf("Percent = %f\n", (float)(tnum*index)/(float)(activeVoxelNumber));
#else
        if(floor(((float)(index)/1000.0f))==((float)(index)/1000.0f))
            printf("Percent = %f\n", (float)(index)/(float)(activeVoxelNumber));

#endif // _OPENMP
    }


    delete []activeVoxel;
    nifti_image_write(ResultPmap);
    nifti_image_free(ResultPmap);
    nifti_image_free(SVBM_image);
    nifti_image_free(VBM_image_grp0);
    nifti_image_free(VBM_image_grp1);

    time(&end);
    int minutes = (int)floorf(float(end-start)/60.0f);
    int seconds = (int)(end-start - 60*minutes);
    if(verbose>0)
        cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;
    return 0;
}

float PermuteTest(float & SVBM_Vector,int SVBM_image_size,float & VBM_Vector,int VBM_Vector_grp0_size,int VBM_Vector_grp1_size,  float * VBM_PVAL_Pos, float * VBM_PVAL_Neg,float * sVBM_PVAL_Pos, float * sVBM_PVAL_Neg)
{

    unsigned int seed=time(NULL);
    float RandomSampling[MSVBM_Permutations]= {0.0f};
    RandomSampling[0]=RandomSampling[1];
    int partition=RAND_MAX/2;

    float * VBM_Vector_ptr=&VBM_Vector;
    float * SVBM_Vector_prt=&SVBM_Vector;

    float DiffVBM=0.0f;
    float mean_grp0=0.0f;
    for(int grp0_index=0; grp0_index<(VBM_Vector_grp0_size); grp0_index++ )
    {
        mean_grp0+=VBM_Vector_ptr[grp0_index];
    }
    float mean_grp1=0.0f;
    for(int grp1_index=0; grp1_index<(VBM_Vector_grp1_size); grp1_index++ )
    {
        mean_grp1+=VBM_Vector_ptr[VBM_Vector_grp0_size+grp1_index];
    }
    DiffVBM=(mean_grp0/VBM_Vector_grp0_size)-(mean_grp1/VBM_Vector_grp1_size);
    int Diffsup_count=0;
    int Diffinf_count=0;

    float grp0_mean, grp1_mean, grp0_samp, grp1_samp;
    for(int sample=0; sample<MSVBM_Permutations; sample++)
    {
        grp0_mean=0;
        grp1_mean=0;
        grp0_samp=0;
        grp1_samp=0;
        for(int subsample=0; subsample<(VBM_Vector_grp0_size+VBM_Vector_grp1_size); subsample++ )
        {

            if(rand_r(&
                      seed)>partition)
            {
                grp0_mean+=VBM_Vector_ptr[subsample];
                grp0_samp++;
            }
            else
            {
                grp1_mean+=VBM_Vector_ptr[subsample];
                grp1_samp++;
            }
        }
        grp0_mean/=grp0_samp;
        grp1_mean/=grp1_samp;
        RandomSampling[sample]=grp0_mean-grp1_mean;


        if(RandomSampling[sample]>DiffVBM)
        {
            Diffsup_count++;
        }

        if(RandomSampling[sample]<DiffVBM)
        {
            Diffinf_count++;
        }



    }


    float histogram_svbm[100]= {0.0f};
    float histogram_null[100]= {0.0f};

    for(int sample=0; sample<MSVBM_Permutations; sample++)
    {
        // 100 samples between -1 and 1
        float location=(RandomSampling[sample]>1)?99.0f:(RandomSampling[sample]<(-1)?0.0f:((RandomSampling[sample]+1.0f)*49.5f));
        histogram_null[(int)ceil(location)]+=(location-floor(location));
        histogram_null[(int)floor(location)]+=1-(location-floor(location));
    }

    float meanSVBM=0.0f;
    for(int sample=0; sample<SVBM_image_size; sample++)
    {
        // 100 samples between -1 and 1
        float location=(SVBM_Vector_prt[sample]>1)?99.0f:(SVBM_Vector_prt[sample]<(-1)?0.0f:((SVBM_Vector_prt[sample]+1.0f)*49.5f));
        histogram_svbm[(int)ceil(location)]+=(location-floor(location));
        histogram_svbm[(int)floor(location)]+=1-(location-floor(location));
        meanSVBM+=SVBM_Vector_prt[sample];
    }
    meanSVBM/=SVBM_image_size;


    float normfactor_null=0;
    float normfactor_svbm=0;
    for(int sample=0; sample<100; sample++)
    {
        // 100 samples between -1 and 1
        normfactor_null+=histogram_null[sample];
        normfactor_svbm+=histogram_svbm[sample];
    }

    for(int sample=0; sample<100; sample++)
    {
        // 100 samples between -1 and 1
        histogram_null[sample]/=normfactor_null;
        histogram_svbm[sample]/=normfactor_svbm;
    }

    sVBM_PVAL_Pos[0]=0.0f;
    sVBM_PVAL_Neg[0]=0.0f;
    if(meanSVBM>0)
    {
        sVBM_PVAL_Neg[0]=1;
        for(int sample=0; sample<100; sample++)
        {
            // 100 samples between -1 and 1
            sVBM_PVAL_Pos[0]+=histogram_null[sample]*histogram_svbm[sample];
        }
    }
    else
    {
        sVBM_PVAL_Pos[0]=1;
        for(int sample=0; sample<100; sample++)
        {
            // 100 samples between -1 and 1
            sVBM_PVAL_Neg[0]+=histogram_null[sample]*histogram_svbm[sample];
        }
    }




    VBM_PVAL_Pos[0]=1-(float)Diffsup_count/(float)MSVBM_Permutations;
    VBM_PVAL_Neg[0]=1-(float)Diffinf_count/(float)MSVBM_Permutations;

    return 0;

}


