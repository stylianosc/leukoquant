#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
// [leukoquant patch] See _seg_GIF.h for why NiftyReg's headers need to come
// before NiftySeg's here (shared nifti1_io.h struct-collision guard --
// whichever copy is #included first is the one whose RNifti setup runs).
#include "_reg_resampling.h"
#include "Maths.hpp"
#include "_reg_tools.h"

#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_EM.h"

using namespace std;
#define SegPrecisionTYPE float


void Usage(char *exec)
{
    printf("\n  Usage ->\t%s -in <filename> -mask <filename> -out <filename>\n\n",exec);
    printf("  * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-in <filename>\t\t| Filename of the input image\n");
    printf("\t-sulc_out <filename>\t\t| Filename of the Sulcal image\n");
    printf("\t-seg_out <filename>\t\t| Filename of the Seg image\n");
    printf("\t-priors <filename>\t| Filename of the brain-mask of the input image\n");
    printf("\t-mask <filename>\t| Filename of the brain-mask of the input image\n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}


void no_memory ()
{
    cout << "Failed to allocate memory!\n";
    exit (1);
}



int main(int argc, char **argv)
{
    try
    {

        set_new_handler(no_memory);

        if (argc <=1)
        {
            Usage(argv[0]);
            return 0;
        }

        seg_EM_Params * segment_param = new seg_EM_Params [1]();

        // Set defaults (SEG_PARAM constructor sets everything to zero)
        segment_param->MRF_strength=0.4f;
        segment_param->bias_order=3;
        segment_param->flag_Bias=1;
        segment_param->flag_MRF=1;
        segment_param->MRF_strength=0.2;
        segment_param->flag_PV_model=1;
        segment_param->flag_SG_deli=1;
        segment_param->maxIteration=100;
        segment_param->minIteration=4;
        segment_param->OutliernessRatio=0.01;
        segment_param->verbose_level=1;

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
            else if((strcmp(argv[i], "-in") == 0 || strcmp(argv[i], "--in") == 0) && (i+1)<argc)
            {
                segment_param->filename_T1 = argv[++i];
                segment_param->flag_T1=1;
            }
            else if((strcmp(argv[i], "-seg_out") == 0 || strcmp(argv[i], "--seg_out") == 0) && (i+1)<argc)
            {
                segment_param->filename_out = argv[++i];
                segment_param->flag_out=1;
            }
            else if((strcmp(argv[i], "-sulc_out") == 0 || strcmp(argv[i], "--sulc_out") == 0) && (i+1)<argc)
            {
                segment_param->filename_out_outlier = argv[++i];
                segment_param->flag_out=1;
            }
            else if((strcmp(argv[i], "-mask") == 0 || strcmp(argv[i], "--mask") == 0) && (i+1)<argc)
            {
                segment_param->filename_mask=argv[++i];
                segment_param->flag_mask=1;
            }
            else if((strcmp(argv[i], "-priors") == 0 || strcmp(argv[i], "--priors") == 0) && (i+1)<argc)
            {
                segment_param->filename_priors= (char **) calloc(1,sizeof(char *));
                segment_param->filename_priors[0]=argv[++i];
                segment_param->flag_manual_priors=1;
                nifti_image * tmpread = nifti_image_read(segment_param->filename_priors[0],false);
                segment_param->numb_classes=tmpread->nt;
                if(segment_param->numb_classes<2)
                {
                    cout<<"Number of classes has to be bigger than 1";
                    return 0;
                }
                nifti_image_free(tmpread);
                segment_param->flag_priors4D=true;
            }
            else
            {
                fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
                Usage(argv[0]);
                return 1;
            }
        }

        if(!segment_param->flag_T1)
        {
            fprintf(stderr,"Err:\tThe T1 image name has to be defined.\n");
            Usage(argv[0]);
            return 1;
        }

        if(segment_param->flag_out==0)
        {
            fprintf(stderr,"Err:\tThe output image name has to be defined.\n");
            Usage(argv[0]);
            return 1;
        }

        // READING T1
        nifti_image * InputImage=nifti_image_read(segment_param->filename_T1,true);
        if(InputImage == NULL)
        {
            fprintf(stderr,"* Error when reading the T1 image: %s\n",segment_param->filename_T1);
            return 1;
        }
        if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
            seg_changeDatatype<SegPrecisionTYPE>(InputImage);

        InputImage->dim[4]=InputImage->nt=(InputImage->nt<1)?1:InputImage->nt;
        InputImage->dim[5]=InputImage->nu=(InputImage->nu<1)?1:InputImage->nu;
        if(InputImage->nu>1)
        {
            InputImage->dim[5]=InputImage->nu;
            InputImage->dim[4]=1;
        }
        else if(InputImage->nt>1)
        {
            InputImage->dim[5]=InputImage->nt;
            InputImage->dim[4]=1;
        }
        else if(InputImage->nt>1)
        {
            InputImage->dim[5]=1;
            InputImage->dim[4]=1;
        }
        InputImage->dim[0]=5;
        nifti_update_dims_from_array(InputImage);
        ImageSize * CurrSize = new ImageSize [1]();
        CurrSize->numel=(int)(InputImage->nx*InputImage->ny*InputImage->nz);
        CurrSize->xsize=InputImage->nx;
        CurrSize->ysize=InputImage->ny;
        CurrSize->zsize=InputImage->nz;
        CurrSize->usize=(InputImage->nu>1)?InputImage->nu:1;
        CurrSize->tsize=(InputImage->nt>1)?InputImage->nt:1;

        nifti_image * Mask=NULL;
        if(segment_param->flag_mask)
        {
            Mask = nifti_image_read(segment_param->filename_mask,true);
            if(Mask->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<SegPrecisionTYPE>(Mask);

            if(Mask == NULL)
            {
                fprintf(stderr,"* Error when reading the mask image: %s\n",segment_param->filename_mask);
                return 1;
            }
            if( (Mask->nx != InputImage->nx) || (Mask->ny != InputImage->ny) || (Mask->nz != InputImage->nz) )
            {
                fprintf(stderr,"* Error: Mask image not the same size as input\n");
                return 1;
            }
            Mask->dim[4]=Mask->nt=(Mask->nt<1)?1:Mask->nt;
            Mask->dim[5]=Mask->nu=(Mask->nu<1)?1:Mask->nu;
            nifti_update_dims_from_array(Mask);
        }


        nifti_image *  Priors = nifti_image_read(segment_param->filename_priors[0],true);
        if(Priors == NULL)
        {
            fprintf(stderr,"* Error when reading the Prior image: %s\n",segment_param->filename_priors[0]);
            return 1;
        }
        if( (Priors->nx != InputImage->nx) || (Priors->ny != InputImage->ny) || (Priors->nz != InputImage->nz) )
        {
            fprintf(stderr,"* Error: Prior image ( %s ) not the same size as input\n",segment_param->filename_priors[0]);
            return 1;
        }
        seg_changeDatatype<SegPrecisionTYPE>(Priors);


        seg_EM SEG(Priors->nt,InputImage->dim[5],InputImage->dim[4]);
        SEG.SetInputImage(InputImage);
        if(segment_param->flag_mask)
            SEG.SetMaskImage(Mask);
        if(segment_param->flag_manual_priors)
            SEG.SetPriorImage(Priors);

        SEG.SetVerbose(segment_param->verbose_level);
        SEG.SetFilenameOut(segment_param->filename_out);
        SEG.SetMaximalIterationNumber(segment_param->maxIteration);
        SEG.SetMinIterationNumber(segment_param->minIteration);
        SEG.SetRelaxation(0.5,1);

        if(segment_param->flag_Bias)
            SEG.SetBiasField(segment_param->bias_order,segment_param->Bias_threshold);
        if(segment_param->flag_MRF)
            SEG.SetMRF(segment_param->MRF_strength);
        SEG.Run_EM();

        if(segment_param->verbose_level>0)
        {
            cout << "Runing Sulci Detection."<<endl;
        }
        nifti_image * Result = SEG.GetResult();
        if(segment_param->verbose_level>0)
        {
            cout << "Saving Seg."<<endl;
        }
        nifti_image_write(Result);

        SegPrecisionTYPE * ResPtr = static_cast<SegPrecisionTYPE *>(Result->data);
        bool * MaskPtr = static_cast<bool *>(Mask->data);

        bool * NewImg=new bool [Result->nx*Result->ny*Result->nz];
        SegPrecisionTYPE * NewImgFloat=new SegPrecisionTYPE [Result->nx*Result->ny*Result->nz];
        ;
        for(int i=0; i<(Result->nx*Result->ny*Result->nz); i++)
        {
            NewImgFloat[i]=1;//+5*(ResPtr[i+2*(Result->nx*Result->ny*Result->nz)]+(MaskPtr[i]==0));
            NewImg[i]=(((ResPtr[i+2*(Result->nx*Result->ny*Result->nz)]+ResPtr[i+(Result->nx*Result->ny*Result->nz)])+(MaskPtr[i]==0)))>0.5;


        }
        GaussianFilter4D_cArray(NewImgFloat,1.5f,CurrSize);
        if(segment_param->verbose_level>0)
        {
            cout << "Calculating GeoDist."<<endl;
        }
        float * GeoDist=DoubleEuclideanDistance_3D(NewImg,NewImgFloat,CurrSize);
        GaussianFilter4D_cArray(GeoDist,0.7f,CurrSize);

        for(int i=0; i<(Result->nx*Result->ny*Result->nz); i++)
        {
            NewImgFloat[i]=0;
        }
        if(segment_param->verbose_level>0)
        {
            cout << "Finding Max."<<endl;
        }
        for(int nx=1; nx<(Result->nx-1); nx++)
        {
            for(int ny=1; ny<(Result->ny-1); ny++)
            {
                for(int nz=1; nz<(Result->nz-1); nz++)
                {


                    float curvx=pow(GeoDist[(nx-1)+Result->nx*(ny+nz*Result->ny)]-2*GeoDist[(nx)+Result->nx*(ny+nz*Result->ny)]+GeoDist[(nx+1)+Result->nx*(ny+nz*Result->ny)],1);
                    float curvy=pow(GeoDist[(nx)+Result->nx*((ny-1)+nz*Result->ny)]-2*GeoDist[(nx)+Result->nx*(ny+nz*Result->ny)]+GeoDist[nx+Result->nx*((ny+1)+nz*Result->ny)],1);
                    float curvz=pow(GeoDist[(nx)+Result->nx*(ny+(nz-1)*Result->ny)]-2*GeoDist[(nx)+Result->nx*(ny+nz*Result->ny)]+GeoDist[nx+Result->nx*(ny+(nz+1)*Result->ny)],1);

                    float gradx=(GeoDist[(nx-1)+Result->nx*(ny+nz*Result->ny)]-GeoDist[(nx+1)+Result->nx*(ny+nz*Result->ny)])/2;
                    float grady=(GeoDist[(nx)+Result->nx*((ny-1)+nz*Result->ny)]-GeoDist[(nx)+Result->nx*((ny+1)+nz*Result->ny)])/2;
                    float gradz=(GeoDist[(nx)+Result->nx*(ny+(nz-1)*Result->ny)]-GeoDist[(nx)+Result->nx*(ny+(nz+1)*Result->ny)])/2;
                    float gradnorm=sqrt(pow(gradx,2)+pow(grady,2)+pow(gradz,2));


                    NewImgFloat[(nx)+Result->nx*(ny+nz*Result->ny)]=-1*(curvx+curvy+curvz)*(1-gradnorm)*(double)(((curvx+curvy+curvz)*(1-gradnorm))<0)*(double)(MaskPtr[(nx)+Result->nx*(ny+nz*Result->ny)]);
                    NewImgFloat[(nx)+Result->nx*(ny+nz*Result->ny)]=(NewImgFloat[(nx)+Result->nx*(ny+nz*Result->ny)]>0.1)?NewImgFloat[(nx)+Result->nx*(ny+nz*Result->ny)]:0;
                }
            }
        }

        for(int i=0; i<(Result->nx*Result->ny*Result->nz); i++)
        {
            ResPtr[i]=NewImgFloat[i];
        }

        if(segment_param->verbose_level>0)
        {
            cout << "Saving Sulci."<<endl;
        }
        Result->dim[4]=1;
        nifti_set_filenames(Result, segment_param->filename_out_outlier,0,0);
        nifti_image_write(Result);
        nifti_image_free(Result);
        nifti_image_free(InputImage);
        nifti_image_free(Mask);

        free(segment_param->filename_priors);

        delete [] segment_param;
    }

    catch(std::exception & e)
    {
        std::cerr << "Standard exception: " << e.what() << std::endl;
    }

    catch(...)
    {
        std::cerr << "Unhandled Exception: Something went wrong! Please report the error to mjorgecardoso"<<(char) 64<<"gmail.com" << std::endl;
    }
    return 0;
}


