#include "_seg_CT_colon.h"
#include <iostream>
#include <time.h>


void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -in <filename> -out <filename>.\n\n",exec);
    printf("\t* * Mandatory * *\n");
    printf("\t-in <filename>\t\tFilename of the input image segmentation\n");
    printf("\t-out <filename>\t\tFilename of the brainmask of the input image\n");

    printf("\t* * Options * *\n");
    printf("\t-max_iter <int>\t\t\tMaximum number of iterations (default = 100)\n");
    printf("\t-v <int>\t\t\tVerbose level [0 = off, 1 = verbose, 2 = debug] (default = 0)\n");
    printf("\t-mrf_beta <float>\t\tMRF prior strength [off = 0] (default = 0.2) \n");
    printf("\t-minvol <int>\t\t\tMinimum size of the structure to keep in liters (default = 0.25)\n");
    printf("\t-dilDir <int>\t\t\tInverse gravity direction (UP) [x=1,-x=-1,y=2,-y=-2,+z=3,-z=-3] (default=1)\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}


int main(int argc, char **argv)
{


    SEG_PARAM * segment_param = new SEG_PARAM [1]();
    segment_param->maxIteration=100;
    segment_param->flag_T1=0;
    segment_param->flag_out=0;
    segment_param->flag_mask=0;
    segment_param->flag_MRF=1;
    segment_param->flag_Bias=1;
    segment_param->flag_SG_deli=1;
    segment_param->flag_bc_out=0;
    segment_param->relax_factor=1;
    segment_param->flag_PV_model=0;
    segment_param->verbose_level=0;
    segment_param->flag_manual_priors=0;
    segment_param->bias_order=4;
    segment_param->MRF_strength=0.2f;
    segment_param->numb_classes=3;
    //stores the minimum volume of the object
    segment_param->relax_gauss_kernel=0.25;
    //stores the direction of dilation (minus the gravity direction - UP)
    segment_param->Bias_threshold=1;




    char * filename_out=NULL;
    char * filename_in=NULL;
    for(int i=1; i<argc; i++)
    {
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0)
        {
            filename_in=argv[++i];
        }
        else if(strcmp(argv[i], "-out") == 0)
        {
            filename_out = argv[++i];
        }
        else if(strcmp(argv[i], "-minvol") == 0)
        {
            //stores the minimum volume of the object
            segment_param->relax_gauss_kernel=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-dilDir") == 0)
        {
            //stores the direction of dilation (minus the gravity direction - UP)
            segment_param->Bias_threshold=(int)(atof(argv[++i]));
        }
        else if(strcmp(argv[i], "-v") == 0)
        {
            segment_param->verbose_level=(int)atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-mrf_beta") == 0)
        {
            segment_param->MRF_strength=(SegPrecisionTYPE)atof(argv[++i]);
            if(segment_param->bias_order==0)
            {
                segment_param->flag_MRF=0;
            }
        }
        else if(strcmp(argv[i], "-max_iter") == 0)
        {
            segment_param->maxIteration=atoi(argv[++i]);
        }
        else
        {
            fprintf(stderr,"Err:\tParameter %s unknown->\n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }


    if(filename_in == NULL)
    {
        fprintf(stderr,"* Error: No input defined\n");
        return 1;
    }

    nifti_image * CT=nifti_image_read(filename_in,true);
    if(CT == NULL)
    {
        fprintf(stderr,"* Error when reading the input Segmentation image\n");
        return 1;
    }
    seg_changeDatatype<SegPrecisionTYPE>(CT);


    if(filename_out == NULL)
    {
        fprintf(stderr,"* Error: No output defined\n");
        return 1;
    }

    nifti_image * Result = CT_colon(CT, filename_out, segment_param);

    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(CT);
    return 0;
}


