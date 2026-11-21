#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_FMM.h"
#include "_seg_EM.h"

using namespace std;
#define SegPrecisionTYPE float


void SmallUsage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf(" PBSI:\n  Usage .\t%s -img <filename_tp1> <filename_tp2> -seg <segmentation_tp1> <segmentation_tp1> [options]\n\n",exec);
    printf("\tSee the help for more details (-h).\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void Usage(char *exec)
{

    printf("\n PBSI:\n  Usage .\tseg_PBSI -img <filename_tp1> <filename_tp2> -seg <segmentation_tp1> <segmentation_tp1> [options]\n\n");
    printf("  * * * * * * * * * * * Mutually exclusive options * * * * * * * * * * * *\n\n");
    printf("  -clip <low_clip_tp1> <up_clip_tp1> <low_clip_tp2> <up_clip_tp2> \t\t The lower and upper clip intensities\n");
    printf("  -brain \t\t\t For brain BSI\n");
    printf("  -brainJoint \t\t\t For brain BSI\n");
    printf("  -hippo \t\t\t Hippo BSI (still not implemented)\n\n");
    printf("  * * * * * * * * * * * * * * * * * * Options * * * * * * * * * * * * * * * * *\n\n");
    printf("  -ker_size <int> \t\t The kernel size for the dilation and erosion operations. [default=1]\n");
    printf("  -dil <int> \t\t Dilation previous segmentation. [default=2]\n");
    printf("  -ero <int> \t\t Erotion previous segmentation. [default=0]\n");
    printf("  -window <float> \t\t Window cutoff percentage [between 0 and 0.4, default=0.2]. \n\n");
    printf("  -v <int>\t\t Verbose level [0 = off, 1 = on] (default = 0)\n");
    printf("  -save_diff <filename> \t\t The file name of the difference image. \n\n");
    printf("  -save_xor <filename> \t\t The file name of the xor image. \n\n");
    printf("  -save_seg <filename> \t\t The file name of the seg image. \n\n");
    printf("  -save_outlier <filename> \t\t The file name of the outlierness image. \n\n");
    printf("  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    return;
}


void volume(char * filename_Seg_TP1,char * filename_Seg_TP2,double BSI_VALUE)
{
    float VOLUME_VALUE_TP1=0;
    float VOLUME_VALUE_TP2=0;
    float PROB_DIFF=0;
    int i;


    nifti_image * Seg_TP1=nifti_image_read(filename_Seg_TP1,true);
    seg_changeDatatype<float>(Seg_TP1);

    nifti_image * Seg_TP2=nifti_image_read(filename_Seg_TP2,true);
    seg_changeDatatype<float>(Seg_TP2);

    float * Seg_TP1_ptr=static_cast<float *>(Seg_TP1->data);
    float * Seg_TP2_ptr=static_cast<float *>(Seg_TP2->data);

    for(i=0; i<(Seg_TP1->nx*Seg_TP1->ny*Seg_TP1->nz); i++)
    {
        VOLUME_VALUE_TP1 +=Seg_TP1_ptr[i];
        VOLUME_VALUE_TP2 +=Seg_TP2_ptr[i];
        PROB_DIFF +=Seg_TP1_ptr[i]-Seg_TP2_ptr[i];
    }
    VOLUME_VALUE_TP1 *=(Seg_TP1->dx*Seg_TP1->dy*Seg_TP1->dz);
    VOLUME_VALUE_TP2 *=(Seg_TP2->dx*Seg_TP2->dy*Seg_TP2->dz);
    PROB_DIFF *=(Seg_TP2->dx*Seg_TP2->dy*Seg_TP2->dz);
    cout << "PBVC=[" << (BSI_VALUE*100/VOLUME_VALUE_TP1) << "," << (BSI_VALUE*100/VOLUME_VALUE_TP2) <<"]%"<<endl;
    cout << "VOLUME TP1=" << VOLUME_VALUE_TP1 << "mm3 " << endl;
    cout << "VOLUME TP2=" << VOLUME_VALUE_TP2 << "mm3"<<endl;
    cout << "DIFF VOL=" << VOLUME_VALUE_TP1-VOLUME_VALUE_TP2 << "mm3"<<endl;
    cout << "PROB DIFF=" << PROB_DIFF << "mm3"<<endl;
    cout << "VOXEL VOLUME=" << (Seg_TP1->dx*Seg_TP1->dy*Seg_TP1->dz) << "mm3"<<endl;
    cout << "TOTAL VOXELS=" << i << endl;
}


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        Usage(argv[0]);
        return 0;
    }

    char * filename_Image_TP1=NULL;
    char * filename_Image_TP2=NULL;
    char * filename_Seg_TP1=NULL;
    char * filename_Seg_TP2=NULL;
    char * filename_DIFF=NULL;
    char * filename_XOR=NULL;
    char * filename_SEG_OUT=NULL;
    char * filename_Outlier_OUT=NULL;
    float window=0.2;
    float clip_low=0;
    float clip_high=0;
    float clip_low_tp1=0;
    float clip_high_tp1=0;
    float clip_low_tp2=0;
    float clip_high_tp2=0;
    int ker_size=1;
    int mode=-1;
    int verbose=0;
    int dil=4;
    int ero=2;

    for(int i=1; i<argc; i++)
    {
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0 ||
                strcmp(argv[i], "--HELP")==0)
        {
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-img") == 0 && (i+2)<argc)
        {
            filename_Image_TP1 = argv[++i];
            filename_Image_TP2 = argv[++i];
        }
        else if(strcmp(argv[i], "-v") == 0 && (i+1)<argc)
        {
            verbose = atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-seg") == 0  && (i+2)<argc)
        {
            filename_Seg_TP1 = argv[++i];
            filename_Seg_TP2 = argv[++i];
        }
        else if(strcmp(argv[i], "-clip") == 0)
        {
            mode=0;
            clip_low = atof(argv[++i]);
            clip_high = atof(argv[++i]);
            clip_low_tp2 = atof(argv[++i]);
            clip_high_tp2 = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-window") == 0)
        {
            window = atof(argv[++i]);
            window=(window>0.4)?0.4:((window<0)?0:window);
        }
        else if(strcmp(argv[i], "-dil") == 0)
        {
            dil = atoi(argv[++i]);
            dil=dil<0?0:dil;
        }
        else if(strcmp(argv[i], "-ero") == 0)
        {
            ero = atoi(argv[++i]);
            ero=ero<0?0:ero;
        }
        else if(strcmp(argv[i], "-brain") == 0)
        {
            mode=1;
        }
        else if(strcmp(argv[i], "-brainJoint") == 0)
        {
            mode=2;
        }
        else if(strcmp(argv[i], "-hippo") == 0)
        {
            mode=3;
        }
        else if(strcmp(argv[i], "-ker_size") == 0)
        {
            ker_size = atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-save_diff") == 0)
        {
            filename_DIFF = argv[++i];
        }
        else if(strcmp(argv[i], "-save_xor") == 0)
        {
            filename_XOR = argv[++i];
        }
        else if(strcmp(argv[i], "-save_seg") == 0)
        {
            filename_SEG_OUT = argv[++i];
        }
        else if(strcmp(argv[i], "-save_outlier") == 0)
        {
            filename_Outlier_OUT = argv[++i];
        }
        else
        {
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }


    if(filename_Image_TP1!=NULL && filename_Image_TP2!=NULL && filename_Seg_TP1!=NULL && filename_Seg_TP2!=NULL)
    {

        // READING Images
        nifti_image * Image_TP1=nifti_image_read(filename_Image_TP1,true);
        if(Image_TP1 == NULL)
        {
            fprintf(stderr,"* Error when reading the input image: %s\n", filename_Image_TP1);
            return 1;
        }
        seg_changeDatatype<float>(Image_TP1);

        nifti_image * Image_TP2=nifti_image_read(filename_Image_TP2,true);
        if(Image_TP2 == NULL)
        {
            fprintf(stderr,"* Error when reading the input image: %s\n", filename_Image_TP2);
            return 1;
        }
        seg_changeDatatype<float>(Image_TP2);
        if(Image_TP1->nx != Image_TP2->nx && Image_TP1->ny != Image_TP2->ny && Image_TP1->nz != Image_TP2->nz)
        {
            fprintf(stderr,"* Error: the input images are not the same size");
            return 1;
        }
        //*********************************
        //***** READING Probabilities *****
        //*********************************
        nifti_image * Seg_TP1=nifti_image_read(filename_Seg_TP1,true);
        if(Seg_TP1 == NULL)
        {
            fprintf(stderr,"* Error when reading the input image: %s\n", filename_Seg_TP1);
            return 1;
        }
        seg_changeDatatype<float>(Seg_TP1);
        if(Image_TP1->nx != Seg_TP1->nx && Image_TP1->ny != Seg_TP1->ny && Image_TP1->nz != Seg_TP1->nz)
        {
            fprintf(stderr,"* Error: the input images are not the same size");
            return 1;
        }

        nifti_image * Seg_TP2=nifti_image_read(filename_Seg_TP2,true);
        if(Seg_TP2 == NULL)
        {
            fprintf(stderr,"* Error when reading the input image: %s\n", filename_Seg_TP2);
            return 1;
        }
        seg_changeDatatype<float>(Seg_TP2);
        if(Image_TP1->nx != Seg_TP2->nx && Image_TP1->ny != Seg_TP2->ny && Image_TP1->nz != Seg_TP2->nz)
        {
            fprintf(stderr,"* Error: the input images are not the same size");
            return 1;
        }
        //*********************************
        // BSI starts
        // ********************************

        ImageSize * CurrSize = new ImageSize [1]();
        CurrSize->numel=(int)(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz);
        CurrSize->xsize=Image_TP1->nx;
        CurrSize->ysize=Image_TP1->ny;
        CurrSize->zsize=Image_TP1->nz;
        CurrSize->usize=(Image_TP1->nu>1)?Image_TP1->nu:1;
        CurrSize->tsize=(Image_TP1->nt>1)?Image_TP1->nt:1;

        float * Image_TP1_ptr=static_cast<float *>(Image_TP1->data);
        float * Image_TP2_ptr=static_cast<float *>(Image_TP2->data);
        float * Seg_TP1_ptr=static_cast<float *>(Seg_TP1->data);
        float * Seg_TP2_ptr=static_cast<float *>(Seg_TP2->data);

        float * Seg_TP1_union_TP2=new float [Image_TP1->nx*Image_TP1->ny*Image_TP1->nz];
        float * Seg_TP1_intersection_TP2=new float [Image_TP1->nx*Image_TP1->ny*Image_TP1->nz];

        for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
        {
            Seg_TP1_union_TP2[i]=max(Seg_TP1_ptr[i],Seg_TP2_ptr[i]);
            Seg_TP1_intersection_TP2[i]=min(Seg_TP1_ptr[i],Seg_TP2_ptr[i]);
        }

        // SEGMENTATION SECTION

        Dillate(Seg_TP1_union_TP2, ker_size, CurrSize);
        Erosion(Seg_TP1_intersection_TP2, ker_size, CurrSize);

        float * Seg_XOR=new float [Image_TP1->nx*Image_TP1->ny*Image_TP1->nz];

        for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
        {
            float a=Seg_TP1_ptr[i];//Seg_TP1_union_TP2[i];
            float b=Seg_TP2_ptr[i];//Seg_TP1_intersection_TP2[i];
            float p1=(a*(1-b));
            float p2=((1-a)*(b));
            Seg_XOR[i]=(p1+p2-p1*p2);
        }


        // VARIOUS MODES

        // IF IT IS A BRAIN
        if(mode==1)
        {
            if(verbose>0)
            {
                cout << "Segmenting TP1";
            }
            float * means=NULL;
            float * stds=NULL;
            nifti_image * Tmp_IMG_TP1=nifti_image_read(filename_Image_TP1,true);
            seg_changeDatatype<float>(Tmp_IMG_TP1);
            nifti_image * Tmp_Seg_TP1=nifti_image_read(filename_Seg_TP1,true);
            seg_changeDatatype<float>(Tmp_Seg_TP1);
            float * Tmp_Seg_TP1_ptr=static_cast<float *>(Tmp_Seg_TP1->data);
            for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
            {
                Seg_TP1_ptr[i]=Seg_TP1_ptr[i]>0.5;
            }

            if(dil>0) Dillate(Seg_TP1_ptr, dil, CurrSize);
            Close_Forground_ConnectComp<float,float>(Seg_TP1_ptr,Tmp_Seg_TP1_ptr,CurrSize);
            if (ero>0) Erosion(Seg_TP1_ptr, ero, CurrSize);

            seg_EM * SEG= new seg_EM(3,1,1);
            SEG->SetInputImage(Tmp_IMG_TP1);
            SEG->SetMaskImage(Tmp_Seg_TP1);
            SEG->SetVerbose(verbose);
            SEG->SetMinIterationNumber(5);
            SEG->SetMaximalIterationNumber(20);
            SEG->SetRegValue(1);
            SEG->SetOutlierness(2.5f, 0.01f);
            SEG->SetMRF(0.4f);
            if(filename_SEG_OUT!=NULL)
            {
                SEG->SetFilenameOut(filename_SEG_OUT);
            }
            SEG->Run_EM();
            nifti_image * OUTLIER = NULL;
            if(filename_Outlier_OUT!=NULL)
            {
                string tp1_name=filename_Outlier_OUT;
                tp1_name="tp1_"+tp1_name;
                OUTLIER = SEG->GetOutlierness((char *)tp1_name.c_str());
                if(verbose>0)
                {
                    cout << "Saving outlierness" <<endl;
                }
                nifti_image_write(OUTLIER);
                nifti_image_free(OUTLIER);
            }
            OUTLIER = NULL;
            nifti_image * Result = NULL;
            if(filename_SEG_OUT!=NULL)
            {
                Result = SEG->GetResult();
                if(verbose>0)
                {
                    cout << "Saving segmentation" <<endl;
                }
                string tp1_name=filename_SEG_OUT;
                tp1_name="tp1_"+tp1_name;
                nifti_set_filenames(Result, tp1_name.c_str(),0,0);
                nifti_image_write(Result);
                nifti_image_free(Result);
            }
            Result = NULL;
            if(verbose>0)
            {
                cout << " -> Done" <<endl;
            }
            means=SEG->GetMeans();
            stds=SEG->GetSTD();
            SEG->~seg_EM();
            nifti_image_free(Tmp_IMG_TP1);
            nifti_image_free(Tmp_Seg_TP1);
            if(verbose>0)
            {
                cout << "Means TP1=["<<means[0]<<","<<means[1]<<","<<means[2]<<"]"<<endl;
            }
            if(verbose>0)
            {
                cout << "STD TP1=["<<stds[0]<<","<<stds[1]<<","<<stds[2]<<"]"<<endl;
            }
            clip_low=means[0];
            clip_high=means[1];
            clip_low_tp1=clip_low+window*(clip_high-clip_low);
            clip_high_tp1=clip_high-window*(clip_high-clip_low);


            if(verbose>0)
            {
                cout << endl<< "Segmenting TP2";
            }
            float * mean_TP2=NULL;
            float * std_TP2=NULL;
            nifti_image * Tmp_IMG_TP2=nifti_image_read(filename_Image_TP2,true);
            seg_changeDatatype<float>(Tmp_IMG_TP2);
            nifti_image * Tmp_Seg_TP2=nifti_image_read(filename_Seg_TP2,true);
            seg_changeDatatype<float>(Tmp_Seg_TP2);
            float * Tmp_Seg_TP2_ptr=static_cast<float *>(Tmp_Seg_TP2->data);
            for(int i=0; i<(Image_TP2->nx*Image_TP2->ny*Image_TP2->nz); i++)
            {
                Seg_TP2_ptr[i]=Seg_TP2_ptr[i]>0.5;
            }

            if(dil>0) Dillate(Seg_TP2_ptr,dil, CurrSize);
            Close_Forground_ConnectComp<float,float>(Seg_TP2_ptr,Tmp_Seg_TP2_ptr,CurrSize);
            if (ero>0) Erosion(Seg_TP1_ptr, ero, CurrSize);

            seg_EM * SEG2= new seg_EM(3,1,1);
            SEG2->SetInputImage(Tmp_IMG_TP2);
            SEG2->SetMaskImage(Tmp_Seg_TP2);
            SEG2->SetVerbose(verbose);
            SEG2->SetMinIterationNumber(5);
            SEG2->SetMaximalIterationNumber(20);
            SEG2->SetRegValue(1);
            SEG2->SetOutlierness(2.5f, 0.01f);
            SEG2->SetMRF(0.4f);
            SEG2->Run_EM();
            if(verbose>0)
            {
                cout << " -> Done" <<endl;
            }
            if(filename_Outlier_OUT!=NULL)
            {
                string tp2_name=filename_Outlier_OUT;
                tp2_name="tp2_"+tp2_name;
                OUTLIER = SEG2->GetOutlierness((char *)tp2_name.c_str());
                if(verbose>0)
                {
                    cout << "Saving outlierness" <<endl;
                }
                nifti_image_write(OUTLIER);
                nifti_image_free(OUTLIER);
            }
            OUTLIER = NULL;
            if(filename_SEG_OUT!=NULL)
            {
                Result = SEG2->GetResult();
                if(verbose>0)
                {
                    cout << "Saving segmentation" <<endl;
                }
                string tp2_name=filename_SEG_OUT;
                tp2_name="tp2_"+tp2_name;
                nifti_set_filenames(Result, tp2_name.c_str(),0,0);
                nifti_image_write(Result);
                nifti_image_free(Result);
            }
            Result = NULL;

            mean_TP2=SEG2->GetMeans();
            std_TP2=SEG2->GetSTD();

            SEG2->~seg_EM();
            nifti_image_free(Tmp_IMG_TP2);
            nifti_image_free(Tmp_Seg_TP2);
            if(verbose>0)
            {
                cout << "Means TP2=["<<mean_TP2[0]<<","<<mean_TP2[1]<<","<<mean_TP2[2]<<"]"<<endl;
            }
            if(verbose>0)
            {
                cout << "STD TP2=["<<std_TP2[0]<<","<<std_TP2[1]<<","<<std_TP2[2]<<"]"<<endl;
            }



            clip_low=mean_TP2[0];
            clip_high=mean_TP2[1];
            clip_low_tp2=clip_low+window*(clip_high-clip_low);
            clip_high_tp2=clip_high-window*(clip_high-clip_low);
        }
        if(mode==2)
        {
            if(verbose>0)
            {
                cout << "Segmenting TP1 and TP2"<<endl;
            }
            float * means=NULL;
            float * stds=NULL;
            nifti_image * Tmp_IMG_Joint=nifti_copy_nim_info(Image_TP1);
            //seg_changeDatatype<float>(Tmp_IMG_Joint);
            Tmp_IMG_Joint->dim[0]=5;
            Tmp_IMG_Joint->dim[4]=1;
            Tmp_IMG_Joint->dim[5]=2;
            Tmp_IMG_Joint->datatype=DT_FLOAT32;
            nifti_update_dims_from_array(Tmp_IMG_Joint);
            nifti_datatype_sizes(Tmp_IMG_Joint->datatype,&Tmp_IMG_Joint->nbyper,&Tmp_IMG_Joint->swapsize);
            //free(Tmp_IMG_Joint->data);
            Tmp_IMG_Joint->data = (void *) calloc(Tmp_IMG_Joint->nvox, sizeof(SegPrecisionTYPE));
            float * Tmp_Image_Joint_ptr = static_cast<float *>(Tmp_IMG_Joint->data);
            for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
            {
                Tmp_Image_Joint_ptr[i]=(Image_TP1_ptr[i]);
                Tmp_Image_Joint_ptr[i+(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz)]=(Image_TP2_ptr[i]);
            }

            nifti_image * Tmp_Seg_Joint=nifti_image_read(filename_Seg_TP1,true);
            seg_changeDatatype<float>(Tmp_Seg_Joint);
            float * Tmp_Seg_Joint_ptr=static_cast<float *>(Tmp_Seg_Joint->data);


            for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
            {
                Seg_TP1_ptr[i]=(Seg_TP1_ptr[i]+Seg_TP2_ptr[i])/2;
            }
            if(dil>0) Dillate(Seg_TP1_ptr, dil, CurrSize);
            for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
            {
                Seg_TP1_ptr[i]=Seg_TP1_ptr[i]>0.5;
            }

            Close_Forground_ConnectComp<float,float>(Seg_TP1_ptr,Tmp_Seg_Joint_ptr,CurrSize);
            if (ero>0) Erosion(Seg_TP1_ptr, ero, CurrSize);

            seg_EM * SEG= new seg_EM(3,1,1);

            SEG->SetInputImage(Tmp_IMG_Joint);
            SEG->SetMaskImage(Tmp_Seg_Joint);
            SEG->SetVerbose(verbose);
            SEG->SetMinIterationNumber(5);
            SEG->SetMaximalIterationNumber(20);
            SEG->SetOutlierness(2.5f, 0.01f);
            SEG->SetMRF(0.4f);

            SEG->Run_EM();

            if(filename_SEG_OUT!=NULL)
            {
                SEG->SetFilenameOut(filename_SEG_OUT);
            }
            if(verbose>0)
            {
                cout << "Finished Segmentation" <<endl;
            }
            means=SEG->GetMeans();
            stds=SEG->GetSTD();


            nifti_image * OUTLIER = NULL;
            if(filename_Outlier_OUT!=NULL)
            {
                OUTLIER = SEG->GetOutlierness(filename_Outlier_OUT);
                if(verbose>0)
                {
                    cout << "Saving outlierness" <<endl;
                }
                nifti_image_write(OUTLIER);
                nifti_image_free(OUTLIER);
            }
            nifti_image * Result = NULL;
            if(filename_SEG_OUT!=NULL)
            {
                Result = SEG->GetResult();
                if(verbose>0)
                {
                    cout << "Saving segmentation" <<endl;
                }
                nifti_image_write(Result);
                nifti_image_free(Result);
            }

            SEG->~seg_EM();
            nifti_image_free(Tmp_IMG_Joint);
            nifti_image_free(Tmp_Seg_Joint);
            if(verbose>0)
            {
                cout << "Means TP1=["<<means[0]<<","<<means[1]<<","<<means[2]<<"]"<<endl;
            }
            if(verbose>0)
            {
                cout << "Means TP2=["<<means[3]<<","<<means[4]<<","<<means[5]<<"]"<<endl;
            }
            if(verbose>0)
            {
                cout << "STD TP1=["<<stds[0]<<","<<stds[1]<<","<<stds[2]<<"]"<<endl;
            }
            if(verbose>0)
            {
                cout << "STD TP2=["<<stds[3]<<","<<stds[4]<<","<<stds[5]<<"]"<<endl;
            }

            clip_low=means[0];
            clip_high=means[1];
            clip_low_tp1=clip_low+window*(clip_high-clip_low);
            clip_high_tp1=clip_high-window*(clip_high-clip_low);
            clip_low=means[3];
            clip_high=means[4];
            clip_low_tp2=clip_low+window*(clip_high-clip_low);
            clip_high_tp2=clip_high-window*(clip_high-clip_low);
        }
        else if(mode==-1 || mode>2)
        {
            // If no mode is given, then find the min and max values within the current XOR area

            //find_minimum_of_two_images
            float currmin_TP1=1.0e30;
            float currmin_TP2=1.0e30;
            for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
            {
                if(Seg_XOR[i]>0.05)
                {
                    currmin_TP1=(Image_TP1_ptr[i]<currmin_TP1)?Image_TP1_ptr[i]:currmin_TP1;
                    currmin_TP2=(Image_TP2_ptr[i]<currmin_TP2)?Image_TP2_ptr[i]:currmin_TP2;
                }
            }
            //find_maximum_of_two_images
            float currmax_TP1=-1.0e30;
            float currmax_TP2=-1.0e30;
            for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
            {
                if(Seg_XOR[i]>0.05)
                {
                    currmax_TP1=(Image_TP1_ptr[i]>currmax_TP1)?Image_TP1_ptr[i]:currmax_TP1;
                    currmax_TP2=(Image_TP2_ptr[i]>currmax_TP2)?Image_TP2_ptr[i]:currmax_TP2;
                }
            }
            clip_low_tp1=currmin_TP1;
            clip_high_tp1=currmax_TP1;
            clip_low_tp2=currmin_TP2;
            clip_high_tp2=currmax_TP2;

        }


        if(verbose>0)
        {
            cout << "BSI Calculation:"<<endl<<"Clip values = [" <<clip_low_tp1<<","<<clip_high_tp1 <<"]"<<" and [" <<clip_low_tp2 << ","<<clip_high_tp2<<"]\n";
        }
        // ***********************        CALCULATE BSI VALUE       *******************
        double BSI_VALUE=0;
        float * Seg_DIFF=new float [Image_TP1->nx*Image_TP1->ny*Image_TP1->nz];

        float normalisation_tp1=(clip_high_tp1-clip_low_tp1);
        float normalisation_tp2=(clip_high_tp2-clip_low_tp2);

        for(int i=0; i<(Image_TP1->nx*Image_TP1->ny*Image_TP1->nz); i++)
        {
            float cur_t1_intensity=Image_TP1_ptr[i]>clip_high_tp1?clip_high_tp1:(Image_TP1_ptr[i]<clip_low_tp1?clip_low_tp1:Image_TP1_ptr[i]);
            cur_t1_intensity=(cur_t1_intensity-clip_low_tp1)/normalisation_tp1;

            float cur_t2_intensity=Image_TP2_ptr[i]>clip_high_tp2?clip_high_tp2:(Image_TP2_ptr[i]<clip_low_tp2?clip_low_tp2:Image_TP2_ptr[i]);
            cur_t2_intensity=(cur_t2_intensity-clip_low_tp2)/normalisation_tp2;

            float xorval=Seg_XOR[i];

            BSI_VALUE +=xorval*(cur_t1_intensity-cur_t2_intensity);
            Seg_DIFF[i]=xorval*(cur_t1_intensity-cur_t2_intensity);
        }
        BSI_VALUE *=(Image_TP1->dx*Image_TP1->dy*Image_TP1->dz);
        cout << fixed;
        cout << "BSI=" << BSI_VALUE << "mm3"<<endl;
        volume(filename_Seg_TP1,filename_Seg_TP2,BSI_VALUE);

        if(filename_DIFF!=NULL)
        {
            nifti_image * DIFF_NII = nifti_copy_nim_info(Image_TP1);
            DIFF_NII->datatype=DT_FLOAT32;
            nifti_set_filenames(DIFF_NII,filename_DIFF,0,0);
            nifti_update_dims_from_array(DIFF_NII);
            nifti_datatype_sizes(DIFF_NII->datatype,&DIFF_NII->nbyper,&DIFF_NII->swapsize);
            DIFF_NII->data = Seg_DIFF;
            Seg_DIFF=NULL;
            nifti_image_write(DIFF_NII);
            nifti_image_free(DIFF_NII);
        }
        else
        {
            delete [] Seg_DIFF;
        }

        if(filename_XOR!=NULL)
        {
            nifti_image * XOR_NII = nifti_copy_nim_info(Image_TP1);
            XOR_NII->datatype=DT_FLOAT32;
            nifti_set_filenames(XOR_NII,filename_XOR,0,0);
            nifti_update_dims_from_array(XOR_NII);
            nifti_datatype_sizes(XOR_NII->datatype,&XOR_NII->nbyper,&XOR_NII->swapsize);
            XOR_NII->data = Seg_XOR;
            Seg_XOR=NULL;
            nifti_image_write(XOR_NII);
            nifti_image_free(XOR_NII);
        }
        else
        {
            delete [] Seg_XOR;
        }

        Image_TP1_ptr=NULL;
        nifti_image_free(Image_TP1);
        Image_TP2_ptr=NULL;
        nifti_image_free(Image_TP2);
        Seg_TP1_ptr=NULL;
        nifti_image_free(Seg_TP1);
        Seg_TP2_ptr=NULL;
        nifti_image_free(Seg_TP2);
        delete [] Seg_TP1_union_TP2;
        delete [] Seg_TP1_intersection_TP2;
    }



    return 0;
}



