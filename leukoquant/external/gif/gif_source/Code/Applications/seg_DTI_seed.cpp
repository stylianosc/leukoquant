#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s <input> <[0,1]=[seed,waypoint]> <output>.\n\n",exec);
    printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}
int main(int argc, char **argv)
{

    if (argc != 5)
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


    char * filename_in=argv[1];
    nifti_image * InputImage=nifti_image_read(filename_in,true);
    if(InputImage == NULL)
    {
        fprintf(stderr,"* Error when reading the input Segmentation image\n");
        return 1;
    }
    if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
    {
        seg_changeDatatype<SegPrecisionTYPE>(InputImage);
    }
    SegPrecisionTYPE * InputImagePtr = static_cast<SegPrecisionTYPE *>(InputImage->data);
    float * imageslice=new float [(InputImage->nx)*(InputImage->nz)];
    int planeindex=0;
    float maxsum=-10000.0f;

    bool mode=atoi(argv[2]);
    int count_in=atoi(argv[3]);

    // Getting best plane
    for(int iy=0; iy<InputImage->ny; iy++)
    {
        for(int ix=0; ix<InputImage->nx; ix++)
        {
            for(int iz=0; iz<InputImage->nz; iz++)
            {
                imageslice[ix+iz*InputImage->nx]=InputImagePtr[ix+(iy+iz*InputImage->ny)*InputImage->nx];
            }
        }

        int * order=quickSort_order(imageslice,(InputImage->nx)*(InputImage->nz));
        delete[]order;
        float cursum=0;


        for(int i=0; i<count_in; i++)
        {
            cursum+=imageslice[(InputImage->nx)*(InputImage->nz)-1-i];
        }
        if(cursum>maxsum)
        {
            planeindex=iy;
            maxsum=cursum;
        }
    }


    char * filename_out=argv[argc-1];
    nifti_image * OutputImage = nifti_copy_nim_info(InputImage);
    OutputImage->datatype=NIFTI_TYPE_UINT8;
    nifti_set_filenames(OutputImage,filename_out,0,0);
    OutputImage->dim[4]=OutputImage->nt=1;
    OutputImage->dim[5]=OutputImage->nu=1;
    OutputImage->dim[6]=OutputImage->nv=1;
    OutputImage->dim[7]=OutputImage->nw=1;
    OutputImage->dim[0]=(int)(OutputImage->dim[1]>1)+(int)(OutputImage->dim[2]>1)+(int)(OutputImage->dim[3]>1)+(int)(OutputImage->dim[4]>1)+(int)(OutputImage->dim[5]>1)+(int)(OutputImage->dim[6]>1)+(int)(OutputImage->dim[7]>1);
    nifti_update_dims_from_array(OutputImage);
    nifti_datatype_sizes(OutputImage->datatype,&OutputImage->nbyper,&OutputImage->swapsize);

    OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned char));
    unsigned char * OutputImagePtr = static_cast<unsigned char *>(OutputImage->data);
    for(unsigned int i=0; i<(OutputImage->nvox); i++)
    {
        OutputImagePtr[i]=0;
    }


    if(mode)
    {

        for(int ix=0; ix<InputImage->nx; ix++)
        {
            for(int iz=0; iz<InputImage->nz; iz++)
            {
                imageslice[ix+iz*InputImage->nx]=InputImagePtr[ix+(planeindex+iz*InputImage->ny)*InputImage->nx];
            }
        }
        int * order=quickSort_order(imageslice,InputImage->nx*InputImage->nz);

        for(int besti=0; besti<count_in; besti++)
        {
            int index_xz=order[InputImage->nx*InputImage->nz-1-besti];
            int indexz=floor((float)(index_xz)/(float)(InputImage->nx));
            int indexx=index_xz-floor((float)(index_xz)/(float)(InputImage->nx))*InputImage->nx;
            //cout <<index_xz<<"->("<< indexx<<","<< indexz<<")"<<endl;
            OutputImagePtr[indexx+(planeindex+indexz*InputImage->ny)*InputImage->nx]=1;

        }

        delete[]order;

    }
    else
    {

        //getting best index

        int bestindex_x=0;
        int bestindex_z=0;
        float bestindex_val=-100000.0f;
        for(int ix=0; ix<InputImage->nx; ix++)
        {
            for(int iz=0; iz<InputImage->nz; iz++)
            {
                imageslice[ix+iz*InputImage->nx]=InputImagePtr[ix+(planeindex+iz*InputImage->ny)*InputImage->nx];
                if(imageslice[ix+iz*InputImage->nx]>bestindex_val)
                {
                    bestindex_x=ix;
                    bestindex_z=iz;
                    bestindex_val=imageslice[ix+iz*InputImage->nx];
                }
            }
        }

        //int * order=quickSort_order(imageslice,InputImage->nx*InputImage->nz);
        //int bestindex=order[(InputImage->nz*InputImage->nx)-1];
        //delete[]order;

        cout << "best plane = "<<planeindex<<" with a sum distance of "<<maxsum<<endl;
        cout << "best index in this plane is = ("<<bestindex_x<<","<<bestindex_z<<")"<<endl;



        float * imageline=new float [(InputImage->nx)];
        float * lineranking=new float [(InputImage->nz)];
        for(int iz=0; iz<InputImage->nz; iz++)
        {
            for(int ix=0; ix<InputImage->nx; ix++)
            {
                imageline[ix]=InputImagePtr[ix+(planeindex+iz*InputImage->ny)*InputImage->nx];
            }
            int * order=quickSort_order(imageline,InputImage->nx);
            delete[]order;
            float cursum=0;
            for(int i=0; i<5; i++)
            {
                cursum+=imageline[(InputImage->nx)-1-i];
            }
            lineranking[iz]=cursum;

        }

        int * order2=quickSort_order(lineranking,(InputImage->nz));

        int bestindex_to_use[3]= {0};
        bestindex_to_use[0]=bestindex_z;

        int bestindexall[10]= {0};

        for(int i=0; i<10; i++)
        {
            bestindexall[i]=order2[(InputImage->nz)-1-i]-bestindex_z;
            //cout<<bestindexall[i]<<","<<order2[(InputImage->nz)-1-i]<<endl;
        }

        cout<<endl;
        for(int i=0; i<10; i++)
        {
            if(bestindexall[i]!=0 && bestindex_to_use[1]==0)
            {
                if(fabs(bestindexall[i])==1.0f)
                {
                    bestindex_to_use[1]=floor(order2[(InputImage->nz)-1-i]);
                }
                else
                {
                    if(fabs(bestindexall[i])==2.0f)
                    {
                        bestindex_to_use[1]=floor(order2[(InputImage->nz)-1-i]);
                        bestindex_to_use[2]=((float)(bestindex_to_use[0])+(float)(bestindex_to_use[1]))/2;
                    }
                }
            }
            else if(bestindexall[i]!=0 && bestindex_to_use[1]>0 && bestindex_to_use[2]==0)
            {
                if(fabs(bestindexall[i])==1.0f || fabs(order2[(InputImage->nz)-1-i]-bestindex_to_use[1])==1 )
                {
                    bestindex_to_use[2]=floor(order2[(InputImage->nz)-1-i]);
                }
            }
        }
        //for(int i=0; i<3; i++){
        //cout<<bestindex_to_use[i]<<endl;
        //}
        delete[] order2;
        for(int besti=0; besti<3; besti++)
        {
            int iz=bestindex_to_use[besti];
            for(int ix=0; ix<InputImage->nx; ix++)
            {
                imageline[ix]=InputImagePtr[ix+(planeindex+iz*InputImage->ny)*InputImage->nx];
            }
            int * order=quickSort_order(imageline,InputImage->nx);
            for(int i=0; i<5; i++)
            {
                OutputImagePtr[order[InputImage->nx-1-i]+(planeindex+iz*InputImage->ny)*InputImage->nx]=1;
            }
            delete [] order;
        }
        delete [] imageline;
    }






    nifti_image_write(OutputImage);
    nifti_image_free(OutputImage);

    delete [] imageslice;

    nifti_image_free(InputImage);
    return 0;
}


