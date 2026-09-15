#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
//#include "_seg_Topo.h"

using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s <input> <gm_thr> <output>.\n\n",exec);
    printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}
float CalcL1(long int index, long int *dim_array, float *L1, float *Gradient,float *VoxelSize, float scaling, int imagesize);
float CalcL0(long int index, long int *dim_array, float *L0, float *Gradient,float *VoxelSize, float scaling, int imagesize);
bool isInTheBorderOf(long int index, long int *dim_array, unsigned char *LabelImage, unsigned char lab);

int main(int argc, char **argv)
{

    if (argc != 4)
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

    float threshold=atof(argv[2]);
    char * filename_out=argv[3];


    // LAPLACE

    int index6n[6];
    int index6i=0;
    int kx,ky,kz;
    kx=ky=kz=0;
    for (int kx=-1; kx<=1; kx+=2) index6n[index6i++]=kx+InputImage->nx*(ky+InputImage->ny*kz);
    kx=ky=kz=0;
    for (int ky=-1; ky<=1; ky+=2) index6n[index6i++]=kx+InputImage->nx*(ky+InputImage->ny*kz);
    kx=ky=kz=0;
    for (int kz=-1; kz<=1; kz+=2) index6n[index6i++]=kx+InputImage->nx*(ky+InputImage->ny*kz);


    int index=0;
    int imgsize=InputImage->nz*InputImage->ny*InputImage->nx;
    float * Laplace = new float [imgsize];
    float * Laplace2 = new float [imgsize];
    unsigned char * LaplaceMask = new unsigned char [imgsize];
    for (int iz=1; iz<(InputImage->nz); iz++)
    {
        for (int iy=1; iy<(InputImage->ny); iy++)
        {
            for (int ix=1; ix<(InputImage->nx); ix++)
            {
                LaplaceMask[index]=0;
                if(InputImagePtr[index+imgsize]>(1-threshold))
                {
                    LaplaceMask[index]=1;
                    Laplace[index]=0.5;
                    Laplace2[index]=0.5;
                }
                else
                {
                    Laplace[index]=1;
                    Laplace2[index]=1;
                    if((InputImagePtr[index]+InputImagePtr[index+imgsize]+InputImagePtr[index+2*imgsize])==0) // then we are outside the mask
                    {
                        LaplaceMask[index]=0;
                        Laplace[index]=0;
                        Laplace2[index]=0;
                    }
                    else if(InputImagePtr[index]>InputImagePtr[index+2*imgsize]) // then we are in the WM
                    {
                        LaplaceMask[index]=2;
                        Laplace[index]=1;
                        Laplace2[index]=1;
                    }
                    else  // then we are in the CSF
                    {
                        LaplaceMask[index]=0;
                        Laplace[index]=0;
                        Laplace2[index]=0;
                    }
                }
                index++;
            }
        }
    }

    cout << "==>> Solving the Laplace Equation (MaxIter=1000)" << endl;
    flush(cout);
    int iternumb=0;
    float epsilon=10000.0f;
    float epsilonThres=1e-5f;
    float CurrentEnergy=0.0f;
    float LastEnergy=0.0f;

    while(epsilon>epsilonThres && iternumb<500)
    {
        if((float(iternumb)/10.0f)==round(float(iternumb)/10.0f))
        {
            cout <<"iter="<< iternumb/10.0f<<"  eps="<<epsilon<<endl;
        }
        for(int i=0; i<imgsize; i++)
        {
            if(LaplaceMask[i]==1)
            {
                int curcount=0;
                Laplace2[i]=0;
                for(int neighb=0; neighb<6; neighb++)
                {
                    curcount++;
                    Laplace2[i]+=Laplace[i+index6n[neighb]];
                }
                Laplace2[i]=Laplace2[i]/curcount;
            }
        }
        CurrentEnergy=0;
        for(int i=0; i<imgsize; i++)
        {
            if(LaplaceMask[i]==1)
            {
                Laplace[i]=Laplace2[i];
                CurrentEnergy+=Laplace2[i];
            }
        }
        iternumb++;
        if (iternumb>1)
        {
            epsilon=fabs((LastEnergy-CurrentEnergy)/LastEnergy);
        }
        LastEnergy=CurrentEnergy;
    }
    cout << "==>> Finished solving the Laplace Equation" << endl;

    delete [] Laplace2;



    float * Gradient=new float [imgsize*3];
    float GradSum;

    for(int index=0; index<imgsize; index++)
    {
        if(LaplaceMask[index]==1)
        {

            Gradient[index] = (Laplace[index+index6n[0]]-Laplace[index+index6n[1]])/2;
            Gradient[index+1*imgsize]=(Laplace[index+index6n[2]]-Laplace[index+index6n[3]])/2;
            Gradient[index+2*imgsize]=(Laplace[index+index6n[4]]-Laplace[index+index6n[5]])/2;
            GradSum=sqrt(Gradient[index]*Gradient[index]+Gradient[index+1*imgsize]*Gradient[index+1*imgsize]+Gradient[index+2*imgsize]*Gradient[index+2*imgsize]);
            if(GradSum!=0)
            {
                Gradient[index]= Gradient[index]/GradSum;
                Gradient[index+1*imgsize]= Gradient[index+1*imgsize]/GradSum;
                Gradient[index+2*imgsize]= Gradient[index+2*imgsize]/GradSum;
            }
            else
            {
                Gradient[index]=1;
                Gradient[index+1*imgsize]=1;
                Gradient[index+2*imgsize]=1;

            }
        }
    }
    delete [] Laplace;


    //***************** GET L0 ******************
    cout << "==>> Calc L0" << endl;
    flush(cout);
    typedef std::multimap<float, long int> MapLengthType;
    typedef std::pair<float, long  int> PairType;
    typedef std::map<long int, MapLengthType::iterator > MapIndexType;
    long int indexcalc;

    MapLengthType MapLength;
    MapIndexType MapIndex;
    float * VoxelSize= (float *) calloc(3, sizeof(float));
    VoxelSize[0]=InputImage->dx;
    VoxelSize[1]=InputImage->dy;
    VoxelSize[2]=InputImage->dz;

    /* CRETING THE TAGED IMAGE AND FILLING THE SORTED HEAP FOR THE ORDERED
             * TRAVERSAL METHOD IN THE PURE GM AREA
             *
             * If Solved  => SolvedVisited = 2
             * If Visited => SolvedVisited = 1
             */

    long int * dim_array = (long int *) calloc(3, sizeof(long int));
    dim_array[0]=InputImage->nx;
    dim_array[1]=InputImage->ny;
    dim_array[2]=InputImage->nz;

    unsigned char * SolvedVisited = (unsigned char *) calloc(imgsize, sizeof(unsigned char));
    float * L0 = (float *) calloc(imgsize, sizeof(float));

    /* SOLVING THE EUCLIDEAN PDE BY THE ORDERED TRAVERSAL METHOD FOR THE OUTPV AREA */
    iternumb=1;
    LastEnergy=2;
    CurrentEnergy=1;

    for(int i=0; i<imgsize; i++)
    {
        L0[i]=0;
    }

    while(fabs((CurrentEnergy-LastEnergy)/CurrentEnergy)>epsilonThres && iternumb<6)
    {
        //while(iternumb<4){
        LastEnergy=CurrentEnergy;
        for(int i=0; i<imgsize; i++)
        {
            if(LaplaceMask[i]!=1)
            {
                SolvedVisited[i]=(unsigned char)2;
            }
            else
            {
                if(isInTheBorderOf(i, dim_array, LaplaceMask, 0)==1)
                {
                    L0[i]=CalcL0(i, dim_array, L0, Gradient,VoxelSize,1,imgsize);
                    MapIndex[i] = MapLength.insert(PairType(Laplace[i], i));
                    SolvedVisited[i]=1;
                }
                else
                {
                    SolvedVisited[i]=0;
                }
            }
        }



        MapLengthType::iterator tempminimum=MapLength.begin();
        while(MapLength.empty()==false)
        {
            tempminimum=MapLength.begin();
            index=tempminimum->second;
            SolvedVisited[index]=2;
            MapLength.erase(MapLength.begin());
            MapIndex.erase(index);


            for(long int delta=-1; delta<2; delta+=2)
            {
                /*     X Direction     */
                indexcalc=(long int)(index+delta);
                if(SolvedVisited[indexcalc]!=2)
                {
                    if(SolvedVisited[indexcalc]==1)
                    {
                        MapLength.erase(MapIndex[indexcalc]);
                        MapIndex.erase(indexcalc);
                    }
                    L0[indexcalc]=CalcL0(indexcalc, dim_array, L0, Gradient,VoxelSize,1,imgsize);
                    MapIndex[indexcalc] = MapLength.insert(PairType(Laplace[indexcalc], indexcalc));
                    SolvedVisited[indexcalc]=1;
                }
                /*     Y Direction     */
                indexcalc=(long int)(index+delta*dim_array[0]);
                if( SolvedVisited[indexcalc]!=2)
                {
                    if(SolvedVisited[indexcalc]==1)
                    {
                        MapLength.erase(MapIndex[indexcalc]);
                        MapIndex.erase(indexcalc);
                    }
                    L0[indexcalc]=CalcL0(indexcalc, dim_array, L0, Gradient,VoxelSize,1,imgsize);
                    MapIndex[indexcalc] = MapLength.insert(PairType(Laplace[indexcalc], indexcalc));
                    SolvedVisited[indexcalc]=1;
                }
                /*     Z Direction     */
                indexcalc=(long int)(index+delta*dim_array[0]*dim_array[1]);
                if( SolvedVisited[indexcalc]!=2)
                {
                    if(SolvedVisited[indexcalc]==1)
                    {
                        MapLength.erase(MapIndex[indexcalc]);
                        MapIndex.erase(indexcalc);
                    }
                    L0[indexcalc]=CalcL0(indexcalc, dim_array, L0, Gradient,VoxelSize,1,imgsize);
                    MapIndex[indexcalc] = MapLength.insert(PairType(Laplace[indexcalc], indexcalc));
                    SolvedVisited[indexcalc]=1;

                }
            }
        }
        CurrentEnergy=0;
        for(int i=0; i<imgsize; i++)
        {
            CurrentEnergy+=L0[i]!=L0[i]?0:L0[i];
        }
        cout << "Iteration = "<<iternumb<<"  -> FreeEnergy = " << fabs((CurrentEnergy-LastEnergy)/(CurrentEnergy)) << endl;

        iternumb++;
    }

    //***************** GET L1 ******************
    cout << "==>> Calc L1" << endl;
    flush(cout);

    //unsigned char * SolvedVisited = (unsigned char *) calloc(imgsize, sizeof(unsigned char));
    float * L1 = (float *) calloc(imgsize, sizeof(float));

    /* SOLVING THE EUCLIDEAN PDE BY THE ORDERED TRAVERSAL METHOD FOR THE OUTPV AREA */
    iternumb=1;
    LastEnergy=2;
    CurrentEnergy=1;

    for(int i=0; i<imgsize; i++)
    {
        L1[i]=0;
    }

    while(fabs((CurrentEnergy-LastEnergy)/CurrentEnergy)>epsilonThres && iternumb<6)
    {
        //while(iternumb<4){
        LastEnergy=CurrentEnergy;
        for(int i=0; i<imgsize; i++)
        {
            if(LaplaceMask[i]!=1)
            {
                SolvedVisited[i]=(unsigned char)2;
            }
            else
            {
                if(isInTheBorderOf(i, dim_array, LaplaceMask, 2)==1)
                {
                    L1[i]=CalcL1(i, dim_array, L1, Gradient,VoxelSize,1,imgsize);
                    MapIndex[i] = MapLength.insert(PairType(Laplace[i], i));
                    SolvedVisited[i]=1;
                }
                else
                {
                    SolvedVisited[i]=0;
                }
            }
        }



        MapLengthType::iterator tempminimum=MapLength.begin();
        while(MapLength.empty()==false)
        {
            tempminimum=MapLength.begin();
            index=tempminimum->second;
            SolvedVisited[index]=2;
            MapLength.erase(MapLength.begin());
            MapIndex.erase(index);


            for(long int delta=-1; delta<2; delta+=2)
            {
                /*     X Direction     */
                indexcalc=(long int)(index+delta);
                if(SolvedVisited[indexcalc]!=2)
                {
                    if(SolvedVisited[indexcalc]==1)
                    {
                        MapLength.erase(MapIndex[indexcalc]);
                        MapIndex.erase(indexcalc);
                    }
                    L1[indexcalc]=CalcL1(indexcalc, dim_array, L1, Gradient,VoxelSize,1,imgsize);
                    MapIndex[indexcalc] = MapLength.insert(PairType(Laplace[indexcalc], indexcalc));
                    SolvedVisited[indexcalc]=1;
                }
                /*     Y Direction     */
                indexcalc=(long int)(index+delta*dim_array[0]);
                if( SolvedVisited[indexcalc]!=2)
                {
                    if(SolvedVisited[indexcalc]==1)
                    {
                        MapLength.erase(MapIndex[indexcalc]);
                        MapIndex.erase(indexcalc);
                    }
                    L1[indexcalc]=CalcL1(indexcalc, dim_array, L1, Gradient,VoxelSize,1,imgsize);
                    MapIndex[indexcalc] = MapLength.insert(PairType(Laplace[indexcalc], indexcalc));
                    SolvedVisited[indexcalc]=1;
                }
                /*     Z Direction     */
                indexcalc=(long int)(index+delta*dim_array[0]*dim_array[1]);
                if( SolvedVisited[indexcalc]!=2)
                {
                    if(SolvedVisited[indexcalc]==1)
                    {
                        MapLength.erase(MapIndex[indexcalc]);
                        MapIndex.erase(indexcalc);
                    }
                    L1[indexcalc]=CalcL1(indexcalc, dim_array, L1, Gradient,VoxelSize,1,imgsize);
                    MapIndex[indexcalc] = MapLength.insert(PairType(Laplace[indexcalc], indexcalc));
                    SolvedVisited[indexcalc]=1;

                }
            }
        }
        CurrentEnergy=0;
        for(int i=0; i<imgsize; i++)
        {
            CurrentEnergy+=L1[i]!=L1[i]?0:L1[i];
        }
        cout << "Iteration = "<<iternumb<<"  -> FreeEnergy = " << fabs((CurrentEnergy-LastEnergy)/(CurrentEnergy)) << endl;

        iternumb++;
    }



    // FINAL

    cout << "Finished ==>> Saving results" << endl;
    nifti_image * OutputImage = nifti_copy_nim_info(InputImage);
    OutputImage->datatype=NIFTI_TYPE_FLOAT32;
    nifti_set_filenames(OutputImage,filename_out,0,0);
    OutputImage->dim[4]=2;
    OutputImage->dim[5]=0;
    OutputImage->dim[6]=OutputImage->nv=0;
    OutputImage->dim[7]=OutputImage->nw=0;
    OutputImage->dim[0]=(int)(OutputImage->dim[1]>=1)+(int)(OutputImage->dim[2]>=1)+(int)(OutputImage->dim[3]>=1)+(int)(OutputImage->dim[4]>=1)+(int)(OutputImage->dim[5]>=1)+(int)(OutputImage->dim[6]>=1)+(int)(OutputImage->dim[7]>=1);
    nifti_update_dims_from_array(OutputImage);
    nifti_datatype_sizes(OutputImage->datatype,&OutputImage->nbyper,&OutputImage->swapsize);
    OutputImage->data = (void *) calloc(2*imgsize, sizeof(float));
    float * OutputImagePtr = static_cast<float *>(OutputImage->data);
    //for(int i=0; i<(imgsize); i++){OutputImagePtr[i]=(float)L0[i];}
    //for(int i=0; i<(imgsize); i++){OutputImagePtr[i+imgsize]=(float)L1[i];}
    //for(int i=0; i<(imgsize); i++){OutputImagePtr[i+2*imgsize]=(float)L0[i]+(float)L1[i];}
    for(int i=0; i<(imgsize); i++)
    {
        OutputImagePtr[i]=((((float)L0[i]/((float)L0[i]+(float)L1[i]))<0.5|| LaplaceMask[i]==0)*InputImagePtr[i+imgsize]);
    }
    for(int i=0; i<(imgsize); i++)
    {
        OutputImagePtr[i+imgsize]=((((float)L0[i]/((float)L0[i]+(float)L1[i]))>0.5 || LaplaceMask[i]==2)*InputImagePtr[i+imgsize]);
    }
    nifti_image_write(OutputImage);
    nifti_image_free(OutputImage);

    nifti_image_free(InputImage);


    return 0;
}


float CalcL0(long int index, long int *dim_array, float *L0, float *Gradient,float *VoxelSize, float scaling,  int imagesize)
{

    float a=VoxelSize[0];
    float b=VoxelSize[1];
    float c=VoxelSize[2];
    long int IndexPlusDeltaX;
    long int IndexPlusDeltaY;
    long int IndexPlusDeltaZ;

    if(Gradient[index]<=0)
    {
        IndexPlusDeltaX=index-1;
    }
    else
    {
        IndexPlusDeltaX=index+1;
    }
    if(Gradient[index+imagesize]<=0)
    {
        IndexPlusDeltaY=index-dim_array[0];
    }
    else
    {
        IndexPlusDeltaY=index+dim_array[0];
    }
    if(Gradient[index+2*imagesize]<=0)
    {
        IndexPlusDeltaZ=index-dim_array[0]*dim_array[1];
    }
    else
    {
        IndexPlusDeltaZ=index+dim_array[0]*dim_array[1];
    }
    float divis=(b*c*fabsf(Gradient[index])+a*c*fabsf(Gradient[index+1*imagesize])+a*b*fabsf(Gradient[index+2*imagesize]));
    return divis>0?((scaling*a*b*c+b*c*fabsf(Gradient[index])*L0[IndexPlusDeltaX]+a*c*fabsf(Gradient[index+1*imagesize])*L0[IndexPlusDeltaY]+a*b*fabsf(Gradient[index+2*imagesize])*L0[IndexPlusDeltaZ])/divis):0;
}

float CalcL1(long int index, long int *dim_array, float *L1, float *Gradient,float *VoxelSize, float scaling, int imagesize)
{

    long int IndexPlusDeltaX;
    long int IndexPlusDeltaY;
    long int IndexPlusDeltaZ;
    float a=VoxelSize[0];
    float b=VoxelSize[1];
    float c=VoxelSize[2];

    if(Gradient[index]>0)
    {
        IndexPlusDeltaX=index-1;
    }
    else
    {
        IndexPlusDeltaX=index+1;
    }
    if(Gradient[index+1*imagesize]>0)
    {
        IndexPlusDeltaY=index-dim_array[0];
    }
    else
    {
        IndexPlusDeltaY=index+dim_array[0];
    }
    if(Gradient[index+2*imagesize]>0)
    {
        IndexPlusDeltaZ=index-dim_array[0]*dim_array[1];
    }
    else
    {
        IndexPlusDeltaZ=index+dim_array[0]*dim_array[1];
    }
    float divis=(b*c*fabsf(Gradient[index])+a*c*fabsf(Gradient[index+1*imagesize])+a*b*fabsf(Gradient[index+2*imagesize]));
    return divis>0?((scaling*a*b*c+b*c*fabsf(Gradient[index])*L1[IndexPlusDeltaX]+a*c*fabsf(Gradient[index+1*imagesize])*L1[IndexPlusDeltaY]+a*b*fabsf(Gradient[index+2*imagesize])*L1[IndexPlusDeltaZ])/divis):0;
}

bool isInTheBorderOf(long int index, long int *dim_array, unsigned char *LabelImage, unsigned char lab)
{
    if ((index+dim_array[0]*dim_array[1])<(dim_array[2]*dim_array[0]*dim_array[1]) && (index-dim_array[0]*dim_array[1])>0  )
    {
        if(LabelImage[index+1]==lab ||
                LabelImage[index-1]==lab ||
                LabelImage[index+dim_array[0]]==lab ||
                LabelImage[index-dim_array[0]]==lab ||
                LabelImage[index+dim_array[0]*dim_array[1]]==lab ||
                LabelImage[index-dim_array[0]*dim_array[1]]==lab)
        {
            return 1;
        }
    }
    return 0;
}
