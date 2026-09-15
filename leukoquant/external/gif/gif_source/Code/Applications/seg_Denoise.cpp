#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_FMM.h"

#define kernel_size_predefined 5
#define kernel_size_predefined_square kernel_size_predefined*kernel_size_predefined
#define kernel_size_predefined_cube kernel_size_predefined_square*kernel_size_predefined
using namespace std;
#define SegPrecisionTYPE float

void dct_forward(float DCT_Patch [kernel_size_predefined*kernel_size_predefined*kernel_size_predefined],float DCT_1D_kernel [kernel_size_predefined*kernel_size_predefined], float DCT_Coefficients [kernel_size_predefined*kernel_size_predefined*kernel_size_predefined])
{

    float DCT_temp[kernel_size_predefined*kernel_size_predefined*kernel_size_predefined];
    float DCT_weight[kernel_size_predefined]= {sqrtf(1.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined)};

    // xdir
    int iz, iy, ix, coef_index;
    for(iz=0; iz<kernel_size_predefined; iz++)
    {
        for(iy=0; iy<kernel_size_predefined; iy++)
        {
            int index=(kernel_size_predefined*(iy+kernel_size_predefined*iz));
            for(coef_index=0; coef_index<kernel_size_predefined; coef_index++)
            {
                int cur_index=index;
                float curr_coef_value=0;
                for(ix=0; ix<kernel_size_predefined; ix++, cur_index++)
                {
                    curr_coef_value+=DCT_Patch[cur_index]*DCT_1D_kernel[ix+coef_index*kernel_size_predefined];
                }
                DCT_Coefficients[index+coef_index]=DCT_weight[coef_index]*curr_coef_value;
            }
        }
    }

    //ydir
    for(iz=0; iz<kernel_size_predefined; iz++)
    {
        for(ix=0; ix<kernel_size_predefined; ix++)
        {
            int index=((kernel_size_predefined*kernel_size_predefined*iz)+ix);
            for(coef_index=0; coef_index<kernel_size_predefined; coef_index++)
            {
                int cur_index=index;
                float curr_coef_value=0;
                for(iy=0; iy<kernel_size_predefined; iy++, cur_index+=(kernel_size_predefined))
                {
                    curr_coef_value+=DCT_Coefficients[cur_index]*DCT_1D_kernel[iy+coef_index*kernel_size_predefined];
                }
                DCT_temp[index+coef_index*kernel_size_predefined]=DCT_weight[coef_index]*curr_coef_value;
            }
        }
    }

    // zdir
    for(iy=0; iy<kernel_size_predefined; iy++)
    {
        for(ix=0; ix<kernel_size_predefined; ix++)
        {
            int index=(kernel_size_predefined*(iy)+ix);
            for(coef_index=0; coef_index<kernel_size_predefined; coef_index++)
            {
                int cur_index=index;
                float curr_coef_value=0;
                for(iz=0; iz<kernel_size_predefined; iz++, cur_index+=(kernel_size_predefined_square))
                {
                    curr_coef_value+=DCT_temp[cur_index]*DCT_1D_kernel[iz+coef_index*kernel_size_predefined];
                }
                DCT_Coefficients[index+coef_index*kernel_size_predefined_square]=DCT_weight[coef_index]*curr_coef_value;
            }
        }
    }
    return;
}




void dct_inverse(float DCT_Coefficients [kernel_size_predefined*kernel_size_predefined*kernel_size_predefined],float DCT_1D_kernel [kernel_size_predefined*kernel_size_predefined], float DCT_Patch [kernel_size_predefined*kernel_size_predefined*kernel_size_predefined])
{
    float DCT_weight[kernel_size_predefined]= {sqrtf(1.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined)};
    // xdir
    int iz, iy, ix;
    for(int i=0; i<pow(kernel_size_predefined,3); i++)
    {
        DCT_Patch[i]=0;
    }

    for(int kiz=0; kiz<kernel_size_predefined; kiz++)
    {
        for(int kiy=0; kiy<kernel_size_predefined; kiy++)
        {
            for(int kix=0; kix<kernel_size_predefined; kix++)
            {
                if(DCT_Coefficients[kix+(kiy+kiz*kernel_size_predefined)*kernel_size_predefined]!=0)
                {
                    float DCT_Coef_Precompute=DCT_Coefficients[kix+(kiy+kiz*kernel_size_predefined)*kernel_size_predefined]*DCT_weight[kix]*DCT_weight[kiy]*DCT_weight[kiz];
                    for(iz=0; iz<kernel_size_predefined; iz++)
                    {
                        for(iy=0; iy<kernel_size_predefined; iy++)
                        {
                            for(ix=0; ix<kernel_size_predefined; ix++)
                            {
                                int index=(kernel_size_predefined*(iy+kernel_size_predefined*iz))+ix;
                                DCT_Patch[index]+=DCT_Coef_Precompute*DCT_1D_kernel[ix+kix*kernel_size_predefined]*DCT_1D_kernel[iy+kiy*kernel_size_predefined]*DCT_1D_kernel[iz+kiz*kernel_size_predefined];
                            }
                        }
                    }
                }
            }
        }
    }


    return;
}

float dct_inverse_center(float DCT_Coefficients [kernel_size_predefined*kernel_size_predefined*kernel_size_predefined],float DCT_1D_kernel [kernel_size_predefined*kernel_size_predefined])
{
    float DCT_weight[kernel_size_predefined]= {sqrtf(1.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined),sqrtf(2.0f/kernel_size_predefined)};
    // xdir
    float curr_value=0;
    for(int kiz=0; kiz<kernel_size_predefined; kiz++)
    {
        for(int kiy=0; kiy<kernel_size_predefined; kiy++)
        {
            for(int kix=0; kix<kernel_size_predefined; kix++)
            {
                curr_value+=DCT_Coefficients[kix+(kiy+kiz*kernel_size_predefined)*kernel_size_predefined]*DCT_weight[kix]*DCT_weight[kiy]*DCT_weight[kiz]*DCT_1D_kernel[(kernel_size_predefined/2)+kix*kernel_size_predefined]*DCT_1D_kernel[(kernel_size_predefined/2)+kiy*kernel_size_predefined]*DCT_1D_kernel[(kernel_size_predefined/2)+kiz*kernel_size_predefined];
            }
        }
    }

    return curr_value;
}


void SmallUsage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf(" Denoise:\n  Usage ->\t%s -in <filename> -<Type of Denoising> [OPTIONS]\n\n",exec);
    printf("\tSee the help for more details (-h).\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void Usage(char *exec)
{

    printf("\n Label Fusion:\n  Usage ->\t%s -in <filename> -<Type of Label Fusion> [OPTIONS]\n\n",exec);
    printf("  * * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * *\n\n");
    printf("  -in <filename>\t\t| Filename of the 4D integer label image\n\n");
    printf("  \t\t- Type of Classifier Fusion (mutually exclusive) -\n\n");
    printf("  -NLMF <k> <n> <img> <tmpl> \t| STEPS algorithm\n");
    printf("  \t\t\t\t| Size of the kernel (k), number of local labels to use (n),\n\n");


    printf("  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    return;
}



int main(int argc, char **argv)
{
    if (argc < 2)
    {
        Usage(argv[0]);
        return 0;
    }

    time_t start,end;
    time(&start);
    char * filename_Image=NULL;
    char * filename_OUT=NULL;
    float DCT_threshold=1;

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

        else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc)
        {
            filename_Image = argv[++i];
        }
        else if(strcmp(argv[i], "-out") == 0)
        {
            filename_OUT = argv[++i];
        }
        else if(strcmp(argv[i], "-NLMF") == 0)
        {
            DCT_threshold = atof(argv[++i]);
        }
        else
        {
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }

    // READING T1
    nifti_image * Image=nifti_image_read(filename_Image,true);
    if(Image == NULL)
    {
        fprintf(stderr,"* Error when reading the image: %s\n",filename_Image);
        return 1;
    }
    seg_changeDatatype<float>(Image);


    SegPrecisionTYPE * Image_PTR= static_cast<float *>(Image->data);
    float image_mean=0;
    int image_number=0;

    for(unsigned int i=0; i<Image->nvox; i++,Image_PTR++)
    {
        image_mean+=*Image_PTR;
        image_number++;
    }
    image_mean=image_mean/image_number;

    float image_std=0;
    Image_PTR= static_cast<float *>(Image->data);
    for(unsigned int i=0; i<Image->nvox; i++,Image_PTR++)
    {
        image_std+=(*Image_PTR-image_mean)*(*Image_PTR-image_mean);
    }
    image_std=sqrt(image_std/image_number);
    Image_PTR= static_cast<float *>(Image->data);

    DCT_threshold= DCT_threshold*image_std;

    float DCT_1D_kernel[kernel_size_predefined*kernel_size_predefined];
    float DCT_kernel_normalization=3.14159265/kernel_size_predefined;
    for(int j=0; j<kernel_size_predefined; j++)
    {
        for(int i=0; i<kernel_size_predefined; i++)
        {
            DCT_1D_kernel[i+j*kernel_size_predefined]=cosf((DCT_kernel_normalization*((float)i+0.5))*j);
        }
    }



    int kernel_shift=(int)floor(kernel_size_predefined/2);
    int Direction_Shift[3]= {1,Image->nx,Image->nx*Image->ny};

    nifti_image * Result = nifti_copy_nim_info(Image);
    Result->datatype=DT_FLOAT32;
    nifti_set_filenames(Result,filename_OUT,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(SegPrecisionTYPE));
    float * Resultdata = static_cast<float *>(Result->data);
    float * Buffer=new float [Result->nvox];
    float * Coefficients_image=new float [Result->nvox];

    for(unsigned int i=0; i<Image->nvox; i++)
    {
        Resultdata[i]=0;
        Coefficients_image[i]=0;
        Buffer[i]=0;
    }



    cout<< "Find Oracle"<< endl;
    //#ifdef _OPENMP
    //#pragma omp parallel for shared(Buffer,Result,Image_PTR,DCT_1D_kernel)
    //#endif
    float DCT_Coefficients[kernel_size_predefined*kernel_size_predefined*kernel_size_predefined];
    float DCT_Patch[kernel_size_predefined*kernel_size_predefined*kernel_size_predefined];
    for (int iz=kernel_shift; iz<(Image->nz-kernel_shift); iz++)
    {
        cout<<iz<<endl;
        for (int iy=kernel_shift; iy<(Image->ny-kernel_shift); iy++)
        {
            for (int ix=kernel_shift; ix<(Image->nx-kernel_shift); ix++)
            {
                // Copy patch to DCT_Patch array
                int indexCentre=((Image->nx*(iy+Image->ny*iz))+ix);
                for (int diz=(-kernel_shift); diz<(kernel_shift); diz++)
                {
                    for (int diy=(-kernel_shift); diy<(kernel_shift); diy++)
                    {
                        for (int dix=(-kernel_shift); dix<(kernel_shift); dix++)
                        {
                            DCT_Patch[(dix+kernel_shift)+(diy+kernel_shift)*kernel_size_predefined+(diz+kernel_shift)*kernel_size_predefined_square]=
                                Image_PTR[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]];
                        }
                    }
                }

                dct_forward(DCT_Patch,DCT_1D_kernel,DCT_Coefficients);
                float number_of_coefficients=0;
                for(int i=0; i<(kernel_size_predefined_cube); i++)
                {
                    DCT_Coefficients[i]=(fabs(DCT_Coefficients[i])<DCT_threshold)?0:DCT_Coefficients[i];
                    if(fabs(DCT_Coefficients[i])>0)
                    {
                        number_of_coefficients+=1;
                    }
                }
                number_of_coefficients=(1)/(1+number_of_coefficients);

                dct_inverse(DCT_Coefficients,DCT_1D_kernel,DCT_Patch);

                for (int diz=(-kernel_shift); diz<(kernel_shift); diz++)
                {
                    for (int diy=(-kernel_shift); diy<(kernel_shift); diy++)
                    {
                        for (int dix=(-kernel_shift); dix<(kernel_shift); dix++)
                        {
                            Buffer[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]]+=
                                number_of_coefficients*DCT_Patch[(dix+kernel_shift)+(diy+kernel_shift)*kernel_size_predefined+(diz+kernel_shift)*kernel_size_predefined_square];
                            Coefficients_image[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]]+=number_of_coefficients;
                        }
                    }
                }
            }
        }
    }


    for(unsigned int i=0; i<Image->nvox; i++)
    {
        Buffer[i]=Buffer[i]/Coefficients_image[i];
        Coefficients_image[i]=0;
    }
    cout<< endl;

    //      #ifdef _OPENMP
    //      #pragma omp parallel for shared(Buffer, Result,Image_PTR,DCT_1D_kernel)
    //      #endif

    cout<< "Solve for oracle based coefficients"<< endl;
    float DCT_Coefficients_oracle[kernel_size_predefined*kernel_size_predefined*kernel_size_predefined];
    for (int iz=kernel_shift; iz<(Image->nz-kernel_shift); iz++)
    {
        cout<<iz<<endl;
        for (int iy=kernel_shift; iy<(Image->ny-kernel_shift); iy++)
        {
            for (int ix=kernel_shift; ix<(Image->nx-kernel_shift); ix++)
            {
                // Copy patch to DCT_Patch array
                int indexCentre=((Image->nx*(iy+Image->ny*iz))+ix);
                for (int diz=(-kernel_shift); diz<(kernel_shift); diz++)
                {
                    for (int diy=(-kernel_shift); diy<(kernel_shift); diy++)
                    {
                        for (int dix=(-kernel_shift); dix<(kernel_shift); dix++)
                        {
                            DCT_Patch[(dix+kernel_shift)+(diy+kernel_shift)*kernel_size_predefined+(diz+kernel_shift)*kernel_size_predefined_square]=Buffer[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]];
                        }
                    }
                }
                dct_forward(DCT_Patch,DCT_1D_kernel,DCT_Coefficients_oracle);

                for (int diz=(-kernel_shift); diz<(kernel_shift); diz++)
                {
                    for (int diy=(-kernel_shift); diy<(kernel_shift); diy++)
                    {
                        for (int dix=(-kernel_shift); dix<(kernel_shift); dix++)
                        {
                            DCT_Patch[(dix+kernel_shift)+(diy+kernel_shift)*kernel_size_predefined+(diz+kernel_shift)*kernel_size_predefined_square]=Image_PTR[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]];
                        }
                    }
                }
                dct_forward(DCT_Patch,DCT_1D_kernel,DCT_Coefficients);

                float number_of_coefficients=0;
                for(int i=0; i<(kernel_size_predefined_cube); i++)
                {
                    DCT_Coefficients[i]=(fabs(DCT_Coefficients_oracle[i])<(DCT_threshold))?0:DCT_Coefficients[i];
                    if(fabs(DCT_Coefficients[i])>0)
                    {
                        number_of_coefficients+=1;
                    }
                }
                number_of_coefficients=(1)/(1+number_of_coefficients);

                dct_inverse(DCT_Coefficients,DCT_1D_kernel,DCT_Patch);
                for (int diz=(-kernel_shift); diz<(kernel_shift); diz++)
                {
                    for (int diy=(-kernel_shift); diy<(kernel_shift); diy++)
                    {
                        for (int dix=(-kernel_shift); dix<(kernel_shift); dix++)
                        {
                            Resultdata[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]]+=number_of_coefficients*DCT_Patch[(dix+kernel_shift)+(diy+kernel_shift)*kernel_size_predefined+(diz+kernel_shift)*kernel_size_predefined_square];
                            Coefficients_image[indexCentre+dix*Direction_Shift[0]+diy*Direction_Shift[1]+diz*Direction_Shift[2]]+=number_of_coefficients;
                        }
                    }
                }
            }
        }
    }

    for(unsigned int i=0; i<Image->nvox; i++)
    {
        Resultdata[i]=Resultdata[i]/Coefficients_image[i];
    }

    cout<< endl;

    delete [] Buffer;
    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(Image);

    time(&end);

    int minutes = (int)floorf(float(end-start)/60.0f);
    int seconds = (int)(end-start - 60*minutes);
    cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;
    return 0;
}


