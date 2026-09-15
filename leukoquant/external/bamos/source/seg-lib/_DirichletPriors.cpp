
#include "_DirichletPriors.h"

using namespace std;
static int KernelSize=8;
inline float pow_int1(const float base,
                             int exp);

// Returns the value of the digamma function applied to x
float Digamma(float x){
    if(x<0){
        return 0;
    }
    else if(x>0 && x<1E-5){
        return -gamma-1.0/x;
    }
    else if(x> 40){
        return logf(x)-1.0/(2.0*x)-1.0/(12*pow_int1(x,2))+1.0/(120*pow_int1(x,4))-1.0/(252*pow_int1(x,6));
    }
    else{

        int t=(int)(40-x)+1;
        float z=x+t;
        float tmpDigamma=Digamma(z);
        float tmpCorr=0;
        for(int i=1;i<=t;i++){
            tmpCorr+=1.0/(x+i-1);
        }
        return tmpDigamma-tmpCorr;
    }
}

// Returns the value taken by the derivative of the digamma function in x
float DigammaDer(float x){
    if(x<0){
        return 0;
    }
    else if(x>0 && x<1E-5){
        return 1.0/(x*x);
    }
    else if(x> 40){
        return 1.0/(x)+1.0/(2.0*pow_int1(x,2))+1.0/(6*pow_int1(x,3))-1.0/(30*pow_int1(x,5))+1.0/(42*pow_int1(x,7));
    }
    else{

        int t=(int)(40-x)+1;
        float z=x+t;
        float tmpDigamma=DigammaDer(z);
        float tmpCorr=0;
        for(int i=1;i<=t;i++){
            tmpCorr+=1.0/((x+i-1)*(x+i-1));
        }
        return tmpDigamma+tmpCorr;
    }
}

// Returns the value of x such that psi(x) = y using the Newton method to invert the digamma function
float NewtonDigamma(float y){
    float InitX;
    switch (y>=-2.22){
    case 0 :
        InitX=-1.0/(y+gamma);
        break;
    default :
        InitX=expf(y)+0.5;
        break;
    }
    float OldX=InitX;
    float NewX=InitX;
    int iter=0;
    while(fabs((OldX-NewX))>1E-15 || iter<=1){
        OldX=NewX;
        NewX=OldX-(Digamma(OldX)-y)/DigammaDer(OldX);
        iter++;
    }
    cout<< "Number of iterations for the inversion is "<<iter<<endl;
    return NewX;
}

// Returns the array of
float * AlphaCalculation(float * p){

    // Initialisation and Normalisation of the frequencies used for the determination of the alpha parameters
    float Proportions[MaxSupport];
    for(int i=0;i<MaxSupport;i++){
        if(p[i]>1||p[i]<=0){
            Proportions[i]=DefaultProp;
        }
        else{
            Proportions[i]=p[i];
        }
    }
    float NormalisationFactor=0;
    for(int i=0;i<MaxSupport;i++){
        NormalisationFactor+=Proportions[i];
    }
    for(int i=0;i<MaxSupport;i++){
        Proportions[i]/=NormalisationFactor;
    }

    // Initialisation of the values for alpha
    float Alpha_tmp[MaxSupport];
    float Alpha_new[MaxSupport];
    for(int i=0;i<MaxSupport;i++){
        Alpha_tmp[i]=1.0/1000;//Proportions[i];
        Alpha_new[i]=Alpha_tmp[i];
    }
    float RelativeChange=0;
    int iter=0;
//    while(RelativeChange>0.000001 || iter<=1){
    for(int l=0;l<10;l++){
        float SumAlpha=0;
        for(int i=0;i<MaxSupport;i++){
            Alpha_tmp[i]=Alpha_new[i];
            SumAlpha+=Alpha_tmp[i];
        }
        float DigammaSum=Digamma(SumAlpha);
        float RelativeChange_tmp=0;
        for(int i=0;i<MaxSupport;i++){
            Alpha_new[i]=NewtonDigamma(DigammaSum+logf(Proportions[i]));
//            float test1=Digamma(Alpha_new[i]);
//            float test2=DigammaSum+logf(Proportions[i]);
            RelativeChange_tmp+=fabs((Alpha_new[i]-Alpha_tmp[i]));///Alpha_tmp[i]);
            cout<<((Alpha_new[i]-Alpha_tmp[i]))<<" ";
        }
        cout<<endl;
        RelativeChange=RelativeChange_tmp/(float)MaxSupport;
        iter++;
    }
    float * AlphaResult=new float[MaxSupport];
    for(int i=0;i<MaxSupport;i++){
        AlphaResult[i]=Alpha_new[i];
    }
    return AlphaResult;
}

float BandwidthCalculation(int * Counts,SEG_PARAMETERS * segment_param){
    float Mean=0;
    float PartVar=0;
    float NormFact=0;
    for(int i=0;i<MaxSupport;i++){
        Mean+=(float)i*Counts[i]/MaxSupport;
        PartVar+=(float)i/MaxSupport*(float)i/MaxSupport*Counts[i];
        NormFact+=Counts[i];
    }
    float StdCount= sqrt(PartVar/NormFact-(Mean/NormFact)*(Mean/NormFact));
    return segment_param->BWfactor*powf((4*pow_int1(StdCount,5)/(3*NormFact)),0.2);
}

float GaussianKernel(float x1, float x2, float BW){
    return 1.0/sqrt(2*M_PI)*expf(-0.5*((x1-x2)/BW)*((x1-x2)/BW));
}

float * DirectKernelSmoothing(int * Counts, SEG_PARAMETERS * segment_param){
    int TotCounts=0;
    for(int i=0;i<MaxSupport;i++){
        TotCounts+=Counts[i];
    }
    float BWKernel=BandwidthCalculation(Counts,segment_param);
    float * KSResult=new float[MaxSupport];
    float tmpGK=0;
    for(int i=0;i<MaxSupport;i++){
        float Numerator=0;
        float Denominator=0;
        for(int j=0;j<MaxSupport;j++){
            tmpGK=GaussianKernel((float)i/MaxSupport,(float)j/MaxSupport,BWKernel);
            Numerator+=(float)Counts[j]/TotCounts*tmpGK;
            Denominator+=tmpGK;
        }
        KSResult[i]=Numerator/Denominator;
    }
    return KSResult;
}

float * DerivativeLocalRegression(int i, int * Counts, int order, float BW, float * Coeffs){
    float * DerivativeValue= new float[order];
    for(int o=0;o<order;o++){
        DerivativeValue[o]=0;
    }
    int TotCounts=0;
    for(int j=0;j<MaxSupport;j++){
        TotCounts+=Counts[i];
    }
    for(int j=0;j<MaxSupport;j++){
        float GK=GaussianKernel((float)i/MaxSupport,(float)j/MaxSupport,BW);
        float Diff=(float)i/MaxSupport-(float)j/MaxSupport;
        float ExpValue=0;
        for(int o=0;o<order;o++){
            ExpValue+=Coeffs[0]*pow_int1(Diff,o);
        }
        ExpValue=expf(ExpValue);
        for(int o=0;o<order;o++){
            DerivativeValue[o]+=GK*(pow_int1(Diff,o)*ExpValue-Counts[j]*pow_int1(Diff,o));
        }

    }
    return DerivativeValue;
}

float GradientDescentPoissonRegression(int i, int * Counts, int order,SEG_PARAMETERS* segment_param){
    float BW=BandwidthCalculation(Counts,segment_param);
    float RelativeChange=0;
    int iter=0;
    float * NewX=new float[order];
    float * OldX=new float[order];
    for(int o=0;o<order;o++){
        NewX[o]=-1;
        OldX[o]=-1;
    }
    while((RelativeChange>0.0001 || iter<=1) && iter <5000){
        float * Grad=DerivativeLocalRegression(i,Counts,order,BW,NewX);
        float step=1;
        for(int o=0;o<order;o++){
            OldX[o]=NewX[o];
            step=fabs(Grad[o])/10.0;
            NewX[o]=OldX[o]-step*Grad[o];
        }
        if(Grad!=NULL){
            delete [] Grad;
            Grad=NULL;
        }
        RelativeChange=fabs(OldX[0]-NewX[0]);
        iter++;
    }
    cout<<iter<<" ";
    float ValueResult=expf(NewX[0]);
    if(ValueResult<=1E-6){
        ValueResult=1E-6;
    }
    delete [] NewX;
    delete [] OldX;
    return ValueResult;
}

float * PoissonRegressionKernelSmoothing(int * Counts, int order,SEG_PARAMETERS * segment_param){
    float * KSResult=new float[MaxSupport];
    float NormFactor=0;
    for(int i=0;i<MaxSupport;i++){
        KSResult[i]=GradientDescentPoissonRegression(i, Counts, order,segment_param);
        NormFactor+=KSResult[i];
    }
    for(int i=0;i<MaxSupport;i++){
        KSResult[i]/=NormFactor;
    }
    return KSResult;
}

float * MultivariateJointHistogram(vector<int *> Counts,vector<int*> CountDistribution){
    int numbClasses=Counts.size();
    int sumCheck=0;
    // First check that the total count is equal in all of the arrays
    for(int c=0;c<numbClasses;c++){
        int sumCounts=0;
        for(int l=0;l<MaxSupport;l++){
            sumCounts+=Counts[c][l];
        }
        if(sumCheck==0){
            sumCheck=sumCounts;
        }
        if(sumCheck !=sumCounts){
            cout << "Pb of counts between files"<<endl;
            return NULL;
        }
    }
    for(int m=0;m<numbClasses;m++){
//    for(int i=0;i<sumCheck;i++){
//        cout<< CountDistribution[m][i]<<" ";
//    }
//    cout<<endl;
    }
    // Initialise the joint histogram
    int SizeHistogram=(int)pow_int1(MaxSupport,numbClasses);
    float * CountHistogram=new float[SizeHistogram];
    for(int i=0;i<SizeHistogram;i++){
        CountHistogram[i]=0;
    }
    int * Shift=new int[numbClasses];
    for(int c=0;c<numbClasses;c++){
        Shift[c]=(int)pow_int1(MaxSupport,c);
    }
    for(int d=0;d<sumCheck;d++){
        int indHist=0;
        for(int c=0;c<numbClasses;c++){
            indHist+=CountDistribution[c][d]*Shift[c];
//            cout<< c<<" "<<Shift[c]<<"    ";
        }
        CountHistogram[indHist]+=1.0/(float)sumCheck;
    }
    int CountNonZero=0;
    for(int i=0;i<SizeHistogram;i++){
        if(CountHistogram[i]>0){
            CountNonZero++;
        }
    }
    return CountHistogram;
}

float * GaussianBlurring(float * CountHistogram, float gauss_std, vector<int> dim,bool NormalisedSum){
    int numbDim=dim.size();
    int kernelsizemin=(int)floorf(gauss_std*KernelSize);
    int kernelsizemax=(int)ceilf(gauss_std*KernelSize);
    int kernelsize=0;
    (kernelsizemin/2)*2==kernelsizemin?kernelsize=kernelsizemax:kernelsize=kernelsizemin;

    // Construction of the kernel
    float * Kernel=new float[kernelsize];
    int kernelradius=kernelsize/2;
    for(int i =0;i<kernelsize;i++){
        Kernel[i]=1.0/(sqrt(2.0*M_PI)*gauss_std)*expf(-0.5*pow((i-kernelradius)/gauss_std,2));
    }

    // Construction of the Shift array
    int * Shift=new int[numbDim];
    for(int c=0;c<numbDim;c++){
        if(c>0){
        Shift[c]=Shift[c-1]*dim[c-1];
        }
        else{
            Shift[c]=1;
        }
    }

    // Copying the initial array to blur
    int SizeHistogram=Shift[numbDim-1]*dim[numbDim-1];
    float * Gaussian_tmp=new float[SizeHistogram];
    float * Gaussian_tmp2=new float[SizeHistogram];
    int CountNonZero=0;
    int CountZeroB=0;
    for(int i=0;i<SizeHistogram;i++){
        Gaussian_tmp[i]=CountHistogram[i];
        if(Gaussian_tmp[i]>0){
            CountNonZero++;
        }
    }

    for(int c=0;c<numbDim;c++){ // To do the blurring in each direction
        for(int j=0;j<SizeHistogram/(dim[c]);j++){
            int i=0;
            //            int RemainDim[numbDim-1];
            //            int IndRemainDim[numbDim-1];
            //            int RemainShift[numbDim-1];
            //            int FinalShift[numbDim-1];
                        int * RemainDim=new int[numbDim-1];
                        int * IndRemainDim=new int[numbDim-1];
                        int * RemainShift=new int[numbDim-1];
                        int * FinalShift=new int[numbDim-1];
            for(int i=0;i<numbDim-1;i++){
                RemainShift[i]=1;
                FinalShift[i]=1;
            }
            for(int d=0;d<numbDim;d++){
                if(d!=c){
                    RemainDim[i]=dim[d];
                    IndRemainDim[i]=d;
                    i++;
                }
            }
            for(int i=0;i<numbDim-1;i++){
                if(i>0){
                    RemainShift[i]=RemainShift[i-1]*RemainDim[i-1];
                }
            }
            int t=j;
            int tmp=0;
            int * Index=new int[numbDim-1];
            for(int i=numbDim-2;i>=0;i--){
                tmp=t/RemainShift[i];
                t=t-tmp*RemainShift[i];
                Index[i]=tmp;
            }
            for(int i=0;i<numbDim-1;i++){
                if(c<IndRemainDim[i]){
                    FinalShift[i]=dim[c]*RemainShift[i];
                }
                else{
                    FinalShift[i]=RemainShift[i];
                }
            }
            int FinalIndex=0;
            for(int i=0;i<numbDim-1;i++){
                FinalIndex+=Index[i]*FinalShift[i];
            }
            for(int i=0;i<dim[c];i++){
                float tmpkernel=0;
                 int FinalIndexB=FinalIndex+i*Shift[c];
                for(int k=0;k<kernelsize;k++){
                    if(i-kernelsize/2+k>=0 && i-kernelsize/2+k<dim[c]){
                        tmpkernel+=Kernel[k]*Gaussian_tmp[FinalIndexB+(-kernelsize/2+k)*Shift[c]];
                    }
                }
                Gaussian_tmp2[FinalIndexB]=tmpkernel;
            }
            if(RemainShift!=NULL){
                delete [] RemainShift;
                RemainShift=NULL;
            }
            if(RemainDim!=NULL){
                delete [] RemainDim;
                RemainDim=NULL;
            }
            if(IndRemainDim!=NULL){
                delete [] IndRemainDim;
                IndRemainDim=NULL;
            }
            if(FinalShift!=NULL){
                delete [] FinalShift;
                FinalShift=NULL;
            }
            if(Index!=NULL){
                delete [] Index;
                Index=NULL;
            }
        }
        float sumNorm=0;

        for(int l=0;l<SizeHistogram;l++){
            if(Gaussian_tmp2[l]>0){
                CountZeroB++;
            }
            sumNorm+=Gaussian_tmp2[l];
            Gaussian_tmp[l]=Gaussian_tmp2[l];   
        }
        if (NormalisedSum) {
            for(int l=0;l<SizeHistogram;l++){
                Gaussian_tmp[l]/=sumNorm;
            }
        }

    }
    if (Shift!=NULL) {
        delete [] Shift;
        Shift=NULL;
    }
    if (Kernel!=NULL) {
        delete [] Kernel;
        Kernel=NULL;
    }
    delete [] Gaussian_tmp2;
    Gaussian_tmp2=NULL;
    return Gaussian_tmp;
}



