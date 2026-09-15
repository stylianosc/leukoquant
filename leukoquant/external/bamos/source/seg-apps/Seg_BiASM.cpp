
//  main.cpp
//  TreeSeg
//
//  Created by Carole Sudre on 08/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//  Modified by Stylianos Charalampous on 08/05/2026.
//

#include "_Seg_InitAndRun.h"


void Usage(char *exec)
{
    printf("\n  BaMoS Statistical Segmentation:\n  Usage ->\t%s -in <int number of images> <filename1> <filename2> [OPTIONS]\n\n",exec);
     printf("  * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
     printf("\t-in <number of images to segment> <filename1> <filename2>\t\t| Filename of the input image\n");
     printf("\t-out <number of outputs> <filename1> <filename2>\t\t| Filename1 with all leaves of the segmented image and Filename2 with only general classes\n");
     printf("\t\t\t\t| The input images should be 3D images.\n");
     //  printf("\t\t\t\t| Images should be coregistered\n\n");
     printf("  \t\t- Select one of the following (mutually exclusive) -\n\n");
     printf("\t-priors <n> <fnames>\t| The number of priors (n>0) and their filenames. Priors should be registered to the input images\n");
     printf("\t\t\t\t| Priors are mandatory for now. Won't be in the next release.\n");
     printf("\t-nopriors <n>\t\t| The number of classes (n>0) \n\n");
     printf("\t-txt_in <filename> \t| Name of the file containing all informations to recreate Tree. Usually used in combination with DP options to rework on current model obtained \n.");
     printf("\t-txt_out <filename> \t| Name of the file under which to store information on the result once everything performed \n");

     printf("  * * * * * * * * * * * * * * * * GENERAL OPTIONS * * * * * * * * * * * * * * * * *\n\n");

     printf("\t-mask <filename>\t| Filename of the brain-mask of the input image\n");
     printf("\t-quantMax <float between 0 and 1>\t| Percentage of voxels to consider as high outliers (default=0.99) \n");
     printf("\t-max_iter <int>\t\t| Maximum number of iterations (default = 100)\n");
     printf("\t-min_iter <int>\t\t| Minimum number of iterations (default = 6)\n");
     printf( "\t-norm_thresh <float>\t| Threshold used for stop for EM convergence when not applying BF correction\n");
     printf("\t-BiASM <bool> \t\t BiASM run ON [1] or OFF [0]. Default is ON\n");
     printf("\t-quantMin <float> : \t | Value between 0 and 1 indicating the lower threshold for the intensities quantile cut off (default is 0) \n");
     printf("\t-quantMax <float> : \t | Value between 0 and 1 indicating the upper threshold for the intensities quantile cut off. -quantMax must be higher than -quantMin (default is 1) \n");
     printf("\t-splitAccept <int> \t\t | Value deciding about the behavior to adopt when wanting to split a small class \n");
     printf("\t\t 0 : if the class weight is under weightMini/5 relative to main Node, the split operation cannot be performed and the EMf is not continued \n");
     printf("\t\t 1 : if the class weight is under weightMini/5 relative to main Node, the split operation cannot be performed but the EM is still conducted and the operation checked for acceptance \n");
     printf("\t\t other : the split operation is conducted whatever the weight for the class to split \n");
     printf("\t-in_DC <filename>\t| Filename of the data already corrected and log-transformed to be used for the analysis\n");
     printf("\t-priorDGM <filename>\t| Filename of the priors to belong to the DGM. Required if doing the correction for juxtacortical lesions\n");
     printf("\t-juxtaCorr <bool>\t| Flag for juxtacortical correction (default=0)\n");
     printf("\t-progMod <number of steps> <number of modalities to add per step>\t| Indicates whether the modalities should be considered all together or progressively. Recommanded: 0 for 2 modalities; 2 2 1 for 3 modalities\n");
     printf("\t-averageIO <weight> <Filename atlas inlier> <Filename atlas outlier>\t| Indicates filenames to use as I/O atlases. Weight give the importance attributed to inlier atlas\n");
     printf("\t-averageGC <weight> <Filename segmentation1> <Filename segmentation2>\t| Indicates filenames to use as I/O atlases. Weight give the importance attributed to inlier atlas\n");


     printf(" * * * * BIAS FIELD (BF) CORRECTION OPTIONS * * * * \n\n");

     printf("\t-max_iterBF <int>\t\t| Maximum number of iterations used for the BF correction (default = 100)\n");
     printf("\t-bc_order <int>\t\t| Polynomial order for the bias field [off = 0, max = 5] (default = 3) \n");
     printf("\t-bc_thresh <float>\t| Bias field correction will run only if the ratio of improvement is below bc_thresh (default=0 [OFF]) \n");
     printf("\t-bf_out <filename>\t| Output the bias field image\n");
     printf("\t-bc_out <filename>\t| Output the bias corrected image. If not provided automatically given and saved anyway under DataCorrected_${nameoutpu}\n\n");
     printf("\t-BFP <bool>\t| Boolean flag to indicate if bias field should be progressive\n\n");

     printf(" * * * * Options on the BiASM part * * * * \n\n");
     printf("\t-init_split <int> \t| Choice of type of initialisation wanted for the split operation. NOT IN USE CURRENTLY \n.");
     printf("\t-CEM <bool> 0 for OFF, 1 for ON \t| Flag on application of classification EM for the BiASM steps : \n EM embedded for each general class is only applied on the voxels classified in general class by current hard segmentation \n");
     printf("\t-AtlasWeight <number of inputs> <list float between 0 and 1> \t| Weight attributed to the current smoothed segmentation when readapting the priors. Choose 0 to avoid any change of the statistical atlases (default = 0). A different weighting can be chose at each level. If only one value is given, the same value is applied at all levels\n");
     printf("\t-AtlasSmoothing <number of inputs > <list of floats> \t | Standard deviation attributed to the Gaussian filter used to smooth current segmentation when modifying atlases (default = 0.3) \n");
     printf("\t-KernelSize <int> \t| Size of the kernel to use for the gaussian filter used for the modification of the atlases (default = 3) \n");
     printf("\t-PriorsKept <int> \t|Choice on the behavior to adopt concerning the statistical atlases\n");
     printf("\t\t 0 : The statistical atlases are not kept after the first of the EM algorithm \n");
     printf("\t\t 1 : Default choice. The statistical atlases are not modified throughout the algorithm \n");
     printf("\t\t 2 : The statistical atlases are replaced by a smoothed version of the segmentation result from EM first. To be used only with an outlier model in order to use a robust enough initial segmentation \n");
     printf("\t\t 3 : The statistical atlases are replaced by a smoothed version of the segmentation result each type a new model is accepted \n");
     printf("\t\t 4 : To be used only in conjonction with uniform type change case 5. The statistical atlases at the first level are replaced only after EM first by the smoothed version of the segmentation result. The atlases related to the outliers subclasses (Level 2) are replaced by their segmentation result before normalisation each time a new model is accepted \n");
     printf("\t-SMOrder <int> \t| Type of ordering of the split and merge operations list 0 for the diagonal order, 1 for the vertical one (default = 1)\n");
     printf("\t-CommonChanges <bool> 0 for OFF, 1 for ON \t| Flag on authorising operations on more than one general class at a time. Option only available if the vertical order has been chosen (default = 0) \n\n");
     printf("\t-acceptType <int> : \t| Under what conditions a change in the BIC induces a change in the model\n");
     printf("\t\t 0 : The change is authorised as soon as BICnew > BICold\n");
     printf("\t\t 1 : And above : default behavior: the change is accepted only if relative change is above acceptance threshold set up with -accept_thresh\n");
     printf("\t-accept_thresh <float> : \t| Value of the threshold above which the relative change in BIC will lead to the acceptation of the new model if the -acceptType 1 is chosen. Default value is 1E-5 \n");
     printf("\t-deleteUW <bool> : \t | Behavior adopted with subclasses with very low weight. 1 to put it ON, 0 for OFF. By default ON. Delete the subclasses whose weight relative to their general class is below minimum weight after accepting a model");
     printf("\t-miniW <float> : \t | Value of the minimal weight for a subclass relative to its general class (Level 1). The model modification is not accepted if one of the gaussian distribution stemming from the model change is underweight. Same value used when the option -deleteUW is ON. (default = 0.01)\n");
     printf("\t-BICFP <bool> : \t | Choice for the definition of the number of free parameters when calculating the BIC. 0 without mixing coefficients, 1 with (default=1)\n");

     printf("* * * * Options when using known knowledge or influencing current classes number with known distribution * * * *\n\n");
     printf("\t-DP <int> <filename1> ,<filename2>... int with number of filenames to count and then the names of the corresponding filenames \t| Filenames with the number of subclasses in the used distributions \n.");
     printf("\t-DistClassInd <bool> 0 if NO, 1 if YES \t| If subclass distribution count provided, should we considered the general classes independent or not. Non independence aspect only possible if only 1 DP file is provided. Otherwise independent by default \n");
     printf("If only 1 file is used, distributions have to be read as lines of the file (1 line per general class) and the file is to be used as a description of the population. When more than 1 file is available, means we are using Count Files and not population description\n");
     printf("Therefore with more than 1 file, only the independent DP option is available\n");
     printf("\t WARNING With -DP, the order of the filenames (or of the lines if only one file is used) must be consistent with the order of the general classes set in the priors or elsewhere \n");
     printf("\t-Count <int> <List of filenames> \t| int number of general classes (must be coherent with the rest of the options) and the subsequent list of filenames where to store the count of the different subclasses \n");
     printf("Used to record for each run of BiASM on a population the number of subclasses obtained for each general class. Each file correspond to a general class \n ");
     printf("\t-smoothing_order <int> \t| Type of smoothing imposed to the categorical distribution brought by the -DP option when considered independent\n");
     printf("\t\t 0 : Simple Gaussian filtering \n");
     printf("\t\t >1 : Smoothing applying Simonoff method with Poisson regularisation \n");
     printf("\t-BWfactor <float> \t| Factor to use to modify the bandwidth calculation when smoothing the categorical distribution indepedently \n");
     printf("\t-DPGauss_std <float> \t| Standard deviation to use for the Gaussian blurring when using the non independent option for the distribution knowledge \n");
     printf("\t-Countmod <bool> \n\n");

     printf("* * * * OPTIONS ON POSSIBLE OUTLIER MODEL * * * * \n\n");
     printf("\t-outliersM <int> \t|Choice among the different outliers model :\n");
     printf("\t\t 0 : By default, no outliers model considered");
     printf("\t\t 1 : Model where a classical mixture between outliers and inliers is considered at the first level and the priors for the general classes are set up at the second level");
     printf("\t\t 2 : Model where a class of outliers is considered at the same level as the other general classes (GM WM CSF) more easily used");
     printf("\t-outliersW <float> \t| Value of the lambda chosen to define the probability to belong to the outlier class");
     printf("\t-uniformTC <int> \t|Choice of the behavior to adopt when in presence of a uniform distribution to split\n");
     printf("\t\t 0 : The uniform distribution is transformed into a gaussian distribution\n");
     printf("\t\t 2 : The uniform distribution is transformed into 2 gaussian distributions  as in the classical gaussian case\n");
     printf("\t\t 3 : The uniform distribution is transformed into 1 gaussian and 1 uniform distribution and a new node at level 1 is created for each gaussian newly formed. A statistical atlas is given and comes from the smoothed segmentation of the uniform class\n");
     printf("\t\t 4 : The uniform distribution is transformed into 1 gaussian and 1 uniform at level 2 under the outlier node and a classical mixture is considered under the outlier node \n");
     printf("\t\t 5 : The uniform distribution is transformed into 1 gaussian and 1 uniform at level 2 but following atlases originating from the parameters initialisation of the newly formed gaussian.\n");
     printf("\t-unifSplitW <float> \t| Value of the weight given to the uniform distribution when split into 1 gaussian distribution and a new uniform distribution. Only needed in case of -uniformTC=2 or 4");
     printf("\t-varInitUnif <int> : \t|Choice of the type of variance initialisation to adopt for the gaussian distribution when splitting a uniform distribution for -uniformTC >= 3 \n");
     printf("\t\t 0 : The value of the initial mean corresponds to the first maximum of the blurred histogram corresponding to the uniform distribution. All values of the variance calculated for the class under the uniform distribution divided by 10 \n");
     printf("\t\t 1 : The parameters of the newly formed gaussian distribution are chosen among the results of a kmeans classification with k=2 for the uniform distribution. The choice of the parameters among the results is defined with -init_splitUnif \n");
     printf("\t\t >1 : By default, the initial mean is the maximum of the blurred histogram and the variance is set isotropic taking as value the minimum between the size of one bin of the histogram and the calculated variance divided by two\n");
     printf("-init_splitUnif <int> \t| Choice of the parameters to use to initialise the newly formed gaussian distribution when using the kmeans classification to define the parameters (choice 1 for -varInitUnif)\n");
     printf("\t\t 0 : The parameters with the smallest determinant for the variance are chosen for the new gaussian distribution\n");
     printf("\t\t 1 : Default. The parameters corresponding to the heavier subclass is chosen even if it corresponds to a large variance");
     printf("-miniW <float> \t| Minimal weight for which a subclass is allowed to be created. Has to be between 0.001 and 0.1 (default=0.01)\n");
     printf("-CovPriors <int> \t| Choice of behavior to adopt for the constraint overt the covariance matrix\n");
//     TO DO LIST OF COV PRIORS OPTIONS
     printf("-TypicalityAtlas <bool> \t| Indicates if a typicality atlas should be created to enhance the sensitivity to outliers\n");
     printf("-VLkappa <float> \t| Mahalanobis threshold used to build the typicality atlas (default=3)\n");


     printf("* * * * MRF OPTIONS * * * * \n\n");
     printf("\t-MRF <int> \t| Choice among 4 different options about types of MRF we might use  \n");
     printf("\t\t 0 : By default no MRF applied. Way of putting MRF OFF \n");
     printf("\t\t 1 : MRF applied during the obtention of the model\n");
     printf("\t\t 2 : MRF applied a posteriori on a last EM run after convergence of the model\n");
     printf("\t\t 3 : MRF applied both during the obtention of the model and once afterwards (With G Matrix or not, possibly different)\n");
     printf("\t-GMRF <filename> \t| Setting of the G matrix to use with the MRF application. Informations for the neighborhood relationships and the value in the G matrix when MRF applied in the obtention of the model\n");
     printf("\t-GMRFPost <filename> \t| Gives the filename to use to obtain information for construction of GMatrix used a posteriori for an MRF \n");


     printf("\t-scratchDir <PathString> \t| Path to scratch directory for debug output files (default: /scratch0)\n");
     printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
     return;
 }

// argc is the number of arguments and argv the corresponding array of parameters
int main(int argc, char **argv)
{

    if (argc < 1)
    {
        Usage(argv[0]);
        return 0;
    }

    // Set all default parameters to be used in the segmentation process
    SEG_PARAMETERS * segment_param = new SEG_PARAMETERS [1]();
    segment_param->maxIteration=100;
    segment_param->maxIterationBF=100;
    segment_param->minIteration=6;
    segment_param->BWfactor=1;
    segment_param->choiceInitSplit=1;
    segment_param->maxNumbLeaves=40;
    segment_param->DPGauss_std=1;
    segment_param->ConvThreshold=1E-5;
    segment_param->Bias_threshold=1E-5;
    segment_param->flag_progressiveBFC=0;
    segment_param->AtlasWeight.push_back(0);
    segment_param->AtlasSmoothing.push_back(1); // Allow for same default for nipype Change on 20/11/2017
    segment_param->KernelSize=3;
    segment_param->quantMin=0.01;
    segment_param->quantMax=0.99;
    segment_param->OutliersWeight=0;
    segment_param->UnifSplitWeight=0.5;
    segment_param->WeightMini=0.01;
    segment_param->VarianceInitUnif=1;
    segment_param->OutliersMod=0;
    segment_param->DistInitUnif=4;
    segment_param->flag_Outliers=0;
    segment_param->flag_input=0;
    segment_param->numbmodal=0;
    segment_param->flag_out=0;
    segment_param->flag_DP=0;
    segment_param->flag_DistClassInd=1;
    segment_param->flag_Count=0;
    segment_param->flag_mask=0;
    segment_param->flag_Bias=1;
    segment_param->flag_BiASM=1;
    segment_param->flag_OutBrain=1;
    segment_param->flag_OutBrain=0;
    segment_param->flag_CommonChanges=0;
    segment_param->flag_EMfirst_out=0;
    segment_param->flag_Countmod=0;
    segment_param->flag_bf_out=0;
    segment_param->flag_bc_out=0;
    segment_param->flag_data_out=0;
    segment_param->flag_manual_priors=0;
    segment_param->numb_classes=0;
    segment_param->flag_intxt=0;
    segment_param->flag_CEM=0;
    segment_param->flag_MRF=0;
    segment_param->flag_MRFIn=0;
    segment_param->flag_MRFPost=0;
    segment_param->flag_MRFOut=0;
    segment_param->flag_GMatrix=0;
    segment_param->flag_GMatrixPost=0;
    segment_param->flag_DeleteUnderWeight=1;
    segment_param->flag_bnbComb=0;
    segment_param->flag_unifTot=0;
    segment_param->flag_RefinedSOrdering=1;
    segment_param->smoothing_order=2;
    segment_param->bias_order=3;
    segment_param->PriorsKept=1;
    segment_param->SMOrder=0; // Change on 20/11/2017 To allow for same default in nipype
    segment_param->flag_PriorsAdaptedOut=0;
//    segment_param->uniformTypeChange=5;
    segment_param->uniformTypeChange=4;
    segment_param->choiceInitSplitUnifKMeans=0; // Change on 20/11/2017 Same default as in nipype
    segment_param->flag_KMeansModif=1;
    segment_param->AcceptanceThreshold=1E-5;
    segment_param->AcceptanceType=1;
    segment_param->SplitAccept=0; // Change default for same as nipype // Change on 20/11/2017
    segment_param->BICFP=1;
    segment_param->flag_JuxtaCorrection=1;
    segment_param->CovPriorsMerge=0;
    segment_param->CovPriorsSplit=0;
    segment_param->flag_CovPriors=0;
    segment_param->flag_optMRFOut=1;
    segment_param->flag_OutlierAtlas=-1;
    segment_param->class_keptpriorsmax=-1;
    segment_param->Mahal=3;
    segment_param->ProgressivePriors=3;
    segment_param->MaxRunEM=250; // Change default for same as nipype
    segment_param->flag_savePriors=0;
    segment_param->IndexGM=0;
    segment_param->IndexWM=1;
    segment_param->IndexCSF=2;
    segment_param->IndexOut=3;
    segment_param->MeanPriors=10;
    segment_param->VLkappa=3;
    segment_param->flag_TypicalityAtlas=0;
    segment_param->flag_BoostAtlas=0;

    /* read the input parameters */
    cout << "Number of arguments is "<<argc<<endl;
    cout<< "Arguments are "<< argv<<endl;
    for(int i=1;i<argc;i++){

        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc){

            segment_param->numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if(segment_param->numbmodal<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+segment_param->numbmodal)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<segment_param->numbmodal; m++){
                    segment_param->filename_input.push_back(argv[++i]);
                    segment_param->NumberAddedModalities.push_back(1); // By default, consider that one modality is added at a time
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_input=1;

        }
        else if(strcmp(argv[i], "-DP") == 0 && (i+1)<argc){

            segment_param->numbDP=atoi(argv[++i]); // atoi convert string to integer

            if((i+segment_param->numbDP)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<segment_param->numbDP; m++){
                    segment_param->filename_DP.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_DP=1;

        }

        else if(strcmp(argv[i], "-Count") == 0 && (i+1)<argc){

            segment_param->numbCount=atoi(argv[++i]); // atoi convert string to integer

            if((i+segment_param->numbCount)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<segment_param->numbCount; m++){
                    segment_param->filename_Count.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_Count=1;

        }

        else if(strcmp(argv[i], "-out") == 0 && (i+1)<argc){
            //          segment_param->filename_out = argv[++i];
            //          segment_param->flag_out=1;
            int numberFilesOut=atoi(argv[++i]);

            if((i+numberFilesOut)<argc && numberFilesOut<=2){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numberFilesOut; m++){
                    segment_param->filename_out.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -output are incomplete or not correct\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_out=1;

        }
        else if(strcmp(argv[i], "-mask") == 0 && (i+1)<argc){
            segment_param->filename_mask=argv[++i];
            segment_param->flag_mask=1;
        }
        else if(strcmp(argv[i], "-inDC_flag") == 0 && (i+1)<argc){
            segment_param->flag_inDC=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-inDC") == 0 && (i+1)<argc){
            segment_param->filename_inDC=argv[++i];
            segment_param->flag_inDC=1;
            segment_param->flag_inDCFile=1;
        }
        else if(strcmp(argv[i], "-priorDGM") == 0 && (i+1)<argc){
            segment_param->filename_DGMPrior=argv[++i];
            segment_param->flag_DGMPrior=1;
        }
        else if(strcmp(argv[i], "-priorKC") == 0 && (i+1)<argc){
            segment_param->class_keptpriorsmax=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-juxtaCorr") == 0 && (i+1)<argc){
            segment_param->flag_JuxtaCorrection=atoi(argv[++i]);
//            segment_param->flag_DGMPrior=1;
        }
        else if(strcmp(argv[i], "-meanPriors") == 0 && (i+1)<argc){
            segment_param->MeanPriors=atof(argv[++i]);
            segment_param->flag_meanPriors=1;
            //            segment_param->flag_DGMPrior=1;
        }
        else if(strcmp(argv[i], "-NormMask") == 0 && (i+1)<argc){
            segment_param->flag_NormMask=(bool) atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-priors_out") == 0 && (i+1)<argc){ // In case Priors are used for an outlier model
            segment_param->numb_classes_out=atoi(argv[++i]);
            if((i+segment_param->numb_classes)<argc){// Read filenames of the priors. There must be as many priors filenames
                // as classes defined in segment_param->numb_classes
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int classnum=0; classnum<segment_param->numb_classes_out; classnum++){
                    segment_param->filename_priors_out.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -priors are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_manual_priors_out=1;
        }
        else if(strcmp(argv[i], "-progMod") == 0 && (i+1)<argc){ // In case we add progressively or not the modalities, determines the scheme that is wanted
            int numbProgmodaSteps=atoi(argv[++i]);
            segment_param->NumberAddedModalities.clear();
//            segment_param->NumberAddedModalities.shrink_to_fit();
            cout<<"new size of NumberAddedModalities is"<<segment_param->NumberAddedModalities.size();
            if((i+numbProgmodaSteps)<argc){// Reads and pushback in numberAddedModalities
                for(int pms=0; pms<numbProgmodaSteps; pms++){
                    segment_param->NumberAddedModalities.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -progMod are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
        }
        
        else if(strcmp(argv[i], "-priors") == 0 && (i+1)<argc){ // In case some priors are used
            segment_param->numb_classes=atoi(argv[++i]); // atoi convert string to integer
            if(segment_param->numb_classes<2){
                cout<<"Number of classes has to be bigger than 1";
                return 0;
            }
            if((i+segment_param->numb_classes)<argc){// Read filenames of the priors. There must be as many priors filenames
                // as classes defined in segment_param->numb_classes
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int classnum=0; classnum<segment_param->numb_classes; classnum++){
                    segment_param->filename_priors.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -priors are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_manual_priors=1;
        }
        
        else if(strcmp(argv[i], "-averageIO") == 0 && (i+3)<argc){ // In case we use two other images to make input of IO priors
            segment_param->weightIOAverage=atof(argv[++i]);
            segment_param->filename_IOAverage.push_back(argv[++i]);
            segment_param->filename_IOAverage.push_back(argv[++i]);
            segment_param->flag_IOAverage=1;
        }
        else if(strcmp(argv[i], "-averageGC") == 0 && (i+3)<argc){ // In case we use two other segmentations to do average for priors of general classes
            segment_param->weightGCAverage=atof(argv[++i]);
            segment_param->filename_GCAverage.push_back(argv[++i]);
            segment_param->filename_GCAverage.push_back(argv[++i]);
            segment_param->flag_GCAverage=1;
        }
        else if(strcmp(argv[i], "-AdaptTransform") == 0 && (i+1)<argc){ // In case some priors are used
            int numbAdaptFiles=atoi(argv[++i]);

            if((i+numbAdaptFiles)<argc){// Read filenames of the priors. There must be as many priors filenames
                // as classes defined in segment_param->numb_classes
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int nf=0; nf<numbAdaptFiles; nf++){
                    segment_param->filename_AdaptTransform.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -AdaptTransform are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_AdaptTransform=1;
            
        }
        else if(strcmp(argv[i], "-bnb") == 0 && (i+1)<argc){ // // List of the anatomical priors to be considered under brain. Needed with OM8 when model is 3 layered from the beginning
            int numbBrainPriors=atoi(argv[++i]);
            if((i+numbBrainPriors)<argc){// Read filenames of the priors. There must be as many priors filenames
                // as classes defined in segment_param->numb_classes
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int classnum=0; classnum<numbBrainPriors; classnum++){
                    segment_param->vecBP.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter brain priors combination are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_bnbComb=1;
        }
        else if(strcmp(argv[i], "-priors_bnb") == 0 && (i+1)<argc){
            if((i+2)<argc){// Read filenames of the priors. There must be as many priors filenames
                // as classes defined in segment_param->numb_classes
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int classnum=0; classnum<2; classnum++){
                    segment_param->filename_priors_bnb.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -priors_bnb are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_manual_priors_bnb=1;
        }
        

        else if(strcmp(argv[i], "-PriorsAdaptedOut") == 0 && (i+1)<argc){ // In case some priors are used
            segment_param->numb_classes=atoi(argv[++i]); // atoi convert string to integer
            if(segment_param->numb_classes<2){
                cout<<"Number of classes has to be bigger than 1";

            }
            else if((i+segment_param->numb_classes)<argc){// Read filenames of the priors. There must be as many priors filenames
                // as classes defined in segment_param->numb_classes
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int classnum=0; classnum<segment_param->numb_classes; classnum++){
                    segment_param->filename_PriorsAdaptedOut.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -priors are incomplete\n");
                Usage(argv[0]);
                return 1;
            }
            segment_param->flag_PriorsAdaptedOut=1;
        }

        else if(strcmp(argv[i], "-nopriors") == 0 && (i+1)<argc){
            segment_param->flag_manual_priors=0;
            segment_param->numb_classes=atoi(argv[++i]);// Read the number of classes fixed for the segmentation
        }

        else if(strcmp(argv[i], "-bc_order") == 0 && (i+1)<argc){
            segment_param->bias_order=(int)(atof(argv[++i])); // ????? why not atoi
            if(segment_param->bias_order==0){
                segment_param->flag_Bias=0;
            }
            else{
                segment_param->flag_Bias=1;
            }
        }
        else if(strcmp(argv[i], "-BFP") == 0 && (i+1)<argc){ // Indicates if we want progressive bias field correction or not
            int BFP=atoi(argv[++i]);
            segment_param->flag_progressiveBFC=(bool)BFP<=0?0:1; // ????? why not atoi
        }
        else if(strcmp(argv[i], "-outliersW") == 0 && (i+1)<argc){
            segment_param->OutliersWeight=atof(argv[++i]); // the outliers parameters set to 1 the flag and gives the corresponding uniform weight initially attributed to the outlier class
            segment_param->OutliersWeight=segment_param->OutliersWeight>0.5?0.5:segment_param->OutliersWeight;
            segment_param->OutliersWeight=segment_param->OutliersWeight<=0?0:segment_param->OutliersWeight;
            if(segment_param->OutliersWeight>0){
                segment_param->flag_Outliers=1;
            }
            if (segment_param->OutliersWeight>=0.04) {
                segment_param->VLkappa=2;
            }
        }
        else if(strcmp(argv[i], "-outliersM") == 0 && (i+1)<argc){
            segment_param->OutliersMod=atoi(argv[++i]); // the outliers mode can be 0, 1 or 2, 0 without outliers, 1 with the 2 level model of outliers, and 2 with constant priors over outliers at same level
            segment_param->OutliersMod=segment_param->OutliersMod>8?3:segment_param->OutliersMod;
            segment_param->OutliersMod=segment_param->OutliersMod<=0?0:segment_param->OutliersMod;
            if (segment_param->OutliersMod>0) {
                segment_param->quantMin=0;
                segment_param->quantMax=1;
            }
        }
        else if (strcmp(argv[i], "-unifSplitW") == 0 && (i+1)<argc){ // gives the value of the weight given to the uniform distribution when splitting it into one gaussian and one uniform
            segment_param->UnifSplitWeight=atof(argv[++i]);
            segment_param->UnifSplitWeight=segment_param->UnifSplitWeight>0.999?0.999:segment_param->UnifSplitWeight;
            segment_param->UnifSplitWeight=segment_param->UnifSplitWeight<=0.001?0.001:segment_param->UnifSplitWeight;
        }
        else if (strcmp(argv[i], "-miniW") == 0 && (i+1)<argc){ // gives the value of the minimal weight for which we allow a subclass to be created
            segment_param->WeightMini=atof(argv[++i]);
            segment_param->WeightMini=segment_param->WeightMini>0.1?0.1:segment_param->WeightMini;
            segment_param->WeightMini=segment_param->WeightMini<=0.001?0.001:segment_param->WeightMini;
        }
        else if (strcmp(argv[i], "-deleteUW") == 0 && (i+1)<argc){ // gives the choice for the deletion of underweight subclasses.
            int DeleteUW=atoi(argv[++i]);
            segment_param->flag_DeleteUnderWeight=DeleteUW>0?1:0;
        }
        else if (strcmp(argv[i], "-CovPriors") == 0 && (i+1)<argc){ // gives the choice for the deletion of underweight subclasses.
            int CovPriors=atoi(argv[++i]);
            if (CovPriors>0) {
                segment_param->flag_CovPriors=1;
                segment_param->CovPriorsType=CovPriors;
                segment_param->CovPriorsType=segment_param->CovPriorsType>10?5:segment_param->CovPriorsType;
            }
            else{
                segment_param->flag_CovPriors=0;
                segment_param->CovPriorsType=0;
                segment_param->CovTest=0;
            }
        }
        else if (strcmp(argv[i], "-covTest") == 0 && (i+1)<argc){ // gives the choice for the deletion of underweight subclasses.
            int covTest=atoi(argv[++i]);
            segment_param->CovTest=covTest<=0?0:covTest;
            segment_param->CovTest=covTest>3?3:covTest;
        }
        else if (strcmp(argv[i], "-covPriorsSplit") == 0 && (i+1)<argc){ // gives the choice for the deletion of underweight subclasses.
            int covPriorsSplit=atoi(argv[++i]);
            segment_param->CovPriorsSplit=covPriorsSplit<=0?0:covPriorsSplit;
            segment_param->CovPriorsSplit=covPriorsSplit>3?3:covPriorsSplit;
        }
        else if (strcmp(argv[i], "-outliersCombined") == 0 && (i+1)<argc){ // gives the choice for the deletion of underweight subclasses.
            int numbOutliersPriors=atoi(argv[++i]);
            for (int o=0; o<numbOutliersPriors; o++) {
                vector<int> OutlierPriorsCombination;
                while (atoi(argv[++i])!=-1) {
                    OutlierPriorsCombination.push_back(atoi(argv[i]));
                }
                segment_param->OutliersCombined.push_back(OutlierPriorsCombination);
            }
            cout<<"Combination of outlier priors"<<endl;
        }
        else if (strcmp(argv[i], "-covPriorsMerge") == 0 && (i+1)<argc){ // gives the choice for the deletion of underweight subclasses.
            int covPriorsMerge=atoi(argv[++i]);
            segment_param->CovPriorsMerge=covPriorsMerge<=0?0:covPriorsMerge;
            segment_param->CovPriorsMerge=covPriorsMerge>=1?1:covPriorsMerge;
        }
        else if (strcmp(argv[i], "-varInitUnif") == 0 && (i+1)<argc){ // gives the choice for the variance initialisation of the gaussian distribution when splitting the uniform distribution
            segment_param->VarianceInitUnif=atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-unifTot") == 0 && (i+1)<argc){ // when using the variance kMeans for the initialisation of the gaussian parameters for the uniform split, allows to consider at end of list of split operation the side of uniform distribution with higher variance after kmeans separation
            int UnifTot=atoi(argv[++i]);
            segment_param->flag_unifTot=(bool)UnifTot>0?1:0;
        }
        else if (strcmp(argv[i], "-distInitUnif") == 0 && (i+1)<argc){ // gives the distance at which to consider mean when initialising mean for kmeans in gaussian issued from uniform split initialisation
            segment_param->DistInitUnif=atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-splitAccept") == 0 && (i+1)<argc){ // gives the choice for the condition on which a split is accepted
            segment_param->SplitAccept=atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-TypicalityAtlas") == 0 && (i+1)<argc){ // decides if typicality maps will be used as atlases
            segment_param->flag_TypicalityAtlas=atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-BoostAtlas") == 0 && (i+1)<argc){ // decides if typicality maps will be used as atlases
            segment_param->flag_BoostAtlas=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-bc_thresh") == 0 && (i+1)<argc){
            segment_param->Bias_threshold=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-OutlierAtlas") == 0 && (i+1)<argc){
            segment_param->flag_OutlierAtlas=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-ProgressivePriors") == 0 && (i+1)<argc){
            segment_param->ProgressivePriors=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-MahalOutlier") == 0 && (i+1)<argc){
            segment_param->Mahal=atof(argv[++i]);
            segment_param->flag_OutlierAtlas=1;
        }
        else if(strcmp(argv[i], "-VLkappa") == 0 && (i+1)<argc){
            segment_param->VLkappa=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-norm_thresh") == 0 && (i+1)<argc){
            segment_param->ConvThreshold=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-accept_thresh") == 0 && (i+1)<argc){
            segment_param->AcceptanceThreshold=atof(argv[++i]);
        }
        else if(strcmp(argv[i],"-smoothing_order")==0 && (i+1)<argc){
            segment_param->smoothing_order=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-init_split")==0 && (i+1)<argc){
            segment_param->choiceInitSplit=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-init_splitUnif")==0 && (i+1)<argc){
            segment_param->choiceInitSplitUnifKMeans=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-BiASM")==0 &&(i+1)<argc){
            segment_param->flag_BiASM=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-Saving")==0 &&(i+1)<argc){
            segment_param->flag_savePriors=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-OutBrain")==0 &&(i+1)<argc){
            segment_param->flag_OutBrain=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-optMRFOut")==0 &&(i+1)<argc){
            segment_param->flag_optMRFOut=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-CEM")==0 &&(i+1)<argc){
            segment_param->flag_CEM=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-PriorsKept")==0 &&(i+1)<argc){
            segment_param->PriorsKept=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-DistClassInd")==0 &&(i+1)<argc){
            segment_param->flag_DistClassInd=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-SMOrder")==0 &&(i+1)<argc){
            segment_param->SMOrder=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-CommonChanges")==0 &&(i+1)<argc){
            segment_param->flag_CommonChanges=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-Countmod")==0 &&(i+1)<argc){
            segment_param->flag_Countmod=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-BWfactor")==0 &&(i+1)<argc){
            segment_param->BWfactor=atof(argv[++i]);
            cout<<"BW factor is "<<segment_param->BWfactor<<endl;
        }
        else if(strcmp(argv[i],"-AtlasWeight")==0 &&(i+1)<argc){
            int numbWeights=atoi(argv[++i]);
            if (i+numbWeights>argc) {
                fprintf(stderr,"Err:\t incompatible number of weights.\n");
                return 1;
            }
            if (numbWeights>0) {
                segment_param->AtlasWeight.clear();
            }
            cout << "AtlasWeight are ";
            for (int w=0; w<numbWeights; w++) {
                float WeightToPush=atof(argv[++i]);
                WeightToPush=WeightToPush>1?1:WeightToPush;
                WeightToPush=WeightToPush<0?0:WeightToPush;
                segment_param->AtlasWeight.push_back(WeightToPush);
                cout<<WeightToPush<<" ";
            }
            cout<<endl;
//            segment_param->AtlasWeight=atof(argv[++i]);
//            if(segment_param->AtlasWeight<0){
//                cout<<"Weight has to be between 0 and 1"<<endl;
//                segment_param->AtlasWeight=0;
//            }
//            else if(segment_param->AtlasWeight>1){
//                cout<<"Weight has to be between 0 and 1"<<endl;
//                segment_param->AtlasWeight=1;
//            }
//            cout<<"AtlasWeight is "<<segment_param->AtlasWeight<<endl;
        }
        else if(strcmp(argv[i],"-AtlasSmoothing")==0 &&(i+1)<argc){
            int numbSmoothing=atoi(argv[++i]);
            if (i+numbSmoothing>argc) {
                fprintf(stderr,"Err:\t incompatible number of smoothing.\n");
                return 1;
            }
            cout << "AtlasSmoothing are ";
            if (numbSmoothing>0) {
                cout<<numbSmoothing<<" ";
                segment_param->AtlasSmoothing.clear();
            }
            for (int w=0; w<numbSmoothing; w++) {
                float SmoothingToPush=atof(argv[++i]);
                segment_param->AtlasSmoothing.push_back(SmoothingToPush);
                cout<<SmoothingToPush<<" ";
            }
            cout<<endl;
            
//            segment_param->AtlasSmoothing=atof(argv[++i]);
//            cout<<"AtlasSmoothing is "<<segment_param->AtlasSmoothing<<endl;
        }
        else if(strcmp(argv[i],"-AtlasAveraging")==0 &&(i+1)<argc){
            int numbSmoothing=atoi(argv[++i]);
            if (i+numbSmoothing>argc) {
                fprintf(stderr,"Err:\t incompatible number of smoothing.\n");
                return 1;
            }
            cout << "AtlasAveraging are ";
            if (numbSmoothing>0) {
                cout<<numbSmoothing<<" ";
                segment_param->AtlasAveraging.clear();
            }
            for (int w=0; w<numbSmoothing; w++) {
                float SmoothingToPush=atof(argv[++i]);
                segment_param->AtlasAveraging.push_back(SmoothingToPush);
                cout<<SmoothingToPush<<" ";
            }
            cout<<endl;

//            segment_param->AtlasSmoothing=atof(argv[++i]);
//            cout<<"AtlasSmoothing is "<<segment_param->AtlasSmoothing<<endl;
        }

        else if(strcmp(argv[i],"-KernelSize")==0 &&(i+1)<argc){
            segment_param->KernelSize=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-quantMin")==0 &&(i+1)<argc){
            segment_param->quantMin=atof(argv[++i]);
            if (segment_param->OutliersMod>0) {
                segment_param->quantMin=0;
            }
            cout<<"quantMin is "<<segment_param->quantMin<<endl;
        }
        else if(strcmp(argv[i],"-quantMax")==0 &&(i+1)<argc){
            segment_param->quantMax=atof(argv[++i]);
            if (segment_param->OutliersMod>0) {
                segment_param->quantMin=0;
            }
            cout<<"quantMax is "<<segment_param->quantMax<<endl;
        }
        else if(strcmp(argv[i],"-EM_out")==0 &&(i+1)<argc){
            segment_param->flag_EMfirst_out=1;
            segment_param->filename_EMfirst=argv[++i];
        }


        else if(strcmp(argv[i], "-bc_out") == 0 && (i+1)<argc){
            segment_param->filename_datacorrected = argv[++i];
            segment_param->flag_bc_out=1;
        }
        else if(strcmp(argv[i], "-bf_out")==0 && (i+1)<argc){
            segment_param->filename_correction=argv[++i];
            segment_param->flag_bf_out=1;
        }
        else if(strcmp(argv[i], "-max_iter") == 0 && (i+1)<argc){
            segment_param->maxIteration=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-min_iter") == 0 && (i+1)<argc){
            segment_param->minIteration=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-max_iterBF") == 0 && (i+1)<argc){
            segment_param->maxIterationBF=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-maxRunEM") == 0 && (i+1)<argc){
            segment_param->MaxRunEM=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-txt_out")==0 && (i+1)<argc){
            segment_param->filename_datatxt=argv[++i];
            segment_param->flag_data_out=1;
        }
        else if(strcmp(argv[i], "-MRF_out")==0 && (i+1)<argc){
            segment_param->filename_MRFOut=argv[++i];
            segment_param->flag_MRFOut=1;
        }
        else if(strcmp(argv[i], "-MRF")==0 && (i+1)<argc){
            switch(atoi(argv[++i])){
            case 0:{ // No MRF applied
                segment_param->flag_MRF=0;
                segment_param->flag_MRFPost=0;
                segment_param->flag_MRFIn=0;
                break;
            }
            case 1 :{ // MRF applied during the processing of the model
                segment_param->flag_MRF=1;
                segment_param->flag_MRFPost=0;
                segment_param->flag_MRFIn=1;
                break;
            }
            case 2:{ // MRF applied only after the processing of the model
                segment_param->flag_MRF=0;
                segment_param->flag_MRFPost=1;
                segment_param->flag_MRFIn=0;
                break;
            }
            case 3:{ // MRF applied both during and after the obtention of the model
                segment_param->flag_MRF=1;
                segment_param->flag_MRFPost=1;
                segment_param->flag_MRFIn=1;
                break;
            }
            default:{ // in case of any other key chosen : by default no MRF chosen
                segment_param->flag_MRF=0;
                segment_param->flag_MRFPost=0;
                segment_param->flag_MRFIn=0;
                break;
            }
            }
        }
        else if(strcmp(argv[i], "-GMRFPost")==0 && (i+1)<argc){
            segment_param->filename_GMatrixPost=argv[++i];
            segment_param->flag_GMatrixPost=1;
            segment_param->flag_MRFPost=1;
        }
        else if(strcmp(argv[i], "-GMRF")==0 && (i+1)<argc){
            segment_param->filename_GMatrix=argv[++i];
            segment_param->flag_GMatrix=1;
            segment_param->flag_GMatrixIn=1;
            segment_param->flag_MRF=1;
            cout<<"Check seg param"<<endl;
        }

        else if(strcmp(argv[i], "-txt_in")==0 && (i+1)<argc){
            segment_param->filename_intxt=argv[++i];
            segment_param->flag_intxt=1;
            segment_param->flag_input=1;
        }
        else if(strcmp(argv[i], "-DPGauss_std")==0 && (i+1)<argc){
            segment_param->DPGauss_std=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-uniformTC")==0 && (i+1)<argc){
            segment_param->uniformTypeChange=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-acceptType")==0 && (i+1)<argc){
            segment_param->AcceptanceType=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-BICFP")==0 && (i+1)<argc){
            segment_param->BICFP=atoi(argv[++i]);
        }

        else if(strcmp(argv[i], "-scratchDir") == 0 && (i+1)<argc){
            g_scratch_dir = argv[++i];
        }

        else{
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }

    if (segment_param->flag_NormMask && !segment_param->flag_mask) {
        cout<< "No mask normalisation if no mask available"<<endl;
        segment_param->flag_NormMask=0;
    }
    if(segment_param->flag_out==0){
        fprintf(stderr,"Err:\tThe output image name has to be defined.\n");
        return 1;
    }

    if(!segment_param->flag_input){
        fprintf(stderr,"Err:\tThe T1 image name has to be defined.\n");
        Usage(argv[0]);
        return 1;
    }
    
    if (segment_param->flag_manual_priors_out) {
        segment_param->OutliersMod=6;
    }
    if (segment_param->flag_manual_priors_bnb) {
        if(segment_param->flag_bnbComb){
            segment_param->OutliersMod=8;
        }
        else{
            segment_param->OutliersMod=7;
        }
    }
    
    if (segment_param->NumberAddedModalities.size()==0) { // If nothing in numberAddedModalities, means that we consider all modalities together from first step
        segment_param->NumberAddedModalities.push_back(segment_param->numbmodal);
    }
    // Check that NumberAddedModalities is compatible with total number of modalities presented
    int sumModAdded=0;
    int numbmodal=segment_param->numbmodal;
    int IndexTotal=-1;
    int sizeNAM=segment_param->NumberAddedModalities.size();
    for (int s=0; s<sizeNAM; s++) {
        sumModAdded+=segment_param->NumberAddedModalities[s];
        if (sumModAdded>numbmodal) {
            IndexTotal=s;
            int Corrected=numbmodal-(sumModAdded-segment_param->NumberAddedModalities[s]);
            segment_param->NumberAddedModalities[s]=Corrected;
            break;
        }
        if (sumModAdded==numbmodal) {
            IndexTotal=s;
            break;
        }
    }
    // In case we are trying to add more modalities than that truly exist
    if (IndexTotal<sizeNAM-1 && IndexTotal>=0) {
        cout<<"More modprog that what can be done"<<endl;
        for (int i=sizeNAM-1; i>IndexTotal; i--) {
            segment_param->NumberAddedModalities.pop_back();
        }
    }
    // In case we have not reached the number of modalities expected, we add all of the remaining ones in the last step.
    if (IndexTotal==-1) {
        cout<<"Not enough modprog compared to available mod"<<endl;
        int ModToAdd=numbmodal-sumModAdded;
        segment_param->NumberAddedModalities.push_back(ModToAdd);
    }
    
    

    //  vector<nifti_image *> VectorPriors=ReadFromFilenamesVector(segment_param->filename_priors);
    //  nifti_image* MaskToUse=ReadFromFilename(segment_param->filename_mask);
    //  vector<nifti_image* > VectorImageToSegment=ReadFromFilenamesVector(segment_param->filename_input);
    //  TreeEM * TreeReconstruct=ReadTreeFromTextFile(g_scratch_dir + "/TestTreeText.txt",segment_param);

    //  int * TestHierarch=TreeReconstruct->GetChild(3)->GetChild(4)->GetHierarchy();
    //  int DepthTree=TreeReconstruct->GetNumberLevels();
    //  vector<TreeEM*> Level2Elements=TreeReconstruct->GetAllTreesFromLevel(2);
    //  TreeReconstruct->SaveTreeInTextFile(g_scratch_dir + "/BW04TestSaving.txt",segment_param);
    
    
    // Correcting / adapting / Changing Outliers Combined to be directly useable afterwards
    if (segment_param->OutliersCombined.size()>0) {
        segment_param->OutliersMod=5;
        // Check for each combination that the prior is indeed available
        vector<vector<int> > OutliersCombinedToCheck=segment_param->OutliersCombined;
        int numbOutliersComb=OutliersCombinedToCheck.size();
        int numbpriors=segment_param->numb_classes;
        int * UsedPriors=new int[numbpriors]; // To count if how many times priors are considered in combinations : allows to check for mistakes of redundancy and add missing priors
        for (int p=0; p<numbpriors; p++) {
            UsedPriors[p]=0;
        }
        for (int o=0; o<numbOutliersComb; o++) {
            int comb=OutliersCombinedToCheck[o].size();
            for (int c=0; c<comb; c++) {
                if (OutliersCombinedToCheck[o][c]>=numbpriors) {
                    segment_param->OutliersMod=3;
                    break;
                }
                cout<<"Prior "<<OutliersCombinedToCheck[o][c]<<" is used"<<endl;
                UsedPriors[OutliersCombinedToCheck[o][c]]++;
            }
        }
        for (int p=0; p<numbpriors; p++) {
            if (UsedPriors[p]>1) { // Impossible to use same atlas twice
                segment_param->OutliersMod=3;
                break;
            }
            if (UsedPriors[p]==0) { // Have to consider other priors as well if not combined together.
                vector<int> PriorsForgotten;
                PriorsForgotten.push_back(p);
                segment_param->OutliersCombined.push_back(PriorsForgotten);
            }
        }
    }
    
    
    TreeEM * TreeTest;
    if(segment_param->flag_intxt){
        if (segment_param->flag_AdaptTransform || segment_param->flag_inDC) {
            vector<nifti_image * > ImageIOVector=ReadFromFilenamesVector(segment_param->filename_IOAverage);
            vector<nifti_image *> PriorsIOVector=SmoothAndAverageImagesForPriors(ImageIOVector[0], ImageIOVector[1], segment_param->weightIOAverage, segment_param,0);
            cout<<"Averaging of IO is done"<<endl;
            int sizeAverage=ImageIOVector.size();
            for (int n=0; n<sizeAverage; n++) {
                nifti_image_free(ImageIOVector[n]);
                ImageIOVector[n]=NULL;
            }
//            nifti_set_filenames(PriorsIOVector[1], g_scratch_dir + "/TestSIO.nii.gz", 0, 0);
//            nifti_image_write(PriorsIOVector[1]);
            vector<nifti_image * > ImageGCVector=ReadFromFilenamesVector(segment_param->filename_GCAverage);

            vector<nifti_image *> PriorsGCVector=SmoothAndAverageImagesForPriors(ImageGCVector[0], ImageGCVector[1], segment_param->weightGCAverage, segment_param,1);
//            nifti_set_filenames(PriorsGCVector[1], g_scratch_dir + "/TestSGC.nii.gz", 0, 0);
//            nifti_image_write(PriorsGCVector[1]);
            sizeAverage=ImageGCVector.size();
            for (int n=0; n<sizeAverage; n++) {
                nifti_image_free(ImageGCVector[n]);
                ImageGCVector[n]=NULL;
            }
            cout<<"Averaging of GC is done"<<endl;
            vector<float> AtlasWeightSaved;
            vector<float> AtlasSmoothSaved;
            int sizeS=segment_param->AtlasSmoothing.size();
            int sizeW=segment_param->AtlasWeight.size();
            for (int w=0; w<sizeW; w++) {
                AtlasWeightSaved.push_back(segment_param->AtlasWeight[w]);
            }
            for (int s=0; s<sizeS; s++) {
                AtlasSmoothSaved.push_back(segment_param->AtlasSmoothing[s]);
            }
            TreeTest=ReadTreeFromTextFileWithAdapt(segment_param->filename_intxt, segment_param, PriorsGCVector, PriorsIOVector, segment_param->filename_AdaptTransform);
            segment_param->AtlasWeight.clear();
            segment_param->AtlasSmoothing.clear();
            if (sizeS>0) {
                for (int s=0; s<sizeS; s++) {
                    segment_param->AtlasSmoothing.push_back(AtlasSmoothSaved[s]);
                }
            }
            if (sizeW>0) {
                for (int w=0; w<sizeW; w++) {
                    segment_param->AtlasWeight.push_back(AtlasWeightSaved[w]);
                }
            }
            nifti_set_filenames(PriorsIOVector[0], (g_scratch_dir + "/TestPriors.nii.gz").c_str(), 0, 0);
            nifti_image_write(PriorsIOVector[0]);
            cout<<"Tree with adaptation rebuilt"<<endl;
//            Priors normalisation
//            int * ModalitiesTest=TreeTest->GetModalities(segment_param);
//            TreeTest->SavePriorsAdaptedHierarchy(segment_param);
            
//            TreeTest->NormaliseDataImage();
            
            TreeTest->NormalisePriors();
//            Adapted priors normalisation
            TreeTest->NormalisePriorsAdapted();
            TreeTest->SavePriorsAdaptedHierarchy(segment_param);
            TreeTest->UpdatePartPriorsAdapted();
            TreeTest->UpdateNonNormResp(segment_param);
            TreeTest->UpdateNormResp();
            TreeTest->UpdateNormRespRoot();
//            TreeTest->ZeroRootRespInUniform();
            TreeTest->SaveAllClasses(g_scratch_dir + "/TestClasses.nii.gz", segment_param);
            if(segment_param->flag_MRF){
                float * GMatrixToSet=NULL;
                if(segment_param->flag_GMatrixIn){
                    GMatrixToSet=TreeTest->PrepareGMatrixFromFile(segment_param->filename_GMatrix, segment_param->flag_optMRFOut);
//                    TreeTest->PrepareGInfoFromFile(segment_param->filename_GMatrix);
//                    GMatrixToSet=TreeTest->CreateGMatrixFromInfo( segment_param->flag_optMRFOut,segment_param->filename_GMatrix);
                    //                        GMatrixToSet =this->PrepareGMatrixFromFile(segment_param->filename_GMatrix,segment_param->flag_optMRFOut);
                    if(GMatrixToSet==NULL){
                        segment_param->flag_GMatrix=0;
                    }
                    TreeTest->SetGMatrix(GMatrixToSet);
                    if (GMatrixToSet!=NULL) {
                        delete [] GMatrixToSet;
                        GMatrixToSet=NULL;
                    }
                }
                else {
                    GMatrixToSet=TreeTest->MRFOptSolveLS(segment_param);
                    if(GMatrixToSet!=NULL){
                        delete [] GMatrixToSet;
                        GMatrixToSet=NULL;
                    }
                }
            }
            if(segment_param->flag_MRF){
                TreeTest->UpdateMRF(segment_param);
            }
            TreeTest->UpdateNonNormWeights();
            //
            TreeTest->UpdateNormWeights();
            
            if (segment_param->flag_CovPriors) {
                TreeTest->ModifyCovPriors(segment_param);
            }
            TreeTest->UpdateMRF(segment_param);
            segment_param->bias_order=0;
            TreeTest->SaveTreeInTextFile(g_scratch_dir + "/TestFirstSave.txt", segment_param);
        }
        else{
        TreeTest=ReadTreeFromTextFile(segment_param->filename_intxt,segment_param,NULL);
        cout<<"Tree reconstructed"<<endl;
        }
    }
    else{
        TreeTest=CreateTreeSegmentationInit(segment_param);
    }
    

    
    TreeTest->ClearSMChecks();
//    TreeTest->SaveTreeInTextFile(g_scratch_dir + "/TestFirstSaveTestb.txt", segment_param);
    
    // Compare numbmodal in segment_param and in TreeTest to see what has to be modified in
    if (segment_param->numbmodal!=TreeTest->GetNumberModalities()) {
        cout<<"Will have to progressively add modalities to Tree from TxtFile"<<endl;
    }
    
    
    if(! segment_param->flag_data_out){
        string FilenameOut=nifti_makebasename(segment_param->filename_out[0].c_str());
        stringstream ss;
        ss << TreeTest->GetNumberModalities();
        string cb=ss.str();
        FilenameOut+="_"+cb+".txt";
        segment_param->filename_datatxt=FilenameOut;
    }
    
    // Modify needed segment_param aspects if not given in order to make all saving operations possible
    // Taking care of PriorsAdapted
    if( !segment_param->flag_PriorsAdaptedOut || segment_param->flag_intxt){
        vector<string> FilenamePAVector;
//        int numbPA=TreeTest->GetPriorsAdaptedVector().size();
        int numbPA=TreeTest->GetPriorsAdaptedVectorParent().size();
        for(int c=0;c<numbPA;c++){
            string FilenamePA=nifti_makebasename(segment_param->filename_out[0].c_str());
            int Index=FilenamePA.find_last_of('/');
            string FilenamePA_b=FilenamePA.substr(0,Index+1);
            string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            stringstream ss;
            ss << c;
            string cb=ss.str();
            FilenamePA=FilenamePA_b+"PriorsAdapted_"+cb+FilenamePA_e+".nii.gz";
            FilenamePAVector.push_back(FilenamePA);
        }
        
        segment_param->filename_PriorsAdaptedOut=FilenamePAVector;
    }
    
    if(segment_param->flag_savePriors){
        TreeTest->SavePriorsAdapted(segment_param);
    }
    
    // Data corrected
    if(! segment_param->flag_bc_out){
        string FilenameBC=nifti_makebasename(segment_param->filename_out[0].c_str());
        int Index=FilenameBC.find_last_of('/');
        string FilenameBC_b=FilenameBC.substr(0,Index+1);
        string FilenameBC_e=FilenameBC.substr(Index+1,FilenameBC.length());
        FilenameBC=FilenameBC_b+"DataCorrected_"+FilenameBC_e+".nii.gz";
        segment_param->filename_datacorrected=FilenameBC;
    }
    
    // Bias field if bc_order >0
    if(! segment_param->flag_bf_out && segment_param->bias_order>0){
        string FilenameBC=nifti_makebasename(segment_param->filename_out[0].c_str());
        int Index=FilenameBC.find_last_of('/');
        string FilenameBC_b=FilenameBC.substr(0,Index+1);
        string FilenameBC_e=FilenameBC.substr(Index+1,FilenameBC.length());
        FilenameBC=FilenameBC_b+"BF_"+FilenameBC_e+".nii.gz";
        segment_param->filename_correction=FilenameBC;
        segment_param->flag_bf_out=1;
    }
    
    // MRF Out
    if(! segment_param->flag_MRFOut){
        string FilenameBC=nifti_makebasename(segment_param->filename_out[0].c_str());
        int Index=FilenameBC.find_last_of('/');
        string FilenameBC_b=FilenameBC.substr(0,Index+1);
        string FilenameBC_e=FilenameBC.substr(Index+1,FilenameBC.length());
        FilenameBC=FilenameBC_b+"MRFOut_"+FilenameBC_e+".nii.gz";
        segment_param->filename_MRFOut=FilenameBC;
        if(segment_param->flag_GMatrix){
            segment_param->flag_MRFOut=1;
        }
    }
    

    if(segment_param->numb_classes<2){
        fprintf(stderr,"Please use either the -priors or -nopriors options with at least 2 classes\n");
        return 1;
    }

    TreeEM * TreeResult;
    if (segment_param->flag_intxt && !segment_param->flag_inDC) {
        TreeTest->SetFlagCovPriors(segment_param->CovPriorsType);
        TreeTest->SetFlagMeanPriors(segment_param->MeanPriors);
        if (segment_param->flag_manual_priors) {
            vector<TreeEM *> GeneralClasses=TreeTest->GetGeneralClassesVector();
            vector<TreeEM *> GeneralOutlierClasses;
            if (TreeTest->GetFlagOutliers()==3) {
                GeneralOutlierClasses=TreeTest->GetNodeOutlier()->GetChildren();
                nifti_image * PriorsOutliers=TreeTest->BuildConstantPriors(TreeTest->GetNodeOutlier()->GetNormWeight());
                nifti_image * PriorsInliers=TreeTest->CreateNormaliseOppositeImage(PriorsOutliers);
                TreeTest->GetChild(0)->SetPriors(PriorsInliers);
                TreeTest->GetChild(1)->SetPriors(PriorsOutliers);
            }
            int numbPriorsFiles=segment_param->filename_priors.size();
            int numbGC=GeneralClasses.size();
            if (numbPriorsFiles==numbGC) {
                for (int c=0; c<numbPriorsFiles; c++) {
                    nifti_image * Priors=ReadFromFilename(segment_param->filename_priors[c].c_str());
                    GeneralClasses[c]->SetPriors(Priors);
                    if (GeneralOutlierClasses.size()>0) {
                        nifti_image * PriorsCopy=GeneralClasses[c]->CopyPriors();
                        GeneralOutlierClasses[c]->SetPriors(PriorsCopy);
                    }
                }
            }
        }
        TreeTest->UpdateNonNormResp(segment_param);
        TreeTest->UpdateNormResp();
        TreeTest->UpdateNormRespRoot();
        TreeTest->NormalisePriors();
        TreeEM * TreeTmp=TreeTest->BuildTreeWithAddedModalityFromExistingModel(TreeTest, segment_param);
        delete TreeTest;
        TreeTest=NULL;
        vector<float> AtlasWeightSaved;
        vector<float> AtlasSmoothSaved;
        int sizeS=segment_param->AtlasSmoothing.size();
        int sizeW=segment_param->AtlasWeight.size();
        for (int w=0; w<sizeW; w++) {
            AtlasWeightSaved.push_back(segment_param->AtlasWeight[w]);
        }
        for (int s=0; s<sizeS; s++) {
            AtlasSmoothSaved.push_back(segment_param->AtlasSmoothing[s]);
        }
        segment_param->AtlasWeight.clear();
        segment_param->AtlasSmoothing.clear();
        TreeResult=TreeTmp->RunFullBiASM_ter(segment_param);
//        TreeResult=TreeTmp->RunCompleteBaMoS(segment_param);
    }
    else if(segment_param->flag_intxt && segment_param->flag_inDC){
//        TreeResult=TreeTest->RunFullBiASM_ter(segment_param);
        float CLL=0;
        int Iter=0;
        float OCLL=0;
//        TreeTest->FullAdaptPriors(segment_param);
        
        
        vector<float> AtlasWeightSaved;
        vector<float> AtlasSmoothSaved;
        int sizeS=segment_param->AtlasSmoothing.size();
        int sizeW=segment_param->AtlasWeight.size();
        for (int w=0; w<sizeW; w++) {
            AtlasWeightSaved.push_back(segment_param->AtlasWeight[w]);
        }
        for (int s=0; s<sizeS; s++) {
            AtlasSmoothSaved.push_back(segment_param->AtlasSmoothing[s]);
        }
        segment_param->AtlasWeight.clear();
        segment_param->AtlasSmoothing.clear();
        
        
        
        
        TreeTest->RunFullEM(CLL, OCLL, Iter, segment_param);
        int sizeInput=segment_param->filename_input.size();
        if(segment_param->flag_FurtherAddedModa || TreeTest->GetNumberModalities()<=sizeInput){
            cout<<"Need to try adding modality in existing scheme"<<endl;
           TreeEM* TreeNewToUse=TreeTest->BuildTreeWithAddedModalityFromExistingModel(TreeTest, segment_param);
            TreeResult=TreeNewToUse->RunBaMoSLoop(segment_param);
            segment_param->flag_BiASM=0;
        }
        else{
            if (segment_param->flag_BiASM) {
                if (sizeS>0) {
                    for (int s=0; s<sizeS; s++) {
                        segment_param->AtlasSmoothing.push_back(AtlasSmoothSaved[s]);
                    }
                }
                if (sizeW>0) {
                    for (int w=0; w<sizeW; w++) {
                        segment_param->AtlasWeight.push_back(AtlasWeightSaved[w]);
                    }
                }
                TreeTest->FullAdaptPriors(segment_param);
                segment_param->AtlasWeight.clear();
                segment_param->AtlasSmoothing.clear();
                TreeResult=TreeTest->RunBaMoSLoop(segment_param);
            }
            else {
                TreeResult=TreeTest;
            }
        }
    }
    else if(!segment_param->flag_intxt || segment_param->flag_BiASM){
        if (segment_param->NumberAddedModalities.size()>0) {
            TreeResult=TreeTest->RunCompleteBaMoS(segment_param);
        }
        else{
        switch(segment_param->SMOrder){
        case 0 :
            //          TreeTest->SaveTmpResultMasked(TreeTest->GetChild(0)->GetNormResp(),g_scratch_dir + "/NormRespInit.nii.gz");
            //          TreeResult=TreeTest->RunFullBiASM(segment_param);
                TreeTest->UpdateNonNormResp(segment_param);
                TreeTest->UpdateNonNormWeights();
                TreeTest->UpdateNormResp();
                TreeTest->UpdateNormRespRoot();
                TreeResult=TreeTest->RunFullBiASM_ter(segment_param);
            break;
        default :
            TreeResult=TreeTest->RunFullBiASM_ter(segment_param);
                
        }
        }
    }

    // Case we apply an MRF step to the obtained model
    if(segment_param->flag_MRFPost){
        float * GMatrixPostToUse;
        segment_param->flag_MRF=1;
        if(segment_param->flag_GMatrixPost){
            GMatrixPostToUse=TreeResult->PrepareGMatrixFromFile(segment_param->filename_GMatrixPost,segment_param->flag_optMRFOut);
            if(GMatrixPostToUse!=NULL){
                segment_param->flag_GMatrix=1;
            }
        }
        else{
            GMatrixPostToUse=TreeResult->MRFOptSolveLS(segment_param);
        }
        TreeResult->SetGMatrix(GMatrixPostToUse);
        delete[] GMatrixPostToUse;
        GMatrixPostToUse=NULL;
        float LL=0;
        float OldLL=0;
        int Iteration=0;
        TreeResult->RunFullEM(LL,OldLL,Iteration,segment_param);
    }




    //      TreeEM * TreeResult=TreeTest->RunFullBiASM_bis(segment_param);
    int SavingFiles=segment_param->filename_out.size();
    switch (SavingFiles){

            case 1 :{
            cout <<"There is only one filename out"<< endl;
        TreeResult->SaveAllClasses(segment_param->filename_out[0],segment_param);
            }
        break;

    default :{
        TreeResult->SaveAllClasses(segment_param->filename_out[0],segment_param);
        TreeResult->SaveGeneralClasses(segment_param->filename_out[1],segment_param);
    }
    }




    if(! segment_param->flag_data_out){
        string FilenameOut=nifti_makebasename(segment_param->filename_out[0].c_str());
        FilenameOut+=".txt";
        segment_param->filename_datatxt=FilenameOut;
    }
    // SaveMRFImage must run before SaveTreeInTextFile: if no MRF was actually
    // computed (GetMRF()==NULL), SaveMRFImage resets flag_MRFOut to 0 so that
    // SaveTreeInTextFile correctly omits the "MRFImage" entry.  Writing the txt
    // first left a stale MRFImage reference pointing to a file that was never
    // created, causing nifti_image_read errors in all downstream Seg_Analysis calls.
    if(segment_param->flag_MRF){
        cout << "saving MRF image" << endl;
        TreeResult->SaveMRFImage(segment_param);
    }

    TreeResult->SaveTreeInTextFile(segment_param->filename_datatxt,segment_param);
    TreeResult->SaveBFCorrectedData(segment_param->filename_datacorrected);

    if (segment_param->AtlasWeight.size()>0) {
        if(segment_param->AtlasWeight[0]>0){
            TreeResult->SavePriorsAdapted(segment_param);
        }
    }




    if(segment_param->flag_Count){
        TreeResult->ModifyCountFiles(segment_param->filename_Count);
    }
    if(segment_param->flag_bf_out){
        cout << "Saving BF correction "<< endl;
        TreeResult->SaveBFCorrection(segment_param->filename_correction);
    }

    delete TreeResult;
    TreeResult=NULL;
    std::cout << "Hello, World!\n";
    delete [] segment_param;
    return EXIT_SUCCESS;
}




