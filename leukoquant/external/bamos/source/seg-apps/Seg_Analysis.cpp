
//  Modified by Stylianos Charalampous on 08/05/2026.
#include "_Seg_Analysis.h"
#include <iostream>
#include<sstream>
#include<string>
#include <algorithm>
using namespace std;
static int TypePackage=1;
static int Default=1;
static vector<int> DefaultElIn;
static vector<int> DefaultElOut;
static vector<string> ModalitiesCode;
static int Side[] = {1, 0, 1, 0, 1,0 ,1,0};



void Usage(char *exec)
{
    printf("\n  EM Segmentation Analysis:\n  Usage ->\t%s -in <filename1> [OPTIONS]\n\n",exec);
    printf("\n  BaMoS Segmentation Analysis:\n  Usage ->\t%s -in <filename1> [OPTIONS]\n\n",exec);
    printf("  * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-in <filename1> \t\t| Filename of the segmentation result to analyse\n");
    printf("\t\t\t\t| The input image should be 2D, 3D or 4D images. 2D images should be on the XY plane.\n");



    printf("* * * * OPTIONS FOR RECONSTRUCTION OF THE MODEL * * * * \n\n");
    printf("\t-inTxt2 <filename1> <filename2>\t\t| Filenames containing the information on the tree structure (filename1) and the corresponding 4D image (filename2) for reconstruction of the model\n");
    printf("\t-inChangePath <PathString> \t\t| Path to the folder where all model files are stored (Tree structure, Tree classes, Mask, Corrected Data\n");
    printf("\t-scratchDir <PathString> \t\t| Path to scratch directory for debug output files (default: /scratch0)\n");
    printf("\t-output_dir <PathString> \t\t| Path to the folder where all output files should be saved\n");
    printf("\t-inModa <int> <list of ints> \t\t| Indicates the modalities to be considered in the model, first the number of modalities then there encoding\n");
    printf("\t\t 1 : T1\n");
    printf("\t\t 2 : T2\n");
    printf("\t\t 3 : FLAIR\n");
    printf("\t\t 4 : PD\n\n");

    printf("* * * * NEEDED DATA FOR NEW LESION SEGMENTATION * * * * \n\n");
    printf("\t-inPriorsCGM <filename> \t\t| Filename of atlas CGM\n");
    printf("\t-inPriorsDGM <filename> \t\t| Filename of atlas DGM\n");
    printf("\t-inPriorsECSF <filename> \t\t| Filename of atlas ECSF\n");
    printf("\t-inPriorsICSF <filename> \t\t| Filename of atlas ICSF\n");
    printf("\t-inPriorsWM <filename> \t\t| Filename of atlas WM\n");
    printf("\t-inPriorsOut <filename> \t\t| Filename of atlas Out (NB)\n\n");

    printf("* * * * OPTIONS FOR CORRECTION OF LESION SEGMENTATION / NEW LESION * * * * \n\n");
    printf("\t-inCST <filename> \t\t| Filename of atlas of corticospinal tracts\n");
    printf("\t-typeVentrSeg <bool> \t\t| Flag indicating that the ventricles will be segmented using the parcellation results \n");
    printf("\t-inVentrSeg <filename> \t\t| Filename for an already created ventricle segmentation \n");
    printf("\t-ParcellationIn <filename> \t\t| Filename of parcellation file (required when typeVentrSeg=1) \n");
    printf("\t-inArtefact <filename> \t\t| Filename of artefact image to avoid some areas to be considered as lesions \n");
    printf("\t-WeightedSeg <int1> <float> <int3> \t\t| The first input concerns the type of Mahalanobis weighting adopted, the second input is the Mahalanobis threshold (default=3), and the third the class used for the comparison (default =1)\n");
    // TO DO TYPE OF WEIGHTED SEG TO BE USED
    printf("\t-Edginess float \t\t| Edginess threshold to be used (default = 0)\n");
    printf("\t-SegType <int>  \t\t| Type of segmentation strategy to adopt\n");
    // TO DO TYPE of SEGTYPE TO BE USED
    printf("\t-SP <bool>  \t\t| Flag for initial septum pellucidum correction (default=0)\n");
    printf("\t-CorrIV <bool>  \t\t| Flag for correction of inside ventricles artefacts (default=1) \n");
    printf("\t-TO <bool> \t\t | Flag to account for all outliers and do the post-processing at the voxel level (default=0)");
    printf("\t-juxtaCorr <bool> \t\t | Flag to try and correct for island of GM \n");
    printf("\t-correct \t\t | if indicated signals that the correction of lesion should be performed \n");
    printf("\t-connect \t\t | if indicated signals that the connected elements should be obtained \n");
    printf("\t-Simple \t\t | if indicated signals that a 3D image should be simply determined for potential lesion (and not 4D one for each outlier class as originally)\n");



    printf("* * * * OPTIONS FOR SPECIFIC OUTPUTS / SOME HAVE TO BE RUN INDEPENDENTLY OF LESION SEGMENTATION * * * * \n\n");
    printf("\t-WMMaps <bool> \t\t | Flag indicating that the Mahalanobis distance should be calculated \n");
    printf("\t-IO <bool> \t\t | Flag indicating if the inlier/outlier segmentation should be output \n");
    printf("\t-Euc <bool> \t\t | Flag indicating that the euclidean distance from given mask (through inLes or inLesCorr) must be output \n");



    printf("* * * * OPTIONS FOR LAPLACE LAYERS BUILDING AND LOCAL SUMMARY * * * * \n\n");
    printf("*** Creation of the Distance maps from the lobar segmentation (requires -inLobes -inLesCorr (can be the full T1 image but must be float datatype) and -mask\n");
    printf("\t-inLobes <integer> <list filenames> \t\t | Number of lobar masks and corresponding files to build the corresponding distance map. \n");

    printf("*** Creation of Local summary \n ");
    printf("\t-LS <bool> \t\t | Flag indicating that the local summary should be performed \n");
    printf("\t-inQuadrant <filename> \t\t | Filename where the zonal separation is stored (in lobes but historically in quadrants) \n");
    printf("\t-inLes <filename> \t\t | Filename with the lesion segmentation to assess locally (alternatively can use -inLesCorr) \n");
    printf("\t-LaplaceNormSol <filename> \t\t | Filename of the layer discretisation \n");


    printf("*** Creation of Laplace based layers \n");
    printf("\t-LapIO <infile> <outfile> \t\t | Names of the two binary segmentation used to delineate the borders of the volume on which to build the Laplace solution. \n");
    printf("\t-nameLap <string> \t\t | name to include in the Laplace layer formulation of file saving\n");
    printf("\t-numbLaminae <integer> <list of integers> \t\t | number of layer numbers to calculate followed by the list of number of layers to calculate over the obtained normalised distance \n\n");


    printf("* * * * POSSIBLE RUN THAT DO NOT REQUIRE BAMOS MODEL * * * * \n\n");
    printf("*** Getting Seg vs Ref statistics for a set of segmentations. The process provides a text file with for each comparison a line with :Name,VolRef,VolSeg,LabelRef,LabelSeg,FP,FN,TP,DSC,AvDist,DE,DEFP,DEFN,OEFP,OEFN,OER,VD,TPR,FPR,TPRc,FPRc,FNRc\n" );
    printf("\t -inMaskVec <number of masks> <list of mask filenames> \t\t | indicates the number of masks to use and the list of the associated names (note that the masks must be in binary format)\n");
    printf("\t -inLesVec <number of Seg> <list of seg filenames> \t\t | indicates the number of segmentations to evaluate and the names of the corresponding files. Must be same number as number of masks \n");
    printf("\t -inRefVec <number of Ref> <list of reference filenames> \t\t | indicates the number of reference images and the corresponding filenames. Must be the same number as inLesVec and inMaskVec \n" );
    printf("\t -inNamesVec <number of Names> <list of names to give> \t\t | indicates the number of evaluations and the name attributed to each\n");
    printf("\t -outFile <name of file> \t\t | Name of the file in which the results of the evaluation will be printed \n\n");

    printf("*** Creating the intensity matching output automatically a text file with the polynomial coefficient for each matching pair and the matched transformed images. \n");
    printf("\t -maskMatch <mask filename> \t\t| Name of mask file (binary) over which the intensity matching will be made \n");
    printf("\t -matchFloat <number of images to match> <filenames> \t\t | Number of files to match and corresponding list of filenames \n");
    printf("\t -matchRef <number of reference images> <filenames> \t\t | Number of reference images for the matching and corresponding list of files (should be same number in ref and float)\n");
    printf("\t -orderFit <integer> \t\t | Order of the polynomial fit performed \n\n");

    printf("*** Averaging images together \n");
    printf("\t-outMean <filename> \t\t | indicates the name of the output file containing the result of the averaging \n");
    printf("\t-meanImages <integer> <filenames> | number of images to average and corresponding filenames \n\n");



    printf("* * * * * NEEDED OPTIONS FOR SIMULATION GENERATION * * * * \n\n");

    printf("*** Generation of Bias field images \n");
    printf("\t-BFGen <float1> <float2> <float3> \t\t | 3 float values to generate random bias field with first value the maximum order, the second the number of modalities and the third the max of range of coefficients \n" );
    printf("\t-BFGenText <textfile> \t\t | Filename for textfile where coefficients of BF are stored \n\n");

    printf("*** Generation of random affine transformation\n");
    printf("\t-AffGen <float> <float> \t\t | 2 float values with the range of rotation in degrees for random choice and range of translation \n");
    printf("\t-OutFile <filename> \t\t | Filename for the output (text file for the rigid transformation) \n");



    printf("\t-HS <filename> \t| Requirement for hard segmentation with filename to store it\n");
    printf("\t-DCh\t| Allows for calculation of the hard Dice coefficient \n");
    printf("\t-DCs\t| Allows for calculation of the soft Dice coefficient \n");
    printf("\t-TP\t| Allows for calculation of the number of true positives \n");
    printf("\t-TN\t| Allows for calculation of the number of true negatives \n");
    printf("\t-FP\t| Allows for calculation of the number of false positives \n");
    printf("\t-FN\t| Allows for calculation of the number of false negatives \n");
    printf("\t-mask <filename>\t| Filename of the brain-mask of the input image\n");
    printf("\t-compJoint <filename>\t| Filename of the segmentation to compare it to. Must contain the same number of classes\n");
    printf("\t-compMult <numberOfPartialSeg> <filename> ... <filename>\t| Number of files to form the segmentation result with their filenames. Must contain the same number of classes\n");

    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}


int main(int argc, char **argv)
{
    
    ModalitiesCode.push_back("T1");
    ModalitiesCode.push_back("T2");
    ModalitiesCode.push_back("FLAIR");
    ModalitiesCode.push_back("PD");

    DefaultElIn.push_back(50);
    DefaultElIn.push_back(51);
    DefaultElIn.push_back(52);
    DefaultElIn.push_back(53);
    
    DefaultElOut.push_back(45);
    DefaultElOut.push_back(46);
    DefaultElOut.push_back(48);
    DefaultElOut.push_back(49);
    DefaultElOut.push_back(37);
    DefaultElOut.push_back(38);
    DefaultElOut.push_back(58);
    DefaultElOut.push_back(59);
    DefaultElOut.push_back(60);
    DefaultElOut.push_back(61);
    
    if (argc < 1)
    {
        Usage(argv[0]);
        return 0;
    }

    // Set all default parameters to be used in the segmentation process
    SEG_ANALYSIS * segment_analysis = new SEG_ANALYSIS [1]();
    segment_analysis->flag_correctCerr=1;
    segment_analysis->flag_acceptJCMS=0;
    segment_analysis->flag_DCh=0;
    segment_analysis->flag_DCs=0;
    segment_analysis->flag_TP=0;
    segment_analysis->flag_TN=0;
    segment_analysis->flag_FP=0;
    segment_analysis->flag_FN=0;
    segment_analysis->flag_HS=0;
    segment_analysis->flag_inVentricleSeg=0;
    segment_analysis->flag_fileOut=0;
    segment_analysis->flag_numbActive=0;
    segment_analysis->flag_removespurious=1;
    segment_analysis->flag_LapAnalysis=0;
    segment_analysis->EdginessThreshold=0;
    segment_analysis->flag_VentrSeg=0;
    segment_analysis->flag_Parc=0;
    segment_analysis->flag_in=0;
    segment_analysis->flag_PV=0;
    segment_analysis->flag_GIFPriors=0;
    segment_analysis->filename_GIFPriors="";
    segment_analysis->flag_mask=0;
    segment_analysis->flag_Les=0;
    segment_analysis->flag_data=0;
    segment_analysis->flag_KLD=0;
    segment_analysis->flag_outliers=0;
    segment_analysis->flag_TextFile=0;
    segment_analysis->flag_SegTot=0;
    segment_analysis->LesionUniformRuleType=3;
    segment_analysis->LesionRuleType=3;
    segment_analysis->flag_connect=0;
    segment_analysis->flag_compJoint=0;
    segment_analysis->flag_compMult=0;
    segment_analysis->flag_Distance=0;
    segment_analysis->flag_NW=0;
    segment_analysis->ExploRadius=2;
    segment_analysis->flag_Secondary=0;
    segment_analysis->IndexWM=1;
    segment_analysis->IndexGM=0;
    segment_analysis->IndexCSF=2;
    segment_analysis->IndexOut=3;
    segment_analysis->IndexECSF=-1;
    segment_analysis->IndexICSF=-1;
    segment_analysis->IndexDGM=-1;
    segment_analysis->IndexCGM=-1;
    segment_analysis->Neigh=26;
    segment_analysis->thresh_quant=0.9;
    segment_analysis->ThresholdSP=6;
    segment_analysis->flag_SPCorrection=1;
    segment_analysis->MiniSize=3;
    segment_analysis->numbLaminae.push_back(4);
    segment_analysis->flag_segType=3;
    segment_analysis->flag_segWeighted=0;
    segment_analysis->weightThreshold=3;
    segment_analysis->weightCompClass=0;
    segment_analysis->flag_inConnect=0;
    segment_analysis->flag_inLesCorr=0;
    segment_analysis->flag_Analysis=0;
    segment_analysis->flag_refLesConnect=0;
    segment_analysis->flag_refLes=0;
    segment_analysis->flag_AcceptedGM=0;
    segment_analysis->flag_correctionLevel=2;
    segment_analysis->flag_CorrectIV=1;
    segment_analysis->flag_IO=1;
    segment_analysis->flag_checkInliers=0;
    segment_analysis->flag_TO=0;
    segment_analysis->vecLeavesToAdd.clear();
    segment_analysis->flag_Saving=0;
    segment_analysis->flag_CodedOutliers=0;
    segment_analysis->vecCheckModalities.push_back(1);
    segment_analysis->vecCheckModalities.push_back(3);
    segment_analysis->flag_oldLesion=0;
    segment_analysis->flag_Gaussian=0;
    segment_analysis->index_Gaussian=1;
    segment_analysis->vecThresh;
    segment_analysis->GTAnalysisType=1;
    segment_analysis->flag_simple=0;
    segment_analysis->flag_WMCheck=0;
    segment_analysis->flag_ITCheck=0;
    segment_analysis->flag_CSFCheck=0;
    segment_analysis->flag_CorExt=0;
    segment_analysis->flag_LapIn=0;
    segment_analysis->val_ext = 2;
    segment_analysis->thresh_cortex=100;
    segment_analysis->verbose_level=LOG_NOTHING;
    /* read the input parameters */
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-inChangePath") == 0 && (i+1)<argc){
            segment_analysis->flag_changePath=1;
            segment_analysis->filename_changePath=argv[++i];
        }

        else if(strcmp(argv[i], "-scratchDir") == 0 && (i+1)<argc){
            g_scratch_dir = argv[++i];
        }

        else if(strcmp(argv[i], "-inMNI2") == 0 && (i+2)<=argc){
            char * MNITextFile=argv[++i];
            segment_analysis->flag_MNITransform=1;
            segment_analysis->filename_MNITransform=MNITextFile;
            segment_analysis->flag_Quadrant=1;
            char * MNITemplateFile=argv[++i];
            segment_analysis->flag_MNITemplate=1;
            segment_analysis->filename_MNITemplate=MNITemplateFile;
        }
        else if(strcmp(argv[i], "-inTxt2") == 0 && (i+2)<=argc){
            char * TreeTextFile=argv[++i];
            segment_analysis->flag_TextFile=1;
            segment_analysis->filename_TextFile=TreeTextFile;
            segment_analysis->flag_Les=1;
            char * TreeSegFile=argv[++i];
            segment_analysis->flag_SegTot=1;
            segment_analysis->filename_SegTot=TreeSegFile;
        }
        else if (strcmp(argv[i],"-SpinalCord") == 0 &&(i+1)<argc){
            segment_analysis->flag_OtherSeg =1;
            segment_analysis->flag_sc=1;
            segment_analysis->dist_ring=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-VesselSeg") == 0 && (i+1)<argc){
            segment_analysis->flag_VesselRuleTextFile=1;
            segment_analysis->flag_OtherSeg=1;
            segment_analysis->filename_VesselRuleTextFile=argv[++i];
        }
        else if(strcmp(argv[i], "-Prion") == 0 && (i+1)<argc){
            segment_analysis->flag_OtherSeg=1;
            segment_analysis->flag_prion=1;
        }
        else if(strcmp(argv[i], "-GetBB") == 0 && (i+1)<argc){
            segment_analysis->flag_getBB=1;
        }
        else if(strcmp(argv[i], "-GetCH") == 0 && (i+1)<argc){
            segment_analysis->flag_getCH=1;
        }
        else if(strcmp(argv[i], "-indexCSF") == 0 && (i+1)<argc){
            segment_analysis->IndexCSF=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-indexGM") == 0 && (i+1)<argc){
            segment_analysis->IndexGM=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-indexWM") == 0 && (i+1)<argc){
            segment_analysis->IndexWM=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-indexOut") == 0 && (i+1)<argc){
            segment_analysis->IndexOut=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-GetEntropy") == 0 && (i+2)<argc){
            segment_analysis->flag_getEntropy=1;
            segment_analysis->RadEntropy=atoi(argv[++i]);
            segment_analysis->MaxValueEntr=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-GetGrad") == 0 && (i+1)<argc){
            segment_analysis->flag_Grad=1;
        }
        else if(strcmp(argv[i], "-inRuleTxt") == 0 && (i+1)<argc){
            segment_analysis->flag_RuleTextFile=1;
            segment_analysis->filename_RuleTextFile=argv[++i];
        }
        else if(strcmp(argv[i], "-output_dir") == 0 && (i+1)<argc){
            segment_analysis->flag_outputDir=1;
            segment_analysis->name_outputDir=argv[++i];
        }
        else if(strcmp(argv[i], "-inRuleCorr") == 0 && (i+1)<argc){
            segment_analysis->flag_RuleCorr=1;
            segment_analysis->filename_RuleCorr=argv[++i];
        }
        else if(strcmp(argv[i], "-inLCRuleTxt") == 0 && (i+1)<argc){
            segment_analysis->flag_LCRuleTextFile=1;
            segment_analysis->filename_LCRuleTextFile=argv[++i];
        }
        else if(strcmp(argv[i], "-inGMPre") == 0 && (i+1)<argc){
            segment_analysis->flag_GMPre=1;
            segment_analysis->filename_GMPre=argv[++i];
        }
        else if(strcmp(argv[i], "-inCSFPre") == 0 && (i+1)<argc){
            segment_analysis->flag_CSFPre=1;
            segment_analysis->filename_CSFPre=argv[++i];
        }
        else if(strcmp(argv[i], "-TestLaplace") == 0 && (i+1)<argc){
            segment_analysis->flag_test=1;
            segment_analysis->filename_test=argv[++i];
        }
        else if(strcmp(argv[i], "-inData") == 0 && (i+1)<argc){

            int numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if(numbmodal<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+numbmodal)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbmodal; m++){
                    segment_analysis->filename_Data.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inData are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_data=1;

        }

        else if(strcmp(argv[i], "-inModa") == 0 && (i+1)<argc){

            int numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if(numbmodal<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+numbmodal)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbmodal; m++){
                    segment_analysis->vecModalities.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inModa are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            if (numbmodal==1 && segment_analysis->vecModalities[0]==3){
                segment_analysis->flag_FLAIRonly = 1;
            }
            segment_analysis->flag_GivenModa=1;

        }

        else if(strcmp(argv[i], "-inGaussian") == 0 && (i+1)<argc){

            int numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if((i+2*numbmodal)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbmodal; m++){
                    segment_analysis->filename_ImagesGaussian.push_back(argv[++i]);
                    segment_analysis->Gaussian_Zscore.push_back(atof(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inModa are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Gaussian=1;

        }
        else if(strcmp(argv[i], "-inOutVec") == 0 && (i+1)<argc){

            int numbOut=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbOut)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbOut; m++){
                    segment_analysis->filename_vectorOut.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inOutVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }
        else if(strcmp(argv[i], "-inOutVec") == 0 && (i+1)<argc){

            int numbOut=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbOut)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbOut; m++){
                    segment_analysis->filename_vectorOut.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inOutVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }
        else if(strcmp(argv[i], "-inLobes") == 0 && (i+1)<argc){

            int numbLobes=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbLobes)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbLobes; m++){
                    segment_analysis->filename_inLobes.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inVentrVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_inLobes=1;

        }
        else if(strcmp(argv[i], "-inMaskVec") == 0 && (i+1)<argc){

            int numbMask=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbMask)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbMask; m++){
                    segment_analysis->filename_vectorMask.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inMaskVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }

        else if(strcmp(argv[i], "-inVentrVec") == 0 && (i+1)<argc){

            int numbMask=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbMask)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbMask; m++){
                    segment_analysis->filename_vectorVentr.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inMaskVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }

        else if(strcmp(argv[i], "-inNamesVec") == 0 && (i+1)<argc){

            int numbNames=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbNames)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbNames; m++){
                    segment_analysis->filename_vectorNames.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inNamesVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }


        else if(strcmp(argv[i], "-inRefVec") == 0 && (i+1)<argc){

            int numbRef=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbRef)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbRef; m++){
                    segment_analysis->filename_vectorRef.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inRefVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }


        else if(strcmp(argv[i], "-inLesVec") == 0 && (i+1)<argc){

            int numbLes=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbLes)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbLes; m++){
                    segment_analysis->filename_vectorLes.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inLesVec are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }

        else if(strcmp(argv[i], "-inVecThresh") == 0 && (i+4)<argc){
            segment_analysis->vecThresh.push_back(atof(argv[++i]));
            segment_analysis->vecThresh.push_back(atof(argv[++i]));
            segment_analysis->vecThresh.push_back(atof(argv[++i]));
            segment_analysis->flag_Evaluation=1;

        }
        else if(strcmp(argv[i], "-Dot") == 0 && (i+1)<argc){
            segment_analysis->filename_dot1=(argv[++i]);
            segment_analysis->flag_Dot=1;

        }

        else if(strcmp(argv[i], "-inVecLesRedNumbers") == 0 && (i+1)<argc){

            int numbRed=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbRed)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int m=0; m<numbRed; m++){
                    segment_analysis->vecLesionReduction.push_back(atof(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inVecLesRedNumbers are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_Evaluation=1;

        }

        else if(strcmp(argv[i], "-inVecLesRed") == 0 && (i+4)<argc){
            segment_analysis->vecLesionReductionGeneral.push_back(atof(argv[++i]));
            segment_analysis->vecLesionReductionGeneral.push_back(atof(argv[++i]));
            segment_analysis->vecLesionReductionGeneral.push_back(atof(argv[++i]));
            segment_analysis->vecLesionReductionGeneral.push_back(atof(argv[++i]));
            segment_analysis->flag_LesionReduction=1;

        }

        else if(strcmp(argv[i], "-inRuleMahal") == 0 && (i+1)<argc){

            int numbRules=atoi(argv[++i]); // atoi convert string to integer
            if(numbRules<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+numbRules)<argc){// Read names of Modalities according to code 1 2 3 (Maybe transformed into string afterwards)
                for(int r=0; r<numbRules; r++){
                    segment_analysis->vecRuleMahal.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inRuleMahal are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_GivenMahalRule=1;

        }

        else if(strcmp(argv[i], "-inBF") == 0 && (i+1)<argc){

            int numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if(numbmodal<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+numbmodal)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbmodal; m++){
                    segment_analysis->filename_BF.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -inBF are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_BF=1;

        }

        else if(strcmp(argv[i], "-inOutliers") == 0 && (i+1)<argc){

            int numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if(numbmodal<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+numbmodal)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbmodal; m++){
                    segment_analysis->filename_Outliers.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_outliers=1;

        }

        else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc){

            int numbmodal=atoi(argv[++i]); // atoi convert string to integer
            if(numbmodal<1){
                cout<<"Number of modalities has to be bigger than 1";
                return 0;
            }
            if((i+numbmodal)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbmodal; m++){
                    segment_analysis->filename_In.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_in=1;

        }
        else if(strcmp(argv[i], "-meanImages") == 0 && (i+1)<argc){

            int numbimages=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbimages)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbimages; m++){
                    segment_analysis->filename_ImagesMean.push_back(argv[++i]);
                }
                cout<<"Size of MeanImages is "<<segment_analysis->filename_ImagesMean.size();
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_mean=1;

        }

        else if(strcmp(argv[i], "-LAvec") == 0 && (i+1)<argc){
            segment_analysis->vecLeavesToAdd.clear();
            int numbleaves=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbleaves)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbleaves; m++){
                    segment_analysis->vecLeavesToAdd.push_back(atoi(argv[++i]));
                }
                cout<<"Size of LAvec is "<<segment_analysis->vecLeavesToAdd.size();
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_TO=1;

        }

        else if(strcmp(argv[i], "-vecConnect") == 0 && (i+1)<argc){

            int numbimages=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbimages)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbimages; m++){
                    segment_analysis->filename_vectorRef.push_back(argv[++i]);
                }
                cout<<"Size of RefImages is "<<segment_analysis->filename_vectorRef.size();
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_CorrectConnect=1;

        }
        else if(strcmp(argv[i], "-vecCorrect") == 0 && (i+1)<argc){

            int numbimages=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbimages)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbimages; m++){
                    segment_analysis->filename_vectorLes.push_back(argv[++i]);
                }
                cout<<"Size of RefImages is "<<segment_analysis->filename_vectorLes.size();
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_CorrectConnect=1;

        }

        else if(strcmp(argv[i], "-matchRef") == 0 && (i+1)<argc){

            int numbimages=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbimages)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbimages; m++){
                    segment_analysis->filename_RefImages.push_back(argv[++i]);
                }
                cout<<"Size of RefImages is "<<segment_analysis->filename_RefImages.size();
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_match=1;

        }
        else if(strcmp(argv[i], "-matchFloat") == 0 && (i+1)<argc){

            int numbimages=atoi(argv[++i]); // atoi convert string to integer
            if((i+numbimages)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<numbimages; m++){
                    segment_analysis->filename_FloatImages.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_match=1;

        }
        else if(strcmp(argv[i], "-orderFit") == 0 && (i+1)<argc){
            segment_analysis->orderMatch=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-Neigh") == 0 && (i+1)<argc){
            segment_analysis->Neigh=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-classGaussian") == 0 && (i+1)<argc){
            segment_analysis->index_Gaussian=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-outFit") == 0 && (i+1)<argc){
            segment_analysis->flag_match=1;
            segment_analysis->filename_saveMatch=argv[++i];
        }
        else if(strcmp(argv[i], "-thresh_cor") == 0 && (i+1)<argc){
            segment_analysis->thresh_cortex=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-outMean") == 0 && (i+1)<argc){
            segment_analysis->flag_mean=1;
            segment_analysis->filename_saveMean=argv[++i];
        }
        else if(strcmp(argv[i], "-inLF_roi") == 0 && (i+2)<argc){
            segment_analysis->flag_inROIFiles=1;
            segment_analysis->filename_inROIFiles.push_back(argv[++i]);
            segment_analysis->filename_inROIFiles.push_back(argv[++i]);
        }
        else if(strcmp(argv[i], "-inLF_ex") == 0 && (i+2)<argc){
            segment_analysis->flag_inExFiles=1;
            segment_analysis->filename_inExFiles.push_back(argv[++i]);
            segment_analysis->filename_inExFiles.push_back(argv[++i]);
        }

        //      else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc){
        //          segment_analysis->flag_in=1;
        //          segment_analysis->filename_In=argv[++i];
        //        }
        else if(strcmp(argv[i], "-inLes") == 0 && (i+1)<argc){
            segment_analysis->flag_in=1;
            segment_analysis->flag_Les=1;
            segment_analysis->flag_inLes=1;
            segment_analysis->filename_InLes=argv[++i];
        }
        else if (strcmp(argv[i], "-inToAnalyse") ==0 && (i+1)<argc){
            segment_analysis->flag_inToAnalyse=1;
            segment_analysis->filename_inToAnalyse=argv[++i];
        }
        else if(strcmp(argv[i], "-inLesCorr") == 0 && (i+1)<argc){
            segment_analysis->flag_inLesCorr=1;
            segment_analysis->flag_Les=1;
            segment_analysis->filename_inLesCorr=argv[++i];
        }
        else if(strcmp(argv[i], "-inReclassif") == 0 && (i+1)<argc){
            segment_analysis->flag_inFPTPReclassif=1;
            segment_analysis->filename_inFPTPReclassif=argv[++i];
        }

        else if(strcmp(argv[i], "-inConnect") == 0 && (i+1)<argc){
            segment_analysis->flag_inConnect=1;
            segment_analysis->filename_inConnect=argv[++i];
            segment_analysis->flag_ConnectRefAnalysis=1;
        }
        else if(strcmp(argv[i], "-inOutText") == 0 && (i+1)<argc){
            segment_analysis->flag_inOutText=1;
            segment_analysis->filename_inOutText=argv[++i];
        }
        else if(strcmp(argv[i], "-inWMSeg") == 0 && (i+1)<argc){
            segment_analysis->flag_inWMSeg=1;
            segment_analysis->filename_inWMSeg=argv[++i];
        }
        else if(strcmp(argv[i], "-inMahal") == 0 && (i+1)<argc){
            segment_analysis->flag_inMahal=1;
            segment_analysis->filename_inMahal=argv[++i];
        }
        else if (strcmp(argv[i], "-inHemi") == 0 && (i+1)<argc){
            segment_analysis->filename_Hemi=argv[++i];
            segment_analysis->flag_Hemi=1;
        }
        else if(strcmp(argv[i], "-inCST") == 0 && (i+1)<argc){
            segment_analysis->flag_CST=1;
            segment_analysis->filename_inCST=argv[++i];
        }
        else if(strcmp(argv[i], "-inOptionText") == 0 && (i+1)<argc){
            segment_analysis->flag_inOptionText=1;
            segment_analysis->inOptionText=argv[++i];
        }
        else if(strcmp(argv[i], "-inPriorsICSF") == 0 && (i+1)<argc){
            segment_analysis->flag_inPriorsICSF=1;
            segment_analysis->filename_inPriorsICSF=argv[++i];
        }
        else if(strcmp(argv[i], "-inPriorsECSF") == 0 && (i+1)<argc){
            segment_analysis->flag_inPriorsECSF=1;
            segment_analysis->filename_inPriorsECSF=argv[++i];
        }
        else if(strcmp(argv[i], "-inPriorsOut") == 0 && (i+1)<argc){
            segment_analysis->flag_inPriorsOut=1;
            segment_analysis->filename_inPriorsOut=argv[++i];
        }
        else if(strcmp(argv[i], "-inPriorsDGM") == 0 && (i+1)<argc){
            segment_analysis->flag_inPriorsDGM=1;
            segment_analysis->filename_inPriorsDGM=argv[++i];
        }
        else if(strcmp(argv[i], "-inPriorsCGM") == 0 && (i+1)<argc){
            segment_analysis->flag_inPriorsCGM=1;
            segment_analysis->filename_inPriorsCGM=argv[++i];
        }
        else if(strcmp(argv[i], "-inPriorsWM") == 0 && (i+1)<argc){
            segment_analysis->flag_inPriorsWM=1;
            segment_analysis->filename_inPriorsWM=argv[++i];
        }
        else if(strcmp(argv[i], "-juxtaCorr") == 0 && (i+1)<argc){
            segment_analysis->flag_juxtaCorrection=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-oldLes") == 0 && (i+1)<argc){
            segment_analysis->flag_oldLesion=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-SP") == 0 && (i+1)<argc){
            segment_analysis->flag_SPCorrection=atoi(argv[++i]);
            cout<<"flag_SPCorrection is "<<segment_analysis->flag_SPCorrection<<endl;
        }
        else if(strcmp(argv[i], "-minSize") == 0 && (i+1)<argc){
            segment_analysis->flag_minSize = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-Neigh") == 0 && (i+1)<argc){
            segment_analysis->Neigh = atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-SegGIF") == 0 && (i+1)<argc){
            segment_analysis->flag_segGIF=atoi(argv[++i]);
            if (segment_analysis->flag_segGIF) {
                segment_analysis->IndexOut=0;
                segment_analysis->IndexCSF=1;
                segment_analysis->IndexCGM=2;
                segment_analysis->IndexWM=3;
                segment_analysis->IndexDGM=4;
                segment_analysis->IndexBrainstem=5;
                segment_analysis->IndexGM=-1;
                segment_analysis->IndexICSF=-1;
                segment_analysis->IndexECSF=-1;
            }
        }
        else if(strcmp(argv[i], "-inVentricleSeg") == 0 && (i+1)<argc){
            segment_analysis->flag_inVentricleSeg=1;
            segment_analysis->filename_inVentricleSeg=argv[++i];
        }
        else if(strcmp(argv[i], "-inDataComp") == 0 && (i+1)<argc){
            segment_analysis->flag_inDataComp=1;
            segment_analysis->filename_inDataComp=argv[++i];
        }
        else if(strcmp(argv[i], "-refLes") == 0 && (i+1)<argc){
            segment_analysis->flag_refLes=1;
            segment_analysis->filename_Ref=argv[++i];
            segment_analysis->flag_Les=1;
            segment_analysis->flag_GT=1;
        }
        else if(strcmp(argv[i], "-refConnect") == 0 && (i+1)<argc){
            segment_analysis->flag_refLesConnect=1;
            segment_analysis->filename_RefConnect=argv[++i];
            segment_analysis->flag_Les=1;
            segment_analysis->flag_GT=0;
        }
        else if(strcmp(argv[i], "-HS") == 0 && (i+1)<argc){
            segment_analysis->flag_HS=1;
            segment_analysis->filename_Out=argv[++i];
        }
        else if(strcmp(argv[i], "-OutFile") == 0 && (i+1)<argc){
            segment_analysis->flag_fileOut = 1;
            segment_analysis->filename_Out=argv[++i];
        }
        else if(strcmp(argv[i], "-IO") == 0 && (i+1)<argc){
            segment_analysis->flag_IO=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-TO") == 0 && (i+1)<argc){
            segment_analysis->flag_TO=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-WMMaps") == 0 && (i+1)<argc){
            segment_analysis->flag_WMMaps=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-Vesselness") == 0 && (i+1)<argc){
            segment_analysis->flag_Vesselness=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-LS") == 0 && (i+1)<argc){
            segment_analysis->flag_LS=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-Layers") == 0 && (i+1)<argc){
            segment_analysis->flag_layers=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-create_layers") == 0 && (i+1)<argc){
            segment_analysis->flag_layercreation=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-Euc") == 0 && (i+1)<argc){
            segment_analysis->flag_Euc=atoi(argv[++i]);
        }

        else if(strcmp(argv[i], "-Saving") == 0 && (i+1)<argc){
            segment_analysis->flag_Saving=atoi(argv[++i]);
        }

        else if(strcmp(argv[i], "-checkInliers") == 0 && (i+1)<argc){
            segment_analysis->flag_checkInliers=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-Distance") == 0 && (i+1)<argc){
            segment_analysis->flag_Distance=1;
            segment_analysis->filename_Dist=argv[++i];
        }
        else if(strcmp(argv[i], "-ExploRadius") == 0 && (i+1)<argc){
            segment_analysis->ExploRadius=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-ThresholdSP") == 0 && (i+1)<argc){
            segment_analysis->ThresholdSP=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-LevelCorrection") == 0 && (i+1)<argc){
            segment_analysis->flag_correctionLevel=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-CorrIV") == 0 && (i+1)<argc){
            segment_analysis->flag_CorrectIV=(bool)atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-SegType") == 0 && (i+1)<argc){
            segment_analysis->flag_segType=atoi(argv[++i]);
            cout<<"SegType chosen is "<<segment_analysis->flag_segType<<endl;
        }
        else if(strcmp(argv[i], "-Edginess") == 0 && (i+1)<argc){
            segment_analysis->EdginessThreshold=atof(argv[++i]);
            cout<<"Edginess chosen is "<<segment_analysis->EdginessThreshold<<endl;
        }
        else if(strcmp(argv[i], "-WeightedSeg") == 0 && (i+3)<argc){
            segment_analysis->flag_segWeighted=atoi(argv[++i]);
            float weightSeg=atof(argv[++i]);
            int weightCompClass=atoi(argv[++i]);
            cout<<"Weightseg is "<<segment_analysis->flag_segWeighted<<" ******* "<<endl;
            if (weightSeg>0) {
                segment_analysis->weightThreshold=weightSeg;
                segment_analysis->weightCompClass=weightCompClass;
                cout<<"WeightSeg chosen is "<<segment_analysis->weightThreshold<<endl;
            }
            else{
                segment_analysis->flag_segWeighted=0;
            }

        }

        else if(strcmp(argv[i], "-LesionRuleType") == 0 && (i+1)<argc){
            segment_analysis->LesionRuleType=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-LesionUniformRuleType") == 0 && (i+1)<argc){
            segment_analysis->LesionUniformRuleType=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-out_txt") == 0 && (i+1)<argc){
            segment_analysis->flag_outTxt=1;
            segment_analysis->filename_OutTxt=argv[++i];
        }
        else if(strcmp(argv[i], "-WMCard") == 0 && (i+1)<argc){
            segment_analysis->flag_WMCard=atoi(argv[++i]);
            //segment_analysis->filename_OutTxt=argv[++i];
        }
        else if(strcmp(argv[i], "-compJoint") == 0 && (i+1)<argc){
            segment_analysis->flag_compJoint=1;
            segment_analysis->filename_compJoint=argv[++i];
        }
        else if(strcmp(argv[i], "-inQuadrant") == 0 && (i+1)<argc){
            segment_analysis->flag_inQuadrant=1;
            segment_analysis->filename_inQuadrant=argv[++i];
        }
        else if(strcmp(argv[i], "-LaplaceNormSol") == 0 && (i+1)<argc){
            segment_analysis->flag_LaplaceSolImage=1;
            segment_analysis->flag_LapIn=1;
            segment_analysis->filename_LaplaceSolImage=argv[++i];
        }
        else if(strcmp(argv[i],"-LapAnalysis")==0 && (i+1)<argc){
            segment_analysis->flag_LapAnalysis=(bool)atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-numbCoded")==0 && (i+1)<argc){
            segment_analysis->flag_CodedOutliers=1;
            segment_analysis->numbCoded=atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-outWM")==0 && (i+1)<argc){
            segment_analysis->flag_outWMI=(bool)atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-getCoM")==0 && (i+1)<argc){
            segment_analysis->flag_getCoM=(bool)atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-outConnect")==0 && (i+1)<argc){
            segment_analysis->flag_outConnect=(bool)atoi(argv[++i]);
        }
        else if(strcmp(argv[i],"-LesWMI")==0 && (i+1)<argc){
            segment_analysis->flag_LesWMI=atof(argv[++i]);
        }
        else if(strcmp(argv[i],"-thresh_RG")==0 && (i+1)<argc){
            segment_analysis->thresh_RG=atof(argv[++i]);
            segment_analysis->flag_RG=1;
        }
        else if(strcmp(argv[i],"-thresh_quant")==0 && (i+1)<argc){
            segment_analysis->thresh_quant=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-ParcellationIn") == 0 && (i+1)<argc){
            segment_analysis->flag_Parc=1;
            segment_analysis->filename_Parc=argv[++i];
        }
        else if(strcmp(argv[i], "-AuthorisedIn") == 0 && (i+1)<argc){
            segment_analysis->flag_inAuthorised=1;
            segment_analysis->filename_inAuthorised=argv[++i];
        }
        else if(strcmp(argv[i], "-typeVentrSeg") == 0 && (i+1)<argc){
            segment_analysis->flag_VentrSeg=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-inArtefact") == 0 && (i+1)<argc){
            segment_analysis->flag_inArtefact=1;
            segment_analysis->filename_Artefact=argv[++i];
        }
        else if(strcmp(argv[i], "-SummarisedIn") == 0 && (i+1)<argc){
            segment_analysis->flag_inSum=1;
            segment_analysis->filename_inSum=argv[++i];
        }
        else if(strcmp(argv[i], "-BFGenText") == 0 && (i+1)<argc){
            segment_analysis->flag_BFGen=1;
            segment_analysis->flag_BFText=1;
            segment_analysis->filename_BFText=argv[++i];
        }
        else if(strcmp(argv[i], "-BFGen") == 0 && (i+3)<argc){
            segment_analysis->flag_BFGen=1;
            segment_analysis->BFParam.push_back(atof(argv[++i]));
            segment_analysis->BFParam.push_back(atof(argv[++i]));
            segment_analysis->BFParam.push_back(atof(argv[++i]));
        }
        else if(strcmp(argv[i], "-AffGen") == 0 && (i+1)<argc){
            segment_analysis->flag_AffGen=1;
            segment_analysis->AffParam.push_back(atof(argv[++i]));
            segment_analysis->AffParam.push_back(atof(argv[++i]));

        }
        else if(strcmp(argv[i],"-LapMap") ==0 && (i+1)<argc){
            segment_analysis->flag_MapLap=atoi(argv[++i]);
            cout <<"MapLap is "<< segment_analysis->flag_MapLap<<endl;
        }
        else if(strcmp(argv[i], "-LapIO") == 0 && (i+2)<argc){
            segment_analysis->flag_vecIOLap=1;
            segment_analysis->flag_LapAnalysis=1;
            segment_analysis->filenames_IOLap.push_back(argv[++i]);
            segment_analysis->filenames_IOLap.push_back(argv[++i]);
        }
        else if(strcmp(argv[i], "-nameLap") == 0 && (i+1)<argc){
            segment_analysis->flag_nameLap=1;
            segment_analysis->nameLap=argv[++i];
        }
        else if(strcmp(argv[i], "-numbLaminae") == 0 && (i+1)<argc){
            segment_analysis->numbLaminae.clear();

            int nlsize=atoi(argv[++i]); // atoi convert string to integer

            if((i+nlsize)<argc){
                for(int m=0; m<nlsize; m++){
                    segment_analysis->numbLaminae.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -numbLaminae are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
        }
        else if(strcmp(argv[i], "-compMult") == 0 && (i+1)<argc){

            segment_analysis->numbMult=atoi(argv[++i]); // atoi convert string to integer

            if((i+segment_analysis->numbMult)<argc){// Read filenames of the input images. There must be as many input filenames
                // as modalities defined in segment_param->numbmodal
                //              segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int m=0; m<segment_analysis->numbMult; m++){
                    segment_analysis->filename_compMult.push_back(argv[++i]);
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -input are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
            segment_analysis->flag_compMult=1;

        }
        else if(strcmp(argv[i], "-ElementsAssociateIn") == 0 && (i+1)<argc){

            int numbElementsIn=atoi(argv[++i]); // atoi convert string to integer

            if((i+numbElementsIn)<argc){
                for(int m=0; m<numbElementsIn; m++){
                    segment_analysis->vecElementIn.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -ElementsAssociateIn are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
        }
        else if(strcmp(argv[i], "-ElementsAssociateOut") == 0 && (i+1)<argc){

            int numbElementsOut=atoi(argv[++i]); // atoi convert string to integer

            if((i+numbElementsOut)<argc){
                for(int m=0; m<numbElementsOut; m++){
                    segment_analysis->vecElementOut.push_back(atoi(argv[++i]));
                }
            }
            else{
                fprintf(stderr,"Err:\tParameter -ElementsAssociateOut are incomplete\n");
                Usage(argv[0]);
                return EXIT_SUCCESS;
            }
        }

        else if(strcmp(argv[i], "-Package") == 0 && (i+1)<argc){
            TypePackage=atoi(argv[++i]);
            if (TypePackage>0) {
                cout<<"Going before package"<<endl;
                segment_analysis->flag_Analysis=1;
                segment_analysis->flag_DCh=1;
                segment_analysis->flag_Acc=1;
                segment_analysis->flag_Corr=1;
                segment_analysis->flag_DE=1;
                segment_analysis->flag_FN=1;
                segment_analysis->flag_FP=1;
                segment_analysis->flag_FPR=1;
                segment_analysis->flag_GT=1;
                segment_analysis->flag_TP=1;
                segment_analysis->flag_VD=1;
                segment_analysis->flag_OE=1;
                segment_analysis->flag_TLLHard=1;
                segment_analysis->flag_Spec=1;
                segment_analysis->flag_Sens=1;
                segment_analysis->flag_PSI=1;
                segment_analysis->flag_AvDist=1;
                segment_analysis->flag_correct=1;
                segment_analysis->flag_connect=1;
                segment_analysis->flag_IO=1;
                cout<<"Going after package"<<endl;
            }
            else if (TypePackage==2) {
                segment_analysis->flag_ConnectRefAnalysis=1;
            }
            else if (TypePackage==0){
                segment_analysis->flag_connect=0;
                segment_analysis->flag_correct=0;
            }
        }

        else if(strcmp(argv[i], "-DCh") == 0){
            segment_analysis->flag_DCh=1;
        }
        else if(strcmp(argv[i], "-DCs") == 0){
            segment_analysis->flag_DCs=1;
        }
        else if(strcmp(argv[i], "-NW") == 0){
            segment_analysis->flag_NW=1;
        }
        else if(strcmp(argv[i], "-GTSave") == 0){
            segment_analysis->flag_GT=1;
            segment_analysis->flag_GTSave=1;
        }
        else if(strcmp(argv[i], "-TP") == 0){
            segment_analysis->flag_TP=1;
        }
        else if(strcmp(argv[i], "-TPR") == 0){
            segment_analysis->flag_TPR=1;
        }
        else if(strcmp(argv[i], "-FPR") == 0){
            segment_analysis->flag_FPR=1;
        }
        else if(strcmp(argv[i], "-VD") == 0){
            segment_analysis->flag_VD=1;
        }
        else if(strcmp(argv[i], "-DE") == 0){
            segment_analysis->flag_DE=1;
        }
        else if(strcmp(argv[i], "-OE") == 0){
            segment_analysis->flag_OE=1;
        }
        else if(strcmp(argv[i], "-FP") == 0){
            segment_analysis->flag_FP=1;
        }
        else if(strcmp(argv[i], "-TN") == 0){
            segment_analysis->flag_TN=1;
        }
        else if(strcmp(argv[i], "-FN") == 0){
            segment_analysis->flag_FN=1;
        }
        else if(strcmp(argv[i], "-Sens") == 0){
            segment_analysis->flag_Sens=1;
        }
        else if(strcmp(argv[i], "-Spec") == 0){
            segment_analysis->flag_Spec=1;
        }
        else if(strcmp(argv[i], "-PSI") == 0){
            segment_analysis->flag_PSI=1;
        }
        else if(strcmp(argv[i], "-Acc") == 0){
            segment_analysis->flag_Acc=1;
        }
        else if(strcmp(argv[i], "-Corr") == 0){
            segment_analysis->flag_Corr=1;
        }
        else if(strcmp(argv[i], "-TLLSoft") == 0){
            segment_analysis->flag_TLLSoft=1;
        }
        else if(strcmp(argv[i], "-TLLHard") == 0){
            segment_analysis->flag_TLLHard=1;
        }
        else if(strcmp(argv[i], "-NumbActive") == 0){
            segment_analysis->flag_numbActive=1;
        }
        else if(strcmp(argv[i], "-mask") == 0 && (i+1)<argc){
            segment_analysis->filename_mask=argv[++i];
            segment_analysis->flag_mask=1;
        }
        else if(strcmp(argv[i], "-inLacMask") == 0 && (i+1)<argc){
            segment_analysis->filename_inLacMask=argv[++i];
            segment_analysis->flag_mask=1;
        }
        else if(strcmp(argv[i], "-maskMatch") == 0 && (i+1)<argc){
            segment_analysis->filename_maskMatch=argv[++i];
            segment_analysis->flag_maskMatch=1;
        }
        else if(strcmp(argv[i], "-KLD") == 0){
            segment_analysis->flag_KLD=1;
        }
        else if(strcmp(argv[i], "-SSD") == 0){
            segment_analysis->flag_SSD=1;
        }
        else if(strcmp(argv[i], "-connect") == 0 ){
            segment_analysis->flag_connect=1;
        }
        else if(strcmp(argv[i], "-GMAccept") == 0 ){
            segment_analysis->flag_AcceptedGM=1;
        }
        else if(strcmp(argv[i], "-correct") == 0 ){
            segment_analysis->flag_correct=1;
            segment_analysis->flag_connect=1;
        }
        else if(strcmp(argv[i], "-PrintConnect") == 0 ){
            segment_analysis->flag_PrintConnect=1;
        }
        else if(strcmp(argv[i], "-GetConnect") == 0 ){
            segment_analysis->flag_GetConnect=1;
        }
        else if(strcmp(argv[i], "-AvDist") == 0){
            segment_analysis->flag_AvDist=1;
        }
        else if(strcmp(argv[i], "-PVLes") == 0){
            segment_analysis->flag_PV=1;
        }
        else if(strcmp(argv[i], "-Simple") == 0){
            segment_analysis->flag_simple=1;
        }
        else if(strcmp(argv[i], "-infarcts") == 0){
            segment_analysis->flag_infarcts=1;
        }
        else if(strcmp(argv[i], "-Secondary") == 0 && (i+1)<argc){
            segment_analysis->flag_Secondary=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-WMCheck") == 0){
            segment_analysis->flag_WMCheck=1;
        }
        else if(strcmp(argv[i], "-CSFCheck") == 0){
            segment_analysis->flag_CSFCheck=1;
        }
        else if(strcmp(argv[i], "-CorExt") == 0){
            segment_analysis->flag_CorExt=1;
            segment_analysis->val_ext == atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-ITCheck") == 0){
            segment_analysis->flag_ITCheck=1;
        }
        else if(strcmp(argv[i], "-verbose_level") == 0){
            segment_analysis->verbose_level=GetLevelMap(atoi(argv[++i]));
        }
        else if(strcmp(argv[i], "-VesselSeed") == 0 && (i+1)<argc){
            segment_analysis->flag_VesselSeed=1;
        }
        else if(strcmp(argv[i], "-CorrCerr") == 0 && (i+1)<argc){
            segment_analysis->flag_correctCerr=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-JCMS") == 0 && (i+1)<argc){
            segment_analysis->flag_acceptJCMS=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-GIFPriors") == 0 && (i+1)<argc){
            segment_analysis->filename_GIFPriors=argv[++i];
            segment_analysis->flag_GIFPriors=1;
            segment_analysis->flag_juxtaCorrection=0;
        }



        else{
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return EXIT_SUCCESS;
        }
    }
    
    if (segment_analysis->flag_refLesConnect) {
        segment_analysis->flag_GT=0;
    }
    if (segment_analysis->flag_inPriorsDGM==0) {
        segment_analysis->flag_juxtaCorrection=0;
    }

    if(segment_analysis->flag_VentrSeg==1 && ! segment_analysis->flag_Parc){
        segment_analysis->flag_VentrSeg=0;
    }
    if(segment_analysis->flag_prion){
        segment_analysis->flag_Evaluation=0;
        if (! segment_analysis->flag_Parc && !segment_analysis->flag_maskMatch){
            segment_analysis->flag_prion=0;
        }
    }
    if(segment_analysis->flag_Dot){
        if(!segment_analysis->flag_inToAnalyse){
            segment_analysis->flag_inToAnalyse=1;
            strcpy(segment_analysis->filename_inToAnalyse,segment_analysis->filename_dot1);
        }
    }
    nifti_image * SegToAnalyse=NULL;
    nifti_image * SegToAnalyseNIV=NULL;
    nifti_image * ImageLesionCorr=NULL;
    nifti_image * BinarySeg=NULL;
    nifti_image * ConnectedGTImage=NULL;
    nifti_image * ConnectLabelLesions=NULL;
    nifti_image * CorrectionJuxta=NULL;
    float * ICSF=NULL;
    float * DGM=NULL;
    int * OrderedLabelsGT=NULL;
    int * OrderedLabels=NULL;
    int * OrderedLabelsCorr=NULL;
    nifti_image * VentricleSeg=NULL;
    nifti_image * DGMSeg=NULL;
    nifti_image * PriorsICSF=NULL;
    nifti_image * PriorsDGM=NULL;
    bool * DGMBool=NULL;
    bool * ECSFBool=NULL;
    bool * CGMBool=NULL;
    bool * VentricleBool=NULL;
    bool * SPRegion=NULL;
    TreeEM * TreeToAnalyse=NULL;
    nifti_image * DataCompImage=NULL;
    nifti_image * WMSeg=NULL;
    nifti_image * SummarisedSeg1=NULL;
    float * ParamWM=NULL;
    float * InvertedCovarianceWM=NULL;
    int numel=0;
    int numbmodal=0;
    string FilenamePA;
    string FilenamePA_b;
    string FilenamePA_e;
    vector<int> Modalities;
    if (segment_analysis->vecModalities.size()>0) {
        numbmodal=segment_analysis->vecModalities.size();
        for (int m=0; m<numbmodal; m++) {
            Modalities.push_back(segment_analysis->vecModalities[m]);
        }
    }

    log_level_t GLOBAL_LEVEL = LOG_WARNING;
// log<LOG_DEBUG>(L"TEST %3% %2% %1%") % 5 % 10 % L"privet";
    log_verb(LOG_DEBUG, "Degug testing 5", GLOBAL_LEVEL);
    log_verb(LOG_NOTHING, "Normally nothing", GLOBAL_LEVEL);
    log_verb(LOG_NOTHING, "Normally nothing warning");
    log_verb(LOG_WARNING, "Warning debug 5",GLOBAL_LEVEL);

    vector<TreeEM *> VectorLeaves;
    int numbleavesTree=0;
    //    TreeBuilding if
    
    //    if (segment_analysis->flag_inPriorsDGM) {
    //        nifti_image * PriorsDGM=ReadFromFilename(segment_analysis->filename_inPriorsDGM);
    //        cout<<"Able to see PriorsDGM..."<<endl;
    //    }

    if((segment_analysis->flag_inToAnalyse || segment_analysis->flag_inLes || segment_analysis->flag_inLesCorr) && segment_analysis->flag_minSize>0 && !segment_analysis->flag_TextFile){
        nifti_image * ToCorrect = NULL;

        if (segment_analysis->flag_inToAnalyse){
            ToCorrect=ReadFromFilename(segment_analysis->filename_inToAnalyse);
            FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        }
        else if (segment_analysis->flag_inLes){
            ToCorrect=ReadFromFilename(segment_analysis->filename_inLes);
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLes);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        }
        else{
            ToCorrect=ReadFromFilename(segment_analysis->filename_inLesCorr);
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        }

        if(segment_analysis->Neigh != 6 && segment_analysis->Neigh != 18 && segment_analysis->Neigh != 26){
            segment_analysis->Neigh = 6;
        }
        float * DataToCorr = static_cast<float*>(ToCorrect->data);
        stringstream sG;
        sG << segment_analysis->flag_minSize;
        string cG=sG.str();
        string FilenameSave = FilenamePA_b + "CorrectSize"+cG+"_"+FilenamePA_e + ".nii.gz";
        int numel = ToCorrect->nx*ToCorrect->ny*ToCorrect->nz;
        float VolVox = ToCorrect->pixdim[1]*ToCorrect->pixdim[2]*ToCorrect->pixdim[3];
        int * Dim=new int[3];
        int * Shift=new int[3];
        Shift[0]=1;
        float * PixDim=new float[3];
        for (int d=0; d<3; d++) {
            Dim[d]=ToCorrect->dim[d+1];
            if (d>0) {
                Shift[d]=Dim[d-1]*Shift[d-1];
            }
            PixDim[d]=ToCorrect->pixdim[d+1];
        }
        int * ComponentLabel = ComponentLabeling(DataToCorr,segment_analysis->Neigh,Dim,Shift);
        int * OrderedLabels=OrderedVolumeLabel(ComponentLabel, segment_analysis->flag_minSize,numel,VolVox);
        float * CorrectData = CorrectByMask(DataToCorr, OrderedLabels, numel, 0);
        int numLabels = GetMaxLabel(ComponentLabel, numel);
        int numCorrect = GetMaxLabel(OrderedLabels, numel);
        cout << "Initial number of labels is "<<numLabels<< " and corrected is "<<numCorrect<<endl;
        int NonZero = CountNonZero(OrderedLabels, numel);
        int NonZero_init = CountNonZero(ComponentLabel,numel);
        cout << "Correction of "<< NonZero_init - NonZero << endl;
        nifti_image * CorrectNii = CreateNiiFromArray(CorrectData, ToCorrect, numel);
        nifti_set_filenames(CorrectNii, FilenameSave.c_str(), 0, 0);
        nifti_image_write(CorrectNii);
        nifti_image_free(ToCorrect);
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_inToAnalyse && segment_analysis->flag_layers){
        nifti_image * LesionToAssess = ReadFromFilename(segment_analysis->filename_inToAnalyse);
        nifti_image * Layers = ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSave=FilenamePA_b+"LayersQualified_"+FilenamePA_e+".nii.gz";
        int* Labels = ComponentLabeling(LesionToAssess, 6);
        int numel = LesionToAssess->nvox;
        float * LayerData = static_cast<float*>(Layers->data);
        float * QualifiedLayer=new float[numel];
        int numb_lesion = GetMaxLabel(Labels, numel);
        for (int i=0;i<numel;i++){
            QualifiedLayer[i]=0;
        }
        for(int l=0;l<numb_lesion;l++){
            bool * LesionBool = CreateLesionBool(Labels, l+1, numel);
            float mean_layer = GetMeanData(LayerData, LesionBool, numel);
            for (int i=0;i<numel;i++){
                QualifiedLayer[i]=LesionBool[i]==1?mean_layer:QualifiedLayer[i];

            }
            delete [] LesionBool;
            LesionBool=NULL;
        }
        nifti_image * QualifiedLayerNii = CreateNiiFromArray(QualifiedLayer, LesionToAssess,numel);
        nifti_set_filenames(QualifiedLayerNii, FilenameSave.c_str(),0,0);
        nifti_image_write(QualifiedLayerNii);
        nifti_image_free(Layers);
        nifti_image_free(LesionToAssess);
        nifti_image_free(QualifiedLayerNii);
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_inToAnalyse && segment_analysis->flag_inDataComp & !segment_analysis->flag_RG){
        /* Creation of analysis stats for DataComp over the mask by inToAnalyse*/
        nifti_image * DataCompNii=ReadFromFilename(segment_analysis->filename_inDataComp);
        nifti_image * InToAnalyseNii=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        float * Data=static_cast<float *>(DataCompNii->data);
        int numel = InToAnalyseNii->nvox;
        int numbmodal=DataCompNii->nvox/numel;
        bool * Mask=TranscribeArray<float,bool>(static_cast<float*>(InToAnalyseNii->data),numel);
        float * MinMulti=GetMinDataMulti(Data,Mask,numel,numbmodal);
        float*  MaxMulti=GetMaxDataMulti(Data,Mask,numel,numbmodal);
        float Min=floorf(GetMin(MinMulti,numbmodal));
        float Max=ceilf(GetMax(MaxMulti,numbmodal));
        float * Histogram= GetMaskHistogram(Data,Mask,Min,Max,0.1,numbmodal,numel);


        FilenamePA=nifti_makebasename(segment_analysis->filename_inDataComp);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());

        string FilenamePA2=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index2=FilenamePA2.find_last_of('/');
        string FilenamePA_b2=FilenamePA2.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b2=segment_analysis->name_outputDir;
        }
        string FilenamePA_e2=FilenamePA2.substr(Index2+1,FilenamePA2.length());


        cout << "numbmodal is "<<numbmodal<<endl;
        string FilenameResultsStat=FilenamePA_b+"StatsResults_"+FilenamePA_e + FilenamePA_e2+".csv";
        string FilenameResultsHisto=FilenamePA_b+"HistoResults_"+FilenamePA_e + FilenamePA_e2+".csv";
        ofstream TxtFile(FilenameResultsStat.c_str());
        ofstream HistoFile(FilenameResultsHisto.c_str());
        int maxI=(int)(Max-Min)/0.1;
        float sizeBin=0.1;
        for(int i=0;i<maxI;i++){
            HistoFile<<i*sizeBin+Min<<",";
        }
        HistoFile << endl;
        for(int m=0;m<numbmodal;m++){
            for(int i=0;i<maxI;i++){
                HistoFile<<Histogram[m*maxI+i]<<",";
            }
            HistoFile<<endl;
        }
        TxtFile<<CountNonZero(Mask,numel)<<endl;
        for(int m=0;m<numbmodal;m++){
            float * Statistics=GetStatisticData(&Data[m*numel],Mask,numel);
            TxtFile << "Mask,Data,Mean,Variance,Min,Max,P50,P25,P75,Skewness,Kurtosis";
            TxtFile << endl;
            TxtFile << segment_analysis->filename_inToAnalyse << ",";
            TxtFile << segment_analysis->filename_inDataComp << ",";
            for(int i=0;i<9;i++){
                TxtFile<<Statistics[i]<<",";
            }
            TxtFile << endl;
            delete [] Statistics;
        }

        nifti_image_free(DataCompNii);
        nifti_image_free(InToAnalyseNii);
        delete [] Mask;
        delete [] Histogram;

        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_VesselSeed && (segment_analysis->flag_inMahal|| segment_analysis->flag_inDataComp) && segment_analysis->flag_inToAnalyse){
        nifti_image * HeatNii=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        nifti_image * MahalNii=NULL;



        
        if(! segment_analysis->flag_inMahal){
            float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            float * WMTempOutliers=CreateLong(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            vector<float *> AdditionWM;
            AdditionWM.push_back(WMTempInliers);
            AdditionWM.push_back(WMTempOutliers);
            float * AddedWM=AddArray(AdditionWM, TreeToAnalyse->GetNumberElements());
            nifti_image * AddedWMNii=CreateNiiFromArray(AddedWM, HeatNii, TreeToAnalyse->GetNumberElements());
            nifti_image * WMInliersNii=CreateNiiFromArray(WMTempInliers, HeatNii, TreeToAnalyse->GetNumberElements());
            nifti_image * DataNii = ReadFromFilename(segment_analysis->filename_inDataComp);
            MahalNii = MahalDistMaps(WMInliersNii, AddedWMNii, DataNii );
        }
        else{
            MahalNii=ReadFromFilename(segment_analysis->filename_inMahal);
        }
        
        nifti_image* MaskNii=ReadFromFilename(segment_analysis->filename_mask);
        nifti_image * NewMask=Binarisation(MaskNii);

        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }

        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSave1=FilenamePA_b+"SeedsCheck_"+FilenamePA_e+".nii.gz";
        string FilenameSave2=FilenamePA_b+"SeedsCheckBis_"+FilenamePA_e+".nii.gz";
        string FilenameSave3=FilenamePA_b+"VesselnessT2_"+FilenamePA_e+".nii.gz";
        string FilenameSave4=FilenamePA_b+"VesselnessT1_"+FilenamePA_e+".nii.gz";
        int Dim[3];
        int Shift[3];
        float PixDim[3];
        vector<int>DimVector;
        Shift[0]=1;
        for(int d=0;d<3;d++){
            Dim[d]=MaskNii->dim[d+1];
            PixDim[d]=MaskNii->pixdim[d+1];
            DimVector.push_back(Dim[d]);
            if (d>0){
                Shift[d]=Shift[d-1]*Dim[d-1];
            }

        }
        bool * MaskData=static_cast<bool*>(NewMask->data);
        int numelmaskedWM=0;
        numel=MaskNii->nvox;
        int * S2L=MakeS2L(MaskData,Dim,numelmaskedWM);
        int * L2S=MakeL2S(MaskData,Dim);
        float* WMMapsT2=ExtractTPFromNii<float>(MahalNii, 3);
        // cout << "Max WMMaps T1 " << GetMax(WMMapsT1, numel)<<endl;
        float * BlurredWMMapsT2=GaussianBlurring(WMMapsT2, 0.5, DimVector,0);
        float * SBWMMapsT2=CreateShort(BlurredWMMapsT2, S2L, numelmaskedWM);
        //cout << "Max WMMaps T1 " << GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements())<<endl;
        //MultiplyFloatArrayBy(SBWMMapsT1, GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements()), TreeToAnalyse->GetNumberMaskedElements());
        //cout << "Max WMMaps T1 " << GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements())<<endl;
        float * VesselnessT2=Vesselness(SBWMMapsT2, Dim, Shift, PixDim, L2S, S2L, numelmaskedWM, -1);
        float * VesselLong=CreateLong(VesselnessT2, L2S, numel);
        nifti_image * VesselnessNii=CreateNiiFromArray(VesselLong,HeatNii,numel);
        nifti_set_filenames(VesselnessNii,FilenameSave3.c_str(),0,0);
        nifti_image_write(VesselnessNii);


        float* WMMapsT1=ExtractTPFromNii<float>(MahalNii, 1);
        // cout << "Max WMMaps T1 " << GetMax(WMMapsT1, numel)<<endl;
        float * BlurredWMMapsT1=GaussianBlurring(WMMapsT1, 0.5, DimVector,0);
        float * SBWMMapsT1=CreateShort(BlurredWMMapsT1, S2L, numelmaskedWM);
        //cout << "Max WMMaps T1 " << GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements())<<endl;
        //MultiplyFloatArrayBy(SBWMMapsT1, GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements()), TreeToAnalyse->GetNumberMaskedElements());
        //cout << "Max WMMaps T1 " << GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements())<<endl;
        float * VesselnessT1=Vesselness(SBWMMapsT1, Dim, Shift, PixDim, L2S, S2L, numelmaskedWM, 1);
        float * VesselLongT1=CreateLong(VesselnessT1, L2S, numel);
        nifti_image * VesselnessNiiT1=CreateNiiFromArray(VesselLongT1,HeatNii,numel);
        nifti_set_filenames(VesselnessNiiT1,FilenameSave4.c_str(),0,0);
        nifti_image_write(VesselnessNiiT1);




        float * HeatData=static_cast<float*>(HeatNii->data);
        vector<float *> VectorMult;
        VectorMult.push_back(HeatData);
        VectorMult.push_back(VesselLongT1);
        float * MultHV=MultiplyElementwise(VectorMult,numel);
        nifti_image * MultHVNii=CreateNiiFromArray(MultHV,HeatNii,numel);
        string FilenameSaveHV=FilenamePA_b+"MultHV"+FilenamePA_e+".nii.gz";
        nifti_set_filenames(MultHVNii,FilenameSaveHV.c_str(),0,0);
        nifti_image_write(MultHVNii);
        nifti_image_free(MultHVNii);
        nifti_image * CodeNii=CreateCodePVSFromMahal(MahalNii);
        float * CodeData=static_cast<float*>(CodeNii->data);
        VectorMult.push_back(CodeData);
        float * MultHVC=MultiplyElementwise(VectorMult,numel);
        nifti_image * MultHVCNii=CreateNiiFromArray(MultHVC,HeatNii,numel);
        string FilenameSaveHVC=FilenamePA_b+"MultHVC"+FilenamePA_e+".nii.gz";
        nifti_set_filenames(MultHVCNii,FilenameSaveHVC.c_str(),0,0);
        nifti_image_write(MultHVCNii);
        nifti_image_free(MultHVCNii);
        string FilenameSaveCode=FilenamePA_b+"Code"+FilenamePA_e+".nii.gz";
        nifti_set_filenames(CodeNii,FilenameSaveCode.c_str(),0,0);
        nifti_image_write(CodeNii);
        nifti_image_free(CodeNii);


        bool * Seeds=CreateSeedsForVessel(HeatNii,MahalNii,50);
        int numel = HeatNii->nvox;
        nifti_image * SeedsNii=CreateNiiFromArray(Seeds,HeatNii,numel);
        nifti_set_filenames(SeedsNii,FilenameSave1.c_str(),0,0);
        nifti_image_write(SeedsNii);
        nifti_image * VesselsNii= VesselFromSeeds(Seeds,MahalNii,MaskNii,HeatData,VesselLong);
        string FilenameSave=FilenamePA_b+"EPVSGuessed50"+FilenamePA_e+".nii.gz";
        string FilenameSaveBis=FilenamePA_b+"EPVSGuessed30Bis"+FilenamePA_e+".nii.gz";
        float * VesselsDataFloat=static_cast<float *>(VesselsNii->data);
        bool * VesselsData=TranscribeArray<float,bool>(VesselsDataFloat,numel);
        nifti_set_filenames(VesselsNii,FilenameSave.c_str(),0,0);
        nifti_image_write(VesselsNii);
        nifti_image_free(VesselsNii);

        bool * SeedsBis=CreateSeedsForVessel_bis(HeatNii,VesselsData,MahalNii,25);
        nifti_image * SeedsNiiBis=CreateNiiFromArray(SeedsBis,HeatNii,numel);
        nifti_set_filenames(SeedsNiiBis,FilenameSave2.c_str(),0,0);
        nifti_image_write(SeedsNiiBis);
        nifti_image * VesselsNiiBis= VesselFromSeeds(SeedsBis,MahalNii,MaskNii,HeatData,VesselLong);

        nifti_set_filenames(VesselsNiiBis,FilenameSaveBis.c_str(),0,0);
        nifti_image_write(VesselsNiiBis);
        nifti_image_free(VesselsNiiBis);
        nifti_image_free(MahalNii);
        nifti_image_free(HeatNii);
        delete [] VesselsData;
        delete [] Seeds;
        return EXIT_SUCCESS;
    }
    
    if((segment_analysis->flag_inToAnalyse || segment_analysis->flag_inLes || segment_analysis->flag_inLesCorr) && segment_analysis->flag_Euc){
        nifti_image * LesionImage=NULL;

        if (segment_analysis->flag_inLesCorr) {
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
            LesionImage=ReadFromFilename(segment_analysis->filename_inLesCorr);
        }
        else if(segment_analysis->flag_inToAnalyse){
            FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
            LesionImage=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        }
        else{
            FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
            LesionImage=ReadFromFilename(segment_analysis->filename_InLes);
        }
        float * Data=static_cast<float*>(LesionImage->data);
        numel = LesionImage->nvox;
        bool * SegLesion=TranscribeArray<float, bool>(Data, numel);
        nifti_image * MaskLesion=NULL;
        if(segment_analysis->flag_mask){
            MaskLesion=ReadFromFilename(segment_analysis->filename_mask);
        }
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSave=FilenamePA_b+"DistanceEuc_"+FilenamePA_e+".nii.gz";
        nifti_image *Distance=EuclideanDistanceImage(LesionImage, SegLesion, MaskLesion);
        nifti_set_filenames(Distance, FilenameSave.c_str(), 0, 0);
        nifti_image_write(Distance);
        nifti_image_free(Distance);
        nifti_image_free(LesionImage);
        delete [] SegLesion;
        if (MaskLesion!=NULL) {
            nifti_image_free(MaskLesion);
        }
        if (segment_analysis->flag_inToAnalyse){
            return EXIT_SUCCESS;
        }
    }

    if(segment_analysis->flag_inDataComp && segment_analysis->flag_sc && segment_analysis->flag_inLesCorr){
        nifti_image * Data=ReadFromFilename(segment_analysis->filename_inDataComp);
        nifti_image * Mask=ReadFromFilename(segment_analysis->filename_inLesCorr);
        string FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());

        string FilenameSave=FilenamePA_b+"PotRing_"+FilenamePA_e+".txt";
        int numel=Mask->nvox;
        int Dim[3];
        vector<int> DimVector;
        int Shift[3];
        float PixDim[3];
        float SumPixDim=0;
        float SumDim=0;
        for (int d=0; d<3; d++) {
            Dim[d]=Mask->dim[d+1];
            DimVector.push_back(Dim[d]);
            SumDim+=Dim[d];
            PixDim[d]=Mask->pixdim[d+1];
            SumPixDim+=PixDim[d];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];
        float * MaskData = static_cast<float*>(Mask->data);
        float * DataData = static_cast<float*>(Data->data);
        float * Eroded = ErosionTemplate(MaskData,3,DimVector,1);
        float * Dilated = ErosionTemplate(MaskData, 3, DimVector,0);
        float * AllowedRegion = SubtractArray(Dilated,Eroded,numel);
        float * PotRidge = CreateRidge(DataData,numel,Dim,Shift,PixDim,2);
        MultiplyElementwiseChange(PotRidge,AllowedRegion,numel);
        float * PotRidgeFin = ThresholdArray<float,float>(PotRidge,0.15,numel);
        nifti_image * PotNii = CreateNiiFromArray(PotRidgeFin,Mask,numel);
        nifti_set_filenames(PotNii,FilenameSave.c_str(),0,0);
        nifti_image_write(PotNii);
        nifti_image_free(PotNii);
        delete [] AllowedRegion;
        delete [] Eroded;
        delete [] Dilated;
        delete [] PotRidge;
        delete [] PotRidgeFin;
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_sc && segment_analysis->flag_inLesCorr && segment_analysis->flag_inVentricleSeg){
        cout << "Correction spurious in plane " << endl;
        string FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        nifti_image * VentrSeg = ReadFromFilename(segment_analysis->filename_inVentricleSeg);
        float * VentrSegData = static_cast<float * >(VentrSeg->data);
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSave=FilenamePA_b+"CorrPlanes_"+FilenamePA_e+".nii.gz";
        nifti_image * LesionToCorr = ReadFromFilename(segment_analysis->filename_inLesCorr);
        float * LesionToCorrData = static_cast<float*>(LesionToCorr->data);

        int numel = LesionToCorr->nvox;
        int Dim[3];
        vector<int> DimVector;
        int Shift[3];
        float PixDim[3];
        float SumPixDim=0;
        float SumDim=0;
        for (int d=0; d<3; d++) {
            Dim[d]=LesionToCorr->dim[d+1];
            DimVector.push_back(Dim[d]);
            SumDim+=Dim[d];
            PixDim[d]=LesionToCorr->pixdim[d+1];
            SumPixDim+=PixDim[d];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];

        int BB[6];
        BB[0]=0;
        BB[1]=Dim[0];
        BB[2]=0;
        BB[3]=Dim[1];
        cout << "Dim are "<< Dim[0] << " "<< Dim[1]<< " "<< Dim[2]<< endl;
        vector<float*> CorrectedVector;
        cout << "Initial elements is "<< CountNonZero(LesionToCorrData,numel);
        bool * Remove=new bool[numel];
        for(int i=0;i<numel;i++){
            Remove[i]=0;
        }
        for(int z=0;z<Dim[2];z++){
            cout << "Treating plane "<<z<<endl;
            BB[4]=z;
            BB[5]=z;
            bool * SliceMask = CreateBoolBoundingBox(BB,Dim);
            float * LesionSlice=MaskArray<float,bool,float>(LesionToCorrData,SliceMask,numel);
            bool * VentrSlice = MaskArray<float,bool,bool>(VentrSegData,SliceMask,numel);
            delete [] SliceMask;
            int * Labels=ComponentLabeling(LesionSlice,6,Dim,Shift,0.5);
            int * OrderedLabels = OrderedVolumeLabel(Labels,2,numel,1);
            int max_lab= GetMaxLabel(OrderedLabels,numel);
            for(int l=0;l<max_lab;l++){
                bool * LesionBool = CreateLesionBool(OrderedLabels,l+1,numel);
                float LabelSize=SumOverMask<float,bool>(LesionSlice,LesionBool,0,numel);
                if (LabelSize<segment_analysis->MiniSize){
                    if (LabelSize<segment_analysis->MiniSize-1){
                        AddElementwiseInPlace(Remove,LesionBool,numel);
                    }
                    if(GetDistanceBetweenSeg(VentrSlice,LesionBool,Dim,Shift,PixDim)<sqrt(PixDim[0]*PixDim[0]+PixDim[1]*PixDim[1])){
                        AddElementwiseInPlace(Remove,LesionBool,numel);
                    }
                }
                delete [] LesionBool;
            }

             float * CorrectedLesion = MaskArray<float,int,float>(LesionSlice,OrderedLabels,numel);
             CorrectedVector.push_back(CorrectedLesion);
             delete [] Labels;
             delete [] OrderedLabels;
             delete [] LesionSlice;
        }
        float * CorrectedLesionFin = AddElementwise(CorrectedVector,numel);
        bool * Kept=OpposeBoolArray(Remove,numel);
        float * CorrectedLesionFinFin=MaskArray<float,bool,float>(CorrectedLesionFin,Kept,numel);
        cout << "Annulled due to size / border is "<< CountNonZero(Remove,numel)<<endl;
        cout << "Final lesion is "<< CountNonZero(CorrectedLesionFinFin, numel)<<endl;
        nifti_image * CorrectedLesionNii=CreateNiiFromArray(CorrectedLesionFinFin, LesionToCorr,numel);
        nifti_set_filenames(CorrectedLesionNii,FilenameSave.c_str(),0,0);
        nifti_image_write(CorrectedLesionNii);
        nifti_image_free(CorrectedLesionNii);
        int sizeVec=CorrectedVector.size();
        for(int i=0;i<sizeVec;i++){
            if (CorrectedVector[i]!=NULL){
                delete [] CorrectedVector[i];
                CorrectedVector[i] = NULL;
            }
        }
        delete [] Remove;
        delete [] Kept;
        delete[] CorrectedLesionFin;
        delete [] CorrectedLesionFinFin;
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_inLobes && segment_analysis->flag_sc){
        string FilenamePA=nifti_makebasename(segment_analysis->filename_inLobes[0]);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());

        string FilenameSave=FilenamePA_b+"NumbComponents_"+FilenamePA_e+".txt";
        ofstream Txt(FilenameSave.c_str());
        for(int i=0;i<segment_analysis->filename_inLobes.size();i++){
            nifti_image * Image = ReadFromFilename(segment_analysis->filename_inLobes[i]);
            int * Labels = ComponentLabeling(Image,6,0.5);
            int numel = Image->nvox;
            float VolVox = 1;
            int * Ordered = OrderedVolumeLabel(Labels,3,numel,VolVox);
            int numbLabels=GetMax(Ordered,numel);
            Txt << segment_analysis->filename_inLobes[i] << " "<<numbLabels<<endl;
            delete []Ordered;
            delete[]Labels;
            nifti_image_free(Image);

        }
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_LapIn && segment_analysis->flag_Parc && segment_analysis->flag_CorExt && segment_analysis->flag_mask){
       cout<<" Attempting expansion of cortext following laplace solution " << endl;
//         Create the mask from which to get the distance map according to used parcellation
        nifti_image * ParcNii = ReadFromFilename(segment_analysis->filename_Parc);
        nifti_image * MaskNii = ReadFromFilename(segment_analysis->filename_mask);
        nifti_image * LaplaceNii = ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
        float * LaplaceSolLong = static_cast<float*>(LaplaceNii->data);
        numel = ParcNii->nvox;
        int Dim[3];
        vector<int> DimVector;
        int Shift[3];
        float PixDim[3];
        float SumPixDim=0;
        float SumDim=0;
        for (int d=0; d<3; d++) {
            Dim[d]=ParcNii->dim[d+1];
            DimVector.push_back(Dim[d]);
            SumDim+=Dim[d];
            PixDim[d]=ParcNii->pixdim[d+1];
            SumPixDim+=PixDim[d];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];
        bool * MaskData=static_cast<bool*>(MaskNii->data);
        int * L2S=MakeL2S(MaskData,Dim);
        int numelmasked=0;
        int * S2L=MakeS2L(MaskData,Dim,numelmasked);
        float * LaplaceSol = CreateShort(LaplaceSolLong, S2L, numelmasked);
        float * ParcData = static_cast<float*>(ParcNii->data);
        int * ParcDataInt = TranscribeArray<float, int>(ParcData, numel);
        cout << "Parc Data min max is "<< GetMin(ParcData,numel) << " "<< GetMax(ParcData,numel) << endl;
        bool * MaskCortex = ThresholdArray<float, bool>(ParcData, segment_analysis->thresh_cortex,numel);
        cout << "Mask Cortex contains "<< CountNonZero(MaskCortex, numel) << endl;
//        Create the euclidean distance map and threshold for mask
        nifti_image * EucDist = EuclideanDistanceImage(ParcNii, MaskCortex, MaskNii);
        int numel = EucDist->nvox;
        float * EucDistData = static_cast<float*>(EucDist->data);
        bool * ValidZone = UpperThresholdArray<float, bool>(EucDistData, segment_analysis->val_ext+1,numel);


//        For loops for all labels in segmentation to create the final distance map and the attributed label
        int * AttributedExtension = new int [numel];
        float * CurrentDist = new float[numel];
        for(int i=0;i<numel;i++){
            AttributedExtension[i]=0;
            CurrentDist[i]=segment_analysis->val_ext+2;
        }
        int max_label = GetMax(ParcDataInt, MaskData, numel);
        for(int l=segment_analysis->thresh_cortex; l<=max_label; l++){

            bool * LobeBool = CreateLesionBool(ParcDataInt, l,numel);
            if (CountNonZero(LobeBool, numel)>0){
                cout << "Treating label "<< l << " that contains "<< CountNonZero(LobeBool,numel) << endl;
            bool * LobeBoolShort = CreateShort(LobeBool,S2L,numelmasked);
            float * EDFloat=DistanceFromLobe(LaplaceSol,ValidZone,LobeBoolShort,L2S,S2L,numelmasked,PixDim,Dim,Shift,ParcNii);
            float * EDFloatLong=CreateLong(EDFloat,L2S,numel);
            for (int i = 0;i<numel;i++){
                if(ValidZone[i]){
                    if(EDFloatLong[i] < CurrentDist[i]){
                        CurrentDist[i] = EDFloatLong[i];
                        AttributedExtension[i]=l;
                    }
                }
            }
            delete [] EDFloat;
            delete [] EDFloatLong;

            delete [] LobeBoolShort;
            }
            delete [] LobeBool;
        }
//        Cleaning to make sure we stay under value
        cout << "Number in original Attribution " << CountNonZero(AttributedExtension, numel) << endl;
        for(int i=0;i<numel;i++){
            if(CurrentDist[i]>segment_analysis->val_ext){
                AttributedExtension[i]=0;
            }
        }
        cout << "Number in final Attribution "<< CountNonZero(AttributedExtension, numel) << endl;
        FilenamePA=nifti_makebasename(segment_analysis->filename_Parc);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSave=FilenamePA_b+"AttributionExt"+FilenamePA_e+".nii.gz";
        nifti_image * ResultNii = CreateNiiFromArray(AttributedExtension, ParcNii, numel);
        nifti_set_filenames(ResultNii, FilenameSave.c_str(), 0,0);
        nifti_image_write(ResultNii);
        nifti_image_free(ResultNii);
        nifti_image_free(ParcNii);
        nifti_image_free(MaskNii);
        nifti_image_free(EucDist);
        delete [] CurrentDist;
        delete [] AttributedExtension;
        delete [] ValidZone;
        delete [] MaskCortex;
        delete [] S2L;
        delete [] L2S;
        return EXIT_SUCCESS;
    }
    
    //    Case where we simply build the lesion zones based on distance to lobes
    if(segment_analysis->flag_inLobes && segment_analysis->flag_inLesCorr && segment_analysis->flag_mask){
        cout << "Lobar creation";
        FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        
        
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameSave=FilenamePA_b+"DistanceChoice_"+FilenamePA_e+".nii.gz";
        string FilenameSavePb=FilenamePA_b+"DistanceChoicePb_"+FilenamePA_e+".nii.gz";
        string FilenameSavePb2=FilenamePA_b+"DistanceChoicePb2_"+FilenamePA_e+".nii.gz";
        //        vector<nifti_image*> VectorNiiLobes;
        vector<nifti_image*> VectorNiiDistances;
        cout << "Loading lesion file";
        nifti_image * LesionNii=ReadFromFilename(segment_analysis->filename_inLesCorr);
        int numel=LesionNii->nvox;
        int Dim[3];
        vector<int> DimVector;
        int Shift[3];
        float PixDim[3];
        float SumPixDim=0;
        float SumDim=0;
        for (int d=0; d<3; d++) {
            Dim[d]=LesionNii->dim[d+1];
            DimVector.push_back(Dim[d]);
            SumDim+=Dim[d];
            PixDim[d]=LesionNii->pixdim[d+1];
            SumPixDim+=PixDim[d];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];
        cout << "Loading mask";
        nifti_image * MaskNii=ReadFromFilename(segment_analysis->filename_mask);
        bool * MaskData=static_cast<bool*>(MaskNii->data);
        int * L2S=MakeL2S(MaskData,Dim);
        int numelmasked=0;
        int * S2L=MakeS2L(MaskData,Dim,numelmasked);
        cout << "Mapping S2L L2S created numel is "<< numel<<" "<<MaskNii->nvox <<endl;
        nifti_image * Test1=CreateNiiFromArray(MaskData, LesionNii, numel);
        cout << "Image created" << endl;
        cout << FilenamePA_b <<" " << FilenamePA_e << endl;
        string FilenameSaveMask=FilenamePA_b+"MaskTest_"+FilenamePA_e+".nii.gz";
        cout << "Name to save image mask" << FilenameSaveMask << endl;
        nifti_set_filenames(Test1, FilenameSaveMask.c_str(), 0, 0);
        cout << "Filename set" << endl;
        // nifti_image_write(Test1);
        cout << "Image written" <<endl;
        nifti_image_free(Test1);
        cout << "Memory freed" << endl;
        Test1=NULL;
        int numbLobes=segment_analysis->filename_inLobes.size();
        cout << "Number of lobes is "<<numbLobes;
        vector<float *> Hemisphere;
        nifti_image * ParcNii = NULL;
        float * ParcData = NULL;
        float * HemisphereLeft = new float[numel];
        float * HemisphereRight  = new float[numel];
        nifti_image * HemiNii = NULL;
        float * HemiData = NULL;
        if (segment_analysis->filename_Parc!=NULL){
                ParcNii = ReadFromFilename(segment_analysis->filename_Parc);
                ParcData = static_cast<float*>(ParcNii->data);
        }
        else if(segment_analysis->filename_Hemi!=NULL){
            cout << "Preparation using existing hemispheres";
            HemiNii = ReadFromFilename(segment_analysis->filename_Hemi);
            HemiData = static_cast<float*>(HemiNii->data);
        }
        else {
            ParcData = new float[numel];
            for (int i = 0; i < numel ; i++){
                ParcData[i] = 0;
            }
        }
        if (ParcData==NULL){
            for(int i=0;i<numel;i++){
                if (HemiData[i] == 1){
                    HemisphereLeft[i] = 1;
                    HemisphereRight[i] =0;
                }
                else if (HemiData[i]==2){
                    HemisphereRight[i] = 1;
                    HemisphereLeft[i] = 0;
                }
                else{
                    HemisphereRight[i] = 0;
                    HemisphereRight[i] =0;
                }
            }
            MultiplyFloatArrayBy(HemisphereLeft,1000,numel);
            MultiplyFloatArrayBy(HemisphereRight,1000,numel);
        }
        else{
            float * HemisphereLeft = CreateHemisphereFromParc(ParcData, 0, numel);
            MultiplyFloatArrayBy(HemisphereLeft, 1000, numel);
            float * HemisphereRight = CreateHemisphereFromParc(ParcData, 1, numel);
            cout << "Sum Hemisphere Right "<< GetSum(HemisphereRight, numel) << endl;
            MultiplyFloatArrayBy(HemisphereRight,1000, numel);

        }
                cout << "New sum HR" << GetSum(HemisphereRight, numel) << endl;
        Hemisphere.push_back(HemisphereLeft);
        Hemisphere.push_back(HemisphereRight);
        for(int i=0;i<numbLobes-2;i++){
            nifti_image * LobelNii=ReadFromFilename(segment_analysis->filename_inLobes[i]);
            bool * LobeBool=static_cast<bool*>(LobelNii->data);
            nifti_image * Test=CreateNiiFromArray(LobeBool, LesionNii, numel);
            nifti_image * EDImage = NULL;
            if(segment_analysis->flag_LaplaceSolImage && segment_analysis->flag_inVentricleSeg){
                nifti_image * Laplace=ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
                float * LaplaceSol = static_cast<float*>(Laplace->data);
                nifti_image * VentrSeg = ReadFromFilename(segment_analysis->filename_inVentricleSeg);
                bool * VentrBoolLong=static_cast<bool *> (VentrSeg->data);
                bool * VentrBool=CreateShort(VentrBoolLong,S2L,numelmasked);
                bool * LobeBoolShort = CreateShort(LobeBool,S2L,numelmasked);
                float * EDFloat=DistanceFromLobe(LaplaceSol,NULL,LobeBoolShort,L2S,S2L,numelmasked,PixDim,Dim,Shift,LesionNii);
                float * EDFloatLong=CreateLong(EDFloat,L2S,numel);

                AddElementwiseInPlace(EDFloatLong, Hemisphere[Side[i]], numel);
                EDImage=CreateNiiFromArray(EDFloatLong,Laplace,numel);
                delete [] EDFloat;
                delete [] EDFloatLong;
                nifti_image_free(VentrSeg);
                nifti_image_free(Laplace);
                Laplace=NULL;
                VentrSeg=NULL;
                delete [] VentrBool;
                delete [] LobeBoolShort;
            }
            
            else{
                EDImage=EuclideanDistanceImage(LesionNii, LobeBool, MaskNii);
                float * EDFloatLong = static_cast<float*>(EDImage->data);
                AddElementwiseInPlace(EDFloatLong, Hemisphere[Side[i]], numel);
                cout << " new max is "<< GetMaxArrayMasked(EDFloatLong,MaskData, numel);
            }
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
            Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenameSavePA=FilenamePA_b+"Distance_"+FilenamePA_e+".nii.gz";
            nifti_set_filenames(EDImage, FilenameSavePA.c_str(), 0, 0);
            nifti_image_write(EDImage);
            VectorNiiDistances.push_back(EDImage);
            nifti_image_free(LobelNii);
            nifti_image_free(Test);
            LobelNii=NULL;
        }
        for (int i=numbLobes-2; i<numbLobes; i++) {
            nifti_image * LobeNii=ReadFromFilename(segment_analysis->filename_inLobes[i]);
            bool * LobeBool=static_cast<bool*>(LobeNii->data);
            bool * OppLobe=OpposeBoolArray(LobeBool, numel);
            float * OppLobeFloat=TranscribeArray<bool, float>(OppLobe, numel);
            MultiplyFloatArrayBy(OppLobeFloat, numel, numel);
            nifti_image * EDImage=CreateNiiFromArray(OppLobeFloat, LesionNii, numel);
            VectorNiiDistances.push_back(EDImage);
            delete [] OppLobeFloat;
            delete [] OppLobe;
            OppLobe = NULL;
            OppLobeFloat = NULL;
            nifti_image_free(LobeNii);
            LobeNii=NULL;
        }

        vector<float*>VectorDataDistances;
        for(int i=0;i<numbLobes;i++){
            VectorDataDistances.push_back(static_cast<float*>(VectorNiiDistances[i]->data));
        }
        float * ResultsMinDist=new float[numel];
        int CountMask=0;
        for(int i=0;i<numel;i++){
            float maxDist=numel;
            if(!MaskData[i]){
                ResultsMinDist[i]=0;
            }
            else{
                CountMask++;
                for(int n=0;n<numbLobes;n++){
                    if(VectorDataDistances[n][i]<maxDist){
                        maxDist=VectorDataDistances[n][i];
                        ResultsMinDist[i]=n+1;
                    }
                }
            }
        }
        cout << "Number of elements in Mask is "<<CountMask << endl;
        int * ResultsMinDistInt1 = TranscribeArray<float, int>(ResultsMinDist,numel);
        int max_class = GetMax<int>(ResultsMinDistInt1,numel);
        float * BlurredRes=NULL;
        int * class_decide = new int[max_class];
//        float * CopyResults = TranscribeArray<int, float>(ResultsMinDistInt, numel);
        BlurredRes = GaussianBlurring(ResultsMinDist, 1, DimVector);
//        int * BlurredResInt = TranscribeArray<float, int>(BlurredRes, numel);
        int CountPb0 = 0;
        for (int i=0;i<numel;i++){
            int ListNeigh[18];
            if(MaskData[i] && ResultsMinDist[i]!=BlurredRes[i] ){
                for(int c=0;c<max_class;c++){
                    class_decide[c] = 0;
                }
                GetListNeighbours_bis(ListNeigh,i,Dim,Shift,18);
                for(int n = 0;n<18;n++){
                    if (BlurredRes[i]>0){
                        if (ResultsMinDistInt1[ListNeigh[n]] > 0){
                        class_decide[ResultsMinDistInt1[ListNeigh[n]]-1] +=1;
                        }
                        else{
                            CountPb0++;
                        }
                }
                int new_value = GetIndexMax(class_decide, max_class) + 1;
                BlurredRes[i] = new_value;
            }

        }
            if(!MaskData[i]){
                BlurredRes[i] =0;
            }
    }
    cout << "Count Potential list access issue "<< CountPb0 << endl;
    int * ResultsMinDistInt = TranscribeArray<float, int>(BlurredRes,numel);
//    int * ResultsMinDistInt = TranscribeArray<float, int>(ResultsMinDist,numel);
        int numb_valid = CountNonZero(ResultsMinDistInt, numel);
int * MapProblematic1= new int[numel];
int * MapProblematic2= new int[numel];
int CountPbValue=0;
        if(segment_analysis->flag_removespurious){
            cout << "Removing spurious voxels "<< max_class << endl;
            for(int i=1;i<=max_class;i++){
                cout << "treating class "<< i << " ";
                bool * ClassTemp = CreateLesionBool(ResultsMinDistInt, i, numel);
                float * ClassFloat = TranscribeArray<bool, float>(ClassTemp, numel);
                int *ComponentLabels=ComponentLabeling(ClassFloat,6,Dim,Shift);
                stringstream ws;
                ws << i;
                string FilenameSaveTemp=FilenamePA_b+"DistanceChoice_"+ws.str()+"_"+FilenamePA_e+".nii.gz";
                nifti_image * TempNii = CreateNiiFromArray(ComponentLabels, LesionNii,numel);
                nifti_set_filenames(TempNii, FilenameSaveTemp.c_str(),0,0);
                nifti_image_write(TempNii);

//                int * VolThresh=OrderedVolumeLabel(ComponentLabels,0,numel,1);
                int max_ncc = GetMax<int>(ComponentLabels, numel);
                cout << "number of labels is "<< max_ncc << endl;
                for(int l=1;l<=max_ncc;l++){
                    int size_label = GetLabelSize(l, ComponentLabels, numel);
                    if (size_label<15){
                        cout << "for class "<< i<< " label "<< l<<" size "<< size_label << endl;
                        for (int n=0;n<numel;n++){
                            if (ComponentLabels[n]==l){
                                ResultsMinDistInt[n]=0;
                            }
                        }
                    }

                }
                nifti_image_free(TempNii);
                TempNii=NULL;
//                delete [] VolThresh;
                delete [] ComponentLabels;
                delete [] ClassTemp;
                delete [] ClassFloat;
            }

// Go through element and assign the ones that are missing to surrounding voxel tissue
            for (int i=0;i<numel;i++){
                MapProblematic1[i]=0;
                MapProblematic2[i]=0;
//                if(MaskData[i] && ResultsMinDistInt[i]==0 ){
                    if(MaskData[i]  ){
                int ListNeigh[18];
                GetListNeighbours_bis(ListNeigh, i,Dim,Shift,18);


                    for(int c=0;c<max_class;c++){
                        class_decide[c] = 0;
                    }
                    GetListNeighbours_bis(ListNeigh,i,Dim,Shift,18);
                    for(int n = 0;n<18;n++){
                        if (ResultsMinDistInt[ListNeigh[n]]>0){
                        class_decide[ResultsMinDistInt[ListNeigh[n]]-1] +=1;
                        }
                    }


                    int new_value = GetIndexMax(class_decide, max_class)+1;
                    int second = GetIndexSecondMax(class_decide, max_class, new_value-1);

                    if (ResultsMinDistInt[i]==0){
                    ResultsMinDistInt[i] = new_value ;
                    }
                    else {
                        if (ResultsMinDistInt[i]!=new_value || class_decide[new_value-1]==class_decide[second]){
                            CountPbValue++;
                            if(ResultsMinDistInt[i]!=new_value){
                            MapProblematic1[i]=class_decide[new_value-1]-class_decide[second];
                                                cout << "Repartition of label "  ;
                                                for (int c=0;c<max_class;c++){
                                                    cout << class_decide[c]<<" ";
                                                }
                                                cout << endl;
                            if (MapProblematic1[i]>6){
                                ResultsMinDistInt[i]=new_value;
                            }
                            }
                            if(class_decide[new_value-1]==class_decide[second]){
                            MapProblematic2[i]=new_value;
                            cout << "Undecided between "<< new_value << " and "<< second+1 << endl;
                            }
                        }
                    }
                }
            }
            cout << "Pb neighborhood is "<< CountPbValue << endl;
//            float * CopyResults = TranscribeArray<int, float>(ResultsMinDistInt, numel);
//            BlurredRes = GaussianBlurring(CopyResults, 1, DimVector);
//            for (int i=0;i<numel;i++){
//                int ListNeigh[18];
//                if(MaskData[i] && ResultsMinDistInt[i]!=BlurredRes[i] ){
//                    for(int c=0;c<max_class;c++){
//                        class_decide[c] = 0;
//                    }
//                    GetListNeighbours_bis(ListNeigh,i,Dim,Shift,18);
//                    for(int n = 0;n<18;n++){
//                        if (BlurredRes[i]>0){
//                            class_decide[ResultsMinDistInt[ListNeigh[n]]-1] +=1;
//                    }
//                    int new_value = GetIndexMax(class_decide, max_class) + 1;
//                    BlurredRes[i] = new_value;
//                }

//            }
//                if(!MaskData[i]){
//                    BlurredRes[i] =0;
//                }
//        }
//                    delete [] CopyResults;
//        }
}
        int numb_valid2 = CountNonZero(ResultsMinDistInt, numel);
        cout << "Results are of size "<< numb_valid << " then "<< numb_valid2<<endl;
        nifti_image * MapNii = CreateNiiFromArray(MapProblematic1, LesionNii, numel);
        nifti_set_filenames(MapNii, FilenameSavePb.c_str(),0,0);
        nifti_image_write(MapNii);
        nifti_image_free(MapNii);
        delete [] MapProblematic1;
        MapProblematic1 = NULL;

        nifti_image * MapNii2 = CreateNiiFromArray(MapProblematic2, LesionNii, numel);
        nifti_set_filenames(MapNii2, FilenameSavePb2.c_str(),0,0);
        nifti_image_write(MapNii2);
        nifti_image_free(MapNii2);
        delete [] MapProblematic2;
        MapProblematic2 = NULL;


        nifti_image * ResultsRegion=CreateNiiFromArray(ResultsMinDistInt, LesionNii, numel);
        nifti_set_filenames(ResultsRegion, FilenameSave.c_str(), 0, 0);
        nifti_image_write(ResultsRegion);
        delete [] ResultsMinDist;
        delete [] ResultsMinDistInt;

        delete [] class_decide;
        class_decide = NULL;

        if (BlurredRes!=NULL){
            delete [] BlurredRes;
        }
        for(int i=0;i<numbLobes;i++){
            nifti_image_free(VectorNiiDistances[i]);
            VectorNiiDistances[i]=NULL;
        }
        nifti_image_free(ResultsRegion);
        if (ParcNii !=NULL){
            nifti_image_free(ParcNii);
            ParcData = NULL;
        }
        if (ParcData!=NULL){
            delete [] ParcData;
            ParcData = NULL;
        }
        nifti_image_free(LesionNii);
        LesionNii=NULL;
        nifti_image_free(MaskNii);
        MaskNii=NULL;
        delete [] ResultsMinDistInt1;
        ResultsMinDistInt1 = NULL;
        delete [] Hemisphere[0];
        delete [] Hemisphere[1];
        delete [] L2S;
        delete [] S2L;
        delete [] segment_analysis;
        return EXIT_SUCCESS;
    }
    
    
    
    //    Case where we simply extract the volume information locally based on octants and 4 layers.
    if (segment_analysis->flag_LaplaceSolImage && (segment_analysis->flag_inLes) && segment_analysis->flag_inQuadrant && segment_analysis->flag_LS) {
        cout << "Doing Local summary first"<<endl;
        nifti_image * Lesion=ReadFromFilename(segment_analysis->filename_InLes);
        nifti_image * Quadrant=ReadFromFilename(segment_analysis->filename_inQuadrant);
        nifti_image * Layers=ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
        nifti_image * MahalWM=NULL;
        int numbMahal=0;
        
        
        // Building Tree
        
        if (segment_analysis->flag_TextFile){
            cout<<"Building tree from text file"<<endl;
            if (Modalities.size()==0) {
                Modalities = GetModalitiesFromTextFile(segment_analysis);
            }
            TreeToAnalyse=ReadTreeFromTextFile(segment_analysis->filename_TextFile, NULL,segment_analysis->filename_changePath);
            SEG_PARAMETERS * segment_paramAnalyse=new SEG_PARAMETERS[1]();
            segment_paramAnalyse->IndexCSF=segment_analysis->IndexCSF;
            segment_paramAnalyse->IndexGM=segment_analysis->IndexGM;
            segment_paramAnalyse->IndexOut=segment_analysis->IndexOut;
            segment_paramAnalyse->IndexWM=segment_analysis->IndexWM;
            
            numbmodal=TreeToAnalyse->GetNumberModalities();
            
            //        float * DataTest=static_cast<float *>(TreeToAnalyse->GetDataImage()->data);
            //        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+numel]<<endl;
            nifti_image * SegToUse=ReadFromFilename(segment_analysis->filename_SegTot);
            cout<<"SegToUse read"<<endl;
            // Check if number of elements in SegToUse correspond to number of leaves in TreeToAnalyse
            VectorLeaves=TreeToAnalyse->GetAllLeaves();
            numbleavesTree=VectorLeaves.size();
            int numbleavesSeg=SegToUse->nu*SegToUse->nt;
            if (numbleavesSeg!=numbleavesTree) {
                cout<<"Incompatible text file and segmentation"<<endl;
            }
            else {
                float * SegData=static_cast<float *>(SegToUse->data);
                float * SegData_PTR=SegData;
                int * L2S_PTR=TreeToAnalyse->GetL2S();
                int numel=TreeToAnalyse->GetNumberElements();
                int numelmasked=TreeToAnalyse->GetNumberMaskedElements();
                cout<<"ok for number of elements"<<endl;
                for (int l=0; l<numbleavesTree; l++) {
                    L2S_PTR=TreeToAnalyse->GetL2S();
                    SegData_PTR=&SegData[l*numel];
                    float * NormRespNew=new float[numelmasked];
                    for (int i=0; i<numelmasked; i++) {
                        NormRespNew[i]=0;
                    }
                    float * NormRespNew_PTR=NormRespNew;
                    for (int i=0; i<numel; i++, L2S_PTR++,SegData_PTR++) {
                        if (*L2S_PTR>=0) {
                            *NormRespNew_PTR=*SegData_PTR;
                            NormRespNew_PTR++;
                        }
                    }
                    cout<<"NormResp new calculated for "<<l<<"...";
                    VectorLeaves[l]->SetNormResp(NormRespNew);
                    delete [] NormRespNew;
                    NormRespNew=NULL;
                    cout<<"...and set";
                }
                cout<<"TreeToAnalyse rebuilt and saved..."<<endl;
                vector<int> Modalities2 = GetModalitiesFromTextFile(segment_analysis);
                if (Modalities2.size()<=Modalities.size()) {
                    numbmodal=Modalities2.size();
                }
                if (Modalities.size()==0) {
                    Modalities=GetModalitiesFromTextFile(segment_analysis);
                    numbmodal=Modalities.size();
                }
                if (TreeToAnalyse->GetNumberModalities()>numbmodal) {
                    CorrectDataAccordingToModa(TreeToAnalyse, Modalities);
                }
                
                RebuildNormRespFromLeaves(TreeToAnalyse);
                //            TreeToAnalyse->SaveBFCorrectedData(g_scratch_dir + "/TestCorrData.nii.gz");
                //            TreeToAnalyse->SaveAllClasses(g_scratch_dir + "/TestClassesBegin.nii.gz", segment_paramAnalyse);
//                numel=TreeToAnalyse->GetNumberElements();
                cout<<"Tree completely built"<<endl;
                if (segment_analysis->flag_juxtaCorrection) {
                    segment_paramAnalyse->flag_DGMPrior=1;
                    string CopyDGMPrior=segment_analysis->filename_inPriorsDGM;
                    segment_paramAnalyse->filename_DGMPrior=CopyDGMPrior;
                    cout<<"Trying correction for juxtacortical elements in reconstruction of lesions"<<endl;
                    numbmodal=TreeToAnalyse->GetNumberModalities();
                    for (int m=0; m<numbmodal; m++) {
                        segment_paramAnalyse->Modalities.push_back(Modalities[m]);
                        cout << "Modalities "<< m<<" "<< segment_paramAnalyse->Modalities[m]<< endl;
                    }
                    //                TreeToAnalyse->SaveBFCorrectedData(g_scratch_dir + "/TestCorrData.nii.gz");
                    //                TreeToAnalyse->SaveAllClasses(g_scratch_dir + "/TestClassesBegin.nii.gz", segment_paramAnalyse);
                    segment_paramAnalyse->flag_out=0;
                    segment_paramAnalyse->flag_inDCFile=0;
                    segment_paramAnalyse->flag_inDC=0;
                    CorrectionJuxta=TreeToAnalyse->CorrectInitialEMResultForWrongWMHoles(segment_paramAnalyse);
                    
                    cout<<".... juxtacortical correction performed !"<<endl;
                    //                nifti_image * CorrectionJuxta2=HolesJuxta(TreeToAnalyse, segment_analysis);
                    //                float * CorrJuxtaData=NULL;
                    //                float * CorrJuxta2Data=NULL;
                    //
                    //                if (CorrectionJuxta!=NULL) {
                    //                    CorrJuxtaData=static_cast<float*>(CorrectionJuxta->data);
                    //                    if (CorrectionJuxta2!=NULL) {
                    //                        CorrJuxta2Data=static_cast<float*>(CorrectionJuxta2->data);
                    //                        for (int i=0; i<numel; i++) {
                    //                            CorrJuxtaData[i]+=CorrJuxta2Data[i];
                    //                        }
                    //                    }
                    //                }
                    //                nifti_image_free(CorrectionJuxta2);
                    //                CorrectionJuxta2=NULL;
                    
                }
                
                numel=TreeToAnalyse->GetNumberElements();
                //            SEG_PARAMETERS * segment_param=NULL;
                //            TreeToAnalyse->SaveAllClasses(g_scratch_dir + "/StrangeRebuilt.nii.gz", segment_param);
            }
            nifti_image_free(SegToUse);
            SegToUse=NULL;
            //        if (segment_paramAnalyse!=NULL) {
            //            delete segment_paramAnalyse;
            //            segment_paramAnalyse=NULL;
            //        }
            
        }

        cout<<segment_analysis->flag_LesWMI<<" LesWMI"<<endl;


        
        
        
        
        if (TreeToAnalyse!=NULL  && (segment_analysis->flag_inLesCorr || segment_analysis->flag_inLes)) {
            cout <<"Doing WMI analysis"<<endl;
            nifti_image * LesionCorr=NULL;
            if(segment_analysis->flag_inLesCorr){
                LesionCorr=ReadFromFilename(segment_analysis->filename_inLesCorr);
            }
            else if (segment_analysis->flag_inLes){
                LesionCorr=ReadFromFilename(segment_analysis->filename_InLes);
            }
            cout << "Image lesion Read"<<endl;
            float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            float * WMTempOutliers=CreateLong(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            cout << "Temp seg created "<< endl;
            float * LesionCorrData=static_cast<float*>(LesionCorr->data);
            vector<float *> AdditionWM;
            AdditionWM.push_back(WMTempInliers);
            AdditionWM.push_back(WMTempOutliers);
            AdditionWM.push_back(LesionCorrData);
            float * AddedWM=AddArray(AdditionWM, TreeToAnalyse->GetNumberElements());
            cout << "Addition done"<< endl;
            nifti_image * AddedWMNii=CreateNiiFromArray(AddedWM, LesionCorr, TreeToAnalyse->GetNumberElements());
            nifti_image * WMInliersNii=CreateNiiFromArray(WMTempInliers, LesionCorr, TreeToAnalyse->GetNumberElements());
            numel=TreeToAnalyse->GetNumberElements();
            numbmodal=TreeToAnalyse->GetNumberModalities();
            float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
            cout <<"Data corr obtained"<<endl;
            nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
            cout << "Data Nii"<<DataNii<<endl;
            MahalWM= MahalDistMaps(WMInliersNii, AddedWMNii, TreeToAnalyse->GetDataImage() );
            cout <<"Mahal 1 obtained"<<endl;
            MahalWM= MahalDistMaps(WMInliersNii, AddedWMNii, DataNii );
            nifti_image_free(DataNii);
            delete [] DataCorr;
            numbMahal=MahalWM->dim[4];
            cout <<"Mahal obtained "<<endl;
            
            FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            
            string FilenameSave=FilenamePA_b+"LocalSummaryWMIStatsFin_"+FilenamePA_e+".txt";
            string MahalFilenameSave=FilenamePA_b+"MahalTestCorr_"+FilenamePA_e+".nii.gz";
            
            
            nifti_set_filenames(MahalWM,MahalFilenameSave.c_str(),0,0);
            nifti_image_write(MahalWM);
            vector<float *> LSWMIStats=LocalSummaryStats(WMInliersNii, Quadrant, Layers, MahalWM);

            ofstream TxtFileWMI(FilenameSave.c_str());
            PrintLocalSummaryStats(LSWMIStats, TxtFileWMI);
            cout <<"Printed LS WMI "<<endl;
            for (int m=0; m<numbMahal; m++) {
                delete [] LSWMIStats[m];
                LSWMIStats[m]=NULL;
            }
            nifti_image_free(AddedWMNii);
            nifti_image_free(WMInliersNii);

        }
        Layers=ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
        float * LayersData=static_cast<float*>(Layers->data);
        float * QuadrantData=static_cast<float*>(Quadrant->data);
        numel =Layers->nvox;
        float MaxLayers=GetMax(LayersData, numel)-1;
        float MaxLobes=GetMax(QuadrantData, numel);
        cout << "Max layers is "<< MaxLayers << "Max lobes is "<< MaxLobes<<endl;
        int numbTot=MaxLayers*MaxLobes;
        float * LSResults=LocalSummaryLesion_bis(Lesion, Quadrant, Layers);
        int * LSResultsNumber=LocalSummaryLesionNumber_bis(Lesion, Quadrant, Layers);
        float * LSResultsMahal=LocalSummaryMahal_bis(Lesion, Quadrant, Layers, MahalWM);
        FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
        // int numel=Layers->nvox;

        
        //float * LobesData=static_cast<float*>(Quadrant->data);

        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameSave=FilenamePA_b+"LocalSummaryTotBis_"+FilenamePA_e+".txt";
        cout << "Saving in "<< FilenameSave << " "<<LSResultsMahal<<endl;
        ofstream TxtFileLesion(FilenameSave.c_str());
        PrintLocalSummaryTot_bis(LSResults, LSResultsNumber, LSResultsMahal, numbMahal, numbTot, TxtFileLesion);
        nifti_image_free(Lesion);
        nifti_image_free(Quadrant);
        nifti_image_free(Layers);
        delete [] LSResults;
        LSResults=NULL;
        return EXIT_SUCCESS;
    }
    
    if (segment_analysis->flag_LesionReduction && segment_analysis->filename_vectorLes.size()>0) {
        nifti_image * LesionReduced=NULL;
        nifti_image * BlurredLesionNii=NULL;
        nifti_image * NewBlurredLesionNii=NULL;
        int numbLes=segment_analysis->filename_vectorLes.size();
        for (int l=0; l<numbLes; l++) {
            nifti_image * LesNii=ReadFromFilename(segment_analysis->filename_vectorLes[l]);
            nifti_image * VentrNii=ReadFromFilename(segment_analysis->filename_vectorVentr[l]);
            vector<int> DimVector;
            for(int d=0;d<3;d++){
                DimVector.push_back(LesNii->dim[d+1]);
            }
            
            numel = LesNii->nvox;
            //            float * BlurredLesionInit=GaussianBlurring(LesData, 0, DimVector,0);
            //            nifti_image * BlurredLesInitNii=CreateNiiFromArray(BlurredLesionInit, LesNii, numel);
            //            float *BlurredLesDataInit=static_cast<float*>(BlurredLesInitNii->data);
            nifti_image * MaskNii=ReadFromFilename(segment_analysis->filename_vectorMask[l]);
            bool * MaskData=static_cast<bool*>(MaskNii->data);
            int numbRed=segment_analysis->vecLesionReduction.size();
            int sizeVO = segment_analysis->filename_vectorOut.size();
            if (sizeVO>0 && l<sizeVO){
                FilenamePA=nifti_makebasename(segment_analysis->filename_vectorOut[l]);
            }
            else{
                FilenamePA=nifti_makebasename(segment_analysis->filename_vectorLes[l]);
            }
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            int Dim[3];
            int Shift[3];
            for(int d=0;d<3;d++){
                Dim[d]=LesNii->dim[d+1];
                Shift[d]=d>0?Dim[d-1]*Shift[d-1]:1;
            }
            for (int r=0;r<numbRed;r++){
                float * LesData=static_cast<float*>(LesNii->data);
                float * VentrData=static_cast<float*>(VentrNii->data);
                int numbVentr=GetCountAbove<float, bool>(VentrData, MaskData, 0.5, numel);
                //                First Check that the max is 1 in a probabilistic segmentation
                float MaxLes=GetMax(LesData, numel);
                if (MaxLes>1) {
                    cout << "Normalisation to perform..."<<endl;
                    MultiplyFloatArrayBy(LesData, 1/MaxLes, numel);
                }
                float * LesVentrData=AddElementwise(VentrData, LesData, numel);
                
                nifti_image * LesVentrNii=CreateNiiFromArray(LesVentrData, LesNii, numel);
                stringstream rs;
                rs << r+2;
                string FilenameSaveReduced=FilenamePA_b+"Reduced_"+rs.str()+"_"+FilenamePA_e+".nii.gz";
                string FilenameSaveBlurred=FilenamePA_b+"Lesion_"+FilenamePA_e+"_"+rs.str()+".nii.gz";
                vector<float> LesionReductionVector;
                LesionReductionVector.push_back(segment_analysis->vecLesionReductionGeneral[0]);
                LesionReductionVector.push_back(segment_analysis->vecLesionReductionGeneral[1]);
                LesionReductionVector.push_back(segment_analysis->vecLesionReduction[r]);
                float Threshold=segment_analysis->vecLesionReductionGeneral[0];
                numel=LesNii->nvox;
                int InitialNumber=GetCountAbove(LesVentrData, MaskData, Threshold, numel)-numbVentr;
                int NumberOut;
                if (segment_analysis->vecLesionReductionGeneral[1]) {
                    NumberOut=LesionReductionVector[2];
                }
                else{
                    NumberOut=LesionReductionVector[2]*InitialNumber/100;
                }
                cout <<"number to take out is "<<NumberOut<<" and the threshold is "<<Threshold<<endl;
                int InitialAbove=GetCountAbove(LesVentrData, MaskData, Threshold, numel)-numbVentr;
                if (NumberOut==0) {
                    
                    nifti_set_filenames(LesNii, FilenameSaveBlurred.c_str(), 0, 0);
                    nifti_image_write(LesNii);
                }
                else if(InitialAbove<NumberOut){
                    nifti_set_filenames(LesNii, FilenameSaveBlurred.c_str(), 0, 0);
                    nifti_image_write(LesNii);
                }
                else{
                    LesionReduced=ReductionLesion(LesVentrNii, MaskNii, Threshold, NumberOut);

                    float * ReducedData=static_cast<float*>(LesionReduced->data);
                    cout<<"Initial number above threshold at "<< GetCountAbove(LesVentrData, MaskData, Threshold, numel)-numbVentr;
                    cout << "Initial reduction with new number at "<<GetCountAbove(ReducedData , MaskData, Threshold, numel)-numbVentr<<endl;
                    nifti_image_free(LesVentrNii);
                    float * ReducedMVentr=SubtractArray(ReducedData, VentrData, numel);
                    for(int i=0;i<numel;i++){
                        ReducedMVentr[i]=ReducedMVentr[i]<0?0:ReducedMVentr[i];
                    }
//                    float MinData=GetMin(ReducedData, numel);
//                    float MaxData=GetMax(ReducedData, numel);
                    float * BlurredMVentr=GaussianBlurring(ReducedMVentr, segment_analysis->vecLesionReductionGeneral[2], DimVector,0);
                    //                nifti_image * BlurredMVentrNii= CreateNiiFromArray(BlurredMVentr, LesNii, numel);
                    //                nifti_set_filenames(BlurredMVentrNii, (g_scratch_dir + "/TestBlurredMVentr.nii.gz").c_str(),  0, 0);
                    //                nifti_image_write(BlurredMVentrNii);
                    //                nifti_image_free(BlurredMVentrNii);

                    float * BlurredLesion=GaussianBlurring(ReducedData, segment_analysis->vecLesionReductionGeneral[2], DimVector,0);
                    //                BlurredMVentrNii= CreateNiiFromArray(BlurredLesion, LesNii, numel);
                    //                nifti_set_filenames(BlurredMVentrNii, (g_scratch_dir + "/TestBlurredVentr.nii.gz").c_str(),  0, 0);
                    //                nifti_image_write(BlurredMVentrNii);
                    //                nifti_image_free(BlurredMVentrNii);

                    for (int i=0; i<numel; i++) {
                        BlurredLesion[i]=(BlurredLesion[i]-BlurredMVentr[i])>0.1?BlurredMVentr[i]:BlurredLesion[i];
                        BlurredLesion[i]=VentrData[i]>0.5?0:BlurredLesion[i];
                    }
                    delete [] BlurredMVentr;
                    delete [] ReducedMVentr;
                    //                BlurredMVentrNii= CreateNiiFromArray(BlurredLesion, LesNii, numel);
                    //                nifti_set_filenames(BlurredMVentrNii, (g_scratch_dir + "/TestBlurredVentr2.nii.gz").c_str(),  0, 0);
                    //                nifti_image_write(BlurredMVentrNii);
                    //                nifti_image_free(BlurredMVentrNii);
                    //                float * BlurredLesionPVentr=AddElementwise(VentrData, BlurredLesion, numel);
                    //                for (int i=0; i<numel; i++) {
                    //                    BlurredLesionPVentr[i]=BlurredLesionPVentr[i]>1?1:BlurredLesionPVentr[i];
                    //                }
                    cout << "Intermediate reduction with new number at "<<GetCountAbove(BlurredLesion,MaskData, Threshold, numel)<<endl;
                    BlurredLesionNii=CreateNiiFromArray(BlurredLesion, LesNii, numel);
                    int NumberNeeded=InitialNumber-NumberOut;
                    NumberNeeded=NumberNeeded>0?NumberNeeded:0;
//                    MinData=GetMin(BlurredLesion, numel);
//                    MaxData=GetMax(BlurredLesion, numel);
                    //                Perform the transformation to have the appropriate cut off at the needed threshold
                    float * NewBlurredTransform=NULL;
                    if (NumberNeeded>0) {
                        
                        float TFTV=GetThresholdForVolume(BlurredLesion, MaskData, NumberNeeded,numel);
                        vector<vector<float> > AffineCoeffs;
                        vector<float> Affine1;
                        vector<float> Affine2;
                        //                float c=0;
                        //                float b=2*(TFTV*TFTV-Threshold)/(TFTV*(TFTV*TFTV-1));
                        //                float a=Threshold/TFTV-b/TFTV;
                        //                float B=0;
                        //                float A=(Threshold-1)/(TFTV*TFTV-1);
                        //                float C=1-A;
                        //                float c=0;
                        float b=Threshold/TFTV;
                        float a=0;
                        //                a=Threshold/(TFTV*TFTV);
                        //                b=0;
                        float Test=2*(1-Threshold)/((TFTV-1)*(TFTV-1));

//                        float B=(2*(Threshold-TFTV*TFTV)/(TFTV*(TFTV-1)*(TFTV-1)));
                        float B=Test;
                        float A=(Threshold-1)/(TFTV*TFTV-1)-B/(TFTV+1);
                        float C=1-A-B;
                        B=(1-Threshold)/(1-TFTV);
                        C=1-B;
                        A=0;

                        //                float ValueDecrease=-B*(TFTV+1)/(2*(Threshold-1)/(TFTV-1)-2*B);
                        //
                        //                float Testb=2*(Threshold-TFTV*TFTV)/(TFTV*(TFTV*TFTV-1))+(Test)*(1-TFTV)/(1+TFTV);
                        //                float RangeBMin=(Threshold-1)/(TFTV-1);
                        //                B=Test;
                        //                b=Testb;
                        //                A=(Threshold-1)/(TFTV*TFTV-1)-B/(TFTV+1);


                        //                float c=0;
                        //                float b=0;
                        //                float a=Threshold/(TFTV*TFTV);
                        //                float B=(Threshold-1)/(TFTV-1.0/2-TFTV*TFTV/2);
                        //                float Test1=(TFTV-1.0/2-TFTV*TFTV/2);
                        //                float A=-B/2;
                        //                float C=1-A-B;
                        //                float Tan1=2*a*TFTV+b;
                        //                float Tan2=2*A*TFTV+B;



                        Affine1.push_back(0);
                        Affine1.push_back(b);
                        Affine1.push_back(a);
                        Affine2.push_back(C);
                        Affine2.push_back(B);
                        Affine2.push_back(A);
                        AffineCoeffs.push_back(Affine1);
                        AffineCoeffs.push_back(Affine2);


                        //                AffineCoeffs.push_back(0);
                        //                float c=Threshold/(2*ThresholdForTrueVolume);
                        //                float a=(Threshold/2-ThresholdForTrueVolume)/(ThresholdForTrueVolume*ThresholdForTrueVolume-ThresholdForTrueVolume);
                        //                float b=1-c-a;
                        //                AffineCoeffs.push_back(c);
                        //                AffineCoeffs.push_back(b);
                        //                AffineCoeffs.push_back(a);
                        vector<float> TFTVVec;
                        TFTVVec.push_back(TFTV);
                        NewBlurredTransform=ApplyingPolyfitPiecewise(AffineCoeffs,BlurredLesion,TFTVVec,numel);
                        //                for (int i=0; i<numel; i++) {
                        //                    NewBlurredTransform[i]=VentrData[i]>0.5?0:NewBlurredTransform[i];
                        //                }
                    }
                    else{
                        NewBlurredTransform=TranscribeArray<float, float>(BlurredLesion, numel);
                    }
//                    MinData=GetMin(NewBlurredTransform, numel);
//                    MaxData=GetMax(NewBlurredTransform, numel);
                    NewBlurredLesionNii=CreateNiiFromArray(NewBlurredTransform, LesNii, numel);
                    cout << "Final reduction with new number at "<<GetCountAbove(NewBlurredTransform,MaskData, Threshold, numel)<<endl;

                    nifti_set_filenames(NewBlurredLesionNii, FilenameSaveBlurred.c_str(), 0, 0);
                    nifti_set_filenames(LesionReduced, FilenameSaveReduced.c_str(), 0, 0);
                    nifti_image_write(NewBlurredLesionNii);
                    nifti_image_write(LesionReduced);
                    delete [] NewBlurredTransform;
                    NewBlurredTransform=NULL;
                    if (segment_analysis->vecLesionReductionGeneral[3]) {
                        nifti_image_free(LesNii);
                        LesNii=NewBlurredLesionNii;
                    }
                    else{
                        nifti_image_free(NewBlurredLesionNii);
                    }
                    //                nifti_image_free(LesNii);
                    //                LesNii=BlurredLesionNii;

                    nifti_image_free(LesionReduced);
                    nifti_image_free(BlurredLesionNii);
                    LesionReduced=NULL;
                    BlurredLesionNii=NULL;
                }
            }
            nifti_image_free(MaskNii);
            nifti_image_free(LesNii);
            //            nifti_image_free(BlurredLesionNii);
            //            delete [] BlurredLesionInit;
            
        }
        return EXIT_SUCCESS;
    }
    
    if (segment_analysis->flag_Evaluation && segment_analysis->vecThresh.size()>0) {
        cout<<"vecThresh is "<<segment_analysis->vecThresh[2]<<endl;
        cout <<"Probably doing the range analysis"<<endl;
        int numbRef=segment_analysis->filename_vectorRef.size();
        int  numbNames=segment_analysis->filename_vectorNames.size();
        int numbMask=segment_analysis->filename_vectorMask.size();
        if(numbRef!=numbNames || numbNames!=numbRef || numbMask!=numbRef){
            cout<<"Processing evaluation impossible"<<endl;
            return EXIT_SUCCESS;
        }
        else{
            cout << "Properly doing the analysis"<<endl;
            ofstream TxtFileReportEval(segment_analysis->filename_Out);
            TxtFileReportEval<<"Name,";
            int sizeResults=(int)segment_analysis->vecThresh[2]+1;
            float step_bin=fabs(segment_analysis->vecThresh[1]-segment_analysis->vecThresh[0])/segment_analysis->vecThresh[2];
            for (int s=0; s<sizeResults; s++) {
                TxtFileReportEval<<s*step_bin+segment_analysis->vecThresh[0]<<",";
            }
            TxtFileReportEval<<endl;
            for(int i=0;i<numbRef;i++){
                nifti_image * RefNii=ReadFromFilename(segment_analysis->filename_vectorRef[i]);
                nifti_image * MaskNii=ReadFromFilename(segment_analysis->filename_vectorMask[i]);
                float * VolRef=VolumeRange(RefNii,segment_analysis->vecThresh,MaskNii);
                nifti_image_free(RefNii);
                nifti_image_free(MaskNii);
                TxtFileReportEval<<segment_analysis->filename_vectorNames[i]<<",";
                PrintThreshRange(VolRef,segment_analysis->vecThresh, TxtFileReportEval);
                delete [] VolRef;
                VolRef=NULL;
                RefNii=NULL;
                MaskNii=NULL;
            }
        }
        return EXIT_SUCCESS;
    }
    if(segment_analysis->flag_Evaluation && segment_analysis->filename_vectorLes.size()>0){
        cout<<"Performing evaluation"<<endl;
        int numbRef=segment_analysis->filename_vectorRef.size();
        int numbLes=segment_analysis->filename_vectorLes.size();
        int  numbNames=segment_analysis->filename_vectorNames.size();
        int numbMask=segment_analysis->filename_vectorMask.size();
        if(numbRef!=numbLes || numbLes!=numbNames || numbNames!=numbRef || numbMask!=numbRef){
            cout<<"Processing evaluation impossible"<<endl;
            return EXIT_SUCCESS;
        }
        else{
            ofstream TxtFileReportEval(segment_analysis->filename_Out);
            for(int i=0;i<numbRef;i++){
                nifti_image * LesNii=ReadFromFilename(segment_analysis->filename_vectorLes[i]);
                nifti_image * RefNii=ReadFromFilename(segment_analysis->filename_vectorRef[i]);
                nifti_image * MaskNii=ReadFromFilename(segment_analysis->filename_vectorMask[i]);
                EvaluationReport * Eval=CreateEvaluationReport(LesNii,RefNii,MaskNii);
                nifti_image_free(LesNii);
                nifti_image_free(RefNii);
                nifti_image_free(MaskNii);
                TxtFileReportEval<<segment_analysis->filename_vectorNames[i]<<",";
                PrintEvalReport(Eval,TxtFileReportEval);
                delete Eval;
                LesNii=NULL;
                RefNii=NULL;
                MaskNii=NULL;
            }
        }
        return EXIT_SUCCESS;
    }
    

    if (segment_analysis->flag_TextFile) {
        cout<<"Building tree from text file"<<endl;
        if (Modalities.size()==0) {
            Modalities = GetModalitiesFromTextFile(segment_analysis);
        }
        TreeToAnalyse=ReadTreeFromTextFile(segment_analysis->filename_TextFile, NULL,segment_analysis->filename_changePath);
        SEG_PARAMETERS * segment_paramAnalyse=new SEG_PARAMETERS[1]();
        segment_paramAnalyse->IndexCSF=segment_analysis->IndexCSF;
        segment_paramAnalyse->IndexGM=segment_analysis->IndexGM;
        segment_paramAnalyse->IndexOut=segment_analysis->IndexOut;
        segment_paramAnalyse->IndexWM=segment_analysis->IndexWM;

        numbmodal=TreeToAnalyse->GetNumberModalities();
        
        //        float * DataTest=static_cast<float *>(TreeToAnalyse->GetDataImage()->data);
        //        cout<<"Data is at 148 119 120 "<<DataTest[7894932]<<" "<<DataTest[7894932+numel]<<endl;
        nifti_image * SegToUse=ReadFromFilename(segment_analysis->filename_SegTot);
        cout<<"SegToUse read"<<endl;
        // Check if number of elements in SegToUse correspond to number of leaves in TreeToAnalyse
        VectorLeaves=TreeToAnalyse->GetAllLeaves();
        numbleavesTree=VectorLeaves.size();
        int numbleavesSeg=SegToUse->nu*SegToUse->nt;
        if (numbleavesSeg!=numbleavesTree) {
            cout<<"Incompatible text file and segmentation"<<endl;
        }
        else {
            float * SegData=static_cast<float *>(SegToUse->data);
            float * SegData_PTR=NULL;
            int * L2S_PTR=NULL;
            int numel=TreeToAnalyse->GetNumberElements();
            int numelmasked=TreeToAnalyse->GetNumberMaskedElements();
            cout<<"ok for number of elements"<<endl;
            for (int l=0; l<numbleavesTree; l++) {
                L2S_PTR=TreeToAnalyse->GetL2S();
                SegData_PTR=&SegData[l*numel];
                float * NormRespNew=new float[numelmasked];
                for (int i=0; i<numelmasked; i++) {
                    NormRespNew[i]=0;
                }
                float * NormRespNew_PTR=NormRespNew;
                for (int i=0; i<numel; i++, L2S_PTR++,SegData_PTR++) {
                    if (*L2S_PTR>=0) {
                        *NormRespNew_PTR=*SegData_PTR;
                        NormRespNew_PTR++;
                    }
                }
                cout<<"NormResp new calculated for "<<l<<"...";
                VectorLeaves[l]->SetNormResp(NormRespNew);
                delete [] NormRespNew;
                NormRespNew=NULL;
                cout<<"...and set";
            }
            cout<<"TreeToAnalyse rebuilt and saved..."<<endl;
            vector<int> Modalities2 = GetModalitiesFromTextFile(segment_analysis);
            if (Modalities2.size()<=Modalities.size()) {
                numbmodal=Modalities2.size();
            }
            if (Modalities.size()==0) {
                Modalities=GetModalitiesFromTextFile(segment_analysis);
                numbmodal=Modalities.size();
            }
            if (TreeToAnalyse->GetNumberModalities()>numbmodal) {
                CorrectDataAccordingToModa(TreeToAnalyse, Modalities);
            }
            
            RebuildNormRespFromLeaves(TreeToAnalyse);
            //            TreeToAnalyse->SaveBFCorrectedData(g_scratch_dir + "/TestCorrData.nii.gz");
            //            TreeToAnalyse->SaveAllClasses(g_scratch_dir + "/TestClassesBegin.nii.gz", segment_paramAnalyse);
//            numel=TreeToAnalyse->GetNumberElements();
            cout<<"Tree completely built"<<endl;
            if (segment_analysis->flag_juxtaCorrection) {
                segment_paramAnalyse->flag_DGMPrior=1;
                string CopyDGMPrior=segment_analysis->filename_inPriorsDGM;
                segment_paramAnalyse->filename_DGMPrior=CopyDGMPrior;
                cout<<"Trying correction for juxtacortical elements in reconstruction of lesions"<<endl;
                numbmodal=TreeToAnalyse->GetNumberModalities();
                for (int m=0; m<numbmodal; m++) {
                    segment_paramAnalyse->Modalities.push_back(Modalities[m]);
                    cout << "Modalities "<< m<<" "<< segment_paramAnalyse->Modalities[m]<< endl;
                }
                //                TreeToAnalyse->SaveBFCorrectedData(g_scratch_dir + "/TestCorrData.nii.gz");
                //                TreeToAnalyse->SaveAllClasses(g_scratch_dir + "/TestClassesBegin.nii.gz", segment_paramAnalyse);
                segment_paramAnalyse->flag_out=0;
                segment_paramAnalyse->flag_inDCFile=0;
                segment_paramAnalyse->flag_inDC=0;
                CorrectionJuxta=TreeToAnalyse->CorrectInitialEMResultForWrongWMHoles(segment_paramAnalyse);

                cout<<".... juxtacortical correction performed !"<<endl;
                //                nifti_image * CorrectionJuxta2=HolesJuxta(TreeToAnalyse, segment_analysis);
                //                float * CorrJuxtaData=NULL;
                //                float * CorrJuxta2Data=NULL;
                //
                //                if (CorrectionJuxta!=NULL) {
                //                    CorrJuxtaData=static_cast<float*>(CorrectionJuxta->data);
                //                    if (CorrectionJuxta2!=NULL) {
                //                        CorrJuxta2Data=static_cast<float*>(CorrectionJuxta2->data);
                //                        for (int i=0; i<numel; i++) {
                //                            CorrJuxtaData[i]+=CorrJuxta2Data[i];
                //                        }
                //                    }
                //                }
                //                nifti_image_free(CorrectionJuxta2);
                //                CorrectionJuxta2=NULL;
                
            }

//            numel=TreeToAnalyse->GetNumberElements();
            //            SEG_PARAMETERS * segment_param=NULL;
            //            TreeToAnalyse->SaveAllClasses(g_scratch_dir + "/StrangeRebuilt.nii.gz", segment_param);
        }
        nifti_image_free(SegToUse);
        SegToUse=NULL;
        //        if (segment_paramAnalyse!=NULL) {
        //            delete segment_paramAnalyse;
        //            segment_paramAnalyse=NULL;
        //        }

    }

    nifti_image * FinalLes=NULL;

    if(TreeToAnalyse!=NULL && segment_analysis->flag_LesWMI>0 && (segment_analysis->flag_Parc || segment_analysis->flag_inArtefact) && segment_analysis->flag_GIFPriors){
        Rule * LesionRule;
        if (segment_analysis->flag_RuleTextFile) {
            LesionRule=BuildRuleFromTextFile(TreeToAnalyse, segment_analysis);
        }

        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveBorderWM=FilenamePA_b+"BorderWMDGM_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLesOut=FilenamePA_b+"LesOut_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLesOutb=FilenamePA_b+"LesOutb_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLesOutc=FilenamePA_b+"LesOutc_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLesWMI=FilenamePA_b+"LesWMI_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLesCSFI=FilenamePA_b+"LesCSFI_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLesAdded=FilenamePA_b+"LesTotAdded_"+FilenamePA_e+".nii.gz";
        string FilenameSaveArt=FilenamePA_b+"ArtBuilt_"+FilenamePA_e+".nii.gz";
        string FilenameSaveGM=FilenamePA_b+"GMBuilt_"+FilenamePA_e+".nii.gz";

        string FilenameSaveMahal=FilenamePA_b+"MahalBuilt_"+FilenamePA_e+".nii.gz";
        vector<TreeEM *> LesionClassesNew;
        numel=TreeToAnalyse->GetNumberElements();
        numbmodal=TreeToAnalyse->GetNumberModalities();
        float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
        cout <<"Data corr obtained"<<endl;
        float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        float * GMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        float * CorrWMTemp=GetMeanDataMulti(DataCorr,WMTempInliers,numel,numbmodal);
        float * CorrGMTemp=GetMeanDataMulti(DataCorr,GMTempInliers,numel,numbmodal);
        // Try to find FLAIR index
        int IndexFLAIR;
        std::vector<int>::iterator it;
        it = find (Modalities.begin(), Modalities.end(), 3);
        if (it!=Modalities.end()){
            IndexFLAIR=std::distance(Modalities.begin(), it)+1;
        }
        else{
            IndexFLAIR=2;
        }
        if(CorrWMTemp[IndexFLAIR-1]<CorrGMTemp[IndexFLAIR-1]){
            cout << "No WMI correction needed"<<endl;
            delete [] DataCorr;
            delete [] WMTempInliers;
            delete [] GMTempInliers;
            delete [] CorrWMTemp;
            delete [] CorrGMTemp;
        }
        else{
            cout<<"Inclusion from WMI directly needed "<<endl;
            nifti_image * WMInliersNii=HardSegmentationThresholdFromNormResp(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(),TreeToAnalyse,0.5);
            nifti_image * GMInliersNii=CreateNiiFromArray(GMTempInliers,WMInliersNii,numel);
            nifti_set_filenames(GMInliersNii,FilenameSaveGM.c_str(),0,0);
            nifti_image_write(GMInliersNii);

            //CreateNiiFromArray(WMTempInliers, TreeToAnalyse->GetDataImage(), TreeToAnalyse->GetNumberElements());
            nifti_image * MaskNii=CreateNiiFromArray(static_cast<bool*>(TreeToAnalyse->GetMask()->data),TreeToAnalyse->GetDataImage(),numel);

            nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
            cout << "Data Nii"<<DataNii<<endl;
            // First weighted
            nifti_image * LesionTotOut=LesionVoxelwise(TreeToAnalyse,LesionRule,Modalities,segment_analysis);
            cout << "LesionVoxelwise"<<endl;
            nifti_image * SegToAnalyse=ReconstructLesionImageTotWeighted(TreeToAnalyse, LesionClassesNew, LesionTotOut,Modalities,segment_analysis,NULL);
            cout <<"SegToAnalyse"<<endl;
            nifti_set_filenames(SegToAnalyse,FilenameSaveLesOutc.c_str(),0,0);
            nifti_image_write(SegToAnalyse);
            nifti_set_filenames(LesionTotOut,FilenameSaveLesOutb.c_str(),0,0);
            nifti_image_write(LesionTotOut);
            nifti_image * PriorsGIF=ReadFromFilename(segment_analysis->filename_GIFPriors);
            nifti_image* CSFILesNii=NULL;
            float * DataPriorGIF=static_cast<float *>(PriorsGIF->data);
            float * GIFCSF=&DataPriorGIF[numel];
            float * GIFCGM=&DataPriorGIF[2*numel];
            float * GIFWM=&DataPriorGIF[3*numel];
            float * GIFDGM=&DataPriorGIF[4*numel];
            float * GIFBrainstem=&DataPriorGIF[5*numel];



            cout<<"Correcting WMI seg now !"<<WMInliersNii<<" "<<GMInliersNii<<endl;
            CorrectWMI(TreeToAnalyse,WMInliersNii,GMInliersNii,segment_analysis);

            nifti_image * MahalNii=MahalDistMaps(WMInliersNii,MaskNii,DataNii);
            nifti_set_filenames(MahalNii,FilenameSaveMahal.c_str(),0,0);
            nifti_image_write(MahalNii);
            nifti_image_free(DataNii);
            DataNii=NULL;
            float * MahalData=static_cast<float *>(MahalNii->data);
            float * MahalToUse=&MahalData[numel*IndexFLAIR];
            nifti_image * PriorsDGM=NULL;
            float * PriorsDGMData=NULL;
            if(segment_analysis->flag_inPriorsDGM){
                PriorsDGM=ReadFromFilename(segment_analysis->filename_inPriorsDGM);
                PriorsDGMData=static_cast<float*>(PriorsDGM->data);
            }
            float * PriorsWMI=NULL;
            nifti_image * PriorsWM=NULL;
            if(!segment_analysis->flag_inPriorsWM){
                PriorsWMI=new float[numel];
                for(int i=0;i<numel;i++){
                    PriorsWMI[i]=0;
                }
            }
            else{
                PriorsWM=ReadFromFilename(segment_analysis->filename_inPriorsWM);
                PriorsWMI=static_cast<float*>(PriorsWM->data);
            }
            // Then corrected with Parc or Artefacts
            nifti_image * ArtefactNii=NULL;
            bool * ArtefactBool=NULL;
            if(segment_analysis->flag_inArtefact){
                ArtefactNii=ReadFromFilename(segment_analysis->filename_Artefact);
                float * ArtefactData=static_cast<float *>(ArtefactNii->data);
                ArtefactBool=TranscribeArray<float,bool>(ArtefactData,numel);
            }
            else{
                nifti_image * ParcellationNii=ReadFromFilename(segment_analysis->filename_Parc);
                int * ParcellationData=TranscribeArray<float,int>(static_cast<float*>(ParcellationNii->data),numel);
                static const int ArtefactArray[]={5, 12, 16, 24, 31, 32, 33, 39, 40, 50, 51, 62, 63, 64, 65, 72, 73, 74, 48, 49, 101, 102, 117, 118, 139, 140, 171, 172, 187, 188, 167, 168};
                vector<int> VecArtefact (ArtefactArray, ArtefactArray + sizeof(ArtefactArray) / sizeof(ArtefactArray[0]) );
                ArtefactBool=new bool[numel];

                for(int i=0;i<numel;i++){
                    ArtefactBool[i]=0;
                    if(GIFDGM[i]+GIFWM[i]<0.08){
                        ArtefactBool[i]=1;
                    }
                    if(it==Modalities.end() && Modalities.size()==2 && GIFCGM[i]>0.8){
                        ArtefactBool[i]=1;
                    }
                    if(PriorsDGM!=NULL){
                        if(PriorsDGMData[i]>0.5 && MahalToUse[i]<segment_analysis->weightThreshold){
                            ArtefactBool[i]=1;
                        }
                    }
                    if(PriorsWM!=NULL){
                        if(PriorsWMI[i]<0.1 && GIFCGM[i]>0.5){
                            ArtefactBool[i]=1;
                        }
                    }

                }
                int numArt=VecArtefact.size();
                cout<<numArt<<" numArt"<<endl;
                for(int s=0;s<numArt;s++){
                    bool * ArtBool=CreateLesionBool(ParcellationData,VecArtefact[s],numel);
                    cout<<CountNonZero(ArtBool,numel)<<endl;
                    OROperationBool(ArtBool,ArtefactBool,ArtefactBool,numel);
                    cout<<CountNonZero(ArtefactBool,numel)<<endl;
                    delete [] ArtBool;
                    ArtBool=NULL;
                }
                ArtefactNii=CreateNiiFromArray(ArtefactBool,LesionTotOut,numel);
                nifti_set_filenames(ArtefactNii,FilenameSaveArt.c_str(),0,0);
                nifti_image_write(ArtefactNii);
                nifti_image_free(ParcellationNii);
                delete [] ParcellationData;
            }
            nifti_image * CorrLesTotOut=CorrectionNii(SegToAnalyse,ArtefactBool);
            nifti_set_filenames(CorrLesTotOut,FilenameSaveLesOut.c_str(),0,0);
            nifti_image_write(CorrLesTotOut);

            int Dim[3];
            int Shift[3];
            Shift[0]=1;
            for(int d=0;d<3;d++){
                Dim[d]=CorrLesTotOut->dim[d+1];
                if (d>0){
                    Shift[d]=Dim[d-1]*Shift[d-1];
                }
            }
            // Then from WMI

            float * WeightedWMISeg=new float[numel];
            bool * WMThresh=ThresholdArray<float,bool>(WMTempInliers,0.5,numel);
            bool * DGMThresh=ThresholdArray<float,bool>(GIFDGM,0.4,numel);
            float * LesTotData=static_cast<float *>(CorrLesTotOut->data);
            bool * LesThresh=ThresholdArray<float,bool>(LesTotData,0.5,numel);
            OROperationBool(WMThresh,DGMThresh,WMThresh,numel);
            OROperationBool(WMThresh,LesThresh,WMThresh,numel);
            bool * BorderWM=CreateBorderFromBool(WMThresh,Dim,Shift);
            nifti_image * BorderNii=CreateNiiFromArray(BorderWM,MahalNii,numel);
            nifti_set_filenames(BorderNii,FilenameSaveBorderWM.c_str(),0,0);
            nifti_image_write(BorderNii);
            nifti_image_free(BorderNii);





            BorderNii=NULL;
            for(int i=0;i<numel;i++){
                float Weighted=MahalToUse[i]/segment_analysis->flag_LesWMI;
                if (PriorsDGM!=NULL && PriorsDGMData[i]>0.5){
                    Weighted/=2;
                }
                else if (GIFCGM[i]>0.01){
                    Weighted=MahalToUse[i]/(min(2*segment_analysis->flag_LesWMI*GIFCGM[i]+segment_analysis->flag_LesWMI,(float)2.5*segment_analysis->flag_LesWMI));
                }
                if (GIFDGM[i]>0.01 && Weighted>MahalToUse[i]/(1.5*segment_analysis->flag_LesWMI)){
                    Weighted=MahalToUse[i]/(min(1.5*segment_analysis->flag_LesWMI*GIFDGM[i]+segment_analysis->flag_LesWMI,2.5*segment_analysis->flag_LesWMI));
                }
                WeightedWMISeg[i]=Weighted>1?WMTempInliers[i]:WMTempInliers[i]*Weighted;
                if( GMTempInliers[i]>0.4 && (GIFWM[i]+GIFDGM[i]>0.9 || PriorsWMI[i]>0.5)){
                    WeightedWMISeg[i]+=Weighted>1?GMTempInliers[i]:GMTempInliers[i]*Weighted;
                    WeightedWMISeg[i]=WeightedWMISeg[i]>1?1:WeightedWMISeg[i];
                }
                if(BorderWM[i] && ( MahalToUse[i]<segment_analysis->weightThreshold)){
                    WeightedWMISeg[i]=0;
                }
                if(GIFCGM[i]>0.5 && ( MahalToUse[i]<5)){
                    WeightedWMISeg[i]=0;
                }
                if(it==Modalities.end() && BorderWM[i]&& GIFCGM[i]>0.8){
                    WeightedWMISeg[i]=0;
                }
                if(MahalToUse[i]<0){
                    WeightedWMISeg[i]=0;
                }
                if(2*segment_analysis->flag_LesWMI*GIFCGM[i]+segment_analysis->flag_LesWMI>MahalToUse[i] && GIFCGM[i]>0.1){
                    WeightedWMISeg[i]=0;
                }
                WeightedWMISeg[i]=WeightedWMISeg[i]>1?1:WeightedWMISeg[i];

            }
            nifti_image* WMILesNii=CreateNiiFromArray(WeightedWMISeg,SegToAnalyse,numel);
            nifti_image * CorrWMILesNii=CorrectionNii(WMILesNii,ArtefactBool);
            nifti_set_filenames(CorrWMILesNii,FilenameSaveLesWMI.c_str(),0,0);
            nifti_image_write(CorrWMILesNii);

            // Then from CSFI with higher wm priors but only if no FLAIR images.
            float * WeightedCSFISeg=new float[numel];
            if (it==Modalities.end()){ // only if there is no FLAIR modality
                float * CSFTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexCSF)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
                for(int i=0;i<numel;i++){
                    if (GIFWM[i]+GIFBrainstem[i]>GIFCSF[i] && GIFCGM[i]<0.1){
                        WeightedCSFISeg[i]=MahalToUse[i]/segment_analysis->flag_LesWMI;
                        WeightedCSFISeg[i]=WeightedCSFISeg[i]>1?CSFTempInliers[i]:CSFTempInliers[i]*WeightedCSFISeg[i];
                    }
                    else{
                        WeightedCSFISeg[i]=0;
                    }
                    if(MahalToUse[i]<0){
                        WeightedCSFISeg[i]=0;
                    }
                }
                delete [] CSFTempInliers;
                CSFILesNii=CreateNiiFromArray(WeightedCSFISeg,SegToAnalyse,numel);
            }
            // Then from partial volume CSF GM with low CGM and CSF priors


            nifti_image * CorrCSFILesNii=NULL;
            if(it==Modalities.end()){
                CorrCSFILesNii=CorrectionNii(CSFILesNii,ArtefactBool);
            }
            nifti_set_filenames(CorrCSFILesNii,FilenameSaveLesCSFI.c_str(),0,0);
            nifti_image_write(CorrCSFILesNii);
            vector<nifti_image*> PartsLesVec;
            PartsLesVec.push_back(CorrCSFILesNii);
            PartsLesVec.push_back(CorrLesTotOut);
            PartsLesVec.push_back(CorrWMILesNii);
            nifti_image * AddedLesion=AddNii<float>(PartsLesVec);
            float * AddedFloat = static_cast<float*>(AddedLesion->data);
            for(int i=0;i<numel;i++){
                AddedFloat[i]=AddedFloat[i]>1?1:AddedFloat[i];
            }
            int *ComponentLabels=ComponentLabeling(AddedLesion,6,0.75);
            int * VolThresh=OrderedVolumeLabel(ComponentLabels,3,numel,1);
            bool * VolRemain=TranscribeArray<int,bool>(VolThresh,numel);








            float * RemainLesion=MultiplyElementwiseChoose<float,bool,float>(static_cast<float*>(AddedLesion->data),VolRemain,numel);
            FinalLes=CreateNiiFromArray(RemainLesion,AddedLesion,numel);
            nifti_set_filenames(FinalLes,FilenameSaveLesAdded.c_str(),0,0);
            nifti_image_write(FinalLes);
            nifti_image_free(AddedLesion);
            nifti_image_free(CSFILesNii);
            nifti_image_free(WMILesNii);
            nifti_image_free(CorrLesTotOut);
            nifti_image_free(CorrWMILesNii);
            nifti_image_free(SegToAnalyse);
            nifti_image_free(LesionTotOut);
            nifti_image_free(PriorsGIF);
            nifti_image_free(MahalNii);
            nifti_image_free(MaskNii);
            nifti_image_free(WMInliersNii);
            delete [] ComponentLabels;
            delete [] VolThresh;
            delete [] VolRemain;
            delete [] RemainLesion;
            delete [] WMThresh;
            delete [] DGMThresh;
            delete [] LesThresh;
            delete [] BorderWM;
            delete [] ArtefactBool;
            delete [] WMTempInliers;
            delete [] GMTempInliers;
            nifti_image_free(ArtefactNii);
            if(WeightedCSFISeg!=NULL){
                delete [] WeightedCSFISeg;
            }
            delete [] WeightedWMISeg;
        }

        if(!segment_analysis->flag_correct){
            delete LesionRule;
            delete TreeToAnalyse;
            return EXIT_SUCCESS;
        }
    }
    
    cout<< segment_analysis->flag_LaplaceSolImage << segment_analysis->flag_inLesCorr << segment_analysis->flag_inQuadrant << segment_analysis->flag_LS<<endl;
    
    if(segment_analysis->flag_inLesCorr && segment_analysis->flag_refLes && TreeToAnalyse!=NULL ){
        // First Get Mahal Map
        nifti_image * LesionCorr=ReadFromFilename(segment_analysis->filename_inLesCorr);
        numel=TreeToAnalyse->GetNumberElements();
        numbmodal=TreeToAnalyse->GetNumberModalities();
        nifti_image * RefLes=ReadFromFilename(segment_analysis->filename_Ref);
        bool * GTLes=ThresholdArray<float, bool>(static_cast<float *>(RefLes->data), 0.5, numel);
        bool * BoolLes=ThresholdArray<float, bool>(static_cast<float*>(LesionCorr->data), 0.5, numel);
        bool * LesOppose=OpposeBoolArray(BoolLes, numel);
        bool * GTOppose=OpposeBoolArray(GTLes, numel);
        cout << "Image lesion Read"<<endl;
        float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        float * WMTempOutliers=CreateLong(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        
        cout << "Temp seg created "<< endl;
        float * LesionCorrData=static_cast<float*>(LesionCorr->data);
        vector<float *> AdditionWM;
        AdditionWM.push_back(WMTempInliers);
        AdditionWM.push_back(WMTempOutliers);
        AdditionWM.push_back(LesionCorrData);
        float * AddedWM=AddArray(AdditionWM, TreeToAnalyse->GetNumberElements());
        cout << "Addition done"<< endl;
        nifti_image * AddedWMNii=CreateNiiFromArray(AddedWM, LesionCorr, TreeToAnalyse->GetNumberElements());
        MultiplyElementwiseChange(WMTempInliers, GTOppose, numel);
        nifti_image * WMInliersNii=CreateNiiFromArray(WMTempInliers, LesionCorr, TreeToAnalyse->GetNumberElements());
        
        float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
        cout <<"Data corr obtained"<<endl;
        nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
        cout << "Data Nii"<<DataNii<<endl;
//        nifti_image * MahalWM= MahalDistMaps(WMInliersNii, AddedWMNii, TreeToAnalyse->GetDataImage() );
        cout <<"Mahal 1 obtained"<<endl;
        cout << "MahalWM"<<endl;
        nifti_image * MahalWM = MahalDistMaps(WMInliersNii, AddedWMNii, DataNii );
        float * MahalData=static_cast<float*>(MahalWM->data);
        nifti_image_free(DataNii);
        delete [] DataCorr;
        cout <<"Mahal obtained "<<endl;
        numbmodal=MahalWM->dim[4];
        cout << numbmodal <<endl;
        FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameHistogram=FilenamePA_b+"Histogram_"+FilenamePA_e+".txt";
        string FilenameHistogramGT=FilenamePA_b+"HistogramGT_"+FilenamePA_e+".txt";
        string FilenameStat=FilenamePA_b+"Stats"+FilenamePA_e+".txt";
        string FilenameStatGT=FilenamePA_b+"StatsGT"+FilenamePA_e+".txt";
        string FilenameFP=FilenamePA_b+"FP"+FilenamePA_e+".nii.gz";
        string FilenameFN=FilenamePA_b+"FN"+FilenamePA_e+".nii.gz";
        string FilenameTP=FilenamePA_b+"TP"+FilenamePA_e+".nii.gz";
        ofstream TxtHistogram(FilenameHistogram.c_str());
        ofstream TxtStat(FilenameStat.c_str());
        ofstream TxtHistogramGT(FilenameHistogramGT.c_str());
        ofstream TxtStatGT(FilenameStatGT.c_str());
        //        Create Threshold over GT

        bool * FPMap=new bool[numel];
        bool * FNMap=new bool[numel];
        bool * TPMap=new bool[numel];
        
        // Then Create FP Map and FN Map

        ANDOperationBool(BoolLes, GTLes, TPMap, numel);
        cout << "TPMap "<<endl;
        ANDOperationBool(BoolLes, GTOppose, FPMap, numel);
        cout << "FPMap "<<endl;
        ANDOperationBool(GTLes, LesOppose, FNMap, numel);
        cout << "FNMap "<<endl;
        float * FPVal=MultiplyElementwiseChoose<float, bool, float>(static_cast<float*>(LesionCorr->data), FPMap, numel);
        float * FNVal=MultiplyElementwiseChoose<float, bool, float>(static_cast<float*>(LesionCorr->data), FPMap, numel);
        float * TPVal=MultiplyElementwiseChoose<float, bool, float>(static_cast<float*>(LesionCorr->data), TPMap, numel);
        nifti_image * FPNii=CreateNiiFromArray(FPVal, LesionCorr, numel);
        nifti_image * FNNii=CreateNiiFromArray(FNVal, LesionCorr, numel);
        nifti_image * TPNii=CreateNiiFromArray(TPVal, LesionCorr, numel);
        nifti_set_filenames(FPNii, FilenameFP.c_str(), 0, 0);
        nifti_set_filenames(FNNii, FilenameFN.c_str(), 0, 0);
        nifti_set_filenames(TPNii, FilenameTP.c_str(), 0, 0);
        nifti_image_write(FPNii);
        nifti_image_write(FNNii);
        nifti_image_write(TPNii);
        nifti_image_free(FNNii);
        nifti_image_free(FPNii);
        nifti_image_free(TPNii);
        cout << "Images written"<<endl;
        
        // Then Histograms of FP FN GT TP
        float * GTData=static_cast<float *>(RefLes->data);
        float * HistogramFP=GetMaskHistogram(MahalData, FPMap, -20, 20, 0.01, MahalWM->nt, numel);
        float * HistogramFN=GetMaskHistogram(MahalData, FNMap, -20, 20, 0.01, MahalWM->nt, numel);
        float * HistogramTP=GetMaskHistogram(MahalData, TPMap, -20, 20, 0.01, MahalWM->nt, numel);
        float * HistogramGT=GetMaskHistogram(MahalData, GTLes, -20, 20, 0.01, MahalWM->nt, numel);
        float * HistogramFP_gt=GetMaskHistogram(GTData, FPMap, 0, 1, 0.01, 1, numel);
        float * HistogramFN_gt=GetMaskHistogram(GTData, FNMap, 0, 1, 0.01, 1, numel);
        float * HistogramTP_gt=GetMaskHistogram(GTData, TPMap, 0, 1, 0.01, 1, numel);
        
        float * Stats=GetStatisticData(GTData, FPMap, numel);
        TxtStatGT<<"FP ";
        for (int i=0; i<9; i++) {
            TxtStatGT<< Stats[i]<<" ";
        }
        TxtStatGT << endl;
        delete [] Stats;
        Stats=NULL;
        Stats=GetStatisticData(GTData, FNMap, numel);
        TxtStatGT<<"FN ";
        for (int i=0; i<9; i++) {
            TxtStatGT<< Stats[i]<<" ";
        }
        TxtStatGT << endl;
        delete [] Stats;
        Stats=NULL;
        Stats=GetStatisticData(GTData, TPMap, numel);
        TxtStatGT<<"TP ";
        for (int i=0; i<9; i++) {
            TxtStatGT<< Stats[i]<<" ";
        }
        TxtStatGT << endl;
        delete [] Stats;
        Stats=NULL;
        
        
        
        TxtHistogramGT << "Init Values ";
        for (int i=0; i<100; i++) {
            TxtHistogramGT << i*0.01-0.5*0.01 << " ";
        }
        TxtHistogramGT << "FP GT" ;
        for (int i=0; i<100; i++) {
            TxtHistogramGT << HistogramFP_gt[i] << " ";
        }
        TxtHistogramGT << endl;
        TxtHistogramGT << "FN GT ";
        for (int i=0; i<100; i++) {
            TxtHistogramGT << HistogramFN_gt[i] << " ";
        }
        TxtHistogramGT << endl;
        TxtHistogramGT << "TP GT ";
        for (int i=0; i<100; i++) {
            TxtHistogramGT << HistogramTP_gt[i] << " ";
        }
        TxtHistogramGT << endl;
        
        
        
        
        cout << "Histogram obtained" << endl;
        int numbSteps=(int)40/0.01;
        TxtHistogram << "Init Values ";
        for (int i=0; i<numbSteps; i++) {
            TxtHistogram << i*0.01-20-0.5*0.01 << " ";
        }
        TxtHistogram << endl;
        for (int m=0; m<numbmodal; m++) {
            float * Stats=GetStatisticData(&MahalData[m*numel], FPMap, numel);
            TxtStat<<"FP ";
            for (int i=0; i<9; i++) {
                TxtStat<< Stats[i]<<" ";
            }
            TxtStat << endl;
            delete [] Stats;
            Stats=NULL;

            TxtHistogram << "FP Mod"<<m<<" ";
            for (int i=0; i<numbSteps; i++) {
                TxtHistogram << HistogramFP[i+m*numbSteps] << " ";
            }
            TxtHistogram << endl;
        }
        TxtHistogram << endl;
        for (int m=0; m<numbmodal; m++) {
            float * Stats=GetStatisticData(&MahalData[m*numel], FNMap, numel);
            TxtStat<<"FN ";
            for (int i=0; i<9; i++) {
                TxtStat<< Stats[i]<<" ";
            }
            TxtStat<<endl;
            TxtHistogram << "FN Mod"<<m<<" ";
            for (int i=0; i<numbSteps; i++) {
                TxtHistogram << HistogramFN[i+m*numbSteps] << " ";
            }
            TxtHistogram << endl;
        }
        TxtHistogram << endl;
        for (int m=0; m<numbmodal; m++) {
            float * Stats=GetStatisticData(&MahalData[m*numel], TPMap, numel);
            TxtStat << "TP ";
            for (int i=0; i<9; i++) {
                TxtStat<< Stats[i]<<" ";
            }
            TxtStat<<endl;
            TxtHistogram << "TP Mod"<<m<<" ";
            for (int i=0; i<numbSteps; i++) {
                TxtHistogram << HistogramTP[i+m*numbSteps] << " ";
            }
            TxtHistogram << endl;
        }
        TxtHistogram << endl;
        for (int m=0; m<numbmodal; m++) {
            float * Stats=GetStatisticData(&MahalData[m*numel], GTLes, numel);
            TxtStat<< "GT ";
            for (int i=0; i<9; i++) {
                TxtStat<< Stats[i]<<" ";
            }
            TxtStat<<endl;
            TxtHistogram << "GT Mod"<<m<<" ";
            for (int i=0; i<numbSteps; i++) {
                TxtHistogram << HistogramGT[i+m*numbSteps] << " ";
            }
            TxtHistogram << endl;
        }
        delete [] FPMap;
        delete [] FNMap;
        delete [] TPMap;
        delete [] GTLes;
        delete [] BoolLes;
        delete [] LesOppose;
        delete [] GTOppose;
        delete [] HistogramFN;
        delete [] HistogramFP;
        delete [] HistogramGT;
        delete [] HistogramTP;
        delete TreeToAnalyse;
        return EXIT_SUCCESS;

    }
    
    
    if (segment_analysis->flag_inLesCorr && segment_analysis->flag_inMahal && segment_analysis->flag_LS ){
        nifti_image * LesionNii = ReadFromFilename(segment_analysis->filename_inLesCorr);
        nifti_image * MahalNii = ReadFromFilename(segment_analysis->filename_inMahal);
        numbmodal=MahalNii->nu*MahalNii->nt;
        numel=LesionNii->nvox;
        bool * LesionBool=TranscribeArray<float, bool>(static_cast<float*>(LesionNii->data), numel);
        float * MahalData=static_cast<float*>(MahalNii->data);
        float MaxMahal=20;
        float MinMahal=-20;
        float step=0.1;
        FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameHistogram=FilenamePA_b+"Histogram_"+FilenamePA_e+".txt";
        ofstream TxtHistogram(FilenameHistogram.c_str());
        TxtHistogram << MinMahal << " " << MaxMahal << " "<<step<< endl;
        float * HistogramMahal=GetMaskHistogram(MahalData,LesionBool,MinMahal,MaxMahal,step,numbmodal,numel);
        int numbSteps=(int)(MaxMahal-MinMahal)/step;
        TxtHistogram << "Values ";
        for (int i=0; i<numbSteps; i++) {
            TxtHistogram << i*step+MinMahal-0.5*step << " ";
        }
        TxtHistogram << endl;
        for (int m=0; m<numbmodal; m++) {
            TxtHistogram << "Mod"<<m<<" ";
            for (int i=0; i<numbSteps; i++) {
                TxtHistogram << HistogramMahal[i+m*numbSteps] << " ";
            }
            TxtHistogram << endl;
        }
        for (int m=0; m<numbmodal; m++) {
            float * Stats=GetStatisticData(&MahalData[m*numel], LesionBool, numel);
            TxtHistogram << "Mod"<<m<<" ";
            for (int i=0; i<9; i++) {
                TxtHistogram<< Stats[i]<<" ";
            }
            TxtHistogram << endl;
            delete [] Stats;
            Stats=NULL;
        }
        nifti_image_free(MahalNii);
        nifti_image_free(LesionNii);
        delete [] LesionBool;
        delete [] HistogramMahal;
        return EXIT_SUCCESS;
        
    }
    
    
    if (TreeToAnalyse!=NULL && (segment_analysis->flag_inLesCorr)  && segment_analysis->flag_LS && !segment_analysis->flag_inQuadrant) {
        cout<<"Doing Local summary Regional (LFROI)"<<endl;
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameSave=FilenamePA_b+"RegionLFStatsFin_"+FilenamePA_e+".txt";
        string FilenameSaveMasked=FilenamePA_b+"MaskedFin_"+FilenamePA_e+".txt";
        ofstream TxtFile(FilenameSave.c_str());
        
        
        nifti_image * MahalWM=NULL;
        int numbMahal=0;
        nifti_image * LesionCorr=ReadFromFilename(segment_analysis->filename_inLesCorr);
        int Dim[3];
        Dim[0]=LesionCorr->dim[1];
        Dim[1]=LesionCorr->dim[2];
        Dim[2]=LesionCorr->dim[3];
        int Shift[3];
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Dim[0]*Dim[1];
        numel=TreeToAnalyse->GetNumberElements();
        float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        float * GMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());

        nifti_image * WMInliersNii=CreateNiiFromArray(WMTempInliers, LesionCorr, TreeToAnalyse->GetNumberElements());
        float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
        nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
//        MahalWM= MahalDistMaps(WMInliersNii, LesionCorr, TreeToAnalyse->GetDataImage() );
        MahalWM= MahalDistMaps(WMInliersNii, LesionCorr, DataNii );
        nifti_image_free(DataNii);
        

        numbMahal=MahalWM->dim[4];
        // Mask the CSF and the outliers from the studied regions
        vector<float *> AdditionInliers;
        AdditionInliers.push_back(WMTempInliers);
        AdditionInliers.push_back(GMTempInliers);
        float * AddedInliers=AddArray(AdditionInliers, TreeToAnalyse->GetNumberElements());
        bool * ThresholdAddedInliers=ThresholdArray<float, bool>(AddedInliers, 0.5, numel);
        nifti_image * AddedInliersNii=CreateNiiFromArray(ThresholdAddedInliers, LesionCorr, TreeToAnalyse->GetNumberElements());
        // Look at each of the regions separately and get CoM and local Mahal statistics
        int * RegionLabels=ComponentLabeling(LesionCorr, 26);
        //        float * LesionData=static_cast<float*>(LesionCorr->data);
        int * LesionMasked=MultiplyElementwiseChoose<bool, int, int>(ThresholdAddedInliers, RegionLabels, numel);
        nifti_image * LesionMaskedNii=CreateNiiFromArray(LesionMasked, LesionCorr, numel);
        nifti_set_filenames(LesionMaskedNii, FilenameSaveMasked.c_str(), 0, 0);
        nifti_image_write(LesionMaskedNii);
        
        int numbReg=GetMaxLabel(RegionLabels, numel);
        numbmodal=TreeToAnalyse->GetNumberModalities();
        cout << numbReg << "number of regions "<< endl;
        float * MahalWMData = static_cast<float*>(MahalWM->data);
        for (int r=0; r<numbReg; r++) {
            TxtFile << "Region"<<r+1<<" ";
            bool * RegMaskBool=CreateLesionBool(RegionLabels, r+1, numel);
            int CoM=GetCenterGravity(RegMaskBool, Dim);
            int * CoMCorr=CorrespondingCoordinates(CoM, Dim, Shift);
            for (int d=0; d<3; d++) {
                cout << CoMCorr[d]<<" ";
                TxtFile<<CoMCorr[d]<<" ";
            }
            TxtFile << GetLabelSize(r+1, RegionLabels, numel)<<" ";
            for (int m=0; m<numbMahal; m++) {
                float * Stats=GetStatisticData(&MahalWMData[m*numel], RegMaskBool, numel);
                TxtFile<<"Mod"<<m<<" ";
                for (int s=0; s<9; s++) {
                    cout << Stats[s] <<" ";
                    TxtFile << Stats[s]<<" ";
                }
                delete [] Stats;
                Stats=NULL;
            }
            cout << "Treating Data corr now "<< DataCorr<<endl;
            for (int m=0;m<numbmodal;m++){
                float * Stats=GetStatisticData(&DataCorr[m*numel], RegMaskBool, numel);
                cout << "Gotten statistics"<<endl;
                TxtFile<<ModalitiesCode[segment_analysis->vecModalities[m]-1]<<" ";
                for (int s=0; s<9; s++) {
                    cout << Stats[s] <<" ";
                    TxtFile << Stats[s]<<" ";
                }
                delete [] Stats;
                Stats=NULL;
            }
            TxtFile<< endl;
            cout << endl;
            delete [] CoMCorr;
            CoMCorr=NULL;
        }
        delete [] DataCorr;
        nifti_image_free(LesionCorr);
        nifti_image_free(LesionMaskedNii);
        nifti_image_free(MahalWM);
        nifti_image_free(AddedInliersNii);
        nifti_image_free(WMInliersNii);
        delete [] AddedInliers;
        delete [] LesionMasked;
        delete [] WMTempInliers;
        delete [] GMTempInliers;
        delete [] ThresholdAddedInliers;
        delete TreeToAnalyse;
        return EXIT_SUCCESS;

    }
    
    
    if (segment_analysis->flag_LaplaceSolImage && (segment_analysis->flag_inLesCorr) && segment_analysis->flag_inQuadrant && segment_analysis->flag_LS) {
        cout << "Doing Local summary"<<endl;
        nifti_image * Lesion=ReadFromFilename(segment_analysis->filename_inLesCorr);
        nifti_image * Quadrant=NULL;
        nifti_image * Layers=NULL;
        nifti_image * MahalWM=NULL;
        int numbMahal=0;
        if (TreeToAnalyse!=NULL  && segment_analysis->flag_inLesCorr && segment_analysis->flag_inToAnalyse) {
            cout <<"Doing WMI"<<endl;
            Quadrant=ReadFromFilename(segment_analysis->filename_inQuadrant);
            Layers=ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
            nifti_image * LesionCorr=ReadFromFilename(segment_analysis->filename_inLesCorr);
            float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            float * WMTempOutliers=CreateLong(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            
            float * LesionCorrData=static_cast<float*>(LesionCorr->data);
            vector<float *> AdditionWM;
            AdditionWM.push_back(WMTempInliers);
            AdditionWM.push_back(WMTempOutliers);
            AdditionWM.push_back(LesionCorrData);
            float * AddedWM=AddArray(AdditionWM, TreeToAnalyse->GetNumberElements());
            cout << "Added array obtained"<<endl;
            nifti_image * AddedWMNii=CreateNiiFromArray(AddedWM, LesionCorr, TreeToAnalyse->GetNumberElements());
            nifti_image * WMInliersNii=CreateNiiFromArray(WMTempInliers, LesionCorr, TreeToAnalyse->GetNumberElements());
            float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
            nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
//            MahalWM= MahalDistMaps(WMInliersNii, AddedWMNii, TreeToAnalyse->GetDataImage() );
            MahalWM= MahalDistMaps(WMInliersNii, AddedWMNii, DataNii );
            nifti_image_free(DataNii);
            delete [] DataCorr;
            numbMahal=MahalWM->dim[4];
            vector<float *> LSWMIStats=LocalSummaryStats(WMInliersNii, Quadrant, Layers, MahalWM);
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            
            string FilenameSave=FilenamePA_b+"LocalSummaryWMIStats_"+FilenamePA_e+".txt";
            ofstream TxtFileWMI(FilenameSave.c_str());
            PrintLocalSummaryStats(LSWMIStats, TxtFileWMI);
            int sizeLSWMI=LSWMIStats.size();
            cout << "Printed Local WMI"<<endl;
            for (int m=0; m<sizeLSWMI; m++) {
                if (LSWMIStats[m]!=NULL) {
                    delete [] LSWMIStats[m];
                    LSWMIStats[m]=NULL;
                }
            }
            cout << "Cleared LSWMI"<<endl;
            nifti_image_free(AddedWMNii);
            nifti_image_free(WMInliersNii);
            nifti_image_free(Layers);
            nifti_image_free(Quadrant);
            Layers=NULL;
            Quadrant=NULL;
            
        }
        cout <<"Going to second part of LS"<<endl;
        Layers=ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
        Quadrant=ReadFromFilename(segment_analysis->filename_inQuadrant);
        cout << "Loaded images "<<endl;
        float * LayersData=static_cast<float*>(Layers->data);
        float * QuadrantData=static_cast<float*>(Quadrant->data);
        numel =Layers->nvox;
        cout << "Read Images "<<endl;
        float MaxLayers=GetMax(LayersData, numel)-1;
        float MaxLobes=GetMax(QuadrantData, numel);
        cout << "Max layers is "<< MaxLayers << "Max lobes is "<< MaxLobes<<endl;
        int numbTot=MaxLayers*MaxLobes;
        
        float * LSResults=LocalSummaryLesion_bis(Lesion, Quadrant, Layers);
        cout <<"LSResults done"<<endl;
        int * LSResultsNumber=LocalSummaryLesionNumber_bis(Lesion, Quadrant, Layers);
        cout <<"LSResults Number done"<<endl;
        float * LSResultsMahal=LocalSummaryMahal_bis(Lesion, Quadrant, Layers, MahalWM);
        cout << "LSResults Mahal done"<<endl;
        FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        

        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        
        
        string FilenameSave=FilenamePA_b+"LocalSummaryTotBis2_"+FilenamePA_e+".txt";
        
        ofstream TxtFileLesion(FilenameSave.c_str());
        PrintLocalSummaryTot_bis(LSResults, LSResultsNumber, LSResultsMahal, numbMahal,numbTot, TxtFileLesion);
        nifti_image_free(Lesion);
        nifti_image_free(Quadrant);
        nifti_image_free(Layers);
        delete [] LSResults;
        delete [] LSResultsNumber;
        if (LSResultsMahal!=NULL) {
            delete [] LSResultsMahal;
        }
        LSResults=NULL;
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
    }

    
    
    if (segment_analysis->flag_AffGen) {
        float Angle[3];
        float Translation[3];
        float maxA=fabs(segment_analysis->AffParam[0]);
        float minA=-maxA;
        float maxT=fabs(segment_analysis->AffParam[1]);
        float minT=-maxT;
        srand(time(NULL));
        for (int i=0; i<3; i++) {
            Angle[i]=3.14*randu(minA, maxA)/180;
            Translation[i]=randu(minT, maxT);
        }
        ofstream TxtFileAffine(segment_analysis->filename_Out);
        TxtFileAffine << cos(Angle[0])*cos(Angle[1])<< " ";
        TxtFileAffine << cos(Angle[0])*sin(Angle[1])*sin(Angle[2])-sin(Angle[0])*cos(Angle[2])<< " ";
        TxtFileAffine << cos(Angle[0])*sin(Angle[1])*cos(Angle[2])+sin(Angle[0])*sin(Angle[2])<< " ";
        TxtFileAffine << Translation[0]<< " " << endl;
        TxtFileAffine << cos(Angle[1])*sin(Angle[0])<< " ";
        TxtFileAffine << sin(Angle[0])*sin(Angle[1])*sin(Angle[2])+cos(Angle[0])*cos(Angle[2])<< " ";
        TxtFileAffine << sin(Angle[0])*sin(Angle[1])*cos(Angle[2])-cos(Angle[0])*sin(Angle[2])<< " ";
        TxtFileAffine << Translation[1] << " " << endl;
        TxtFileAffine << -sin(Angle[1])<< " ";
        TxtFileAffine << cos(Angle[1])*sin(Angle[2])<< " ";
        TxtFileAffine << cos(Angle[1])*cos(Angle[2])<< " ";
        TxtFileAffine << Translation[2] << " "<< endl;
        TxtFileAffine << 0<< " ";
        TxtFileAffine << 0<< " ";
        TxtFileAffine << 0<< " ";
        TxtFileAffine << 1<< " ";
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
    }
    
    if (segment_analysis->flag_BFGen) {
        cout<<"Generating bias field"<<endl;

        float * BFCoeffs=NULL;
        int numbBF=0;
        int numbCoeffs=0;
        if (segment_analysis->flag_BFText) {
            ifstream textAdapt (segment_analysis->filename_BFText);
            if(!textAdapt){
                std::cout<<"could not open the text file properly for the params adaptation!"<<endl;
            }
            else{
                int m=0;
                std::string lineAdapt;
                while(getline(textAdapt,lineAdapt)){
                    istringstream in(lineAdapt);
                    std:: string type;
                    in >> type;
                    if(type == "BFOrder"){

                        float Order=0;
                        in >> Order;
                        segment_analysis->BFParam.push_back(Order);
                    }
                    if(type == "NumbModal"){
                        
                        float NumbModa=0;
                        in >> NumbModa;
                        segment_analysis->BFParam.push_back(NumbModa);
                        int numbBF=(segment_analysis->BFParam[0]+1)*(segment_analysis->BFParam[0]+2)*(segment_analysis->BFParam[0]+3)/6;
                        //        Generate random number for coefficients
                        numbCoeffs=numbBF*segment_analysis->BFParam[1];
                        cout << numbCoeffs<<endl;
                        BFCoeffs=new float[numbCoeffs];
                    }
                    else{
                        int numbBF=(segment_analysis->BFParam[0]+1)*(segment_analysis->BFParam[0]+2)*(segment_analysis->BFParam[0]+3)/6;
                        //        Generate random number for coefficients
                        //                        int numbCoeffs=numbBF*segment_analysis->BFParam[1];
                        float Coeff=0;
                        BFCoeffs=new float[numbBF];
                        while (in >> Coeff) {
                            BFCoeffs[m]=Coeff;
                            cout << Coeff << endl;
                            m++;
                        }
                    }
                }
            }
        }
        else{
            int numbBF=(segment_analysis->BFParam[0]+1)*(segment_analysis->BFParam[0]+2)*(segment_analysis->BFParam[0]+3)/6;
            //        Generate random number for coefficients
            int numbCoeffs=numbBF*segment_analysis->BFParam[1];
            //        float * BFCoeffs=new float[numbCoeffs];
            float max=fabs(segment_analysis->BFParam[2]);
            float min=-max;
            cout<<"the BF generated coefficients are ";
            BFCoeffs=randuSeries(min, max, numbCoeffs);
        }
        
        for (int i=0; i<numbCoeffs; i++) {
            //            BFCoeffs[i]=randu(min, max);
            cout << BFCoeffs[i] << " ";
        }
        
        cout << endl;
        cout<<"Obtaining the BF..."<<endl;
        nifti_image * DataRef=NULL;
        int numel = 0;
        if (segment_analysis->flag_inDataComp){
            DataRef=ReadFromFilename(segment_analysis->filename_inDataComp);
            numel = DataRef->nx*DataRef->ny*DataRef->nz;
        }
        else{
            DataRef=TreeToAnalyse->GetDataImage();
            numel = TreeToAnalyse->GetNumberElements();
        }
        float * BFCorrection=TreeToAnalyse->GenerateBFCorrection(segment_analysis->BFParam[0], segment_analysis->BFParam[1], BFCoeffs, DataRef);
        cout << segment_analysis->BFParam[1]<<endl;
        cout << "BFCorrection performed"<< BFCoeffs[0]<<endl;
        cout <<" Performing the exponentiation"<<endl;
        float * ExpBF=ExpArray<float>(BFCorrection,numel*segment_analysis->BFParam[1]);
        cout<<"Exponentiation performed"<<endl;
        //        float * LongBFGen=CreateLongPadding<float>(ExpBF, 1, TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        cout<<"Padding performed"<<endl;
        nifti_image * BFGenNii=CreateNiiFromArray(ExpBF,DataRef,numel*segment_analysis->BFParam[1]);
        cout<<"BFNii created"<<endl;
        if (! segment_analysis->flag_fileOut){
            string FilenameTemp;
            if(segment_analysis->flag_inDataComp){
                FilenameTemp= nifti_makebasename(segment_analysis->filename_inDataComp);

            }
            else{
                FilenameTemp = nifti_makebasename(segment_analysis->filename_SegTot);
            }
                int Index=FilenameTemp.find_last_of('/');
                string FilenameTemp_b=FilenameTemp.substr(0,Index+1);
                if(segment_analysis->flag_outputDir){
                    FilenameTemp_b=segment_analysis->name_outputDir;
                }
                string FilenameTemp_e=FilenameTemp.substr(Index+1,FilenameTemp.length());
                string Fin = FilenameTemp_b + "BFGen_" + FilenameTemp_e + ".nii.gz";
                segment_analysis->filename_Out = strdup(Fin.c_str());
        }
        FilenamePA=nifti_makebasename(segment_analysis->filename_Out);
        if (segment_analysis->flag_inDataComp){
            nifti_image_free(DataRef);
            DataRef=NULL;
        }
        
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameSaveCoeffs=FilenamePA_b+"BFCoeffs_"+FilenamePA_e+".txt";
        ofstream TxtFileCoeffs(FilenameSaveCoeffs.c_str());
        TxtFileCoeffs << "BFOrder " << segment_analysis->BFParam[0];
        TxtFileCoeffs << "NumbModal "<< segment_analysis->BFParam[1];
        for (int m=0; m<segment_analysis->BFParam[1]; m++) {
            for (int c=0; c<numbBF; c++) {
                TxtFileCoeffs << BFCoeffs[c+m*numbBF]<<" ";
            }
            TxtFileCoeffs << endl;
        }
        delete [] BFCorrection;
        //        delete [] LongBFGen;
        delete [] BFCoeffs;
        delete [] ExpBF;
        cout << "Image saved in "<< segment_analysis->filename_Out<< endl;
        nifti_set_filenames(BFGenNii, segment_analysis->filename_Out, 0, 0);
        nifti_image_write(BFGenNii);
        nifti_image_free(BFGenNii);
        cout<<"BF saved"<<endl;
        BFCoeffs=NULL;
        BFCorrection=NULL;
        //        LongBFGen=NULL;
        ExpBF=NULL;
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
    }
    
    //    Creation of gaussian noise intensity lesion like lesions according to the values for the intensity in the inlier (normal images)
    if(segment_analysis->flag_Gaussian){
        cout<<"Doing the Gaussian noise"<<endl;
        //        First get the parameters of the distribution given resp inliers of the class to consider for the different images separately. Means that we assume independence between modalities...
        //        Reading images on which to get parameters.
        vector<nifti_image*> BasisImagesVector;
        int numbImages=segment_analysis->filename_ImagesGaussian.size();
        int numelmasked=TreeToAnalyse->GetNumberMaskedElements();
        int * S2L=TreeToAnalyse->GetS2L();
        for(int m=0;m<numbImages;m++){
            BasisImagesVector.push_back(ReadFromFilename(segment_analysis->filename_ImagesGaussian[m]));
        }
        cout <<"Size of BasisImagesVector is "<< BasisImagesVector.size()<<endl;
        //        Get inlier NormResp for wanted index to build Gaussian noise
        float * NormRespInlier=TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->index_Gaussian)->GetNormResp();
        vector<float> MeanVector;
        vector<float> StdVector;
        
        //        Get parameters for the distributions afterwards
        for (int m=0; m<numbImages; m++) {
            float * DataImage=static_cast<float*>(BasisImagesVector[m]->data);
            float * DataImageShort=CreateShort(DataImage, S2L, numelmasked);
            MeanVector.push_back(GetMeanData(DataImageShort, NormRespInlier, numelmasked));
            cout <<"Mean init is"<<MeanVector[m]<<endl;
            StdVector.push_back(sqrtf(GetVarianceData(DataImageShort, NormRespInlier, numelmasked)));
            MeanVector[m]+=StdVector[m]*segment_analysis->Gaussian_Zscore[m];
            cout<< "and mean fin is "<<MeanVector[m]<<endl;
            delete [] DataImageShort;
            DataImageShort=NULL;
        }
        //        Get the intensity distribution over the mask
        cout << MeanVector.size() << " "<< StdVector.size()<< " " << BasisImagesVector[0]->nvox<<" "<<TreeToAnalyse->GetMask()->nvox<<endl;
        nifti_image * GaussianNoiseVector=RandomGaussianNoiseNii(MeanVector, StdVector, BasisImagesVector[0], TreeToAnalyse->GetMask());
        nifti_set_filenames(GaussianNoiseVector, segment_analysis->filename_Out, 0, 0);
        nifti_image_write(GaussianNoiseVector);
        nifti_image_free(GaussianNoiseVector);
        delete TreeToAnalyse;
        TreeToAnalyse=NULL;
        for (int m=0; m<numbImages; m++) {
            nifti_image_free(BasisImagesVector[m]);
            BasisImagesVector[m]=NULL;
        }
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_inToAnalyse && segment_analysis->flag_Grad){ // Create Image Gradient over 3 directions
        nifti_image * SegToAnalyse=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveNii=FilenamePA_b+"Gradient_"+FilenamePA_e+".nii.gz";
        int Dim[3];
        int Shift[3];
        float PixDim[3];
        int numel=SegToAnalyse->nvox;
        for(int d=0;d<3;d++){
            Dim[d]=SegToAnalyse->dim[d+1];
            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
            PixDim[d]=SegToAnalyse->pixdim[d+1];
        }
        float * Data=static_cast<float*>(SegToAnalyse->data);
        int numelmasked=0;
        bool * Mask=ThresholdArray<float,bool>(Data,0,numel);
        cout << CountNonZero(Mask,numel)<< endl;
        int * L2S=MakeL2S(Mask,Dim);
        int * S2L=MakeS2L(Mask,Dim,numelmasked);
        float * DataToGrad=CreateShort(Data,S2L,numelmasked);
        float * GradX=GradientImage(DataToGrad,0,0,Dim,Shift,PixDim,L2S,S2L,numelmasked);
        float * GradY=GradientImage(DataToGrad,1,0,Dim,Shift,PixDim,L2S,S2L,numelmasked);
        float * GradZ=GradientImage(DataToGrad,2,0,Dim,Shift,PixDim,L2S,S2L,numelmasked);
        cout << "Obtained grad " << numelmasked <<endl;
        float * ResultGradient=new float[3*numel];
        for(int i=0;i<numel;i++){
            ResultGradient[i]=GradX[L2S[i]];
            ResultGradient[i+numel]=GradY[L2S[i]];
            ResultGradient[i+2*numel]=GradZ[L2S[i]];
        }
        cout << "Obtained ResultFin" <<endl;
        nifti_image * GradientNii=CreateNiiFromArray(ResultGradient,SegToAnalyse,3*numel);
        nifti_set_filenames(GradientNii,FilenameSaveNii.c_str(),0,0);
        nifti_image_write(GradientNii);
        nifti_image_free(GradientNii);
        nifti_image_free(SegToAnalyse);
        delete [] S2L;
        delete [] L2S;
        delete [] GradX;
        delete [] GradY;
        delete [] GradZ;
        delete [] Mask;
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_Dot){ // Dot product between two vectorial images;
        nifti_image * ToDot1=ReadFromFilename(segment_analysis->filename_dot1);
//        nifti_image * ToDot2=ReadFromFilename(segment_analysis->filename_dot2);
        nifti_image * ToDot2=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        if(!CheckCompatibleDimensions(ToDot1,ToDot2)){
            return EXIT_FAILURE;
        }
        FilenamePA=nifti_makebasename(segment_analysis->filename_dot1);
        string FilenamePA2=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        int Index2=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenamePA_e2=FilenamePA2.substr(Index2+1,FilenamePA2.length());
        string FilenameSaveNii=FilenamePA_b+"Dot_"+FilenamePA_e+"_"+FilenamePA_e2+".nii.gz";
        float * Data1=static_cast<float*>(ToDot1->data);
        float * Data2=static_cast<float*>(ToDot2->data);
        int numel=ToDot1->nx*ToDot1->ny*ToDot1->nz;
        int sizeVec=ToDot1->nt*ToDot1->nu;
        float * Product=MultiplyElementwise(Data1,Data2,numel*sizeVec);
        float * ResultFin=new float[numel];
        for(int i=0;i<numel;i++){
            ResultFin[i]=0;
            for(int d=0;d<sizeVec;d++){
                ResultFin[i]+=Product[i+numel*d];
            }
        }
        nifti_image * ResultNii=CreateNiiFromArray(ResultFin,ToDot1,numel);
        nifti_set_filenames(ResultNii,FilenameSaveNii.c_str(),0,0);
        nifti_image_write(ResultNii);
        nifti_image_free(ResultNii);
        nifti_image_free(ToDot1);
        nifti_image_free(ToDot2);
        delete [] ResultFin;
        delete [] Product;
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_getEntropy){
        cout << "Starting entropy"<<endl;
        nifti_image * SegToAnalyse=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        nifti_image * EntropyResult=CreateEntropy(SegToAnalyse,segment_analysis);
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveNii=FilenamePA_b+"Entropy_"+FilenamePA_e+".nii.gz";
        nifti_set_filenames(EntropyResult,FilenameSaveNii.c_str(),0,0);
        nifti_image_write(EntropyResult);
        nifti_image_free(EntropyResult);
        nifti_image_free(SegToAnalyse);
        cout << "Entropy done !"<<endl;
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_inToAnalyse && segment_analysis->flag_getCH){
        nifti_image * SegToAnalyse=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveLine=FilenamePA_b+"CharMOBB_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLine2=FilenamePA_b+"MOBB_"+FilenamePA_e+".nii.gz";
        int Dim[3];
        int Shift[3];
        float PixDim[3];
        int numel=SegToAnalyse->nvox;
        for(int d=0;d<3;d++){
            Dim[d]=SegToAnalyse->dim[d+1];
            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
            PixDim[d]=SegToAnalyse->pixdim[d+1];
        }

        int * Components=ComponentLabeling(SegToAnalyse,18);
        int * OrderedLabels=OrderedVolumeLabel(Components,1,numel);
        int numbLabels=GetMaxLabel(OrderedLabels,numel);
        float * FinalMinimal=new float[numel];
        float * FinalChar=new float[4*numel];
        for(int i=0;i<numel;i++){
            FinalMinimal[i]=0;
            for(int d=0;d<4;d++){
                FinalChar[i+numel*d]=0;
            }
        }
        for (int l=0;l<numbLabels;l++){
            bool * LesionBool = CreateLesionBool(OrderedLabels,l+1,numel);
            bool * LesionBorder = CreateBorderFromBool(LesionBool,Dim,Shift);
            int * ListCoordLabel = ListCoordinatesBoolSeg(LesionBorder,Dim,Shift);
            int numbLabelBorder = CountNonZero(LesionBorder,numel);
            vector<Face> CH=ConvexHullFrom3D(ListCoordLabel,numbLabelBorder,Dim,PixDim);
            cout << "CH contains "<<CH.size()<<endl;
            vector <Point> ListBB= MMBBFromConvexHull( CH);
            vector<Point> ListP=ListPointsFromFaces(CH);
            float * PointsCH = CreateSegPoints(ListP,Dim,Shift,PixDim);
            MultiplyFloatArrayBy(PointsCH,100,numel);
            float * CharBB=CreateCharMM(LesionBool,ListBB,Dim,Shift,PixDim);
            float * MBBSeg= CreatePVBoxFromBBPoint( ListBB,  Dim,  Shift,  PixDim);
            AddElementwiseInPlace(FinalMinimal,PointsCH,numel);
            AddElementwiseInPlace(FinalMinimal,MBBSeg,numel);
            AddElementwiseInPlace(FinalChar,CharBB,4*numel);
            delete [] LesionBool;
            delete [] LesionBorder;
            delete [] ListCoordLabel;
            delete [] PointsCH;
            delete [] MBBSeg;
            delete [] CharBB;
            nifti_image * MOBBNii=CreateNiiFromArray(FinalMinimal,SegToAnalyse,numel);
            nifti_image * CharBBNii=CreateNiiFromArray(FinalChar,SegToAnalyse,4*numel);
            nifti_set_filenames(MOBBNii,FilenameSaveLine2.c_str(),0,0);
            nifti_set_filenames(CharBBNii,FilenameSaveLine.c_str(),0,0);
            nifti_image_write(MOBBNii);
            nifti_image_write(CharBBNii);
            nifti_image_free(MOBBNii);
            nifti_image_free(CharBBNii);

        }
        nifti_image_free(SegToAnalyse);
        return EXIT_SUCCESS;

    }


    if(segment_analysis->flag_inToAnalyse && segment_analysis->flag_getBB){
        nifti_image * SegToAnalyse=ReadFromFilename(segment_analysis->filename_inToAnalyse);
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveNii=FilenamePA_b+"BBCoM_"+FilenamePA_e+".nii.gz";
        string FilenameSaveEigen=FilenamePA_b+"Eigen_"+FilenamePA_e+".nii.gz";
        string FilenameSaveEigen2=FilenamePA_b+"Eigen2_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLine=FilenamePA_b+"Line_"+FilenamePA_e+".nii.gz";
        string FilenameSaveLine2=FilenamePA_b+"Line2_"+FilenamePA_e+".nii.gz";
        int Dim[3];
        int Shift[3];
        float PixDim[3];
        int numel=SegToAnalyse->nvox;
        for(int d=0;d<3;d++){
            Dim[d]=SegToAnalyse->dim[d+1];
            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
            PixDim[d]=SegToAnalyse->pixdim[d+1];
        }
        int * ResultsBBCoM=new int[numel];
        float * MainEigenVector=new float[3*numel];
        float * MainEigenVector2=new float[3*numel];
        float * LineDrawn=new float[numel];
        float * LineDrawn2=new float[numel];
        for(int i=0;i<numel;i++){
            ResultsBBCoM[i]=0;
            MainEigenVector[i]=0;
            MainEigenVector[i+numel]=0;
            MainEigenVector[i+2*numel]=0;
            MainEigenVector2[i]=0;
            MainEigenVector2[i+numel]=0;
            MainEigenVector2[i+2*numel]=0;
            LineDrawn[i]=0;
        }
        int * Components=ComponentLabeling(SegToAnalyse,18);
        int * OrderedLabels=OrderedVolumeLabel(Components,1,numel);
        int numLabels=GetMaxLabel(OrderedLabels,numel);

        for(int l=0;l<numLabels;l++){
            bool * LabelToTest=CreateLesionBool(OrderedLabels,l+1,numel);
            int * BBCoord=FindBoundingBox(LabelToTest,Dim,Shift);
            bool * BB=CreateBoolBoundingBox(BBCoord,Dim);
            AddElementwiseInPlace<int,bool>(ResultsBBCoM,BB,numel);
            int CoM=GetCenterGravity(LabelToTest,Dim);
            ResultsBBCoM[CoM]=100;
            int LabelSize=CountNonZero(LabelToTest, numel);
            int * ListCoordinates=ListCoordinatesBoolSeg(LabelToTest, Dim, Shift);
            float * ListCoordF=ListCoordinatesTransformed(ListCoordinates, LabelSize, PixDim);
            //        float * Test=TreeToAnalyse->TransposeMatrix(ListCoordF, LabelSize, 3);
            SVD ResultSVD=SVD(ListCoordF, LabelSize, 3);
            float * V=ResultSVD.getV();
            for(int i=0;i<numel;i++){
                if(LabelSize>2){
                    if(LabelToTest[i]){
                        MainEigenVector2[i]=V[6];
                        MainEigenVector2[i+numel]=V[7];
                        MainEigenVector2[i+2*numel]=V[8];
                    }
                }
                else{
                    if(LabelToTest[i]){
                        MainEigenVector2[i]=V[3];
                        MainEigenVector2[i+numel]=V[4];
                        MainEigenVector2[i+2*numel]=V[5];
                    }
                }
            }
            float * LineTemp;
            if(LabelSize>2){
                LineTemp=CreateLineFromCoMBB(LabelToTest,CoM,BBCoord,MainEigenVector,Dim,Shift);
                //             LineTemp2=CreateLineFromEigenCoMBB(&V[6],CoM,BBCoord,Dim,Shift);
            }
            else{
                LineTemp=CreateLineFromCoMBB(LabelToTest,CoM,BBCoord,MainEigenVector,Dim,Shift);
                //                LineTemp2=CreateLineFromEigenCoMBB(&V[3],CoM,BBCoord,Dim,Shift);
            }
            AddElementwiseInPlace<float,float>(LineDrawn,LineTemp,numel);
            //            AddElementwiseInPlace<float,float>(LineDrawn2,LineTemp2,numel);
            delete [] LineTemp;
            delete [] ListCoordF;
            delete [] ListCoordinates;
            delete [] LabelToTest;
            delete [] BBCoord;
            delete [] BB;
        }
        nifti_image * BBCoMNii=CreateNiiFromArray(ResultsBBCoM,SegToAnalyse,numel);
        nifti_image * EigenNii=CreateNiiFromArray(MainEigenVector,SegToAnalyse,3*numel);
        nifti_image * LineNii=CreateNiiFromArray(LineDrawn,SegToAnalyse,numel);
        nifti_set_filenames(BBCoMNii,FilenameSaveNii.c_str(),0,0);
        nifti_set_filenames(EigenNii,FilenameSaveEigen.c_str(),0,0);
        nifti_set_filenames(LineNii,FilenameSaveLine.c_str(),0,0);
        nifti_image_write(EigenNii);
        nifti_image_write(BBCoMNii);
        nifti_image_write(LineNii);
        nifti_image_free(LineNii);
        nifti_image_free(BBCoMNii);
        nifti_image_free(EigenNii);
        nifti_image_free(SegToAnalyse);
        delete [] Components;
        delete [] OrderedLabels;
        delete [] MainEigenVector;
        delete [] MainEigenVector2;
        delete [] LineDrawn;
        delete [] LineDrawn2;
        delete [] ResultsBBCoM;
        return EXIT_SUCCESS;
    }
    
    if(segment_analysis->flag_inROIFiles){ // Label fusion for PIVOT data

        
        FilenamePA=nifti_makebasename(segment_analysis->filename_inROIFiles[0]);
        
        
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveROI=FilenamePA_b+"NewROIAff_"+FilenamePA_e+".nii.gz";
        string FilenameSaveCoM=FilenamePA_b+"NewROIAffCoM_"+FilenamePA_e+".txt";
        
        if (segment_analysis->flag_getCoM) {
            vector<string> FilenameROIVector=ReadListFilenamesFromTextFile(segment_analysis->filename_inROIFiles[1]);
            vector<nifti_image*> VectorImagesROI=ReadFromFilenamesVector(FilenameROIVector);
            nifti_image * ROIRef=ReadFromFilename(segment_analysis->filename_inROIFiles[0]);
            string FilenameCoMRef=FilenamePA_b+"CoM_"+FilenamePA_e+".txt";
            ExtractWriteCoMComponent(ROIRef,26,FilenameCoMRef);
            int numbROI=FilenameROIVector.size();
            for (int n=0; n<numbROI; n++) {
                FilenamePA=nifti_makebasename(FilenameROIVector[n].c_str());
                
                
                int Index=FilenamePA.find_last_of('/');
                FilenamePA_b=FilenamePA.substr(0,Index+1);
                if(segment_analysis->flag_outputDir){
                    FilenamePA_b=segment_analysis->name_outputDir;
                }
                FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
                FilenameCoMRef=FilenamePA_b+"CoM_"+FilenamePA_e+".txt";
                ExtractWriteCoMComponent(VectorImagesROI[n],26,FilenameCoMRef);
            }
            for(int i=0;i<numbROI;i++){
                nifti_image_free(VectorImagesROI[i]);
                VectorImagesROI[i]=NULL;
            }
            // return EXIT_SUCCESS;
        }
        
        vector<string> FilenameROIVector=ReadListFilenamesFromTextFile(segment_analysis->filename_inROIFiles[1]);
        vector<string> FilenameExVector=ReadListFilenamesFromTextFile(segment_analysis->filename_inExFiles[1]);
        nifti_image * ImageRef=ReadFromFilename(segment_analysis->filename_inExFiles[0]);
        //        nifti_image * TestNii=ReadFromFilename(g_scratch_dir + "/PUA012025_WIPAXFLAIRFASTSENSE.nii.gz");
        nifti_image * ROIRef=ReadFromFilename(segment_analysis->filename_inROIFiles[0]);
        cout << "Size FileNames is " << FilenameExVector.size()<< endl;
        vector<nifti_image *> VectorImagesEx=ReadFromFilenamesVector(FilenameExVector);
        cout << "Size VectorImages is " << VectorImagesEx.size()<<endl;
        vector<nifti_image*> VectorImagesROI=ReadFromFilenamesVector(FilenameROIVector);
        
        
        nifti_image * Mask=ReadFromFilename(segment_analysis->filename_mask);
        nifti_image * NewROI=FinalROIPlacement(VectorImagesEx, VectorImagesROI, ImageRef, ROIRef, Mask,FilenameSaveCoM.c_str());
        int * CLTest=ComponentLabeling(NewROI, 26);
        cout << "NumbLabelFin is "<< GetMaxLabel(CLTest, numel);
        delete [] CLTest;
        CLTest=NULL;
        ExtractWriteCoMComponent(NewROI,26,FilenameSaveCoM);
        nifti_set_filenames(NewROI, FilenameSaveROI.c_str(), 0, 0);
        nifti_image_write(NewROI);
        nifti_image_free(NewROI);
        nifti_image_free(ImageRef);
        nifti_image_free(ROIRef);
        int numbEx=VectorImagesEx.size();
        int numbROI=VectorImagesROI.size();
        for(int i=0;i<numbEx;i++){
            nifti_image_free(VectorImagesEx[i]);
            VectorImagesEx[i]=NULL;
        }
        for(int i=0;i<numbROI;i++){
            nifti_image_free(VectorImagesROI[i]);
            VectorImagesROI[i]=NULL;
        }
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_inConnect && segment_analysis->flag_inOutText){
        FilenamePA=nifti_makebasename(segment_analysis->filename_inConnect);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameClassifLac=FilenamePA_b+"ClassifLacChanged_"+FilenamePA_e+".nii.gz";
        nifti_image * ImageToUse=ReadFromFilename(segment_analysis->filename_inConnect);
        numel=ImageToUse->nvox;
        vector<int> ClassifFinal=ReadClassifFromTextFile(segment_analysis);
        int sizeClassif=ClassifFinal.size();
        float * DataConnected=static_cast<float*>(ImageToUse->data);
        int * DataConnectedInt=TranscribeArray<float,int>(DataConnected,numel);
        int maxLabel=GetMaxLabel(DataConnectedInt,numel);
        cout<<"maxLabel and maxClassif are "<<maxLabel<<" "<<sizeClassif<<endl;
        int * DataReclassif=new int[numel];
        for(int i=0;i<numel;i++){
            DataReclassif[i]=0;
        }
        int Dim[3];
        Dim[0]=ImageToUse->dim[1];
        Dim[1]=ImageToUse->dim[2];
        Dim[2]=ImageToUse->dim[3];
        for(int l=0;l<maxLabel;l++){
            vector<int> VectorIndices=GetListIndicesLabel(DataConnectedInt,Dim,l+1);
            int sizeLabel=VectorIndices.size();
            for(int i=0;i<sizeLabel;i++){
                DataReclassif[VectorIndices[i]]=ClassifFinal[l];
            }
        }
        nifti_image * NewClassif=CreateNiiFromArray(DataReclassif,ImageToUse,numel);
        nifti_set_filenames(NewClassif,FilenameClassifLac.c_str(),0,0);
        nifti_image_write(NewClassif);
        delete [] DataConnectedInt;
        delete [] DataReclassif;
        nifti_image_free(ImageToUse);
        nifti_image_free(NewClassif);
        return EXIT_SUCCESS;
    }
    
    if (segment_analysis->flag_inOutText) {
        if (segment_analysis->flag_inConnect) {
            FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
            
            
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenameClassifLac=FilenamePA_b+"ClassifLacChanged_"+FilenamePA_e+".nii.gz";
            
            string FilenameClassifLacFin=FilenamePA_b+"ClassifLacChangedSub_"+FilenamePA_e+".nii.gz";
            string FilenameSave=FilenamePA_b+"OutlierLacTestChanged_"+FilenamePA_e+".txt";

            ofstream TxtFileLesion(FilenameSave.c_str());
            
            
            nifti_image * ImageToUse=ReadFromFilename(segment_analysis->filename_inConnect);
            numel=ImageToUse->nvox;
            vector<Outlier *> VectorOutliers=ReadOutliersFromTextFile(segment_analysis, TreeToAnalyse);

            int * OrderedOutliers= ExtractTPFromNii<int>(ImageToUse, 0);
            int * OutlierClassif=CopyArray(OrderedOutliers, numel);
            int * LesionSubClassif=CopyArray(OrderedOutliers, numel);
            int * GravImage=NULL;
            int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
            float PixDim[3];
            if (numbFinalLesions>0) {
                GravImage=new int[3];
                for (int d=0;d<3; d++) {
                    PixDim[d]=ImageToUse->pixdim[d+1];
                    GravImage[d]=(VectorOutliers[0]->CentreGravity[d]-VectorOutliers[0]->VectorDiffGrav[d]/PixDim[d]);
                }
            }
            int CountModified=0;
            for (int l=0; l<numbFinalLesions; l++) {
                OutlierType newOT=OutlierTypePreliminary(VectorOutliers[l], segment_analysis);
                
                
                if(newOT != VectorOutliers[l]->OutlierClass){
                    CountModified++;
                    cout << "Label  "<<l+1 << "reclassified from " << VectorOutliers[l]->OutlierClass << "to " << newOT<<endl;
                    VectorOutliers[l]->OutlierClass=newOT;
                }
                LesionSubcategories LesionType=OutlierRefinement(VectorOutliers[l], segment_analysis);
                VectorOutliers[l]->LesionType=LesionType;
                cout <<"Label "<<l+1;
                TxtFileLesion<<endl;
                TxtFileLesion<<"Label "<<l+1<<endl;
                PrintOutlierCharacteristics(VectorOutliers[l],TxtFileLesion,numbmodal,GravImage,PixDim,segment_analysis);
                cout<<"...printed"<<endl;
            }
            cout << "Modified classif is "<<CountModified<<endl;
            
            for (int i=0; i<numel; i++) {
                if (OrderedOutliers[i]>0) {
                    OutlierClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->OutlierClass;
                    LesionSubClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->LesionType;
                }
            }
            nifti_image * OutlierClassifNii=CreateNiiFromArray(OutlierClassif, ImageToUse, numel);
            nifti_set_filenames(OutlierClassifNii, FilenameClassifLac.c_str(), 0, 0);
            nifti_image_write(OutlierClassifNii);
            nifti_image_free(OutlierClassifNii);
            
            nifti_image * LesionClassifNii=CreateNiiFromArray(LesionSubClassif, ImageToUse, numel);
            nifti_set_filenames(LesionClassifNii, FilenameClassifLacFin.c_str(), 0, 0);
            nifti_image_write(LesionClassifNii);
            nifti_image_free(LesionClassifNii);
            
            OutlierClassifNii=NULL;
            for (int l=0; l<numbFinalLesions; l++) {
                delete VectorOutliers[l];
                VectorOutliers[l]=NULL;
            }
            delete [] OutlierClassif;
            OutlierClassif=NULL;
            delete [] GravImage;
            GravImage=NULL;
            delete [] OrderedOutliers;
            return EXIT_SUCCESS;
        }
    }
    
    if (segment_analysis->flag_LCRuleTextFile) {
        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        
        
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameFinalLacOrdered=FilenamePA_b+"FinalLacOrdered_"+FilenamePA_e+".nii.gz";
        string FilenameClassifLac=FilenamePA_b+"ClassifLac_"+FilenamePA_e+".nii.gz";
        string FilenameClassifLacSub=FilenamePA_b+"ClassifLacSub_"+FilenamePA_e+".nii.gz";
        
        SimpleRule * BuiltRule=BuildLCRuleFromTextFile(segment_analysis);
        nifti_image * ImageToUse;
        if (!segment_analysis->flag_inConnect) {
            ImageToUse=OutliersLac(BuiltRule, TreeToAnalyse, segment_analysis);
            //            ImageToUse=OutlierSelected(BuiltRule, TreeToAnalyse, segment_analysis);
        }
        else{
            ImageToUse=ReadFromFilename(segment_analysis->filename_inConnect);
        }
        
        numel=ImageToUse->nvox;
        string FilenameSave=FilenamePA_b+"OutlierLacTestb_"+FilenamePA_e+".txt";
        int * ComponentLabelInit=ComponentLabeling(ImageToUse, 6);
        float VolumeVox=ImageToUse->pixdim[1]*ImageToUse->pixdim[2]*ImageToUse->pixdim[3];
        int * OrderedOutliers=OrderedVolumeLabel(ComponentLabelInit, 3, numel,VolumeVox);
        nifti_image * OrderedOutliersNii=CreateNiiFromArray(OrderedOutliers, ImageToUse, numel);
        nifti_set_filenames(OrderedOutliersNii, FilenameFinalLacOrdered.c_str(), 0, 0);
        nifti_image_write(OrderedOutliersNii);
        vector<Outlier *> VectorOutliers=GetVectorOutliers(OrderedOutliers, TreeToAnalyse, segment_analysis);
        ofstream TxtFileLesion(FilenameSave.c_str());
        int * GravImage=NULL;
        int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
        float PixDim[3];
        if (numbFinalLesions>0) {
            GravImage=new int[3];

            for (int d=0;d<3; d++) {
                PixDim[d]=ImageToUse->pixdim[d+1];
                GravImage[d]=(VectorOutliers[0]->CentreGravity[d]-VectorOutliers[0]->VectorDiffGrav[d]/PixDim[d]);
            }
        }
        for (int l=0; l<numbFinalLesions; l++) {
            cout <<"Label "<<l+1;
            TxtFileLesion<<endl;
            TxtFileLesion<<"Label "<<l+1<<endl;
            PrintOutlierCharacteristics(VectorOutliers[l],TxtFileLesion,numbmodal,GravImage,PixDim,segment_analysis);
            cout<<"...printed"<<endl;
        }
        int * OutlierClassif=CopyArray(OrderedOutliers, numel);
        int * LesionSubClassif=CopyArray(OrderedOutliers, numel);
        for (int i=0; i<numel; i++) {
            if (OrderedOutliers[i]>0) {
                OutlierClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->OutlierClass;
                LesionSubClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->LesionType;
            }
        }
        nifti_image * OutlierClassifNii=CreateNiiFromArray(OutlierClassif, ImageToUse, numel);
        nifti_set_filenames(OutlierClassifNii, FilenameClassifLac.c_str(), 0, 0);
        nifti_image_write(OutlierClassifNii);
        nifti_image_free(OutlierClassifNii);
        OutlierClassifNii=NULL;
        
        nifti_image * LesionClassifNii=CreateNiiFromArray(LesionSubClassif, ImageToUse, numel);
        nifti_set_filenames(LesionClassifNii, FilenameClassifLacSub.c_str(), 0, 0);
        nifti_image_write(LesionClassifNii);
        nifti_image_free(LesionClassifNii);
        for (int l=0; l<numbFinalLesions; l++) {
            delete VectorOutliers[l];
            VectorOutliers[l]=NULL;
        }
        delete [] OutlierClassif;
        OutlierClassif=NULL;
        delete [] LesionSubClassif;
        LesionSubClassif=NULL;
        delete [] GravImage;
        GravImage=NULL;
        nifti_image_free(OrderedOutliersNii);
        delete [] OrderedOutliers;
        return EXIT_SUCCESS;

    }
    
    
    if(segment_analysis->flag_GetConnect && segment_analysis->flag_inLes){
        nifti_image * SegToAnalyse=ReadFromFilename(segment_analysis->filename_InLes);
        numel=SegToAnalyse->nvox;
//        float * SegData=static_cast<float*>(SegToAnalyse->data);
        int * CompLabel_I=ComponentLabeling(SegToAnalyse, 26);
        float VolumeVox=SegToAnalyse->pixdim[1]*SegToAnalyse->pixdim[2]*SegToAnalyse->pixdim[3];
        int * CompLabel=OrderedVolumeLabel(CompLabel_I, 1, numel,VolumeVox);
        nifti_image * LabelsCorrespondence=CreateNiiFromArray(CompLabel, SegToAnalyse, numel);
        FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveTemp=FilenamePA_b+"Connected_"+FilenamePA_e+".nii.gz";
        nifti_set_filenames(LabelsCorrespondence, FilenameSaveTemp.c_str(), 0, 0);
        nifti_image_write(LabelsCorrespondence);
        
    }

    
    
    
    
    if (segment_analysis->flag_CodedOutliers) {
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        
        
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveCard=FilenamePA_b+"OutlierCard_"+FilenamePA_e+".txt";
        string FilenameCodedOutliers=FilenamePA_b+"CodedOutliers_"+FilenamePA_e+".nii.gz";
        string FilenameOriginOutliers=FilenamePA_b+"OriginOutliers_"+FilenamePA_e+".nii.gz";
        string FilenameTreatedCSF=FilenamePA_b+"TreatedCSF_"+FilenamePA_e+".nii.gz";
        string FilenamePotLacunes=FilenamePA_b+"PotLacunes_"+FilenamePA_e+".nii.gz";
        string FilenamePotLacunesB=FilenamePA_b+"PotLacunesB_"+FilenamePA_e+".nii.gz";
        string FilenamePotWMH=FilenamePA_b+"PotWMH_"+FilenamePA_e+".nii.gz";
        string FilenamePotIron=FilenamePA_b+"PotIron_"+FilenamePA_e+".nii.gz";
        string FilenameRejected=FilenamePA_b+"Rejected_"+FilenamePA_e+".nii.gz";
        string FilenameRejectedBorder=FilenamePA_b+"RejectedBorder_"+FilenamePA_e+".nii.gz";
        string FilenamePVCC=FilenamePA_b+"PVCC_"+FilenamePA_e+".nii.gz";
        string FilenameLacunePV=FilenamePA_b+"LacunePV_"+FilenamePA_e+".nii.gz";
        string FilenameFinalLac=FilenamePA_b+"FinalLac_"+FilenamePA_e+".nii.gz";
        string FilenameFinalLacOrdered=FilenamePA_b+"FinalLacOrdered_"+FilenamePA_e+".nii.gz";
        string FilenameClassifLac=FilenamePA_b+"ClassifLac_"+FilenamePA_e+".nii.gz";
        
        
        if (segment_analysis->flag_inConnect) {
            nifti_image * ImageToUse=ReadFromFilename(segment_analysis->filename_inConnect);
            numel=ImageToUse->nvox;
            FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenameSave=FilenamePA_b+"OutlierWMHCardb_"+FilenamePA_e+".txt";
            int * ComponentLabelInit=TranscribeArray<float, int>(static_cast<float*>(ImageToUse->data), numel);
            float VolumeVox=ImageToUse->pixdim[1]*ImageToUse->pixdim[2]*ImageToUse->pixdim[3];
            int * OrderedOutliers=OrderedVolumeLabel(ComponentLabelInit, 3, numel,VolumeVox);
            nifti_image * OrderedOutliersNii=CreateNiiFromArray(OrderedOutliers, ImageToUse, numel);
            nifti_set_filenames(OrderedOutliersNii, FilenameFinalLacOrdered.c_str(), 0, 0);
            nifti_image_write(OrderedOutliersNii);
            
            //
            ////            Temporary to check SP details
            //            if (segment_analysis->flag_MNITransform && segment_analysis->flag_MNITemplate) {
            //                nifti_image * MNITemplate=ReadFromFilename(segment_analysis->filename_MNITemplate);
            //
            //                mat44 MNIMatCorrToInd=MNITemplate->qto_ijk;
            //                float CentreMNI_xyz[4];
            //                for (int d=0; d<3; d++) {
            //                    CentreMNI_xyz[d]=0;
            //                }
            //                CentreMNI_xyz[3]=1;
            //
            //                mat44 ImageMatCorrToWorld=OrderedOutliersNii->qto_xyz;
            //                float * TransfoMatrix=ReadAffineFromTextFile(segment_analysis->filename_MNITransform);
            //                float * MNIMatToInd=new float[16];
            //                TranscribeMat44(MNIMatCorrToInd, MNIMatToInd);
            //                float * ImageMatToWorld=new float[16];
            //                //            int SizeMatVec[4];
            //                //            SizeMatVec[0]=4;
            //                //            SizeMatVec[1]=4;
            //                //            SizeMatVec[2]=4;
            //                //            SizeMatVec[3]=1;
            //                //            float * CentreMNI_ijk=TreeToAnalyse->ProductMatrix(MNIMatToInd, CentreMNI_xyz, SizeMatVec);
            //
            //                int SizeTransfo[4];
            //                SizeTransfo[0]=4;
            //                SizeTransfo[1]=4;
            //                SizeTransfo[2]=4;
            //                SizeTransfo[3]=4;
            //                TranscribeMat44(ImageMatCorrToWorld, ImageMatToWorld);
            //                float * Transfo1=TreeToAnalyse->ProductMatrix(TransfoMatrix, ImageMatToWorld,SizeTransfo);
            //                //            float * TransfoTot=TreeToAnalyse->ProductMatrix(MNIMatToInd, Transfo1, SizeTransfo);
            //                //            Normally the direct multiplication of the coordinates vector will give the corresponding coordinate vector in the MNI template thus allowing for the separation in quadrants
            //
            //                //            Create the quadrant image according to code for quadrants
            //                nifti_image * QuadrantResult=QuadrantTransformation(Transfo1,CentreMNI_xyz, OrderedOutliersNii,TreeToAnalyse->GetL2S());
            //                VentricleBool=VentricleSegmentationDirect(TreeToAnalyse, segment_analysis);
            //                bool * CGMBool=NULL;
            //                if (segment_analysis->flag_inPriorsCGM) {
            //                    nifti_image * CGMPriors=ReadFromFilename(segment_analysis->filename_inPriorsCGM);
            //                    CGMBool=ThresholdArray<float, bool>(static_cast<float*>(CGMPriors->data), 0.5, numel);
            //                    nifti_image_free(CGMPriors);
            //                    CGMPriors=NULL;
            //                }
            //                bool * SPPotential=PotentialSPRegion(QuadrantResult, VentricleBool, CGMBool, 2, TreeToAnalyse->GetMask());
            //                nifti_image * SPNii=CreateNiiFromArray(SPPotential, QuadrantResult, numel);
            //                nifti_set_filenames(SPNii, (g_scratch_dir + "/TestSPPot.nii.gz").c_str(), 0, 0);
            //                nifti_image_write(SPNii);
            //                nifti_image_free(SPNii);
            //                nifti_image_free(QuadrantResult);
            //                delete [] SPPotential;
            //                SPNii=NULL;
            //                QuadrantResult=NULL;
            //                SPPotential=NULL;
            //            }
            
            
            //            //////////
            
            
            
            vector<Outlier *> VectorOutliers=GetVectorOutliers(OrderedOutliers, TreeToAnalyse, segment_analysis);
            ofstream TxtFileLesion(FilenameSave.c_str());
            int * GravImage=NULL;
            int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
            float PixDim[3];
            if (numbFinalLesions>0) {
                GravImage=new int[3];

                for (int d=0;d<3; d++) {
                    PixDim[d]=ImageToUse->pixdim[d+1];
                    GravImage[d]=(VectorOutliers[0]->CentreGravity[d]-VectorOutliers[0]->VectorDiffGrav[d]/PixDim[d]);
                }
            }
            for (int l=0; l<numbFinalLesions; l++) {
                cout <<"Label "<<l+1;
                TxtFileLesion<<endl;
                TxtFileLesion<<"Label "<<l+1<<endl;
                PrintOutlierCharacteristics(VectorOutliers[l],TxtFileLesion,numbmodal,GravImage,PixDim,segment_analysis);
                cout<<"...printed"<<endl;
            }
            int * OutlierClassif=CopyArray(OrderedOutliers, numel);
            for (int i=0; i<numel; i++) {
                if (OrderedOutliers[i]>0) {
                    OutlierClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->OutlierClass;
                }
            }
            nifti_image * OutlierClassifNii=CreateNiiFromArray(OutlierClassif, ImageToUse, numel);
            nifti_set_filenames(OutlierClassifNii, FilenameClassifLac.c_str(), 0, 0);
            nifti_image_write(OutlierClassifNii);
            nifti_image_free(OutlierClassifNii);
            OutlierClassifNii=NULL;
            for (int l=0; l<numbFinalLesions; l++) {
                delete VectorOutliers[l];
                VectorOutliers[l]=NULL;
            }
            delete [] OutlierClassif;
            OutlierClassif=NULL;
            delete [] GravImage;
            GravImage=NULL;
            nifti_image_free(OrderedOutliersNii);
            delete [] OrderedOutliers;
            return EXIT_SUCCESS;
        }
        
        

        
        nifti_image * ResultsCoded = CodedOutliers(TreeToAnalyse, segment_analysis);
        nifti_set_filenames(ResultsCoded, FilenameCodedOutliers.c_str(), 0, 0);
        nifti_image_write(ResultsCoded);
        
        nifti_image * ResultsOrigin=OutliersOrigin(TreeToAnalyse, segment_analysis);
        nifti_set_filenames(ResultsOrigin, FilenameOriginOutliers.c_str(), 0, 0);
        nifti_image_write(ResultsOrigin);
        
        float * OutliersCode=static_cast<float *>(ResultsCoded->data);
        float * OutliersOrigin=static_cast<float *>(ResultsOrigin->data);
        numel=ResultsCoded->nvox;
        float * RejectedZone=new float[numel];
        float * PotentialCSFOut=new float[numel];
        int Dim[3];
        int Shift[3];
        for (int d=0; d<3; d++) {
            Dim[d]=ResultsCoded->dim[d+1];
            if (d!=0) {
                Shift[d]=Dim[d-1]*Shift[d-1];
            }
            else{
                Shift[d]=1;
            }
        }
        int * L2S=TreeToAnalyse->GetL2S();
        for (int i=0; i<numel; i++) {
            RejectedZone[i]=0;
            PotentialCSFOut[i]=0;
            if(L2S[i]<0){
                RejectedZone[i]=1;
                PotentialCSFOut[i]=1;
            }
        }
        
        vector<TreeEM*> LeavesCSF;
        LeavesCSF.push_back(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexCSF));
        LeavesCSF.push_back(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexCSF));
        float * AddedNormRespCSF=AddNormResp(LeavesCSF);
        float * LongAddedNormRespCSF=CreateLong(AddedNormRespCSF, L2S, numel);
        //        nifti_image * CSFNiiFirst=CreateNiiFromArray(LongAddedNormRespCSF, ResultsCoded, numel);
        //        nifti_set_filenames(CSFNiiFirst, (g_scratch_dir + "/TestCSFInit.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(CSFNiiFirst);
        //        nifti_image_free(CSFNiiFirst);
        
        
        vector<TreeEM*> LeavesOut;
        LeavesOut.push_back(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexOut));
        LeavesOut.push_back(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexOut));
        float * AddedNormRespOut=AddNormResp(LeavesOut);
        float * LongAddedNormRespOut=CreateLong(AddedNormRespOut, L2S, numel);
        for (int i=0; i<numel; i++) {
            PotentialCSFOut[i]+=LongAddedNormRespCSF[i]+LongAddedNormRespOut[i];
        }
        if (segment_analysis->flag_inPriorsOut) {
            nifti_image * OutPriors=ReadFromFilename(segment_analysis->filename_inPriorsOut);
            float *OutData=static_cast<float *>(OutPriors->data);
            for (int i=0; i<numel; i++) {
                RejectedZone[i]+=(OutData[i]>0.9);
                PotentialCSFOut[i]+=(OutData[i]>0.9);
            }
            nifti_image_free(OutPriors);
            OutPriors=NULL;
        }
        
        nifti_image * CSFNiiFirst=CreateNiiFromArray(PotentialCSFOut, ResultsCoded, numel);
        //        nifti_set_filenames(CSFNiiFirst, (g_scratch_dir + "/TestCSFInit2.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(CSFNiiFirst);
        
        int * CSFFirst_CC=ComponentLabeling(CSFNiiFirst, 6);
        bool * RejectedZoneBool=TranscribeArray<float, bool>(RejectedZone, numel);
        vector<int> LabelDiscardOutCSF=LabelToDiscard(CSFFirst_CC, RejectedZoneBool, Dim, Shift);
        int * TreatedCSFPre=TreatingLabelFromSetDiscard(CSFFirst_CC, LabelDiscardOutCSF, numel);
        nifti_image_free(CSFNiiFirst);
        CSFNiiFirst=CreateNiiFromArray(TreatedCSFPre, ResultsCoded, numel);
        nifti_set_filenames(CSFNiiFirst, FilenameTreatedCSF.c_str(), 0, 0);
        nifti_image_write(CSFNiiFirst);
        nifti_image_free(CSFNiiFirst);
        CSFNiiFirst=NULL;
        
        bool * RejectedZoneBoolPV=new bool[numel];
        bool * TreatedCSFBool=TranscribeArray<int, bool>(TreatedCSFPre, numel);
        bool * CSFFirstBool=TranscribeArray<int, bool>(CSFFirst_CC, numel);
        XOROperationBool(CSFFirstBool, TreatedCSFBool, RejectedZoneBoolPV, numel);
        
        cout<<"PriorsECSF flag is "<<segment_analysis->flag_inPriorsECSF<<endl;
        if (segment_analysis->flag_inPriorsECSF) {
            nifti_image *ECSFPriors=ReadFromFilename(segment_analysis->filename_inPriorsECSF);
            float *ECSFData=static_cast<float *>(ECSFPriors->data);
            for (int i=0; i<numel; i++) {
                RejectedZone[i]+=(ECSFData[i]>0.1)*LongAddedNormRespCSF[i];
            }
            nifti_image_free(ECSFPriors);
            ECSFPriors=NULL;
        }

        
        bool * RejectedZoneDoubtfulIron=ThresholdArray<float,bool>(RejectedZone,0.15, numel);
        bool * RejectedZoneDoubtfulLacunes=ThresholdArray<float,bool>(RejectedZone, 0.15,numel);
        bool * AcceptedDGM=new bool[numel];
        for (int i=0; i<numel; i++) {
            AcceptedDGM[i]=0;
        }
        if (segment_analysis->flag_inPriorsDGM) {
            nifti_image * DGMPriors=ReadFromFilename(segment_analysis->filename_inPriorsDGM);
            float *DGMData=static_cast<float *>(DGMPriors->data);
            for (int i=0; i<numel; i++) {
                AcceptedDGM[i]=DGMData[i]>0.15?1:0;
            }
            nifti_image_free(DGMPriors);
            DGMPriors=NULL;
        }
        bool * NonAcceptedDGM=OpposeBoolArray(AcceptedDGM, numel);
        bool * RejectedZoneThr=ThresholdArray<float,bool>(RejectedZone,0.15,numel);
        OROperationBool(RejectedZoneBoolPV, NonAcceptedDGM, RejectedZoneDoubtfulIron, numel);
        OROperationBool(RejectedZoneBoolPV, AcceptedDGM, RejectedZoneDoubtfulLacunes, numel);
        
        
        //        OROperationBool(RejectedZoneThr, NonAcceptedDGM, RejectedZoneDoubtfulIron, numel);
        //        OROperationBool(RejectedZoneThr, AcceptedDGM, RejectedZoneDoubtfulLacunes, numel);
        
        //        Creation of the code selection vectors according to intensities
        vector<float> OriginSelection;
        vector<float> CodeSelectionIron;
        vector<float> CodeSelectionLacunes;
        vector<float> CodeSelectionWMH;
        vector<float> CodeSelectionDoubtfulIL;
        vector<float> CodeSelectionLacunesNotEnough;
        OriginSelection.push_back(0);
        OriginSelection.push_back(1);
        OriginSelection.push_back(2);
        
        CodeSelectionDoubtfulIL.push_back(13);
        CodeSelectionDoubtfulIL.push_back(23);
        
        CodeSelectionIron.push_back(14);
        CodeSelectionIron.push_back(24);
        


        CodeSelectionLacunes.push_back(21);
        CodeSelectionLacunes.push_back(12);
        CodeSelectionLacunes.push_back(11);
        
        CodeSelectionLacunesNotEnough.push_back(22);
        
        CodeSelectionWMH.push_back(43);
        CodeSelectionWMH.push_back(42);
        CodeSelectionWMH.push_back(41);
        CodeSelectionWMH.push_back(44);
        CodeSelectionWMH.push_back(33);
        CodeSelectionWMH.push_back(34);
        CodeSelectionWMH.push_back(32);
        CodeSelectionWMH.push_back(31);
        
        vector<TreeEM*> PVVectorGMCSF=SelectionPVClasses( TreeToAnalyse, segment_analysis->IndexGM, segment_analysis->IndexCSF);
        
        float * AddedNormResp=AddNormResp(PVVectorGMCSF);
        float * LongNormResp=CreateLong(AddedNormResp, L2S, numel);
        bool * ThresholdedPV=ThresholdArray<float, bool>(LongNormResp, 0.2, numel);
        
        //        nifti_image * PVNii=CreateNiiFromArray(ThresholdedPV, ResultsCoded, numel);

        
        float * PotentialWMH=OutlierSelection(TreeToAnalyse, OutliersCode, OutliersOrigin, CodeSelectionWMH,  OriginSelection, RejectedZoneThr);
        float * PotentialIron=OutlierSelection(TreeToAnalyse, OutliersCode, OutliersOrigin, CodeSelectionIron, OriginSelection, RejectedZoneDoubtfulIron);
        float * PotentialLacunes=OutlierSelection(TreeToAnalyse, OutliersCode, OutliersOrigin, CodeSelectionLacunes, OriginSelection, RejectedZoneBoolPV);
        float * PotentialLacunesNotEnough=OutlierSelection(TreeToAnalyse, OutliersCode, OutliersOrigin, CodeSelectionLacunesNotEnough, OriginSelection, RejectedZoneBoolPV);
        float * PotentialDoubtfulLacunes=OutlierSelection(TreeToAnalyse, OutliersCode, OutliersOrigin, CodeSelectionDoubtfulIL, OriginSelection, RejectedZoneDoubtfulLacunes);
        float * PotentialDoubtfulIron=OutlierSelection(TreeToAnalyse, OutliersCode, OutliersOrigin, CodeSelectionDoubtfulIL, OriginSelection, RejectedZoneDoubtfulIron);
        float * PotentialLacunesB1=AddElementwise(PotentialDoubtfulLacunes, PotentialLacunes, numel);
        float * PotentialIronB=AddElementwise(PotentialDoubtfulIron, PotentialIron, numel);
        
        nifti_image * LacunesNii=CreateNiiFromArray(PotentialLacunesB1, ResultsOrigin, numel);
        nifti_image * IronNii=CreateNiiFromArray(PotentialIronB, ResultsOrigin, numel);
        nifti_image * WMHNii=CreateNiiFromArray(PotentialWMH, ResultsOrigin, numel);
        
        
        //        Performing the correction of the Iron deposition among the possible lacunes Correction done only on PV
        
        bool * PotentialIronBool=ThresholdArray<float, bool>(PotentialIronB, 0.2, numel);
        //        nifti_image * PotIronRejNii=CreateNiiFromArray(PotentialIronBool, ResultsCoded, numel);
        //        nifti_set_filenames(PotIronRejNii, (g_scratch_dir + "/TestIronRej.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(PotIronRejNii);
        //
        //        nifti_image_free(PotIronRejNii);

        bool * DoubtfulLacunesBool=ThresholdArray<float, bool>(PotentialLacunesNotEnough, 0.2, numel);
        nifti_image * DoubtfulLacunesNii=CreateNiiFromArray(DoubtfulLacunesBool, ResultsCoded, numel);
        int * DoubtfulLacunes_CC=ComponentLabeling(DoubtfulLacunesNii, 6);
        vector<int> LabelsDiscardedIron=LabelToDiscard(DoubtfulLacunes_CC, PotentialIronBool, Dim, Shift);
        int * TreatedPVIronCorrected=TreatingLabelFromSetDiscard(DoubtfulLacunes_CC, LabelsDiscardedIron, numel);
        //        nifti_image * TreatedIronCorrectedNii=CreateNiiFromArray(TreatedPVIronCorrected, ResultsCoded, numel);
        //        nifti_set_filenames(TreatedIronCorrectedNii, (g_scratch_dir + "/TestPVIronCorrected.nii.gz").c_str(), 0, 0);
        float * TreatedPVIronCorrectedFloat=TranscribeArray<int, float>(TreatedPVIronCorrected, numel);
        float * PotentialLacunesB=AddElementwise(TreatedPVIronCorrectedFloat, PotentialLacunesB1, numel);
        nifti_image * LacunesNiiB=CreateNiiFromArray(PotentialLacunesB1, ResultsOrigin, numel);
        
        nifti_set_filenames(LacunesNii, FilenamePotLacunes.c_str(), 0, 0);
        nifti_set_filenames(LacunesNiiB, FilenamePotLacunesB.c_str(), 0, 0);
        nifti_set_filenames(WMHNii, FilenamePotWMH.c_str(), 0, 0);
        nifti_set_filenames(IronNii, FilenamePotIron.c_str(), 0, 0);
        
        nifti_image_write(LacunesNii);
        nifti_image_write(LacunesNiiB);
        nifti_image_write(WMHNii);
        nifti_image_write(IronNii);
        
        //        Building refuse zone for lacunes
        
        bool * RejectedZonePV=CopyArray(RejectedZoneBoolPV, numel);
        //        bool * VentricleBool = VentricleSegmentationDirect(TreeToAnalyse,segment_analysis);
        //        bool * RejectedZonePV=new bool[numel];
        //        OROperationBool(RejectedZoneThr, VentricleBool, RejectedZonePV, numel);
        nifti_image * RejectedPVNii=CreateNiiFromArray(RejectedZonePV, ResultsCoded, numel);
        nifti_set_filenames(RejectedPVNii, FilenameRejected.c_str(), 0, 0);
        nifti_image_write(RejectedPVNii);
        nifti_image_free(RejectedPVNii);
        RejectedPVNii=NULL;
        
        bool * BorderExtRejected=CreateBorderExtFromBool(RejectedZonePV, Dim, Shift, NULL);
        nifti_image * RejectedPVNiiB=CreateNiiFromArray(BorderExtRejected, ResultsCoded, numel);
        nifti_set_filenames(RejectedPVNiiB, FilenameRejectedBorder.c_str(), 0, 0);
        nifti_image_write(RejectedPVNiiB);
        nifti_image_free(RejectedPVNiiB);
        RejectedPVNiiB=NULL;
        
        //        Trying to take into account isolated islands of CSF

        //        bool * CSFSeg=ThresholdArray<float, bool>(LongAddedNormRespCSF, 0.2, numel);
        //        nifti_image * CSFNii=CreateNiiFromArray(CSFSeg, ResultsCoded, numel);
        //        int * CSF_CC=ComponentLabeling(CSFNii, 6);
        //        vector<int> LabelsCSFDiscarded=LabelToDiscard(CSF_CC, RejectedZonePV, Dim, Shift);
        //        int * TreatedCSF=TreatingLabelFromSetDiscard(CSF_CC, LabelsCSFDiscarded, numel);
        //        int * TreatedCSF=TreatingPreLabel(CSF_CC, RejectedZonePV, Dim, Shift);
        
        //        nifti_image * CSFTreatedNii=CreateNiiFromArray(TreatedCSFPre, ResultsCoded, numel);
        //        nifti_set_filenames(CSFTreatedNii, (g_scratch_dir + "/TestCSFTreated.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(CSFTreatedNii);
        //        nifti_image_free(CSFTreatedNii);
        //        CSFTreatedNii=NULL;

        float * LacunePV=AddElementwise(PotentialLacunesB, LongNormResp, numel);
        bool * LacunePVB=ThresholdArray<float, bool>(LacunePV, 0.2, numel);
        nifti_image * LacunePVNii=CreateNiiFromArray(LacunePVB, ResultsCoded, numel);
        int * ComponentPV=ComponentLabeling(LacunePVNii, 6);
        vector<int> LabelsLacPVDiscarded=LabelToDiscard(ComponentPV, RejectedZonePV, Dim, Shift);
        int * TreatedPV=TreatingLabelFromSetDiscard(ComponentPV, LabelsLacPVDiscarded, numel);
        delete [] ComponentPV;
        delete [] LacunePVB;
        delete [] LacunePV;
        //        int * TreatedPV=TreatingPreLabel(ComponentPV,RejectedZonePV,Dim,Shift);
        nifti_image * PVCCNii=CreateNiiFromArray(TreatedPV, ResultsCoded, numel);
        nifti_set_filenames(PVCCNii, FilenamePVCC.c_str(), 0, 0);
        nifti_set_filenames(LacunePVNii, FilenameLacunePV.c_str(), 0, 0);
        nifti_image_write(LacunePVNii);
        nifti_image_write(PVCCNii);
        delete [] LongNormResp;
        delete [] AddedNormResp;
        delete [] ThresholdedPV;
        delete [] PotentialLacunesB;
        nifti_image_free(LacunePVNii);
        nifti_image_free(PVCCNii);
        

        //        nifti_image_write(TreatedIronCorrectedNii);
        //        nifti_image_free(TreatedIronCorrectedNii);
        int * TreatedTotLac=AddElementwiseInt(TreatedCSFPre, TreatedPV, numel);

        nifti_image * TreatedTotNii=CreateNiiFromArray(TreatedTotLac, ResultsCoded, numel);
        //        nifti_set_filenames(TreatedTotNii, (g_scratch_dir + "/TestFinalPre.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(TreatedTotNii);
        int * ComponentLacTot=ComponentLabeling(TreatedTotNii, 6);
        nifti_image * FinalLacNii=CreateNiiFromArray(ComponentLacTot, ResultsCoded, numel);
        nifti_image_free(TreatedTotNii);
        nifti_set_filenames(FinalLacNii, FilenameFinalLac.c_str(), 0, 0);
        nifti_image_write(FinalLacNii);
        
        nifti_image_free(ResultsCoded);
        nifti_image_free(ResultsOrigin);
        nifti_image_free(LacunesNii);
        nifti_image_free(WMHNii);
        nifti_image_free(IronNii);
        delete [] RejectedZoneThr;
        delete [] RejectedZone;
        delete [] ComponentLacTot;
        delete [] TreatedTotLac;
        delete [] PotentialIronBool;
        
        

        int * ComponentLabelInit=TranscribeArray<float, int>(static_cast<float*>(FinalLacNii->data), numel);
        float VolumeVox=FinalLacNii->pixdim[1]*FinalLacNii->pixdim[2]*FinalLacNii->pixdim[3];
        int * OrderedOutliers=OrderedVolumeLabel(ComponentLabelInit, 3, numel,VolumeVox);
        nifti_image * OrderedOutliersNii=CreateNiiFromArray(OrderedOutliers, FinalLacNii, numel);
        nifti_set_filenames(OrderedOutliersNii, FilenameFinalLacOrdered.c_str(), 0, 0);
        nifti_image_write(OrderedOutliersNii);
        vector<Outlier *> VectorOutliers=GetVectorOutliers(OrderedOutliers, TreeToAnalyse, segment_analysis);
        ofstream TxtFileLesion(FilenameSaveCard.c_str());
        int * GravImage=NULL;
        int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
        float PixDim[3];
        if (numbFinalLesions>0) {
            GravImage=new int[3];
            
            for (int d=0;d<3; d++) {
                PixDim[d]=FinalLacNii->pixdim[d+1];
                GravImage[d]=(VectorOutliers[0]->CentreGravity[d]-VectorOutliers[0]->VectorDiffGrav[d]/PixDim[d]);
            }
        }
        for (int l=0; l<numbFinalLesions; l++) {
            cout <<"Label "<<l+1;
            TxtFileLesion<<endl;
            TxtFileLesion<<"Label "<<l+1<<endl;
            PrintOutlierCharacteristics(VectorOutliers[l],TxtFileLesion,numbmodal,GravImage,PixDim,segment_analysis);
            cout<<"...printed"<<endl;
        }
        int * OutlierClassif=CopyArray(OrderedOutliers, numel);
        for (int i=0; i<numel; i++) {
            if (OrderedOutliers[i]>0) {
                OutlierClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->OutlierClass;
            }
        }
        nifti_image * OutlierClassifNii=CreateNiiFromArray(OutlierClassif, FinalLacNii, numel);
        nifti_set_filenames(OutlierClassifNii, FilenameClassifLac.c_str(), 0, 0);
        nifti_image_write(OutlierClassifNii);
        nifti_image_free(OutlierClassifNii);
        OutlierClassifNii=NULL;
        for (int l=0; l<numbFinalLesions; l++) {
            delete VectorOutliers[l];
            VectorOutliers[l]=NULL;
        }
        delete [] OutlierClassif;
        delete [] TreatedPV;
        delete [] TreatedPVIronCorrected;
        delete [] TreatedPVIronCorrectedFloat;
        delete [] DoubtfulLacunes_CC;
        delete [] PotentialIronB;
        delete [] PotentialIron;
        delete [] PotentialLacunesB1;
        delete [] PotentialLacunesB;
        delete [] PotentialLacunes;
        delete [] PotentialWMH;
        delete [] PotentialLacunesNotEnough;
        delete [] LongAddedNormRespCSF;
        delete [] LongAddedNormRespOut;
        delete [] CSFFirst_CC;
        delete [] TreatedCSFPre;
        delete [] PotentialCSFOut;
        OutlierClassif=NULL;
        TreatedPV=NULL;
        delete [] GravImage;
        GravImage=NULL;
        delete [] OrderedOutliers;
        OrderedOutliers=NULL;
        nifti_image_free(OrderedOutliersNii);
        delete [] OrderedOutliers;
        delete [] ComponentLabelInit;
        ComponentLabelInit=NULL;
        nifti_image_free(FinalLacNii);
        FinalLacNii=NULL;
        nifti_image_free(DoubtfulLacunesNii);
        DoubtfulLacunesNii=NULL;
        if (CorrectionJuxta!=NULL) {
            nifti_image_free(CorrectionJuxta);
            CorrectionJuxta=NULL;
        }
        
        
        delete TreeToAnalyse;
        TreeToAnalyse=NULL;
        
        
        
        return EXIT_SUCCESS;
    }
    
    
    if (segment_analysis->flag_inFPTPReclassif && !segment_analysis->flag_TextFile) { // Analysis with reclassification of lesions

        if ( !segment_analysis->flag_inLes || !segment_analysis->flag_inConnect) {
            cout<<"No files for possible change after reclassification"<<endl;
            return EXIT_FAILURE;
        }
        else{
            nifti_image * InitLesSeg=ReadFromFilename(segment_analysis->filename_InLes);
            nifti_image * Connect=ReadFromFilename(segment_analysis->filename_inConnect);
            int Dim[3];
            for (int d=0; d<3; d++) {
                Dim[d]=Connect->dim[d+1];
            }
            int numel=Connect->nx*Connect->ny*Connect->nz;
            float * ConnectData=static_cast<float *>(Connect->data);
            float * InitData=static_cast<float*>(InitLesSeg->data);
            float * InitCorrData=CopyArray(InitData,numel);

            int * ConnectInt=TranscribeArray<float, int>(ConnectData, numel);
            map<int,int> DecisionLesion=ReadReclassifDecisionFromFile(segment_analysis->filename_inFPTPReclassif);

            int nlsize=DecisionLesion.size();
            //        Checking if minimum files are there
            cout << nlsize << " is nlsize" << endl;
            for(int i=0;i<numel;i++){
                if(ConnectInt[i]>0 && ConnectInt[i]<nlsize){
                    std::map<int,int>::const_iterator pos = DecisionLesion.find(ConnectInt[i]);
                    if (pos == DecisionLesion.end()) {
                        //handle the error
                    } else {
                        int value = pos->second;
                        InitCorrData[i]=InitCorrData[i]>=1?value*InitCorrData[i]:value;

                    }
                                    }
            }


//            for (int l=0; l<nlsize; l++) {
//                if (!DecisionLesion[l]) {
//                    vector<int> IndicesToModify=GetListIndicesLabel(ConnectInt, Dim, l+1);
//                    int sizeIdx=IndicesToModify.size();
//                    for (int i=0; i<sizeIdx; i++) {
//                        InitCorrData[IndicesToModify[i]]=0;
//                    }
//                }
//            }


            string FilenamePA=nifti_makebasename(segment_analysis->filename_inConnect);
            int Index=FilenamePA.find_last_of('/');
            string OptText="";
            if (segment_analysis->flag_inOptionText){
                OptText=segment_analysis->inOptionText;
            }
            nifti_image * InitCorrLesSeg = CreateNiiFromArray(InitCorrData,Connect,numel);
            string FilenamePA_b=FilenamePA.substr(0,Index+1);
            string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenameReclassif=FilenamePA_b+"ReclassifCorr"+OptText+"_"+FilenamePA_e+".nii.gz";
            nifti_set_filenames(InitCorrLesSeg, FilenameReclassif.c_str(), 0, 0);
            nifti_image_write(InitCorrLesSeg);
            nifti_image_free(InitCorrLesSeg);
            nifti_image_free(Connect);
            nifti_image_free(InitLesSeg);
            delete [] ConnectInt;
            InitLesSeg=NULL;
            InitCorrLesSeg=NULL;
            ConnectInt=NULL;
            Connect=NULL;
        }
        return EXIT_SUCCESS;
    }

    
    if(segment_analysis->flag_mean){
        vector<string> VectorFiles;
        int sizeVec=segment_analysis->filename_ImagesMean.size();
        for(int m=0;m<sizeVec;m++){
            VectorFiles.push_back(segment_analysis->filename_ImagesMean[m]);
        }
        vector<nifti_image *> VectorImagesToMean=ReadFromFilenamesVector(VectorFiles);
        bool * Mask=NULL;
        nifti_image * MaskNii=NULL;
        if(segment_analysis->flag_mask){
            MaskNii=ReadFromFilename(segment_analysis->filename_mask);
            Binarisation(MaskNii);
            Mask=static_cast<bool *>(MaskNii->data);
        }
        if(segment_analysis->filename_saveMean == NULL){
            if (segment_analysis->name_outputDir != NULL){
                string FilenameSaveMean = string(segment_analysis->name_outputDir) + "/Average.nii.gz";
                segment_analysis->filename_saveMean= const_cast<char*>(FilenameSaveMean.c_str());
            }
            else{
                string FilenamePA = nifti_makebasename(segment_analysis->filename_ImagesMean[0]);
                int Index = FilenamePA.find_last_of('/');
                string FilenamePA_b = FilenamePA.substr(0,Index+1);
                string FilenameSaveMean = FilenamePA_b + "Average.nii.gz";
                segment_analysis->filename_saveMean= const_cast<char*>(FilenameSaveMean.c_str());
            }
        }
        if(CheckCompatibilityVectorNii(VectorImagesToMean)){
            nifti_image * ResultMean=MakeMeanFromVectorImages(VectorImagesToMean,Mask);
            nifti_set_filenames(ResultMean, segment_analysis->filename_saveMean, 0, 0);
            nifti_image_write(ResultMean);
            nifti_image_free(ResultMean);
        }
        else {
            cout<<"Impossible to do mean of incompatible images or void vector"<<endl;
        }
        if(Mask!=NULL){
            nifti_image_free(MaskNii);
            MaskNii=NULL;
            Mask=NULL;
        }
        return EXIT_SUCCESS;
    }

    if(segment_analysis->flag_CorrectConnect){
        vector<nifti_image*> ConnectNiiVect;
        vector<nifti_image*> CorrectNiiVect;
        int numbLes = segment_analysis->filename_vectorLes.size();

        for (int i=0; i<numbLes;i++){
            CorrectNiiVect.push_back(ReadFromFilename(segment_analysis->filename_vectorLes[i]));
            ConnectNiiVect.push_back(ReadFromFilename(segment_analysis->filename_vectorRef[i]));
        }
        int numel  = ConnectNiiVect[0]->nvox;
        // Creation of the correction connected
        vector<float*> CorrectedConnect;
        vector<float*> CorrectLes;
        vector<float*> ConnectLes;
        for (int i=0;i<numbLes;i++){
            float * CorrectTemp = static_cast<float*>(CorrectNiiVect[i]->data);
            float * ConnectTemp = static_cast<float*>(ConnectNiiVect[i]->data);
            float * CorrectConnectTemp = new float[numel];
            for (int n=0;n<numel;n++){
                CorrectConnectTemp[n]=0;
                if(ConnectTemp[n]>0 && CorrectTemp[n]==0){
                    CorrectConnectTemp[n]=float(ConnectTemp[n]);
                }

            }
            CorrectedConnect.push_back(CorrectConnectTemp);
            CorrectLes.push_back(CorrectTemp);
            ConnectLes.push_back(ConnectTemp);
        }
        cout << CorrectLes[0] << " " <<  CorrectLes[1] << numel <<endl;
        float * SumAllCorrect = AddElementwise(CorrectLes, numel, NULL);
        float * ValidSumAll = ThresholdArray<float, float>(SumAllCorrect, numbLes*0.5+0.00001, numel);
        float * NonValidAll = SubtractArray(SumAllCorrect, ValidSumAll, numel);
        float * SumAllCorrectConnect = AddElementwise(CorrectedConnect,numel,NULL);
        // Identification of corrected in one but not in one of the others
        for(int l=0;l<numbLes;l++){
            set<float> SetValuesToCorrect;
            set<float> SetValuesUniqueCorrect;
            set<float> SetValuesAll;
            set<float> SetValuesUnique;
            float * ToFurtherCorrect = new float[numel];
            float * ToRecheck = new float[numel];
            float * SmallCorrect = new float[numel];
            float * SumOthersCorrect = SubtractArray(SumAllCorrectConnect, CorrectedConnect[l],numel);
            for(int i=0;i<numel;i++){
                SmallCorrect[i]=0;

                if(SumOthersCorrect[i]>0 & CorrectLes[l][i]>0){
                    SetValuesToCorrect.insert(ConnectLes[l][i]);
                    SmallCorrect[i] = ConnectLes[l][i];
                }
                if(ValidSumAll[i]>0){
                    SetValuesAll.insert(ConnectLes[l][i]);
                }
                if(NonValidAll[i]>0){
                    SetValuesUnique.insert(ConnectLes[l][i]);
                }
                if(SumOthersCorrect[i]==0 & CorrectedConnect[l][i]>0){
                    SetValuesUniqueCorrect.insert(CorrectedConnect[l][i]);
                }
            }
            for(int i=0;i<numel;i++){
                ToFurtherCorrect[i] = 0;
                ToRecheck[i] = 0;
                if(SetValuesToCorrect.find(ConnectLes[l][i]) != SetValuesToCorrect.end() && SetValuesAll.find(ConnectLes[l][i]) == SetValuesAll.end()){
                    ToFurtherCorrect[i] = ConnectLes[l][i];
                }
                if(SetValuesUnique.find(ConnectLes[l][i]) != SetValuesUnique.end() && CorrectLes[l][i]>0 && SetValuesAll.find(ConnectLes[l][i]) == SetValuesAll.end() ){
                    ToRecheck[i] = ConnectLes[l][i];
                }
            }
            nifti_image * NiiToCheck = CreateNiiFromArray(ToRecheck, ConnectNiiVect[l], numel);
            nifti_image * NiiFurtherCorrect = CreateNiiFromArray(ToFurtherCorrect, ConnectNiiVect[l],numel);
            nifti_image * NiiSmallCorrect = CreateNiiFromArray(SmallCorrect,ConnectNiiVect[l],numel );
            string FilenamePA = nifti_makebasename(segment_analysis->filename_vectorRef[l]);
            int Index = FilenamePA.find_last_of('/');
            string FilenamePA_b = FilenamePA.substr(0,Index+1);
            string FilenamePA_e = FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenameSave = FilenamePA_b+"FurtherConnectCorr_" + FilenamePA_e+".nii.gz";
            string FilenameSave2 = FilenamePA_b+"SmallCorrect_"+FilenamePA_e+".nii.gz";
            string FilenameSave3 = FilenamePA_b+"ToRecheck_"+FilenamePA_e+".nii.gz";
            nifti_set_filenames(NiiSmallCorrect,FilenameSave2.c_str(),0,0);
            nifti_image_write(NiiSmallCorrect);
            nifti_image_free(NiiSmallCorrect);
            nifti_set_filenames(NiiToCheck,FilenameSave3.c_str(),0,0);
            nifti_image_write(NiiToCheck);
            nifti_image_free(NiiToCheck);
            nifti_set_filenames(NiiFurtherCorrect, FilenameSave.c_str(),0,0);
            nifti_image_write(NiiFurtherCorrect);
            nifti_image_free(NiiFurtherCorrect);
            delete [] SumOthersCorrect;
            SumOthersCorrect = NULL;
            delete [] ToRecheck;
            ToRecheck=NULL;
            delete [] ToFurtherCorrect;
            ToFurtherCorrect=NULL;
            delete [] SmallCorrect;
            SmallCorrect = NULL;
        }
        delete [] SumAllCorrectConnect;

        SumAllCorrectConnect = NULL;
        for (int l=0;l<numbLes;l++){
            delete [] CorrectedConnect[l];
            CorrectedConnect[l] = NULL;
            nifti_image_free(CorrectNiiVect[l]);
            nifti_image_free(ConnectNiiVect[l]);
            CorrectNiiVect[l] = NULL;
            ConnectNiiVect[l] = NULL;
        }

      return EXIT_SUCCESS;
    }


    if(segment_analysis->flag_RG && segment_analysis->flag_inToAnalyse && segment_analysis->flag_inDataComp && segment_analysis->flag_mask && segment_analysis->flag_inLes){
        nifti_image * TemporaryLesNii =  ReadFromFilename(segment_analysis->filename_inToAnalyse);
        nifti_image * AuthoLesNii = ReadFromFilename(segment_analysis->filename_InLes);
        nifti_image * DataToUseNii = ReadFromFilename(segment_analysis->filename_inDataComp);
        nifti_image * MaskNii = ReadFromFilename(segment_analysis->filename_mask);
        int Dim[3];
        vector<int> DimVec;
        int Shift[3];
        float PixDim[3];
        for (int d=0; d<3; d++) {
            Dim[d]=TemporaryLesNii->dim[d+1];
            DimVec.push_back(Dim[d]);
            PixDim[d]=TemporaryLesNii->pixdim[d+1];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];
        float * TemporaryLes = static_cast<float*>(TemporaryLesNii->data);
        float * AuthoLes = static_cast<float*> (AuthoLesNii->data);
        float * DataToUse = static_cast<float*>(DataToUseNii->data);
        MaskNii=Binarisation(MaskNii);
        bool * Mask = static_cast<bool*>(MaskNii->data);
        int numel = TemporaryLesNii->nvox;
        int numelmasked = CountNonZero(Mask, numel);
        cout << "numel and numelmasked are "<< numel << " "<< numelmasked;
        numelmasked = 0;
        int * L2S = MakeL2S(Mask, Dim);
        int * S2L = MakeS2L(Mask, Dim, numelmasked);
        float * AbsSumGrad=new float[numel];
        float * SAbsSumGrad=new float[numel];
        for (int i=0; i<numel;i++){
            AbsSumGrad[i] = 0;
            SAbsSumGrad[i] = 0;
        }

        float * ShortData=CreateShort(DataToUse, S2L, numelmasked);
        for (int d=0;d<3;d++){
            float * GradTemp = GradientImage(ShortData,d,0,Dim,Shift,PixDim,L2S,S2L,numelmasked);
            float * SGradTemp = SecondDer(ShortData, d,Dim,Shift,PixDim, L2S, S2L, numelmasked);
            float * AbsGradTemp = Absolute(GradTemp,numelmasked);
            float * SAbsGradTemp = Absolute(SGradTemp,numelmasked);
            cout << GetMin(AbsGradTemp, numelmasked) << "min absgrad for "<< d <<endl;
            float * AbsGradTempLong = CreateLong(AbsGradTemp, L2S, numel);
            float * SAbsGradTempLong = CreateLong(SAbsGradTemp, L2S, numel);
            cout << GetMin(AbsGradTempLong, numel) << "min absgrad for "<< d <<endl;
            AddElementwiseInPlace(AbsSumGrad, AbsGradTempLong,numel);
            AddElementwiseInPlace(SAbsSumGrad, SAbsGradTempLong,numel);
            cout << GetMin(AbsSumGrad, numel) << "min absgrad for "<< d <<endl;
            delete[] GradTemp;
            delete[] AbsGradTemp;
            delete[] AbsGradTempLong;
            delete[] SGradTemp;
            delete[] SAbsGradTemp;
            delete[] SAbsGradTempLong;
        }
        float * PercDiffBlurred = new float[numel];
        delete [] ShortData;
        float * BlurredGradient = GaussianBlurring(AbsSumGrad, 1, DimVec, 0);
        cout << GetMax(BlurredGradient, numel) << "max blurgrad "<<endl;
        for (int i=0;i<numel;i++){
            PercDiffBlurred[i]=(BlurredGradient[i]-AbsSumGrad[i])/(AbsSumGrad[i]+0.000001);
        }
        float * AbsData=Absolute(DataToUse, numel);
        float * MulGrad = new float[numel];
        float * MapToFill = CopyArray(TemporaryLes, numel);
        for (int i=0;i<numel;i++){
            MulGrad[i] = AbsData[i] * AbsSumGrad[i] * AuthoLes[i] ;

        }
        float MeanTest = GetMeanData<float, float>(AbsData, TemporaryLes, numel);
        cout << "Mean is over initial "<< MeanTest<<endl;
        float * Filling = ThresholdArray<float,float>(MulGrad, segment_analysis->thresh_RG, numel);
        MultiplyElementwiseInPlace(AbsSumGrad, TemporaryLes,numel);
        float * ToAddAbsGrad  = ThresholdArray<float,float>(AbsSumGrad, 3, numel );

        AddElementwiseInPlace(Filling, TemporaryLes, numel);
        AddElementwiseInPlace(Filling,ToAddAbsGrad,numel);
        for (int i=0;i<numel;i++){
            if (ToAddAbsGrad[i]>0){
                int ListNeighbours[18];
                GetListNeighbours_bis(ListNeighbours, i, Dim, Shift, 18);
                for (int n=0;n<18;n++){
                    if (AbsData[ListNeighbours[n]] < AbsData[i] & AbsData[ListNeighbours[n]]<0.75 * MeanTest){
                        Filling[ListNeighbours[n]] = 0;
                    }
                }
            }
        }



        RegionGrowing(MapToFill, Dim, Shift, 6,Filling,0.75,TemporaryLesNii);
        cout << GetMin(AbsSumGrad, numel) << "min full" << GetMax(AbsSumGrad,numel)<<endl;
        nifti_image * GradAbsNii = CreateNiiFromArray(AbsSumGrad, TemporaryLesNii, numel);
        nifti_image * SGradAbsNii = CreateNiiFromArray(SAbsSumGrad, TemporaryLesNii, numel);
        nifti_image * GradBlurredNii = CreateNiiFromArray(BlurredGradient, TemporaryLesNii,numel);
        nifti_image * PercDiffBlurredNii = CreateNiiFromArray(PercDiffBlurred, TemporaryLesNii, numel);
        nifti_image * MapToFillNii = CreateNiiFromArray(MapToFill, TemporaryLesNii, numel);
                string FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
                int Index=FilenamePA.find_last_of('/');
                FilenamePA_b=FilenamePA.substr(0,Index+1);
                if(segment_analysis->flag_outputDir){
                    FilenamePA_b=segment_analysis->name_outputDir;
                }
                FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
                string FilenameSaveGrad=FilenamePA_b+"AbsGrad_"+FilenamePA_e+".nii.gz";
                string FilenameSaveSGrad=FilenamePA_b+"SAbsGrad_"+FilenamePA_e+".nii.gz";
                string FilenameSaveGradBlurred=FilenamePA_b+"BlurGrad_"+FilenamePA_e+".nii.gz";
                string FilenameSaveGradPerc=FilenamePA_b+"PercGrad_"+FilenamePA_e+".nii.gz";
                string FilenameSaveRG=FilenamePA_b+"RG_"+FilenamePA_e+".nii.gz";
                nifti_set_filenames(GradAbsNii,FilenameSaveGrad.c_str(),0,0);
                nifti_image_write(GradAbsNii);
                nifti_set_filenames(SGradAbsNii,FilenameSaveSGrad.c_str(),0,0);
                nifti_image_write(SGradAbsNii);
                nifti_set_filenames(GradBlurredNii,FilenameSaveGradBlurred.c_str(),0,0);
                nifti_image_write(GradBlurredNii);
                nifti_set_filenames(PercDiffBlurredNii,FilenameSaveGradPerc.c_str(),0,0);
                nifti_image_write(PercDiffBlurredNii);
                nifti_set_filenames(MapToFillNii,FilenameSaveRG.c_str(),0,0);
                nifti_image_write(MapToFillNii);
                nifti_image_free(MapToFillNii);
                delete [] MulGrad;
                delete [] Filling;
                delete [] AbsData;
                nifti_image_free(TemporaryLesNii);
                nifti_image_free(MaskNii);
                nifti_image_free(DataToUseNii);
                nifti_image_free(GradAbsNii);
                nifti_image_free(SGradAbsNii);
                nifti_image_free(GradBlurredNii);
                nifti_image_free(PercDiffBlurredNii);

                return EXIT_SUCCESS;

    }

    if(segment_analysis->flag_RG && segment_analysis->flag_inToAnalyse && segment_analysis->flag_maskMatch){
        nifti_image * InitialSeeds = ReadFromFilename(segment_analysis->filename_inToAnalyse);
        nifti_image * MappingFilling = ReadFromFilename(segment_analysis->filename_maskMatch);
        float * MapToFill = static_cast<float*>(InitialSeeds->data);
        float * Filling = static_cast<float*>(MappingFilling->data);
        int numel = InitialSeeds->nvox;
//        Correct InitialSeeds according to threshold
        for (int i =0;i<numel;i++){
            if(MapToFill[i]>0 && Filling[i]<segment_analysis->thresh_RG){
                MapToFill[i]=0;
            }
        }
        int Dim[3];
        int Shift[3];
        for (int d=0; d<3; d++) {
            Dim[d]=InitialSeeds->dim[d+1];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];
        RegionGrowing(MapToFill, Dim, Shift, 6,Filling,segment_analysis->thresh_RG,InitialSeeds);
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSave=FilenamePA_b+"RegionGrown_"+FilenamePA_e+".nii.gz";
        nifti_image * GrownNii = CreateNiiFromArray(MapToFill,InitialSeeds,numel);
        nifti_set_filenames(GrownNii,FilenameSave.c_str(),0,0);
        nifti_image_write(GrownNii);
        nifti_image_free(InitialSeeds);
        nifti_image_free(MappingFilling);
        nifti_image_free(GrownNii);
        return EXIT_SUCCESS;
    }

    
    if (segment_analysis->flag_IO && TreeToAnalyse!=NULL) {
        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveIOSeg=FilenamePA_b+"IO_"+FilenamePA_e+".nii.gz";
        string FilenameSaveIOMinMax=FilenamePA_b+"IO_"+FilenamePA_e+".txt";
        string FilenameSaveISeg=FilenamePA_b+"I_"+FilenamePA_e+".nii.gz";
        nifti_image * IOSeg=nifti_copy_nim_info(TreeToAnalyse->GetDataImage());
        nifti_image * ISeg=nifti_copy_nim_info(TreeToAnalyse->GetDataImage());
        IOSeg->dim[0]=4;
        IOSeg->dim[4]=2;
        ISeg->dim[0]=3;
        ISeg->dim[4]=1;
        nifti_update_dims_from_array(IOSeg);
        nifti_update_dims_from_array(ISeg);
        IOSeg->data=(void *) calloc(IOSeg->nvox, sizeof(float));
        ISeg->data=(void *) calloc(ISeg->nvox, sizeof(float));
        float * IOSegData=static_cast<float *>(IOSeg->data);
        float * ISegData=static_cast<float *>(ISeg->data);
        float * NormRespInlier=TreeToAnalyse->GetNodeInlier()->GetNormResp();
        //        SaveTmpResult_short(NormRespInlier, g_scratch_dir + "/TestInlier.nii.gz", ISeg, TreeToAnalyse->GetL2S());
        int numelmasked=TreeToAnalyse->GetNumberMaskedElements();
        float maxInlier=GetMax(NormRespInlier, numelmasked);
        cout<<maxInlier<<" is maxInlier "<<endl;
        float * NormRespOutlier=TreeToAnalyse->GetNodeOutlier()->GetNormResp();
        vector<float *> ResultIOVector;
        ResultIOVector.push_back(NormRespInlier);
        ResultIOVector.push_back(NormRespOutlier);
        numel=TreeToAnalyse->GetNumberElements();
        SaveTmpResult_short(ResultIOVector, FilenameSaveIOSeg, IOSeg, TreeToAnalyse->GetL2S());
        cout<<"Initialising to 0"<<endl;
        for (int i=0; i<numel; i++) {
            IOSegData[i]=0;
            IOSegData[i+numel]=0;
            ISegData[i]=0;
        }
//        numelmasked=TreeToAnalyse->GetNumberMaskedElements();
        //        int * S2L=TreeToAnalyse->GetS2L();
        int * L2S=TreeToAnalyse->GetL2S();
        cout<<"Beginning filling"<<endl;
        int j=0;
        for(int i=0;i<numel;i++){
            if(L2S[i]>=0){
                IOSegData[i]=NormRespInlier[j];
                IOSegData[i+numel]=NormRespOutlier[j];
                ISegData[i]=(NormRespOutlier[j]>0.1)|| NormRespInlier[j]<0.01?0:1;
                j++;
            }
        }
        //        for (int i=0; i<numelmasked; i++) {
        //            IOSegData[S2L[i]]=NormRespInlier[i];
        //            IOSegData[S2L[i]+numel]=NormRespOutlier[i];
        //        }
        nifti_set_filenames(IOSeg, FilenameSaveIOSeg.c_str(), 0, 0);
        nifti_image_write(IOSeg);
        nifti_set_filenames(ISeg, FilenameSaveISeg.c_str(), 0, 0);
        nifti_image_write(ISeg);
        //        Obtaining Min and Max corresponding to modalities for the inlier segmented part.
        nifti_image * HardSegIO=HardSegmentation(IOSeg);
        //        nifti_set_filenames(HardSegIO, (g_scratch_dir + "/HardSegIO.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(HardSegIO);
        nifti_image * MaskIO=HardSegmentationThreshold(HardSegIO, -0.5,2);
        //        nifti_set_filenames(MaskIO, (g_scratch_dir + "/MaskIO.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(MaskIO);
        nifti_image * MaskO=HardSegmentationThreshold(HardSegIO, 0.5,2);
        //        nifti_set_filenames(MaskO, (g_scratch_dir + "/MaskO.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(MaskO);
        float * MaskI=SubtractArray(static_cast<float *>(MaskIO->data), static_cast<float *>(MaskO->data), numel);
        nifti_image_free(MaskIO);
        MaskIO=NULL;
        nifti_image_free(MaskO);
        MaskO=NULL;
        ofstream TxtFile(FilenameSaveIOMinMax.c_str());
        float * ImageData=static_cast<float *>(TreeToAnalyse->GetDataImage()->data);
        for (int m=0; m<numbmodal; m++) {
            TxtFile << ModalitiesCode[Modalities[m]-1]<< " ";
            float MinToPush=GetMin(&ImageData[m*numel], MaskI,numel);
            float MaxToPush=GetMaxArrayMasked(&ImageData[m*numel], MaskI,numel);
            TxtFile << MinToPush<<" ";
            TxtFile<< MaxToPush <<" "<<endl;
        }
        nifti_image_free(IOSeg);
        IOSeg=NULL;
        nifti_image_free(ISeg);
        ISeg=NULL;
        nifti_image_free(HardSegIO);
        HardSegIO=NULL;
        if (MaskI!=NULL) {
            delete [] MaskI;
            MaskI=NULL;
        }
        cout << "IO seg obtained "<<endl;
        if(!segment_analysis->flag_match && !segment_analysis->flag_OtherSeg&& !segment_analysis->flag_inPriorsECSF && !segment_analysis->flag_LapAnalysis && !segment_analysis->flag_connect  && !segment_analysis->flag_WMMaps){
            delete TreeToAnalyse;
            return EXIT_SUCCESS;
        }
    }
    
    if(segment_analysis->flag_match){ // Histogram matching to be performed
        nifti_image * MaskMatch=NULL;
        bool * MaskMatchBool=NULL;
        nifti_image * MaskToApply=NULL;
        bool * MaskToApplyBool=NULL;
        if (segment_analysis->flag_maskMatch) {
            MaskMatch=ReadFromFilename(segment_analysis->filename_maskMatch);
            Binarisation(MaskMatch);
            MaskMatchBool=static_cast<bool*>(MaskMatch->data);
            numel = MaskMatch->nvox;
            cout << CountNonZero(MaskMatchBool,numel) << endl;
        }
        else if(TreeToAnalyse!=NULL){
            if(TreeToAnalyse->GetFlagOutliers()>=3){
                numel=TreeToAnalyse->GetNumberElements();
                nifti_image * BasicImage=nifti_copy_nim_info(TreeToAnalyse->GetDataImage());
                BasicImage->dim[0]=3;
                BasicImage->dim[4]=1;
                BasicImage->dim[5]=1;
                BasicImage->data=(void *) calloc(BasicImage->nvox, sizeof(float));
                nifti_update_dims_from_array(BasicImage);
                float * InlierInit=TreeToAnalyse->GetNodeInlier()->GetNormResp();
                float * LongInlier=CreateLong(InlierInit, TreeToAnalyse->GetL2S(), numel);
                nifti_image * InitInlierImage=CreateNiiFromArray(LongInlier, BasicImage, numel);
                nifti_image * MaskMatch=HardSegmentationThreshold(InitInlierImage, 0.5);
                if(LongInlier!=NULL){
                    delete [] LongInlier;
                    LongInlier=NULL;
                }
                if(InitInlierImage!=NULL){
                    nifti_image_free(InitInlierImage);
                    InitInlierImage=NULL;
                }
                Binarisation(MaskMatch);
                MaskMatchBool=static_cast<bool *>(MaskMatch->data);
                nifti_image_free(BasicImage);
                BasicImage=NULL;
            }

            MaskToApply=TreeToAnalyse->GetMask();
            MaskToApplyBool=static_cast<bool *>(MaskToApply->data);
        }
        

        if(segment_analysis->flag_mask){
            MaskToApply=ReadFromFilename(segment_analysis->filename_mask);
            Binarisation(MaskToApply);
            MaskToApplyBool=static_cast<bool *>(MaskToApply->data);
        }
        else if(MaskMatchBool !=NULL){
            MaskToApplyBool = CopyArray(MaskMatchBool, numel);
            cout << "Using the same inliers and application mask " << endl;
            cout << CountNonZero(MaskToApplyBool,numel) << endl;
        }
        string FilenameSaveMatch;
        string FilenameSave;
        string FilenameSaveTxt;
        if (segment_analysis->filename_saveMatch==NULL){
            cout << "Warning! No name where to save matched image - first image name used" << endl;
            string FilenamePA=nifti_makebasename(segment_analysis->filename_FloatImages[0]);
            int Index=FilenamePA.find_last_of('/');
            string FilenamePA_b=FilenamePA.substr(0,Index+1);
            string FilenamePA_e = FilenamePA.substr(Index+1,FilenamePA.length());
            FilenameSaveMatch = FilenamePA_b + "Match_" + FilenamePA_e + ".nii.gz";
             FilenameSave=nifti_makebasename(FilenameSaveMatch.c_str());
            FilenameSaveTxt = FilenamePA_b + "Match_" + FilenamePA_e + ".txt";
                    
                    
        }
        else{
         FilenameSaveMatch=segment_analysis->filename_saveMatch;
         FilenameSave=nifti_makebasename(segment_analysis->filename_saveMatch);
         FilenameSaveTxt=FilenameSave+".txt";
    }
        vector<float > CoeffMatching=CompleteHistogramMatching(segment_analysis, MaskMatchBool);
        ofstream TxtFileResultMatching(FilenameSaveTxt.c_str());
        PrintingCoeffsHistMatch(CoeffMatching, segment_analysis, TxtFileResultMatching);
        

        
        int numbImages=segment_analysis->filename_FloatImages.size();
        int numbModTot=0;
        int sizeFit=segment_analysis->orderMatch+1;
        for(int i=0;i<numbImages;i++){
            cout << "Treating image for matching "<<i<<endl;
            nifti_image * ImageToModify=ReadFromFilename(segment_analysis->filename_FloatImages[i]);
            int numbmodal=ImageToModify->nu*ImageToModify->nt;
            vector<float> CoeffsToUse;
            stringstream nI;
            nI<< i;
            string numbI=nI.str();
            string FilenameSaveFin=FilenameSave+"_"+numbI+".nii.gz";
            for (int m=0; m<numbmodal; m++) {
                for (int o=0; o<segment_analysis->orderMatch+1; o++) {
                    CoeffsToUse.push_back(CoeffMatching[(m+numbModTot)*sizeFit+o]);
                }
            }
            numbModTot+=numbmodal;
            nifti_image * NewImageTest=ApplyingPolyfitMatching(CoeffsToUse, ImageToModify,MaskToApplyBool);
            nifti_set_filenames(NewImageTest, FilenameSaveFin.c_str(), 0, 0);
            cout<<"Images "<<i<<" saved"<<endl;
            nifti_image_write(NewImageTest);
            nifti_image_free(NewImageTest);
            nifti_image_free(ImageToModify);
            NewImageTest=NULL;
            ImageToModify=NULL;
        }
        cout << "Written images matched" << endl;
        if(MaskMatch!=NULL){
            nifti_image_free(MaskMatch);
            MaskMatch=NULL;
        }
        if(MaskToApply!=NULL && segment_analysis->flag_mask){
            nifti_image_free(MaskToApply);
            MaskToApply=NULL;
        }
        return EXIT_SUCCESS;
    }
    

    if (segment_analysis->flag_refLesConnect) {
        ConnectedGTImage=ReadFromFilename(segment_analysis->filename_RefConnect);
    }
    
    
    if (segment_analysis->flag_inVentricleSeg) {
        VentricleSeg=ReadFromFilename(segment_analysis->filename_inVentricleSeg);
        cout << "Getting the ventricle segmentation "<<endl;
    }

    if(segment_analysis->flag_WMMaps && segment_analysis->flag_inDataComp && segment_analysis->flag_maskMatch){
        cout << "Creating Mahal Map" << endl;
        nifti_image * NiiMask = ReadFromFilename(segment_analysis->filename_mask);

        nifti_image * NiiBasis = ReadFromFilename(segment_analysis->filename_maskMatch);
        nifti_image * DataNii = ReadFromFilename(segment_analysis->filename_inDataComp);
        nifti_image * NiiMaskFloat = CreateNiiFromArray(static_cast<bool*>(NiiMask->data),NiiBasis,NiiMask->nvox);
        nifti_image * MahalWMMapsMask= MahalDistMaps(NiiBasis, NiiMaskFloat, DataNii );
        string FilenamePA=nifti_makebasename(segment_analysis->filename_maskMatch);
        string FilenamePA2=nifti_makebasename(segment_analysis->filename_inDataComp);
        int Index2=FilenamePA2.find_last_of('/');
        string FilenamePA_b2=FilenamePA2.substr(0,Index2+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b2=segment_analysis->name_outputDir;
        }
        string FilenamePA_e2=FilenamePA2.substr(Index2+1,FilenamePA2.length());
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameMahal=FilenamePA_b+"Mahal_"+FilenamePA_e+"_"+FilenamePA_e2+".nii.gz";
        nifti_set_filenames(MahalWMMapsMask,FilenameMahal.c_str(),0,0);
        nifti_image_write(MahalWMMapsMask);
        nifti_image_free(MahalWMMapsMask);
        nifti_image_free(NiiMask);
        nifti_image_free(NiiBasis);
        nifti_image_free(DataNii);
        nifti_image_free(NiiMaskFloat);
        return EXIT_SUCCESS;
    }
    
    //    To get directly the Mahal maps for the WM and the corresponding WM Inlier segmentation
    if (TreeToAnalyse!=NULL && segment_analysis->flag_WMMaps && segment_analysis->flag_inLesCorr) {
        

        //        nifti_image * WMSeg=HardSegmentationThresholdFromNormResp(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM), TreeToAnalyse, 0.5);
        string FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameMahal=FilenamePA_b+"MahalWM_"+FilenamePA_e+".nii.gz";
        string FilenameMahalMask=FilenamePA_b+"MahalWMMask_"+FilenamePA_e+".nii.gz";
        string FilenameWMInliers=FilenamePA_b+"WMInliers_"+FilenamePA_e+".nii.gz";
        nifti_image * LesionCorr=ReadFromFilename(segment_analysis->filename_inLesCorr);
        
        int Dim[3];
        vector<int> DimVector;
        int Shift[3];
        float PixDim[3];
        float SumPixDim=0;
        float SumDim=0;
        for (int d=0; d<3; d++) {
            Dim[d]=LesionCorr->dim[d+1];
            DimVector.push_back(Dim[d]);
            SumDim+=Dim[d];
            PixDim[d]=LesionCorr->pixdim[d+1];
            SumPixDim+=PixDim[d];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Shift[1]*Dim[1];
        
        SEG_PARAMETERS * segment_param=new SEG_PARAMETERS();
        segment_param->VLkappa=3;
        float * WMTyp=TreeToAnalyse->MakeTypicalityFloatClass(1, segment_param);
        float * WMTypLong=CreateLong(WMTyp, TreeToAnalyse->GetL2S(), numel);
        delete [] WMTyp;
        nifti_image * WMNiiTyp=CreateNiiFromArray(WMTypLong, LesionCorr, numel);
        delete [] WMTypLong;
        nifti_image * TypAtlas4DNii=TreeToAnalyse->CreateTypicalityAtlas4D(segment_param);
        nifti_image * TypAtlasNii=TreeToAnalyse->CreateTypicalityAtlas(segment_param);
        
        nifti_set_filenames(WMNiiTyp, (g_scratch_dir + "/TestWMNiiTyp.nii.gz").c_str(), 0, 0);
        nifti_set_filenames(TypAtlas4DNii, (g_scratch_dir + "/TestTypAtlas4D.nii.gz").c_str(), 0, 0);
        nifti_set_filenames(TypAtlasNii, (g_scratch_dir + "/TestTypAtlas.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(WMNiiTyp);
        //        nifti_image_write(TypAtlas4DNii);
        //        nifti_image_write(TypAtlasNii);
        nifti_image_free(WMNiiTyp);
        nifti_image_free(TypAtlasNii);
        nifti_image_free(TypAtlas4DNii);
        
        float * WMTempInliers=CreateLong(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        float * WMTempOutliers=CreateLong(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(), TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
        
        float * LesionCorrData=static_cast<float*>(LesionCorr->data);
        for (int i=0;i<numel;i++){
            if (LesionCorrData[i]>0){
                WMTempInliers[i]=WMTempInliers[i]-LesionCorrData[i]>0?WMTempInliers[i]-LesionCorrData[i]>0:0;
            }
        }
        vector<float *> AdditionWM;
        AdditionWM.push_back(WMTempInliers);
        AdditionWM.push_back(WMTempOutliers);
        AdditionWM.push_back(LesionCorrData);
        float * AddedWM=AddArray(AdditionWM, TreeToAnalyse->GetNumberElements());
        nifti_image * AddedWMNii=CreateNiiFromArray(AddedWM, LesionCorr, TreeToAnalyse->GetNumberElements());
        nifti_image * WMInliersNii=CreateNiiFromArray(WMTempInliers, LesionCorr, TreeToAnalyse->GetNumberElements());
        float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
        nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
        nifti_set_filenames(DataNii, (g_scratch_dir + "/TestDataAtlas.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(DataNii);
        if (segment_analysis->flag_maskMatch){
            nifti_image * MaskMatch=ReadFromFilename(segment_analysis->filename_maskMatch);
            Floatisation(MaskMatch);
            nifti_image * MahalWMMapsMask= MahalDistMaps(WMInliersNii, MaskMatch, DataNii );
            nifti_set_filenames(MahalWMMapsMask, FilenameMahalMask.c_str(), 0, 0);
            nifti_image_write(MahalWMMapsMask);
            nifti_image_free(MahalWMMapsMask);
            //            nifti_image_free(MaskFloat);
            nifti_image_free(MaskMatch);
            
        }
        //        nifti_image * MahalWMMaps= MahalDistMaps(WMInliersNii, AddedWMNii, TreeToAnalyse->GetDataImage() );
        nifti_image * MahalWMMaps= MahalDistMaps(WMInliersNii, AddedWMNii, DataNii );
        nifti_image_free(DataNii);
        delete [] DataCorr;
        //        nifti_image * MahalWMMaps= MahalDistMaps(WMInliersNii, AddedWMNii, TreeToAnalyse->GetDataImage() );
        nifti_set_filenames(MahalWMMaps, FilenameMahal.c_str(), 0, 0);
        nifti_image_write(MahalWMMaps);
        if (segment_analysis->flag_Vesselness){
            int * L2S=TreeToAnalyse->GetL2S();
            int * S2L=TreeToAnalyse->GetS2L();
            numel=TreeToAnalyse->GetNumberElements();
            int numbscales=3;
            float scales[3];
            scales[0]=0.25;
            scales[1]=0.7;
            scales[2]=1;
//            float * MultiVesselness=new float[(numbmodal+1)*numel];
            
            for(int m=0;m<=numbmodal;m++){
                float * MultiVesselness=new float[(numbscales)*numel];
                
                string Code;
                if (m==0) {
                    Code="Tot";
                }
                else{
                    Code=ModalitiesCode[segment_analysis->vecModalities[m-1]];
                }
                string FilenameVesselness=FilenamePA_b+"Vesselness_"+Code+"_"+FilenamePA_e+".nii.gz";
                for (int s=0;s<numbscales;s++){
                    int polarity=1;
                    if( m>0&&segment_analysis->vecModalities[m-1]==2){
                        polarity=-1;
                    }
                    float* WMMapsT1=ExtractTPFromNii<float>(MahalWMMaps, m);
                    cout << "Max WMMaps T1 " << GetMax(WMMapsT1, TreeToAnalyse->GetNumberElements())<<endl;
                    float * BlurredWMMapsT1=GaussianBlurring(WMMapsT1, scales[s], DimVector,0);
                    //                nifti_image * BlurredNii=CreateNiiFromArray(BlurredWMMapsT1, LesionCorr, numel);
                    //                nifti_set_filenames(BlurredNii, (g_scratch_dir + "/TestBluriness.nii.gz").c_str(), 0, 0);
                    //                nifti_image_write(BlurredNii);
                    cout << "Max WMMaps T1 " << GetMax(BlurredWMMapsT1, TreeToAnalyse->GetNumberElements())<<endl;
                    float * SBWMMapsT1=CreateShort(BlurredWMMapsT1, S2L, TreeToAnalyse->GetNumberMaskedElements());
                    cout << "Max WMMaps T1 " << GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements())<<endl;
                    //MultiplyFloatArrayBy(SBWMMapsT1, GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements()), TreeToAnalyse->GetNumberMaskedElements());
                    cout << "Max WMMaps T1 " << GetMax(SBWMMapsT1, TreeToAnalyse->GetNumberMaskedElements())<<endl;
                    float * VesselnessT1=Vesselness(SBWMMapsT1, Dim, Shift, PixDim, L2S, S2L, TreeToAnalyse->GetNumberMaskedElements(), polarity);
                    float * VesselLong=CreateLong(VesselnessT1, L2S, TreeToAnalyse->GetNumberElements());
                    for (int i=0; i<numel; i++) {
                        MultiVesselness[i+s*numel]=VesselLong[i];
                    }
                    delete [] VesselnessT1;
                    delete [] VesselLong;
                    delete [] BlurredWMMapsT1;
                    delete [] WMMapsT1;
                    delete [] SBWMMapsT1;

                }

                nifti_image * VesselnessNii=CreateNiiFromArray(MultiVesselness, LesionCorr, (numbmodal+1)*TreeToAnalyse->GetNumberElements());
                nifti_set_filenames(VesselnessNii, FilenameVesselness.c_str(), 0, 0);
                nifti_image_write(VesselnessNii);
                nifti_image_free(VesselnessNii);
            }
        }
        TreeToAnalyse->SaveTmpResult(WMTempInliers, FilenameWMInliers.c_str());
        delete [] WMTempInliers;
        delete [] WMTempOutliers;
        delete [] AddedWM;
        nifti_image_free(AddedWMNii);
        nifti_image_free(MahalWMMaps);
        nifti_image_free(LesionCorr);
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
        
    }
    
    
    
    
    
    
    //    Initialisation of WMSeg / DataComp for Mahal Distance analysis either for GT only or for ref
    
    if ((segment_analysis->flag_TextFile || segment_analysis->flag_inDataComp) && (segment_analysis->flag_inWMSeg || segment_analysis->flag_TextFile) ) { // Case when we also perform the intensity study of the lesions
        cout<<"Initialisation of WMSeg to be done "<<endl;
        //            Get WM True Segmentation
        if (segment_analysis->flag_inWMSeg) {
            WMSeg=ReadFromFilename(segment_analysis->filename_inWMSeg);
        }
        else{ // If we do not have the true segmentation we only consider the WM inlier in the following
            cout<<"considering inlier for study "<<endl;;
            vector<TreeEM *> GeneralClassesVector=TreeToAnalyse->GetGeneralClassesVector();
            float * ShortWM=GeneralClassesVector[segment_analysis->IndexWM]->GetNormResp();
            cout<<"ShortWM obtained "<<endl;
            float * LongWM=CreateLong(ShortWM, TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
            cout<<"LongWM created "<<endl;
            nifti_image * BasicImage=nifti_copy_nim_info(TreeToAnalyse->GetDataImage());
            BasicImage->dim[0]=3;
            BasicImage->dim[4]=1;
            BasicImage->dim[5]=1;
            BasicImage->data=(void *) calloc(BasicImage->nvox, sizeof(float));
            nifti_update_dims_from_array(BasicImage);
            WMSeg=CreateNiiFromArray(LongWM, BasicImage,TreeToAnalyse->GetNumberElements());
            cout<<"WMSeg niftied "<<endl;
            delete [] LongWM;
            LongWM=NULL;
            nifti_image_free(BasicImage);
            BasicImage=NULL;
        }
        cout<<"Tissue rebuilt"<<endl;
        //            Get the Data we have to compare to (either as part of tree or if filename provided as image under filename. Note that consequently the number of modalities might not be the same in one case and another...)
        if (segment_analysis->flag_inDataComp) {
            DataCompImage=ReadFromFilename(segment_analysis->filename_inDataComp);
        }
        else{ // Extract it from TreeTmp
            DataCompImage=CopyFloatNii(TreeToAnalyse->GetDataImage());
            float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
            DataCompImage=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
            delete[]DataCorr;
        }
        cout<<"Data reread"<<endl;
        ParamWM=GetParamFromImages(WMSeg,DataCompImage,0);
        numbmodal=DataCompImage->nu*DataCompImage->nt;
        InvertedCovarianceWM=TreeToAnalyse->InvertMatrix(&ParamWM[numbmodal], numbmodal);
        numel=WMSeg->nx*WMSeg->ny*WMSeg->nz;
        numbmodal=DataCompImage->nt*DataCompImage->nu;
        cout<<"ParamVectorCalculated "<<endl;
    }
    cout<<"and finishing initialisation "<<endl;
    
    if (!segment_analysis->flag_GivenMahalRule) {
        cout<<"Rule mahal to build "<<endl;
        cout<<"numbmodal is "<<numbmodal<<endl;
        cout<<" and Modalities size is "<<Modalities.size()<<endl;
        for (int m=0; m<numbmodal; m++) {
            if (Modalities[m]==2 || Modalities[m]==3 || Modalities[m]==4) {
                segment_analysis->vecRuleMahal.push_back(10*m+1);
            }
        }
    }
    
    cout<<"Initialisation WM done"<<endl;
    
    if (segment_analysis->flag_refLesConnect) {
        ConnectedGTImage=ReadFromFilename(segment_analysis->filename_RefConnect);
        cout<<"Connected image read at "<<ConnectedGTImage<<endl;
        OrderedLabelsGT=TranscribeArray<float, int>(static_cast<float*>(ConnectedGTImage->data), ConnectedGTImage->nvox);
        cout<<"Image transcribed"<<endl;
    }
    
    if (segment_analysis->flag_GT) {// Take care of analysis or at least preparation of GT if need be.
        if (segment_analysis->flag_refLes) {
            nifti_image * GTImage=ReadFromFilename(segment_analysis->filename_Ref);
            if (GTImage!=NULL) { // If we really have a GT image to treat;
                int numelGT=GTImage->nx*GTImage->ny*GTImage->nz;
                //                First hard segmentation compare to threshold 0.5;
                nifti_image * GTSegHard=HardSegmentationThreshold(GTImage, 0.5);
                numel=GTImage->nvox;
                //                Create Connected Component from GT
                int * ConnectedGT=ComponentLabelingNotRefined(GTSegHard, segment_analysis->Neigh);
                //                int CountNonZeroInit=CountNonZero(static_cast<float *>(GTSegHard->data),numel);
                //                int CountNonZeroConnected=CountNonZero(ConnectedGT, numel);
                float VolumeVox=GTSegHard->pixdim[1]*GTSegHard->pixdim[2]*GTSegHard->pixdim[3];
                OrderedLabelsGT=OrderedVolumeLabel(ConnectedGT, 1, numelGT,VolumeVox);
                //                int CountAfterOrder=CountNonZero(OrderedLabelsGT, numel);
                vector<LesionSimple*> VectorLesionGT=GetVectorLesionSimple(OrderedLabelsGT, GTSegHard);
                nifti_image * ConnectedGTImage=CreateNiiFromArray(OrderedLabelsGT, GTImage,numel);
                
                FilenamePA=nifti_makebasename(segment_analysis->filename_Ref);
                int Index=FilenamePA.find_last_of('/');
                FilenamePA_b=FilenamePA.substr(0,Index+1);
                if(segment_analysis->flag_outputDir){
                    FilenamePA_b=segment_analysis->name_outputDir;
                }
                FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
                
                string OptionText;

                if (segment_analysis->flag_inOptionText) {
                    OptionText=segment_analysis->inOptionText;
                }
                if(segment_analysis->flag_VentrSeg==1){
                    OptionText=OptionText+"Ventr1";
                }
                if (segment_analysis->flag_juxtaCorrection) {
                    OptionText=OptionText+"JC1";
                }
                if (segment_analysis->flag_SPCorrection) {
                    OptionText=OptionText+"SP1";
                }
                if (segment_analysis->flag_oldLesion) {
                    OptionText=OptionText+"OL1";
                }
                if (segment_analysis->flag_GTSave) {
                    
                    //                    Creation of according filenames to GT filename

                    string FilenameGTSaved=FilenamePA_b+"Connected_"+FilenamePA_e+"_"+OptionText+".nii.gz";
                    string FilenameGTTxt=FilenamePA_b+"Connected_"+FilenamePA_e+"_"+OptionText+".txt";
                    nifti_set_filenames(ConnectedGTImage, FilenameGTSaved.c_str(), 0, 0);
                    nifti_image_write(ConnectedGTImage);
                    ofstream TxtFile(FilenameGTTxt.c_str());
                    int numbGTLes=VectorLesionGT.size();
                    for (int l=0; l<numbGTLes; l++) {
                        TxtFile << "Lesion "<<l+1<<endl;
                        PrintLesionSimple(VectorLesionGT[l],TxtFile);
                        TxtFile << endl;
                    }
                }

                //            Case where we also perform the intensity analysis for the ground truth
                if (WMSeg!=NULL && DataCompImage!=NULL) { // Case when we also perform the intensity study of the lesions


                    string FilenameGTWMTxt=FilenamePA_b+"Intensity_"+FilenamePA_e+"_"+OptionText+".txt";
                    int numbmodalGT=DataCompImage->nu*DataCompImage->nt;
                    int sizeOfParam=numbmodalGT*(numbmodalGT+1);
                    nifti_image * MahalGT=MahalDistMaps(WMSeg, ConnectedGTImage, DataCompImage);
                    nifti_image * MahalWM=MahalDistMaps(WMSeg, WMSeg, DataCompImage);

                    if (segment_analysis->flag_GTSave) {
                        string MahalWMFilename=FilenamePA_b+"MahalWM_"+FilenamePA_e+".nii.gz";
                        string MahalGTFilename=FilenamePA_b+"MahalGT_"+FilenamePA_e+".nii.gz";
                        nifti_set_filenames(MahalWM, MahalWMFilename.c_str(), 0, 0);
                        nifti_set_filenames(MahalGT, MahalGTFilename.c_str(), 0, 0);
                        nifti_image_write(MahalWM);
                        nifti_image_write(MahalGT);
                    }
                    vector<LesionIntensity*> VectorLesionGTIntensity=GetVectorLesionIntensity(OrderedLabelsGT, ParamWM,InvertedCovarianceWM,DataCompImage,MahalGT);

                    ofstream TxtFileWM(FilenameGTWMTxt.c_str());
                    int numbGTLes=VectorLesionGTIntensity.size();

                    TxtFileWM<< "Label Volume ";
                    for (int m=0; m<numbmodalGT; m++) {
                        TxtFileWM<< "Mean"<< m<<" ";
                    }
                    for (int m=0; m<numbmodalGT; m++) {
                        TxtFileWM<< "Ratio"<<m<<" ";
                    }
                    TxtFileWM<<"GlobalMahal ";
                    for (int m=0; m<numbmodalGT; m++) {
                        TxtFileWM<<"Mahal"<<m<<" ";
                    }
                    TxtFileWM<<"MinMahal ";
                    TxtFileWM<<"MaxMahal ";
                    TxtFileWM<<"QuantInfMahal ";
                    TxtFileWM<<"MedMahal ";
                    TxtFileWM<<"QuantSupMahal ";
                    TxtFileWM<<endl;
                    TxtFileWM<< "WM ";
                    for (int m=0; m<sizeOfParam; m++) {
                        TxtFileWM<<ParamWM[m]<<" ";
                    }
                    TxtFileWM<<endl;

                    //                Print globally for lesion volume before doing it per connected lesion.
                    TxtFileWM<<"Global ";
                    //                int GlobalVolume=CountNonZero(OrderedLabelsGT, numelGT);
                    //                float * LesionParam=GetParamFromImages(GTSegHard, DataCompImage, 0);
                    bool * BoolLesionGlobal=TranscribeArray<int, bool>(OrderedLabelsGT, numelGT);
                    int * GlobalHardGT=TranscribeArray<bool, int>(BoolLesionGlobal, numelGT);
                    if(CountNonZero(GlobalHardGT, numel)>0){
                        LesionIntensity * LesionGlobal=BuildLesionIntensityFromConnectedLabel(GlobalHardGT, 1, ParamWM, InvertedCovarianceWM, WMSeg,MahalGT);
                        PrintLesionIntensity(LesionGlobal, TxtFileWM);
                    }

                    for (int l=0; l<numbGTLes; l++) {
                        TxtFileWM <<l+1<<" ";
                        PrintLesionIntensity(VectorLesionGTIntensity[l],TxtFileWM);
                        //                    TxtFileWM << endl;
                    }
                    if (MahalGT!=NULL) {
                        nifti_image_free(MahalGT);
                        MahalGT=NULL;
                    }
                    if (MahalWM!=NULL) {
                        nifti_image_free(MahalWM);
                        MahalWM=NULL;
                    }
                    if (BoolLesionGlobal!=NULL) {
                        delete [] BoolLesionGlobal;
                        BoolLesionGlobal=NULL;
                    }
                    if (GlobalHardGT!=NULL) {
                        delete [] GlobalHardGT;
                        GlobalHardGT=NULL;
                    }

                }
            }
            cout<<"Finished GT Analysis"<<endl;
        }
        else{
            cout<<"No possible analysis on GT if no reference image to use..."<<endl;
        }
    }

    //    cout<<segment_analysis->filename_TextFile<<endl;

    if(TypePackage==0 && !segment_analysis->flag_LapAnalysis && !segment_analysis->flag_inToAnalyse && !segment_analysis->flag_RuleCorr && !segment_analysis->flag_PrintConnect){
        cout<<"No construction of lesion to be done"<<endl;
        if(ImageLesionCorr!=NULL){
            nifti_image_free(ImageLesionCorr);
            ImageLesionCorr=NULL;
        }
        if (OrderedLabelsGT!=NULL) {
            delete [] OrderedLabelsGT;
            OrderedLabelsGT=NULL;
        }
        if (OrderedLabels!=NULL) {
            delete [] OrderedLabels;
            OrderedLabels=NULL;
        }
        if (OrderedLabelsCorr!=NULL) {
            delete [] OrderedLabelsCorr;
            OrderedLabelsCorr=NULL;
        }
        if (WMSeg!=NULL) {
            nifti_image_free(WMSeg);
            WMSeg=NULL;
        }
        if(ParamWM!=NULL){
            delete [] ParamWM;
            ParamWM=NULL;
        }
        if(InvertedCovarianceWM!=NULL){
            delete [] InvertedCovarianceWM;
            InvertedCovarianceWM=NULL;
        }
        if(DataCompImage!=NULL){
            nifti_image_free(DataCompImage);
            DataCompImage=NULL;
        }
        if (ConnectedGTImage!=NULL) {
            nifti_image_free(ConnectedGTImage);
            ConnectedGTImage=NULL;
        }
        if (TreeToAnalyse!=NULL) {
            delete TreeToAnalyse;
            TreeToAnalyse=NULL;
        }
        return EXIT_SUCCESS;
    }
    
    if(segment_analysis->flag_in==0 && segment_analysis->flag_TextFile==0 && segment_analysis->flag_test==0 && segment_analysis->flag_inConnect==0 && (segment_analysis->flag_LapAnalysis==0&&segment_analysis->flag_inSum==0) && segment_analysis->flag_RuleCorr){
        fprintf(stderr,"Err:\tThe input image name has to be defined.\n");
        Usage(argv[0]);
        if(ImageLesionCorr!=NULL){
            nifti_image_free(ImageLesionCorr);
            ImageLesionCorr=NULL;
        }
        if (OrderedLabelsGT!=NULL) {
            delete [] OrderedLabelsGT;
            OrderedLabelsGT=NULL;
        }
        if (OrderedLabels!=NULL) {
            delete [] OrderedLabels;
            OrderedLabels=NULL;
        }
        if (OrderedLabelsCorr!=NULL) {
            delete [] OrderedLabelsCorr;
            OrderedLabelsCorr=NULL;
        }
        if (WMSeg!=NULL) {
            nifti_image_free(WMSeg);
            WMSeg=NULL;
        }
        if(DataCompImage!=NULL){
            nifti_image_free(DataCompImage);
            DataCompImage=NULL;
        }
        if (ConnectedGTImage!=NULL) {
            nifti_image_free(ConnectedGTImage);
            ConnectedGTImage=NULL;
        }
        if (TreeToAnalyse!=NULL) {
            delete TreeToAnalyse;
            TreeToAnalyse=NULL;
        }
        return EXIT_SUCCESS;
    }
    
    if (TreeToAnalyse!=NULL && segment_analysis->flag_inToAnalyse ) {
        cout << "Analysis of connected components of external segmentation"<<endl;

        nifti_image * ExtToAnalyse=ReadFromFilename(segment_analysis->filename_inToAnalyse);

        numel=ExtToAnalyse->nvox;

        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);


        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        float * ExtToAnalyseData=NULL;
        if(segment_analysis->flag_connect){
            int *Labels = ComponentLabeling(ExtToAnalyse, segment_analysis->Neigh);
            ExtToAnalyseData=TranscribeArray<int,float>(Labels,numel);
        }
        else{
            ExtToAnalyseData=static_cast<float*>(ExtToAnalyse->data);
        }
        cout <<GetMax(ExtToAnalyseData, numel)<< " "<< GetMin(ExtToAnalyseData, numel)<<endl;
        vector<int> ValuesExt=GetListValues<float,int>(ExtToAnalyseData,numel);
        //        nifti_image * CreateNiiTest=CreateNiiFromArray(ExtToAnalyseData, ExtToAnalyse, numel);

        
        //        string FilenameSaveNiiTest=FilenamePA_b+"Test_"+FilenamePA_e+".nii.gz";
        //        nifti_set_filenames(CreateNiiTest, FilenameSaveNiiTest.c_str(), 0, 0);
        //        nifti_image_write(CreateNiiTest);
        //        bool * EqualTest=EqualArray<float, bool>(ExtToAnalyseData, GetMax(ExtToAnalyseData, numel), numel);
        //        cout << CountNonZero(EqualTest, numel)<< " "<< GetMax(ExtToAnalyseData, numel)<< endl;
        vector<vector<Outlier*> > ExtToAnalyseVector=GetVectorOutliersMultiType(ExtToAnalyseData, TreeToAnalyse, segment_analysis);
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_inToAnalyse);
        
        
        string FilenameSaveDescription=FilenamePA_b+"DescriptionToAnalyse_"+FilenamePA_e+".txt";
        string FilenameSaveDescription2=FilenamePA_b+"DescriptionToAnalyse2_"+FilenamePA_e+".txt";
        ofstream TxtFileLesion(FilenameSaveDescription.c_str());
        ofstream TxtFileLesion2(FilenameSaveDescription2.c_str());
        int * GravImage=NULL;
        //        int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
        float PixDim[3];
        if (ExtToAnalyseVector[0].size()>0) {
            GravImage=new int[3];
            
            for (int d=0;d<3; d++) {
                PixDim[d]=ExtToAnalyse->pixdim[d+1];
                GravImage[d]=(ExtToAnalyseVector[0][0]->CentreGravity[d]-ExtToAnalyseVector[0][0]->VectorDiffGrav[d]/PixDim[d]);
            }
        }
        int numbvalues=ValuesExt.size();
//        numel=ExtToAnalyse->nvox;
        for (int v=0; v<numbvalues; v++) {
            int numbFinalLesions=ExtToAnalyseVector[v].size();
            TxtFileLesion2<< "Type of lesion "<< ValuesExt[v]<<endl;

            for (int l=0; l<numbFinalLesions; l++) {
                cout <<"Label "<<l+1;
                TxtFileLesion2<<endl;
                TxtFileLesion2<<"Label "<<l+1<<endl;
                TxtFileLesion<<l+1+ValuesExt[v]*1000<< " ";
                PrintOutlierCharacteristicsForFile(ExtToAnalyseVector[v][l],TxtFileLesion,numbmodal,GravImage,PixDim,segment_analysis);
                PrintOutlierCharacteristics(ExtToAnalyseVector[v][l],TxtFileLesion2,numbmodal,GravImage,PixDim,segment_analysis);
                TxtFileLesion<<endl;
                cout<<"...printed"<<endl;
            }
        }
        delete [] GravImage;
        nifti_image_free(ExtToAnalyse);
        for (int v=0; v<numbvalues; v++) {
            int numbFinalLesions=ExtToAnalyseVector[v].size();
            for (int l=0; l<numbFinalLesions; l++) {
                if(ExtToAnalyseVector[v][l]!=NULL){
                    delete ExtToAnalyseVector[v][l];
                    ExtToAnalyseVector[v][l]=NULL;
                }
            }
        }
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
    }
    
    
    if(segment_analysis->flag_PrintConnect && (TreeToAnalyse!=NULL || segment_analysis->flag_inDataComp) && segment_analysis->flag_inLes){
        nifti_image * SegToAnalyse=ReadFromFilename(segment_analysis->filename_InLes);
        numel=SegToAnalyse->nvox;
        float * SegData=static_cast<float*>(SegToAnalyse->data);
        int * CompLabel_I=ComponentLabeling(SegToAnalyse, 26);
        float VolumeVox=SegToAnalyse->pixdim[1]*SegToAnalyse->pixdim[2]*SegToAnalyse->pixdim[3];
        int * CompLabel=OrderedVolumeLabel(CompLabel_I, 1, numel,VolumeVox);
        nifti_image * DataImage=NULL;
        int numbmodal;
        if(segment_analysis->flag_inDataComp){
            DataImage=ReadFromFilename(segment_analysis->filename_inDataComp);
            numbmodal=DataImage->nt*DataImage->nu;
        }
        else{
            DataImage=TreeToAnalyse->GetDataImage();
            numbmodal=TreeToAnalyse->GetNumberModalities();
            float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
            DataImage=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
            delete [] DataCorr;
        }
        
        int numbLabel=GetMaxLabel(CompLabel, numel);
        int Dim[3];
        int Shift[3];
        float Pixdim[3];
        for (int d=0; d<3; d++) {
            Dim[d]=SegToAnalyse->dim[d+1];
            Pixdim[d]=SegToAnalyse->pixdim[d+1];
            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
        }
        for (int l=0; l<numbLabel; l++) {
            bool * LesionBool=CreateLesionBool(CompLabel, l+1, numel);
            int MeanLab=floorf(10*GetMeanData(SegData, LesionBool, numel));
            int CentreGrav=GetCenterGravity(LesionBool, Dim);
            
            int * BoundingBox=FindBoundingBox(LesionBool, Dim, Shift);
            int MaxBoundingBox=0;
            for (int d=0; d<3; d++) {
                if (abs(BoundingBox[2*d]-BoundingBox[2*d+1])>MaxBoundingBox) {
                    MaxBoundingBox=abs(BoundingBox[2*d]-BoundingBox[2*d+1]);
                }
            }
            int NBBToPrint[6];
            int * CentreGravCoord=CorrespondingCoordinates(CentreGrav, Dim, Shift);
            CorrectedCentreGrav(LesionBool,CentreGravCoord,Dim,Shift,Pixdim);
            int Box=MaxBoundingBox;
            Box=Box<15?15:MaxBoundingBox;
            int LengthBox=2*Box+1;
            for (int d=0; d<3; d++) {
                NBBToPrint[2*d]=CentreGravCoord[d]-Box;
                NBBToPrint[2*d+1]=CentreGravCoord[d]+Box+1;
            }
            vector<float*> ValuesToPrintVector;
            for (int m=0; m<numbmodal; m++) {
                for (int d=0; d<3; d++) {
                    float * PrintedPlane=PrintData(DataImage,m,NBBToPrint,d);
                    ValuesToPrintVector.push_back(PrintedPlane);
                }
            }
            float * TotalToDisplay=new float[(LengthBox+2)*(LengthBox+2)*numbmodal*3];
            int numbTot=(LengthBox+2)*(LengthBox+2)*numbmodal*3;
            for (int i=0; i<numbTot; i++) {
                TotalToDisplay[i]=0;
            }
            int i2; int j2;
            for (int m=0; m<numbmodal; m++) {
                for (int d=0; d<3; d++) {
                    for (int i=0; i<LengthBox; i++) {
                        for (int j=0; j<LengthBox; j++) {
                            if (d==0){
                                i2=LengthBox-1-i;
                                j2=LengthBox-1-j;
                            }
                            else if (d==1){
                                i2=i;
                                j2=LengthBox-1-j;
                            }
                            else{
                                i2=LengthBox-1-j;
                                j2=i;
                            }
                            TotalToDisplay[m*(LengthBox+2)+d*numbmodal*(LengthBox+2)*(LengthBox+2)+i2*(LengthBox+2)*numbmodal+j2]=ValuesToPrintVector[m*3+d][i*LengthBox+j];
                        }
                    }
                }
            }
            FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            stringstream lab;
            stringstream x;
            stringstream y;
            stringstream z;
            stringstream t;
            lab << l+1;
            x << CentreGravCoord[0];
            y << CentreGravCoord[1];
            z << CentreGravCoord[2];
            t<< MeanLab;
            string FilenameSaveTemp=FilenamePA_b+"Label"+lab.str()+"_T"+t.str()+"_CG_x"+x.str()+"_y"+y.str()+"_z"+z.str()+FilenamePA_e+".txt";
            ofstream TxtFile(FilenameSaveTemp.c_str());
            for (int i=0; i<numbmodal*(LengthBox+2); i++) {
                for (int j=0; j<3*(LengthBox+2); j++) {
                    TxtFile << TotalToDisplay[i*3*(LengthBox+2)+j]<< " ";
                }
                TxtFile<<endl;
            }
            delete [] LesionBool;
            delete [] CentreGravCoord;
            delete [] TotalToDisplay;
            delete [] BoundingBox;
            for (int m=0; m<numbmodal; m++) {
                for (int d=0; d<3; d++) {
                    delete [] ValuesToPrintVector[m*3+d];
                    ValuesToPrintVector[m*3+d]=NULL;
                }
            }
        }
        delete [] CompLabel;
        nifti_image_free(SegToAnalyse);
        if (segment_analysis->flag_TextFile){
            delete TreeToAnalyse;
        }
        else if (segment_analysis->flag_inDataComp){
            nifti_image_free(DataImage);
        }
        return EXIT_SUCCESS;
    }
    
    if (segment_analysis->flag_RuleCorr && TreeToAnalyse!=NULL && segment_analysis->flag_inLes ) {
        //        First build rule that will be used for correction assessment
        RuleCorr * Rule=BuildRuleCorrFromTextFile(segment_analysis);
        SegToAnalyse=ReadFromFilename(segment_analysis->filename_InLes);
        vector<nifti_image *> MahalVec;
        vector<float*> MahalDataVec;
        nifti_image * DataImage=TreeToAnalyse->GetDataImage();
        int numelmasked=TreeToAnalyse->GetNumberMaskedElements();
        int * L2S=TreeToAnalyse->GetL2S();
        nifti_image * PriorsDGM = ReadFromFilename(segment_analysis->filename_inPriorsDGM);
        nifti_image * PriorsWM = ReadFromFilename(segment_analysis->filename_inPriorsWM);
        int numel=TreeToAnalyse->GetNumberElements();
        vector<TreeEM *> WMSegVector;
        WMSegVector.push_back(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM));
        WMSegVector.push_back(TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM));
        float * ArtefactData=NULL;
        nifti_image * Artefact=NULL;
        if (segment_analysis->flag_inArtefact) {
            Artefact=ReadFromFilename(segment_analysis->filename_Artefact);
            ArtefactData=static_cast<float*>(Artefact->data);
        }
        
        float * DGMData= static_cast<float*>(PriorsDGM->data);
        float * WMData=static_cast<float*>(PriorsWM->data);
        
        float * AddWM=AddNormResp(WMSegVector);
        float * WMLong=CreateLong(AddWM, TreeToAnalyse->GetL2S(), numel);
        bool * ZonePotential=new bool[numel];
        for(int i=0;i<numel;i++){
            ZonePotential[i]=0;
            if(WMLong[i]>0.3){
                ZonePotential[i]=1;
            }
            if(DGMData[i]+WMData[i]>0.2){
                ZonePotential[i]=1;
            }
            if (ArtefactData!=NULL){
                if (ArtefactData[i]){
                    ZonePotential[i]=0;
                }
            }
        }
        nifti_image * ZonePotNii=CreateNiiFromArray(ZonePotential, PriorsDGM, numel);
        nifti_image_free(PriorsDGM);
        //        nifti_image_free(PriorsICSF);
        nifti_image_free(PriorsWM);
        
        
        for (int m=0; m<numbmodal; m++) {
            bool * ThresholdedTissueCompShort=ThresholdArray<float, bool>(TreeToAnalyse->GetNodeInlier()->GetChild(Rule->ClassesComparison[m])->GetNormResp(), 0.5, numelmasked);
            bool * TissueCompLong=CreateLong(ThresholdedTissueCompShort, L2S, numel);
            nifti_image * TissueCompNii=CreateNiiFromArray(TissueCompLong, SegToAnalyse, numel);
            nifti_image * MahalToPush=MahalDistMaps(TissueCompNii, ZonePotNii, DataImage);
            //nifti_set_filenames(MahalToPush, (g_scratch_dir + "/Mahal.nii.gz").c_str(), 0, 0);
            //nifti_image_write(MahalToPush);
            nifti_image_free(TissueCompNii);
            delete [] TissueCompLong;
            delete [] ThresholdedTissueCompShort;
            TissueCompLong=NULL;
            TissueCompNii=NULL;
            ThresholdedTissueCompShort=NULL;
            float * MahalDataToPush= ExtractTPFromNii<float>(MahalToPush, m+1);
            //            float * MahalData=static_cast<float*>(MahalToPush->data);
            //            float * MahalDataToPush=&MahalData[(m+1)*numel];
            MahalDataVec.push_back(MahalDataToPush);
            nifti_image_free(MahalToPush);
            //            MahalVec.push_back(MahalToPush);
        }
        nifti_image_free(ZonePotNii);
//        float ValueTest=MahalDataVec[1][5935478]-MahalDataVec[2][5935478];
        int * CorrWeight=CreateCorrWeight(TreeToAnalyse, Rule, SegToAnalyse, MahalDataVec, segment_analysis);
        nifti_image * CorrWeightNii=CreateNiiFromArray(CorrWeight, SegToAnalyse, numel);
        
        nifti_image * PotentialPVS= PVSExtraction(Rule,TreeToAnalyse,segment_analysis,MahalDataVec);
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameSaveCorrWeight=FilenamePA_b+"CorrWeight6_"+FilenamePA_e+".nii.gz";
        string FilenameSaveCorrImprove=FilenamePA_b+"CorrWeightImproved6_"+FilenamePA_e+".nii.gz";
        string FilenameSavePVS=FilenamePA_b+"PVSExtraction6_"+FilenamePA_e+".nii.gz";
        string FilenameSavePVSConnect=FilenamePA_b+"PVSConnected6_"+FilenamePA_e+".nii.gz";
        string FilenameSavePVSTxt=FilenamePA_b+"PVSCard6b_"+FilenamePA_e+".txt";
        string FilenameSavePVSType=FilenamePA_b+"PVSSubtypes26_"+FilenamePA_e+".nii.gz";
        string FilenamePVSProba=FilenamePA_b+"PVSProba_"+FilenamePA_e+".nii.gz";
        string FilenameLacunesProba=FilenamePA_b+"LacunesProba_"+FilenamePA_e+".nii.gz";
        
        int * ComponentLabelInit=ComponentLabeling(PotentialPVS, 6);
        float VolumeVox=PotentialPVS->pixdim[1]*PotentialPVS->pixdim[2]*PotentialPVS->pixdim[3];
        int * OrderedOutliers=OrderedVolumeLabel(ComponentLabelInit, 1, numel,VolumeVox);
        nifti_image * OrderedOutliersNii=CreateNiiFromArray(OrderedOutliers, PotentialPVS, numel);
        nifti_set_filenames(OrderedOutliersNii, FilenameSavePVSConnect.c_str(), 0, 0);
        nifti_image_write(OrderedOutliersNii);
        vector<Outlier *> VectorOutliers=GetVectorOutliers(OrderedOutliers, TreeToAnalyse, segment_analysis);
        ofstream TxtFileLesion(FilenameSavePVSTxt.c_str());
        int * GravImage=NULL;
        int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
        float PixDim[3];
        if (numbFinalLesions>0) {
            GravImage=new int[3];
            
            for (int d=0;d<3; d++) {
                PixDim[d]=PotentialPVS->pixdim[d+1];
                GravImage[d]=(VectorOutliers[0]->CentreGravity[d]-VectorOutliers[0]->VectorDiffGrav[d]/PixDim[d]);
            }
        }
        for (int l=0; l<numbFinalLesions; l++) {
            cout <<"Label "<<l+1;
            TxtFileLesion<<endl;
            TxtFileLesion<<"Label "<<l+1<<endl;
            PrintOutlierCharacteristics(VectorOutliers[l],TxtFileLesion,numbmodal,GravImage,PixDim,segment_analysis);
            cout<<"...printed"<<endl;
        }
        int * OutlierClassif=CopyArray(OrderedOutliers, numel);
        int * LesionSubClassif=CopyArray(OrderedOutliers, numel);
        for (int i=0; i<numel; i++) {
            if (OrderedOutliers[i]>0) {
                OutlierClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->OutlierClass;
                LesionSubClassif[i]=VectorOutliers[OrderedOutliers[i]-1]->LesionType;
            }
        }
        bool * LacunesClassif=ThresholdArray<int,bool>(LesionSubClassif,6,numel);
        bool * PVSClassif=UpperThresholdArray<int,bool>(LesionSubClassif,6,numel);
        nifti_image * LacunesNii=CreateNiiFromArray(LacunesClassif,PotentialPVS,numel);
        nifti_image* PVSNii=CreateNiiFromArray(PVSClassif,PotentialPVS,numel);
        float * ProbaLacunes=ProbaCorrWeight(TreeToAnalyse, Rule,  LacunesNii,  MahalDataVec,  segment_analysis);
        float * ProbaPVS=ProbaCorrWeight(TreeToAnalyse, Rule,  PVSNii,  MahalDataVec,  segment_analysis);
        SaveTmpResult(ProbaPVS,FilenamePVSProba.c_str(),PotentialPVS);
        SaveTmpResult(ProbaLacunes,FilenameLacunesProba.c_str(),PotentialPVS);
        nifti_image * OutlierClassifNii=CreateNiiFromArray(LesionSubClassif, PotentialPVS, numel);
        nifti_set_filenames(OutlierClassifNii, FilenameSavePVSType.c_str(), 0, 0);
        nifti_image_write(OutlierClassifNii);
        nifti_image_free(OutlierClassifNii);
        nifti_image_free(LacunesNii);
        nifti_image_free(PVSNii);
        delete [] LacunesClassif;
        delete [] PVSClassif;
        OutlierClassifNii=NULL;

        
        
        
        
        
        
        int * CorrImp=ImproveCorrWeight(TreeToAnalyse, Rule, SegToAnalyse, MahalDataVec, segment_analysis);
        nifti_image * CorrWeightImpNii=CreateNiiFromArray(CorrImp, SegToAnalyse, numel);
        nifti_set_filenames(CorrWeightNii, FilenameSaveCorrWeight.c_str(), 0, 0);
        nifti_image_write(CorrWeightNii);
        nifti_image_free(CorrWeightNii);
        nifti_set_filenames(CorrWeightImpNii, FilenameSaveCorrImprove.c_str(), 0, 0);
        nifti_image_write(CorrWeightImpNii);
        nifti_image_free(CorrWeightImpNii);
        nifti_set_filenames(PotentialPVS, FilenameSavePVS.c_str(), 0, 0);
        nifti_image_write(PotentialPVS);
        nifti_image_free(PotentialPVS);
        for (int m=0; m<numbmodal; m++) {
            //            nifti_image_free(MahalVec[m]);
            delete [] MahalDataVec[m];
        }
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
    }
    
    
    if (segment_analysis->flag_refLes && segment_analysis->flag_inLes && segment_analysis->GTAnalysisType==1 && TreeToAnalyse!=NULL){
        cout<<"Doing Ref Les comparison for other lesions"<<endl;
        if(SegToAnalyse==NULL){
            SegToAnalyse=ReadFromFilename(segment_analysis->filename_InLes);
        }
        nifti_image * SegToCompare=ReadFromFilename(segment_analysis->filename_Ref);
        int * ComponentSeg=ComponentLabeling(SegToAnalyse, 26);
        int * ComponentRef=ComponentLabeling(SegToCompare, 26);
        numel = SegToCompare->nvox;
        bool * BoolSeg=TranscribeArray<int,bool>(ComponentSeg, numel);
        bool * BoolRef=TranscribeArray<int, bool>(ComponentRef, numel);
        bool * FNConnect=new bool [numel];
        bool * FPConnect=new bool [numel];
        for(int i=0;i<numel;i++){
            FNConnect[i]=0;
            FPConnect[i]=0;
        }
        int numbLabelSeg=GetMaxLabel(ComponentSeg, numel);
        int numbLabelRef=GetMaxLabel(ComponentRef, numel);
        for(int l=0;l<numbLabelSeg;l++){
            bool * LesionBool=CreateLesionBool(ComponentSeg, l+1, numel);
            ANDOperationBool(BoolRef, LesionBool, LesionBool, numel);
            if(CountNonZero(LesionBool, numel)==0){ // Then the lesion is a FP
                for(int i=0;i<numel;i++){
                    if(ComponentSeg[i]==l+1){
                        FPConnect[i]=1;
                    }
                }
            }
            delete [] LesionBool;
        }
        for(int l=0;l<numbLabelRef;l++){
            bool * LesionBool=CreateLesionBool(ComponentRef, l+1, numel);
            ANDOperationBool(BoolSeg, LesionBool, LesionBool, numel);
            if(CountNonZero(LesionBool, numel)==0){ // Then the lesion is a FP
                for(int i=0;i<numel;i++){
                    if(ComponentRef[i]==l+1){
                        FNConnect[i]=1;
                    }
                }
            }
            delete [] LesionBool;
        }
        nifti_image * FNNii=CreateNiiFromArray(FNConnect,SegToCompare,numel);
        nifti_image * FPNii=CreateNiiFromArray(FPConnect, SegToCompare, numel);
        int * FPComponent=ComponentLabeling(FPNii, 26);
        int * FNComponent=ComponentLabeling(FNNii,26);
        float VolumeVox=FPNii->pixdim[1]*FPNii->pixdim[2]*FPNii->pixdim[3];
        int * FPComponentO=OrderedVolumeLabel(FPComponent, 1, numel,VolumeVox);
        int * FNComponentO=OrderedVolumeLabel(FNComponent, 1, numel,VolumeVox);
        delete [] FNComponent;
        delete [] FPComponent;
        delete [] BoolRef;
        delete [] BoolSeg;
        delete [] ComponentRef;
        delete [] ComponentSeg;
        vector<Outlier *> FPVector=GetVectorOutliers(FPComponentO, TreeToAnalyse, segment_analysis);
        vector<Outlier*> FNVector=GetVectorOutliers(FNComponentO, TreeToAnalyse, segment_analysis);
        
        
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        string FilenameSaveFN=FilenamePA_b+"FNAnalyse_"+FilenamePA_e+".txt";
        string FilenameSaveFP=FilenamePA_b+"FPAnalyse_"+FilenamePA_e+".txt";
        string FilenameSaveFNNii=FilenamePA_b+"FNAnalyse_"+FilenamePA_e+".nii.gz";
        string FilenameSaveFPNii=FilenamePA_b+"FPAnalyse_"+FilenamePA_e+".nii.gz";
        ofstream TxtFileLesionFN(FilenameSaveFN.c_str());
        ofstream TxtFileLesionFP(FilenameSaveFP.c_str());
        int * GravImage=NULL;
        //        int numbFinalLesions=GetMaxLabel(OrderedOutliers, numel);
        float PixDim[3];
        if (FNVector.size()>0) {
            GravImage=new int[3];
            
            for (int d=0;d<3; d++) {
                PixDim[d]=SegToCompare->pixdim[d+1];
                GravImage[d]=(FNVector[0]->CentreGravity[d]-FNVector[0]->VectorDiffGrav[d]/PixDim[d]);
            }
        }
        numel=SegToCompare->nvox;
        int numbLabelFP=GetMaxLabel(FPComponentO, numel);
        int numbLabelFN=GetMaxLabel(FNComponentO, numel);
        nifti_set_filenames(FPNii, FilenameSaveFPNii.c_str(), 0, 0);
        nifti_set_filenames(FNNii, FilenameSaveFNNii.c_str(), 0, 0);
        nifti_image_write(FPNii);
        nifti_image_write(FNNii);
        
        for (int l=0; l<numbLabelFP; l++) {
            TxtFileLesionFP<<l+1<< " ";
            PrintOutlierCharacteristicsForFile(FPVector[l],TxtFileLesionFP,numbmodal,GravImage,PixDim,segment_analysis);

            TxtFileLesionFP<<endl;
            cout<<"...printed"<<endl;
        }

        for (int l=0; l<numbLabelFN; l++) {
            TxtFileLesionFN<<l+1<< " ";
            PrintOutlierCharacteristicsForFile(FNVector[l],TxtFileLesionFN,numbmodal,GravImage,PixDim,segment_analysis);
            
            TxtFileLesionFN<<endl;
            cout<<"...printed"<<endl;
        }
        delete [] GravImage;
        nifti_image_free(FPNii);
        nifti_image_free(FNNii);
        for (int l=0; l<numbLabelFP; l++) {
            if(FPVector[l]!=NULL){
                delete FPVector[l];
                FPVector[l]=NULL;
            }
        }
        for (int l=0; l<numbLabelFN; l++) {
            if(FNVector[l]!=NULL){
                delete FNVector[l];
                FNVector[l]=NULL;
            }
        }
        delete TreeToAnalyse;
        return EXIT_SUCCESS;

        
    }
    
    
    
    
    ////    Create card and analysis of ground truth
    //    if(segment_analysis->flag_GT){
    ////        Create connected components for 3 different types of neighborhood for the lesions
    //        int * ComponentLabel=ComponentLabeling(SegToAnalyseHard,segment_analysis->Neigh);
    //        int * OrderedLabels=OrderedVolumeLabel(ComponentLabel, 3, TreeToAnalyse->GetNumberElements());
    //
    //        ImageOrderedLesion=nifti_copy_nim_info(SegToAnalyse);
    //        ImageOrderedLesion->data=(void *) calloc(SegToAnalyse->nvox, sizeof(float));
    //        float * ImageOrderedLesion_PTR=static_cast<float *>(ImageOrderedLesion->data);
    //        for (int i=0; i<numel; i++) {
    //            ImageOrderedLesion_PTR[i]=(float)OrderedLabels[i];
    //        }
    //    }
    

    
    // Used to segment something else than lesion such as vessels
    if (segment_analysis->flag_OtherSeg==1 && TreeToAnalyse!=NULL) {
        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        if(segment_analysis->flag_sc){

            string FilenameSave=FilenamePA_b+"SCTest_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameSaveConnect=FilenamePA_b+"ConnectSCTest_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameSaveInit=FilenamePA_b+"InitSCTest_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameSaveMahal=FilenamePA_b+"MahalGM_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameSaveSCTxt=FilenamePA_b+"TextLesionSC_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".txt";
            string FilenameSaveFinal=FilenamePA_b+"FinalSC_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameSaveRing=FilenamePA_b+"Ring_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameSaveRingPot=FilenamePA_b+"PotRing_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";

            string FilenameSaveBC=FilenamePA_b+"BC_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";
            string FilenameLLWM=FilenamePA_b+"LCWM_"+FilenamePA_e+"-"+segment_analysis->inOptionText+".nii.gz";

            Rule* RuleSC=BuildRuleFromTextFile(TreeToAnalyse,segment_analysis,1);
            numbmodal = TreeToAnalyse->GetNumberModalities();
            numel = TreeToAnalyse->GetNumberElements();
            nifti_image * LesionTotTest=LesionVoxelwise(TreeToAnalyse, RuleSC,Modalities,segment_analysis);
            nifti_set_filenames(LesionTotTest,FilenameSaveInit.c_str(),0,0);
            nifti_image_write(LesionTotTest);
            int Dim[3];
            int Shift[3];
            numel=1;
            vector<int> DimVector;
            float PixDim[3];
            for (int d=0; d<3; d++) {
                Dim[d]=LesionTotTest->dim[d+1];
                PixDim[d]=LesionTotTest->pixdim[d+1];
                if (d==0) {
                    Shift[d]=1;
                }
                else{
                    Shift[d]=Shift[d-1]*Dim[d-1];
                }
                DimVector.push_back(Dim[d]);
                numel*=Dim[d];
            }
            nifti_image * PriorsWM = ReadFromFilename(segment_analysis->filename_inPriorsWM);
            nifti_image * PriorsGM = ReadFromFilename(segment_analysis->filename_inPriorsCGM);
            nifti_image * PriorsCSF = ReadFromFilename(segment_analysis->filename_inPriorsECSF);
            segment_analysis->IndexGM=0;
            segment_analysis->IndexWM=0;
            segment_analysis->IndexCSF=1;
            segment_analysis->IndexOut=2;
            float* WMData = static_cast<float*> (PriorsWM->data);
            float* CGMData = static_cast<float*> (PriorsGM->data);
            float * CSFData = static_cast<float*> (PriorsCSF->data);
            float * LesionData = static_cast<float*> (LesionTotTest->data);
            vector<int> RuleVec;

            RuleVec.push_back(0);
            RuleVec.push_back(1);
            RuleVec.push_back(3);
            RuleVec.push_back(1);
            int numelmasked = TreeToAnalyse->GetNumberMaskedElements();
            float * DataTemp = TreeToAnalyse->GetDataBFCorrected();
            float * DataLesion = &DataTemp[numelmasked];
//            float * NormGMTemp = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();
//            float * NormCSFTemp =TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexCSF)->GetNormResp();
//            float MeanGM=GetMeanData(DataLesion,NormGMTemp,numelmasked);
//            float MeanCSF=GetMeanData(DataLesion,NormCSFTemp,numelmasked);
//            if(MeanGM>MeanCSF){
//                RuleVec[1]=0;
//            }
//            float * NormGM = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();
            vector <TreeEM*> LeavesGM=FindClassesToInlier(TreeToAnalyse,segment_analysis->vecModalities, RuleVec);
            vector <float*> ToAdd;
            int leaves_size = LeavesGM.size();
            float * NormGM=new float[numel];
            bool flag_nomahalcorr=0;
            for(int i=0;i<leaves_size;i++){
                ToAdd.push_back(LeavesGM[i]->GetNormResp());
            }
            if (leaves_size == 0){
                cout << "No possible class to use !!!" << endl;
                flag_nomahalcorr =1 ;
//                segment_analysis->weightThreshold *=2;
                RuleVec[0]=1;
                RuleVec[1]=0;
                 vector <TreeEM*> LeavesCSFValid=FindClassesToInlier(TreeToAnalyse,segment_analysis->vecModalities,RuleVec);
                 leaves_size = LeavesCSFValid.size();
                 cout << "new leaves size is "<< leaves_size<<endl;
                 for(int i=0;i<leaves_size;i++){
                     ToAdd.push_back(LeavesCSFValid[i]->GetNormResp());
                 }

                 NormGM = AddElementwise(ToAdd,TreeToAnalyse->GetNumberMaskedElements(),NULL);
                 float * NormGMTest= TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexCSF)->GetNormResp();


                 float meanGM = GetMeanData(DataLesion,NormGMTest,numelmasked);
                 float varGM = GetVarianceData(DataLesion,NormGMTest,numelmasked);
                 float meanCSF = GetMeanData(DataLesion,NormGM,numelmasked);
                 float varCSF = GetVarianceData(DataLesion,NormGM,numelmasked);
                 cout << varCSF << " " << varGM<<endl;
                    float Dist=sqrt((meanGM-meanCSF)*(meanGM-meanCSF)/varGM);
                    cout << "Dist is "<< Dist<<endl;
                 if (meanCSF-meanGM>0){
                    segment_analysis->weightThreshold += Dist; //*sqrt(varCSF/varGM);
                 }
                 else{
                     if (leaves_size==0){
                         delete [] NormGM;
                         NormGM= TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexCSF)->GetNormResp();
                         segment_analysis->weightThreshold = Dist/2; //*sqrt(varCSF/varGM);
                     }
                     segment_analysis->weightThreshold *=0.8;
                 }

                 cout << "New Threshold is "<< segment_analysis->weightThreshold << endl;
                 }
            else{
                NormGM = AddElementwise(ToAdd,TreeToAnalyse->GetNumberMaskedElements(),NULL);
                cout << NormGM <<endl;
            }
            nifti_image * GMSeg = HardSegmentationThresholdFromNormResp(NormGM,TreeToAnalyse,0.5);
            SaveTmpResult_short(NormGM,g_scratch_dir + "/TestSC.nii.gz",PriorsCSF,TreeToAnalyse->GetL2S());
            float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
            nifti_image * DataNii = CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
            float * GMNorm = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();
            float * WMNorm = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp();
            float * WMOut = TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexWM)->GetNormResp();
            float * Mask = new float[numel];
            bool * CSFBool = new bool[numel];
            bool * WMGMBool = new bool[numel];
            float * WMGMSeg = new float[numel];
            int * L2S = TreeToAnalyse->GetL2S();
            for(int i=0;i<numel;i++){
                Mask[i]=0;
                CSFBool[i] =0;
                WMGMBool[i]=0;
                WMGMSeg[i]=0;
                if (L2S[i]>=0){
                    Mask[i]=1;
                }
                if(CSFData[i]>0){
                    CSFBool[i]=1;
                }
                if (WMData[i]+CGMData[i]>0.05){
                    WMGMBool[i]=1;
                }
                if(WMOut[L2S[i]]+WMNorm[L2S[i]]>0.5){
                    WMGMSeg[i]=1;
                }
            }
//            nifti_image * WMGMSegNii = CreateNiiFromArray(WMGMBoolSeg,PriorsWM,numel);
            int * ComponentWMGMSeg = ComponentLabeling(WMGMSeg,6,Dim,Shift,0.5);
            int * OrderedVolWMGM = OrderedVolumeLabel(ComponentWMGMSeg,3,numel);
            bool * LargestComponentWMGMSeg=CreateLesionBool(OrderedVolWMGM,1,numel);
            SaveTmpResult(LargestComponentWMGMSeg,FilenameLLWM.c_str(),PriorsWM);
            bool * MaskBool=TranscribeArray<float,bool>(Mask,numel);
            bool * BorderMask =CreateBorderFromBool(MaskBool,Dim,Shift);
            nifti_image * FullMaskNii = CreateNiiFromArray(Mask,PriorsWM,numel);
            nifti_image* MahalNii= MahalDistMaps(GMSeg,FullMaskNii,DataNii);
            nifti_image * DistCSF = EuclideanDistanceImage(PriorsCSF, CSFBool,FullMaskNii);
            nifti_image * DistMask = EuclideanDistanceImage(PriorsCSF, BorderMask,FullMaskNii);
            float * Skeleton = ExtractionSkeletonBis(WMGMBool, numel, Dim,Shift,PixDim,PriorsWM);
            bool * SkeletonBool = TranscribeArray<float, bool>(Skeleton,numel);
            nifti_image * DistSkeleton = EuclideanDistanceImage(PriorsCSF,SkeletonBool,FullMaskNii);
            SaveTmpResult(Skeleton,g_scratch_dir + "/SkeletonTest.nii.gz",DistSkeleton);
            SaveTmpResult(WMGMBool,g_scratch_dir + "/SkeletonTestInit.nii.gz",PriorsCSF);
            int CorrespondingInit[3];
            int Corresponding[3];
            bool * PotentialRing=new bool[numel];
            for(int i=0;i<numel;i++){
                PotentialRing[i]=0;
            }
            float * RingPot = CreateRidge(&DataCorr[numel],numel,Dim,Shift,PixDim,2);
            nifti_image * RingPotNii = CreateNiiFromArray(RingPot, PriorsWM,numel);
            nifti_set_filenames(RingPotNii,FilenameSaveRingPot.c_str(),0,0);
            nifti_image_write(RingPotNii);
            nifti_set_filenames(MahalNii,FilenameSaveMahal.c_str(),0,0);
            float * MahalToUse=static_cast<float*> (MahalNii->data);
            float * MahalData = &MahalToUse[2*numel];
            float * MahalDataCorr = &MahalToUse[numel];
            bool * BorderCSF=CreateBorderFromBool(CSFBool,Dim,Shift);
            for (int i=0;i<numel;i++){
                if (BorderCSF[i]==1){
                    float Value=MahalData[i];
                    int Potential=i;
                    CorrespondingCoordinates_bis(CorrespondingInit,i,Dim,Shift);
                    vector<int> IndicesToExplore=GetIndicesToExploreFloatDist(i,segment_analysis->dist_ring,Dim,Shift,PixDim,1);
                    int numbexpl = IndicesToExplore.size();
                    for(int n=0;n<numbexpl;n++){
                        CorrespondingCoordinates_bis(Corresponding,IndicesToExplore[n],Dim,Shift);
                        if(Value>MahalData[IndicesToExplore[n]] && (PotentialRing[IndicesToExplore[n]]==0)&& Corresponding[2]==CorrespondingInit[2]){
                            Potential=IndicesToExplore[n];
                            Value=MahalData[IndicesToExplore[n]];
                        }
                    }
                   PotentialRing[Potential]=1;
                }
            }
//            float * DistSkeletonData=static_cast<float*>(DistSkeleton->data);
//            int * Inside = CreateInsideFromBorder(PotentialRing, SkeletonBool,DistSkeletonData, Dim,Shift,PixDim,1);
            nifti_image * RingNii = CreateNiiFromArray(PotentialRing,PriorsWM,numel);
//            nifti_image * RingNiiInit = CreateNiiFromArray(Inside,PriorsWM,numel);
            nifti_set_filenames(DistMask, (g_scratch_dir + "/MaskTest.nii.gz").c_str(),0,0);
            nifti_set_filenames(RingNii,FilenameSaveRing.c_str(),0,0);
            nifti_image_write(RingNii);
            nifti_image_write(DistMask);
            bool* LesionBool = ThresholdArray<float, bool>(LesionData,0.5,numel);
            float MaxMahal = GetMaxArrayMasked<float, bool>(MahalData,LesionBool,numel);
            if (MaxMahal > 0){
                MaxMahal=0;
            }
            cout << "Max Mahal is "<< MaxMahal << endl;
//            float * DistData=static_cast<float*>(DistCSF->data);
            float * DistMaskData=static_cast<float*>(DistMask->data);
                       nifti_image_write(MahalNii);
            float Added=0;
            for (int i=0;i<numel;i++){
                if (MahalData[i] < MaxMahal && LesionData[i]<0.9 && LargestComponentWMGMSeg[i]==1 ){
                    LesionData[i] += (GMNorm[L2S[i]]+WMNorm[L2S[i]])/2;//-NormGM[L2S[i]];
                    Added+=LesionData[i];
                }
                if (WMData[i]+CGMData[i]<0.2){
                    if(WMOut[L2S[i]]+WMNorm[L2S[i]]<=0.5){
//                        cout << "correction no seg"<<endl;
                        LesionData[i]=0;
                    }
                    if(DistMaskData[i]<2){
//                        cout << "correction dist mask"<<endl;
                        LesionData[i]=0;
                    }
                    if((MahalDataCorr[i]<-1) & !flag_nomahalcorr){
//                        cout << "correction mahal corr"<<endl;
                        LesionData[i]=0;
                    }
                }
                if((MahalDataCorr[i]<-1) & !flag_nomahalcorr){
//                    cout << "correction mahal corr" << endl;
                    LesionData[i]=0;
                }
                if(LargestComponentWMGMSeg[i]==0){
                    LesionData[i]=0;
                }
                if (PotentialRing[i]){
                    vector<int> IndicesClose = GetIndicesToExploreFloatDist(i,PixDim[0]+0.1,Dim,Shift,PixDim,1);

                    int numbclose=IndicesClose.size();
                    int no_correct=0;
                    int no_ring=0;
                    for(int n=0;n<numbclose;n++){
                        if(LesionData[IndicesClose[n]]){
                            no_correct+=1;
                        }
                        if(PotentialRing[IndicesClose[n]]){
                            no_ring+=1;
                        }
                    }

                    no_correct = no_correct>2?1:0;
                    no_ring = no_ring>2?1:0;
                    no_correct *=no_ring;
                    LesionData[i]=no_correct*LesionData[i];
                }

//                if(CGMData[i]>0.275){
//                    LesionData[i]=0;
//                }

                float mult_factor = fabs(MahalData[i]) > segment_analysis->weightThreshold?1:fabs(MahalData[i])/segment_analysis->weightThreshold;
//                cout << "mult_factor is "<<mult_factor<<endl;
                LesionData[i] = LesionData[i]>0.25?mult_factor:0;
            }

            cout << "Added is " << Added<<endl;
            int * LabelData = ComponentLabeling(LesionData,segment_analysis->Neigh,Dim,Shift,0.5);
            int * VolOrdered = OrderedVolumeLabel(LabelData, 6,numel);
            nifti_set_filenames(LesionTotTest,FilenameSaveBC.c_str(),0,0);
            nifti_image_write(LesionTotTest);
            vector <Outlier * > OutlierSC = GetVectorOutliers(VolOrdered, TreeToAnalyse, segment_analysis);
            ofstream SCText(FilenameSaveSCTxt.c_str());
            int * GravImage=NULL;
            int numbFinalLesions=GetMaxLabel(VolOrdered, numel);
            if (numbFinalLesions>0) {
                GravImage=new int[3];

                for (int d=0;d<3; d++) {
                    PixDim[d]=DistCSF->pixdim[d+1];
                    GravImage[d]=(OutlierSC[0]->CentreGravity[d]-OutlierSC[0]->VectorDiffGrav[d]/PixDim[d]);
                }
            }
            for (int l=0; l<numbFinalLesions; l++) {
                cout <<"Label "<<l+1;
                SCText<<endl;
                SCText<<"Label "<<l+1<<endl;
                PrintOutlierCharacteristics(OutlierSC[l],SCText,numbmodal,GravImage,PixDim,segment_analysis);
                cout<<"...printed"<<endl;
            }

            SaveTmpResult(VolOrdered,FilenameSaveConnect,PriorsWM);
            float * FinalLesion = new float[numel];
            for(int i=0;i<numel;i++){
                FinalLesion[i]=0;
            }
            for(int l=0;l<numbFinalLesions;l++){
                cout << "Outlier class for label "<< l+1 << "is "<< OutlierSC[l]->OutlierClass<<endl;
                if (OutlierSC[l]->OutlierClass>0){
                    for(int i=0;i<numel;i++){
                        if (VolOrdered[i]==l+1){
                            FinalLesion[i]=LesionData[i];
                        }
                    }
                }
            }
            nifti_image * FinalLesionNii = CreateNiiFromArray(FinalLesion,PriorsWM,numel);
            nifti_set_filenames(FinalLesionNii, FilenameSaveFinal.c_str(),0,0);
            nifti_set_filenames(LesionTotTest,FilenameSave.c_str(),0,0);
            nifti_image_write(LesionTotTest);
            nifti_image_write(FinalLesionNii);
            nifti_image_free(FinalLesionNii);
            nifti_image_free(PriorsWM);
            nifti_image_free(PriorsGM);
            nifti_image_free(FullMaskNii);
            nifti_image_free(PriorsCSF);
            nifti_image_free(MahalNii);
            nifti_image_free(DistCSF);
            delete [] LesionBool;
            nifti_image_free(LesionTotTest);
            return EXIT_SUCCESS;
        }
        if (segment_analysis->flag_prion ){
            nifti_image * ParcellationNii=ReadFromFilename(segment_analysis->filename_Parc);
            int numel=ParcellationNii->nvox;
            int numbmodal=TreeToAnalyse->GetNumberModalities();
            float * ParcellationData=static_cast<float*>(ParcellationNii->data);
            bool* CaudateMask=new bool[numel];
            if(segment_analysis->flag_maskMatch){
                nifti_image * CaudateTempNii=ReadFromFilename(segment_analysis->filename_maskMatch);
                float * CaudateTempData=static_cast<float*>(CaudateTempNii->data);
                for(int i=0;i<numel;i++){
                    CaudateMask[i]=(bool)CaudateTempData[i];
                }
                nifti_image_free(CaudateTempNii);
                CaudateTempNii=NULL;
            }
            else{
                for (int i=0;i<numel;i++){
                    CaudateMask[i]=(ParcellationData[i]==37||ParcellationData[i]==38)?1:0;
                }
            }
            nifti_image * CaudateNii=CreateNiiFromArray(CaudateMask,ParcellationNii,numel);
            float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
            cout <<"Data corr obtained"<<endl;
            nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
            bool * CGMData=ThresholdArray<float,bool>(ParcellationData,100,numel);
            nifti_image * CGMNii = CreateNiiFromArray(CGMData,ParcellationNii,numel);
            nifti_image * MahalMaps=MahalDistMaps(CaudateNii,CGMNii,DataNii);
            float * MahalToUse=static_cast<float*>(MahalMaps->data);
            vector<float> Thresholds;
            if(segment_analysis->vecThresh.size()==0){
                Thresholds.push_back(3);
                Thresholds.push_back(4);
                Thresholds.push_back(5);
            }
            else{
                for(float t=segment_analysis->vecThresh[0];t<segment_analysis->vecThresh[1];t+=segment_analysis->vecThresh[2]){
                    Thresholds.push_back(t);
                }
            }
            int ThreshSize=Thresholds.size();
            float * SegmentationTest=new float[numel*Thresholds.size()];
            for(int i=0;i<numel;i++){
                for(int t=0;t<ThreshSize;t++){
                    if(ParcellationData[i]>100 || (ParcellationData[i]>23 && ParcellationData[i]<35)
                            || ParcellationData[i]==37 || ParcellationData[i]==38 ||
                            (ParcellationData[i]>54 && ParcellationData[i]<63) ||
                            ParcellationData[i]==76 || ParcellationData[i]==77){
                        SegmentationTest[t*numel+i]=MahalToUse[2*numel+i]>Thresholds[t]/2?MahalToUse[2*numel+i]/Thresholds[t]:0;

                    }
                    else{
                        SegmentationTest[t*numel+i]=0;
                    }
                }
            }
            nifti_image * SegmentationPrion=CreateNiiFromArray(SegmentationTest,MahalMaps,ThreshSize*numel);
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenameSaveMahal=FilenamePA_b+"Mahal_"+FilenamePA_e+".nii.gz";
            string FinalSaveOtherSeg=FilenamePA_b+"OtherSeg_"+FilenamePA_e+".nii.gz";
            nifti_set_filenames(SegmentationPrion,FinalSaveOtherSeg.c_str(),0,0);
            nifti_set_filenames(MahalMaps,FilenameSaveMahal.c_str(),0,0);
            nifti_image_write(MahalMaps);
            nifti_image_write(SegmentationPrion);
            nifti_image_free(MahalMaps);
            nifti_image_free(SegmentationPrion);
            delete [] CaudateMask;
            nifti_image_free(CGMNii);
            nifti_image_free(ParcellationNii);
            return EXIT_SUCCESS;

        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveOtherSeg=FilenamePA_b+"OtherSeg_"+FilenamePA_e+".txt";
        string FinalSaveOtherSeg=FilenamePA_b+"OtherSeg_"+FilenamePA_e+".nii.gz";
        ofstream TxtFile(FilenameSaveOtherSeg.c_str());
        Rule * OtherRule=BuildRuleFromTextFile(TreeToAnalyse, segment_analysis,0);
        vector<TreeEM *> OtherSegClasses=FindLesionClasses(TreeToAnalyse, OtherRule, Modalities,segment_analysis);
        int numbOtherClasses=OtherSegClasses.size();
        if (numbOtherClasses>0) {
            TxtFile << "NumbOtherClasses "<<numbOtherClasses<<endl;
            for (int les=0; les<numbOtherClasses; les++) {
                vector<int> HierarchyOtherClass=OtherSegClasses[les]->GetHierarchyVector();
                int lengthHierarchy=HierarchyOtherClass.size();
                for (int l=0; l<lengthHierarchy; l++) {
                    TxtFile<<HierarchyOtherClass[l]<<" ";
                }
                TxtFile<<endl;
            }
        }
        else{
            TxtFile<<"NumbOtherSegClasses 0"<<endl;
        }
        //        nifti_image * OtherSegClassesImage=ReconstructLesionImage(TreeToAnalyse, OtherRule, Modalities,segment_analysis);
        //        nifti_set_filenames(OtherSegClassesImage, (g_scratch_dir + "/StrangeRebuilt.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(OtherSegClassesImage);
        //        nifti_image * OtherSegPartsFromUniform=LesionFromUniform(TreeToAnalyse, OtherRule, Modalities,OtherSegClasses,segment_analysis);
        nifti_image * OtherSegTot=ReconstructLesionImageTot(TreeToAnalyse, OtherRule, Modalities, segment_analysis);
        nifti_set_filenames(OtherSegTot, FinalSaveOtherSeg.c_str(), 0, 0);
        nifti_image_write(OtherSegTot);
        cout<<"OtherSeg done "<<endl;
    }
    
    
    //    // Used to segment something else than lesion such as vessels
    //    if (segment_analysis->flag_LacuneSeg==1 && TreeToAnalyse!=NULL && ) {
    //        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
    //        int Index=FilenamePA.find_last_of('/');
    //        FilenamePA_b=FilenamePA.substr(0,Index+1);
    //        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
    //        string FilenameSaveOtherSeg=FilenamePA_b+"PotLacSeg_"+FilenamePA_e+".txt";
    //        string FinalSaveOtherSeg=FilenamePA_b+"PotLacSeg_"+FilenamePA_e+".nii.gz";
    //        ofstream TxtFile(FilenameSaveOtherSeg.c_str());
    //
    ////        First build the potential WM + DGM + ICSF Segmentation and WM+DGM Segmentation
    //        WMDGMVentrSeg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,Ventr_priors,Ventr_idx,L2SToUse,ThreshPriors);
    //        WMDGMSeg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,Ventr_priors,Ventr_idx,L2SToUse);
    //
    //
    ////        First give the hypointense outliers
    //
    ////        Then take care of the inliers of GM or CSF that are very close to CSF or GM/CSF PV and should be considered as outliers
    //        Rule * OtherRule=BuildRuleFromTextFile(TreeToAnalyse, segment_analysis,0);
    //        vector<TreeEM *> OtherSegClasses=FindLesionClasses(TreeToAnalyse, OtherRule, Modalities,segment_analysis);
    //        int numbOtherClasses=OtherSegClasses.size();
    //        if (numbOtherClasses>0) {
    //            TxtFile << "NumbOtherClasses "<<numbOtherClasses<<endl;
    //            for (int les=0; les<numbOtherClasses; les++) {
    //                vector<int> HierarchyOtherClass=OtherSegClasses[les]->GetHierarchyVector();
    //                int lengthHierarchy=HierarchyOtherClass.size();
    //                for (int l=0; l<lengthHierarchy; l++) {
    //                    TxtFile<<HierarchyOtherClass[l]<<" ";
    //                }
    //                TxtFile<<endl;
    //            }
    //        }
    //        else{
    //            TxtFile<<"NumbOtherSegClasses 0"<<endl;
    //        }
    //        nifti_image * OtherSegClassesImage=ReconstructLesionImage(TreeToAnalyse, OtherRule, Modalities,segment_analysis);
    //        //        nifti_set_filenames(OtherSegClassesImage, (g_scratch_dir + "/StrangeRebuilt.nii.gz").c_str(), 0, 0);
    //        //        nifti_image_write(OtherSegClassesImage);
    //        //        nifti_image * OtherSegPartsFromUniform=LesionFromUniform(TreeToAnalyse, OtherRule, Modalities,OtherSegClasses,segment_analysis);
    //        nifti_image * OtherSegTot=ReconstructLesionImageTot(TreeToAnalyse, OtherRule, Modalities, segment_analysis);
    //        nifti_set_filenames(OtherSegTot, FinalSaveOtherSeg.c_str(), 0, 0);
    //        nifti_image_write(OtherSegTot);
    //        cout<<"OtherSeg done "<<endl;
    //    }
    
    string FilenameRefConnectAnalysis;
    string FilenameRefConnectCorrAnalysis;
    string FilenameIntensityConnectAnalysis;
    string FilenameIntensityConnectCorrAnalysis;
    string FilenameVentricle;
    if(segment_analysis->flag_outTxt){
        string FilenameTxtOut=segment_analysis->filename_OutTxt;
        int Index=FilenameTxtOut.find_last_of('.');
        string FilenameTmp=FilenameTxtOut.substr(0,Index);
        FilenameRefConnectAnalysis=FilenameTmp+"_refConnect.txt";
        FilenameIntensityConnectAnalysis+"_refIntensity.txt";
        if (segment_analysis->flag_correct || segment_analysis->flag_inLesCorr) {
            FilenameRefConnectCorrAnalysis=FilenameTmp+"_refConnect-corr.txt";
            FilenameIntensityConnectCorrAnalysis=FilenameTmp+"_refIntensity-corr.txt";
        }
    }
    if (segment_analysis->flag_ConnectRefAnalysis) {
        string FilenameTmp=nifti_makebasename(segment_analysis->filename_inConnect);
        FilenameRefConnectAnalysis=FilenameTmp+"_refConnect.txt";
        FilenameIntensityConnectAnalysis+"_refIntensity.txt";
        if (segment_analysis->flag_correct || segment_analysis->flag_inLesCorr) {
            FilenameRefConnectCorrAnalysis=FilenameTmp+"_refConnect-corr.txt";
            FilenameIntensityConnectCorrAnalysis=FilenameTmp+"_refIntensity-corr.txt";
        }
    }
    
    string FilenameSave;
    if (segment_analysis->flag_SegTot) {
        FilenameSave=nifti_makebasename(segment_analysis->filename_SegTot);
    }
    else if(segment_analysis->flag_in){
        FilenameSave=nifti_makebasename(segment_analysis->filename_InLes);
    }
    else if(segment_analysis->flag_inLesCorr){
        FilenameSave=nifti_makebasename(segment_analysis->filename_inLesCorr);
    }
    else if(segment_analysis->flag_inConnect){
        FilenameSave=nifti_makebasename(segment_analysis->filename_inConnect);
    }
    
    
    FilenameIntensityConnectAnalysis=FilenameSave+"_refIntensity.txt";
    if (segment_analysis->flag_correct || segment_analysis->flag_inLesCorr) {
        stringstream CorrLev;
        CorrLev <<segment_analysis->flag_correctionLevel;
        FilenameIntensityConnectCorrAnalysis=FilenameSave+"_refIntensity-corr_CL"+CorrLev.str()+".txt";
    }
    
    if (segment_analysis->flag_TextFile && !segment_analysis->flag_inVentricleSeg) {
        string FilenameTxtIn=segment_analysis->filename_TextFile;
        int Index=FilenameTxtIn.find_last_of('.');
        string FilenameTmp=FilenameTxtIn.substr(0,Index);
        FilenameVentricle=FilenameTmp+"_Ventricles.nii.gz";
    }
    cout << FilenameVentricle<<endl;
    if (segment_analysis->flag_LaplaceSolImage && segment_analysis->flag_layercreation && segment_analysis->flag_Parc){
        nifti_image*  LaplaceSolNii = ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
        nifti_image * MaskNii = ReadFromFilename(segment_analysis->filename_mask);
        nifti_image * ParcNii = ReadFromFilename(segment_analysis->filename_Parc);
        cout << LaplaceSolNii->nvox << " " << ParcNii->nvox <<" " << MaskNii->nvox << endl;
        int Dim[3];
        int Shift[3];
        numel=1;
        vector<int> DimVector;
        float PixDim[3];
        for (int d=0; d<3; d++) {
            Dim[d]=ParcNii->dim[d+1];
            PixDim[d]=ParcNii->pixdim[d+1];
            if (d==0) {
                Shift[d]=1;
            }
            else{
                Shift[d]=Shift[d-1]*Dim[d-1];
            }
            DimVector.push_back(Dim[d]);
            numel*=Dim[d];
        }

        float * ParcData  = static_cast<float*>(ParcNii->data);
        bool * MaskData = static_cast<bool*>(MaskNii->data);
        float * LaplaceSol = static_cast<float*>(LaplaceSolNii->data);
        int * L2S = MakeL2S(MaskData, Dim);
        int numelmasked=0;
        int * S2L = MakeS2L(MaskData, Dim, numelmasked);
        cout << numel << " is numel and numelmasked is "<< numelmasked <<endl;
        float * Parenchyma = UpperThresholdArray<float, float>(ParcData, 99, numel);
        float * ParenchymaCorr = ThresholdArray<float, float>(Parenchyma, 4.5, numel);
        bool * ParenchymaBool = TranscribeArray<float, bool>(ParenchymaCorr, numel);
        bool * ParenchShort = CreateShort(ParenchymaBool,S2L,numelmasked);

        float * Ventricles = CreateVentriclesFromParc(ParcData, numel);
        float * CorrectionVentricles = new float[numel];
        for (int i=0; i<numel ; i++){
            Ventricles[i] = LaplaceSol[i] > 0?0:Ventricles[i];

        }
        bool * VentriclesBool = TranscribeArray<float, bool>(Ventricles, numel);
        bool * VentriclesShort = CreateShort(VentriclesBool, S2L, numelmasked);


        string FilenamePA=nifti_makebasename(segment_analysis->filename_LaplaceSolImage);
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());

        float * NormDist = NormalisedLaplaceLength(LaplaceSol, VentriclesShort, ParenchShort, L2S, S2L, numelmasked, PixDim, Dim, Shift);
        int numblaysize = segment_analysis->numbLaminae.size();
        for (int nl=0; nl<numblaysize; nl++){
            int * ShortLayers=CreateLayersFromNormDistSol(NormDist,segment_analysis->numbLaminae[nl],numelmasked,0.0001);
            int * LongLayers = CreateLong(ShortLayers, L2S, numel);
            nifti_image * LongLayersNii = CreateNiiFromArray(LongLayers, ParcNii, numel);
            stringstream nlName;
            nlName << segment_analysis->numbLaminae[nl];
            string FilenameSave = FilenamePA_b + "Layers_" + nlName.str() + FilenamePA_e + ".nii.gz";
            nifti_set_filenames(LongLayersNii, FilenameSave.c_str(), 0,0);
            nifti_image_write(LongLayersNii);
            nifti_image_free(LongLayersNii);
            LongLayersNii=NULL;
            delete [] LongLayers;
            LongLayers = NULL;
            delete [] ShortLayers;
            ShortLayers = NULL;
        }
        nifti_image_free(MaskNii);
        nifti_image_free(LaplaceSolNii);
        nifti_image_free(ParcNii);
        delete [] S2L;
        delete [] L2S;
        delete [] NormDist;
        delete [] VentriclesBool;
        delete [] VentriclesShort;
        delete [] Ventricles;
        delete [] Parenchyma;
        delete [] ParenchymaCorr;
        delete [] ParenchymaBool;
        delete [] ParenchShort;
        return EXIT_SUCCESS;
    }
    if (segment_analysis->flag_LapAnalysis && segment_analysis->flag_vecIOLap && !segment_analysis->flag_inSum) {
        int sizeVecIO=segment_analysis->filenames_IOLap.size();
        if (sizeVecIO!=2) {
            cout<<"Improper number of elements in IOLapVector"<<endl;
        }
        else {
            FilenamePA=nifti_makebasename(segment_analysis->filenames_IOLap[0]);
            int Index=FilenamePA.find_last_of('/');
            FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            bool * Mask=NULL;
            nifti_image * MaskNii=NULL;
            if (segment_analysis->flag_mask) {
                MaskNii=ReadFromFilename(segment_analysis->filename_mask);
                Mask=static_cast<bool*>(MaskNii->data);
            }
            nifti_image * InImage=ReadFromFilename(segment_analysis->filenames_IOLap[0]);
            nifti_image * OutImage=ReadFromFilename(segment_analysis->filenames_IOLap[1]);
            int numel=1;
            int Dim[3];
            int Shift[3];
            vector<int> DimVector;
            float PixDim[3];
            for (int d=0; d<3; d++) {
                Dim[d]=InImage->dim[d+1];
                PixDim[d]=InImage->pixdim[d+1];
                if (d==0) {
                    Shift[d]=1;
                }
                else{
                    Shift[d]=Shift[d-1]*Dim[d-1];
                }
                DimVector.push_back(Dim[d]);
                numel*=Dim[d];
            }
            string FilenamePA_DilatedOut=FilenamePA_b+"DilatedOut_"+FilenamePA_e+".nii.gz";

            // Case where we are doing the layers within the lesions
            if(strcmp(segment_analysis->filenames_IOLap[0],segment_analysis->filenames_IOLap[1])==0){
                int * Labels=ComponentLabeling(InImage,6);
                int * OrderedLabels=OrderedVolumeLabel(Labels,2,numel);
                numel=InImage->nvox;
                string FilenamePA_Comp=FilenamePA_b+"LesionComp_"+FilenamePA_e+".nii.gz";
                string FilenamePA_CoM=FilenamePA_b+"LesionCoM_"+FilenamePA_e+".nii.gz";
                string FilenamePA_LapTemp=FilenamePA_b+"LaplaceTemp_"+FilenamePA_e+".nii.gz";

                nifti_image * LabelsNii=CreateNiiFromArray(OrderedLabels,InImage,numel);
                nifti_set_filenames(LabelsNii,FilenamePA_Comp.c_str(),0,0);
                nifti_image_write(LabelsNii);

                string NameLap;
                if (segment_analysis->flag_nameLap) {
                    NameLap=segment_analysis->nameLap;
                }
                int maxLabel=GetMaxLabel(OrderedLabels,numel);
                bool * CoMBool=new bool[numel];
                for(int i=0;i<numel;i++){
                    CoMBool[i]=0;
                }

                int nlsize=segment_analysis->numbLaminae.size();
                vector<float *> FinalResult;
                float * FinalNormDist=new float [numel];
                for(int i=0;i<numel;i++){
                    FinalNormDist[i]=0;
                }
                for(int n=0;n<nlsize;n++){
                    float * FinalResultTemp=new float[numel];
                    for(int i=0;i<numel;i++){
                        FinalResultTemp[i]=0;
                    }
                    FinalResult.push_back(FinalResultTemp);
                }

                for(int l=0;l<maxLabel;l++){
                    cout <<" treating label "<<l<<endl;
                    bool * OutBool=CreateLesionBool(OrderedLabels,l+1,numel);
                    bool * LesionCheck=CopyArray(OutBool,numel);
                    float * ExtractedSkel=ExtractionSkeleton(OutBool,numel,Dim,Shift,PixDim,LabelsNii);
                    bool *FinalSkel=NULL;
                    if(ExtractedSkel!=NULL){
                        SaveTmpResult(ExtractedSkel,g_scratch_dir + "/CheckSkel.nii.gz",LabelsNii);
                        float * AbsExtractedSkel=Absolute(ExtractedSkel,numel);
                        FinalSkel=ThresholdArray<float,bool>(AbsExtractedSkel,1.5,numel);
                        delete [] AbsExtractedSkel;
                    }
                    float * NormDist=NULL;
                    if (FinalSkel!=NULL && CountNonZero(FinalSkel,numel)>0){
                        NormDist=NormDistWithinLesionRidge(LesionCheck,FinalSkel,numel,DimVector,Dim,Shift,PixDim);
                    }
                    else{
                        NormDist= NormDistWithinLesion(OutBool,numel,DimVector,Dim, Shift, PixDim, 0, LesionCheck,LabelsNii);
                    }

                    stringstream ls;
                    ls << l;

                    cout<<"Max Norm Dist is "<<GetMax(NormDist,LesionCheck,numel)<< "for label "<<l<<endl;

                    for(int i=0;i<numel;i++){
                        if (FinalNormDist[i]==0){
                            FinalNormDist[i]+=NormDist[i];
                        }
                    }
                    for (int nl=0; nl<nlsize; nl++) {
                        stringstream nlName;
                        nlName << segment_analysis->numbLaminae[nl];
                        //                                int * LayersNormDist=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_s,segment_analysis->numbLaminae[nl],numelmaskedLapOut);
                        //                                int * LongLayers=CreateLong(LayersNormDist,L2S_Out,numel);
                        int * LongLayers=CreateLayersFromNormDistSol(NormDist,segment_analysis->numbLaminae[nl],numel,0.0001);
                        for(int i=0;i<numel;i++){
                            FinalResult[nl][i]+=LongLayers[i];
                            FinalNormDist[i]+=NormDist[i];
                            if (LesionCheck[i]&& LongLayers[i]==0){
                                LongLayers[i]=1;
                                FinalResult[nl][i]=1;
                            }
                        }
                        string FilenamePA_NormDist=FilenamePA_b+"NormDist_Label"+ls.str()+NameLap+"_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";

                        string FilenamePA_LaplaceLayers=FilenamePA_b+"LaplaceLayers_Label"+ls.str()+NameLap+"_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                        SaveTmpResult(LongLayers, FilenamePA_LaplaceLayers, InImage);
                        SaveTmpResult(NormDist,FilenamePA_NormDist.c_str(),InImage);
                        //                                SaveTmpResult_short(LayersNormDist, FilenamePA_LaplaceLayers, InImage, L2S_Out);
                        //                                if (LayersNormDist!=NULL) {
                        //                                    delete [] LayersNormDist;
                        //                                    LayersNormDist=NULL;
                        //                                }
                        if (LongLayers!=NULL) {
                            delete [] LongLayers;
                            LongLayers=NULL;
                        }
                    }

                    //                            float *DistanceNormalisedLaplace  =CreateLong(DistanceNormalisedLaplace_s, L2S_Out, numel);
                    //                            string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplace_Label"+ls.str()+NameLap+"_"+FilenamePA_e+".nii.gz";
                    //                            nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplace, InImage,numel);
                    //                            nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
                    //                            nifti_image_write(NormLaplaceDistImage);
                    //                            nifti_image_free(NormLaplaceDistImage);
                    //                            NormLaplaceDistImage=NULL;
                    //                            delete [] DistanceNormalisedLaplace_s;
                    //                            delete [] Out_short;
                    //                            delete [] Border0;
                    //                            delete [] Border1;
                    //                            delete [] In_short;
                    //                            Out_short=NULL;
                    //                            Border0=NULL;
                    //                            Border1=NULL;
                    //                            In_short=NULL;
                    //                            DistanceNormalisedLaplace_s=NULL;
                    delete [] NormDist;
                    delete [] OutBool;
                    delete [] LesionCheck;

                }
                for(int n=0;n<nlsize;n++){
                    stringstream nlName;
                    nlName << segment_analysis->numbLaminae[n];
                    string LaplaceLayersTemp=FilenamePA_b+"LaplaceLayers"+NameLap+"_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                    SaveTmpResult(FinalResult[n],LaplaceLayersTemp.c_str(),InImage);
                }
                string NormDistName=FilenamePA_b+"NormDist"+NameLap+"_"+FilenamePA_e+".nii.gz";

                SaveTmpResult(FinalNormDist,NormDistName.c_str(),InImage);
                //                nifti_image_free(InImage);
                //                InImage=CreateNiiFromArray(CoMBool,OutImage,numel);
                //                nifti_set_filenames(InImage,FilenamePA_CoM.c_str(),0,0);
                //                nifti_image_write(InImage);
            }
            else{
                bool CompDim=CheckCompatibleDimensions(InImage, OutImage);
                if (CompDim) {

                    int * L2SToUse=NULL;
                    if (Mask!=NULL) {
                        L2SToUse=MakeL2S(Mask, Dim);
                    }
                    else {
                        L2SToUse=new int[numel];
                        for (int i=0; i<numel; i++) {
                            L2SToUse[i]=i;
                        }
                    }
                    float * InData=static_cast<float *>(InImage->data);
                    float * OutData=static_cast<float *>(OutImage->data);
                    bool * InBool=new bool[numel];
                    bool * OutBool=new bool[numel];
                    cout << "numel is "<<numel<<endl;
                    for (int i=0; i<numel; i++) {
                        InBool[i]=InData[i]>0.5?1:0;
                        OutBool[i]=OutData[i]>0.5?1:0;
                    }
                    float * TmpOut=TranscribeArray<bool, float>(OutBool, numel);
                    cout << "Number In InData"<< CountNonZero(InBool,numel) << " "<< CountNonZero(InData, numel)<<endl;

                    float * DilatedOut=Erosion_bis(TmpOut, 3, DimVector, 0);
                    nifti_image * DilatedNii=CreateNiiFromArray(DilatedOut,OutImage,numel);
                    nifti_set_filenames(DilatedNii,FilenamePA_DilatedOut.c_str(),0,0);
                    nifti_image_write(DilatedNii);

                    delete [] TmpOut;
                    TmpOut=NULL;

                    bool * MaskDilOut=TranscribeArray<float, bool>(DilatedOut, numel);

                    int * L2S_Out=MakeL2S(OutBool, Dim);
                    int numelmaskedLapOut=0;
                    int * S2L_Out=MakeS2L(OutBool, Dim, numelmaskedLapOut);

                    int * LabelLaplace=CreateLaplaceLabelling(InBool,OutBool,Dim);
                    SaveTmpResult(LabelLaplace, g_scratch_dir + "/TestLabel.nii.gz", InImage);
                    int numelmaskedbOut=0;
                    int * L2SbOut=MakeL2S(MaskDilOut, Dim);
                    int * S2LbOut=MakeS2L(MaskDilOut, Dim, numelmaskedbOut);
                    if (MaskDilOut!=NULL) {
                        delete [] MaskDilOut;
                        MaskDilOut=NULL;
                    }
                    int * Label_s=CreateShort(LabelLaplace, S2LbOut, numelmaskedbOut);
                    cout<<GetMax<int>(Label_s, numelmaskedbOut)<<endl;
                    //            SaveTmpResult_short(Label_sPar, g_scratch_dir + "/TestLabelB.nii.gz", WMDGMSeg, L2SbPar);
                    if (LabelLaplace!=NULL) {
                        delete [] LabelLaplace;
                        LabelLaplace=NULL;
                    }
                    float * LaplaceSol=NULL;
                    cout <<"Mapping is "<< segment_analysis->flag_MapLap<<endl;
                    if (!segment_analysis->flag_LaplaceSolImage){
                        float * LaplaceSol_short=SolvingLaplaceEquation(Label_s, Dim, Shift, PixDim, 5, L2SbOut, S2LbOut, numelmaskedbOut,InImage);
                         LaplaceSol=CreateLong(LaplaceSol_short,L2SbOut,numel);
                        if (LaplaceSol_short!=NULL) {
                            delete [] LaplaceSol_short;
                        }
                        LaplaceSol_short=NULL;
                       }
                    else{
                        nifti_image * LaplaceNii = ReadFromFilename(segment_analysis->filename_LaplaceSolImage);
                        LaplaceSol = static_cast<float*>(LaplaceNii->data);
                    }
                    
                    cout <<"Laplace solution obtained"<<endl;
                    // Clearing memory for all that won't be needed afterwards
                    delete [] L2SbOut;
                    delete [] S2LbOut;
                    delete [] Label_s;

                    
                    L2SbOut=NULL;
                    S2LbOut=NULL;
                    Label_s=NULL;

                    string FilenamePA_LaplaceSol=FilenamePA_b+"LaplaceSol"+"_"+FilenamePA_e+".nii.gz";
                    SaveTmpResult(LaplaceSol,FilenamePA_LaplaceSol,InImage);
                    bool * In_short=CreateShort(InBool, S2L_Out, numelmaskedLapOut);
                    bool * Out_short=CreateShort(OutBool, S2L_Out, numelmaskedLapOut);
                    bool * Border0=CreateBorderFromBool_L2S(In_short, Dim, Shift, L2S_Out, S2L_Out, numelmaskedLapOut);
                    //            SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", SegToAnalyse, L2S);
                    bool * Border1=CreateBorderFromBool_L2S(Out_short, Dim, Shift, L2S_Out, S2L_Out, numelmaskedLapOut);
                    float * LaplaceSol_s=CreateShort(LaplaceSol, S2L_Out, numelmaskedLapOut);
                    SaveTmpResult_short(LaplaceSol_s,g_scratch_dir + "/LaplaceSol.nii.gz",InImage,L2S_Out);
                    if (segment_analysis->flag_MapLap){
                        cout <<"Do Mapping" << endl;
                        float * MappingLaplace = LaplaceMappingExt( LaplaceSol_s, Out_short, L2S_Out, S2L_Out, numelmaskedLapOut, PixDim, Dim, Shift,InImage);
                        SaveTmpResult_short(MappingLaplace,g_scratch_dir + "/MappingTest.nii.gz",InImage,L2S_Out);
                        nifti_image * CorrespondingNii = ReadFromFilename(segment_analysis->filename_inLesCorr);
                        float * CorrectedMapping = CorrectMapping(MappingLaplace, Border1,L2S_Out,S2L_Out,Dim,Shift,PixDim,numelmaskedLapOut);
                        SaveTmpResult_short(CorrectedMapping,g_scratch_dir + "/MappingTest2.nii.gz",InImage,L2S_Out);

                        float * Mapped = PerformMapping(CorrectedMapping, CorrespondingNii, L2S_Out);
                        SaveTmpResult(Mapped,g_scratch_dir + "/Mapfinal.nii.gz",InImage);
                        cout << "Mapping done!";
                        nifti_image_free(CorrespondingNii);
                        delete [] MappingLaplace;
                        delete [] Mapped;

                    }
                    cout <<"Perform Distance normalisation"<<endl;
                    float * DistanceNormalisedLaplace_s=NormalisedLaplaceLength(LaplaceSol_s, In_short, Out_short, L2S_Out, S2L_Out, numelmaskedLapOut, PixDim, Dim, Shift,InImage);
                    string NameLap;
                    if (segment_analysis->flag_nameLap) {
                        NameLap=segment_analysis->nameLap;
                    }
                    int nlsize=segment_analysis->numbLaminae.size();
                    for (int nl=0; nl<nlsize; nl++) {
                        stringstream nlName;
                        nlName << segment_analysis->numbLaminae[nl];
                        int * LayersNormDist=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_s,segment_analysis->numbLaminae[nl],numelmaskedLapOut);
                        string FilenamePA_LaplaceLayers=FilenamePA_b+"LaplaceLayers"+NameLap+"_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                        SaveTmpResult_short(LayersNormDist, FilenamePA_LaplaceLayers, InImage, L2S_Out);
                        if (LayersNormDist!=NULL) {
                            delete [] LayersNormDist;
                            LayersNormDist=NULL;
                        }
                    }
                    
                    float *DistanceNormalisedLaplace  =CreateLong(DistanceNormalisedLaplace_s, L2S_Out, numel);

                    string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplace"+NameLap+"_"+FilenamePA_e+".nii.gz";
                    nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplace, InImage,numel);
                    nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
                    nifti_image_write(NormLaplaceDistImage);
                    nifti_image_free(NormLaplaceDistImage);
                    NormLaplaceDistImage=NULL;
                    delete [] DistanceNormalisedLaplace_s;
                    delete [] Out_short;
                    delete [] Border0;
                    delete [] Border1;
                    delete [] In_short;
                    Out_short=NULL;
                    Border0=NULL;
                    Border1=NULL;
                    In_short=NULL;
                    DistanceNormalisedLaplace_s=NULL;
                }
            }
            if (MaskNii!=NULL) {
                nifti_image_free(MaskNii);
                MaskNii=NULL;
            }
            if (InImage!=NULL) {
                nifti_image_free(InImage);
                InImage=NULL;
            }
            if (OutImage!=NULL) {
                nifti_image_free(OutImage);
                OutImage=NULL;
            }
        }
        return EXIT_SUCCESS;

    }
    
    if (segment_analysis->flag_LapAnalysis && !segment_analysis->flag_TextFile && segment_analysis->flag_inSum) {
        
        nifti_image * SummarisedSeg=ReadFromFilename(segment_analysis->filename_inSum);
        int Dim[3];
        numel=1;
        for (int d=0; d<3; d++) {
            Dim[d]=SummarisedSeg->dim[d+1];
            numel*=Dim[d];
        }
        
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_inSum);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameSaveTest=FilenamePA_b+"TestSegLaplace_"+FilenamePA_e+".nii.gz";
        nifti_image * MaskTot=NULL;
        bool * Mask=NULL;
        int * L2SToUse=NULL;
        if (segment_analysis->flag_mask) {
            MaskTot=ReadFromFilename(segment_analysis->filename_mask);
            Binarisation(MaskTot);
            Mask=static_cast<bool*>(MaskTot->data);
            L2SToUse=MakeL2S(Mask, Dim);
        }
        
        
        else {
            L2SToUse=new int[numel];
            for (int i=0; i<numel; i++) {
                L2SToUse[i]=i;
            }
        }
        cout<<"Going into Laplace analysis"<<endl;
        int nlsize=segment_analysis->numbLaminae.size();
        vector<int> Ventr_idx;
        if (segment_analysis->IndexICSF >-1) {
            Ventr_idx.push_back(segment_analysis->IndexICSF);
        }
        else if (segment_analysis->IndexCSF >-1){
            Ventr_idx.push_back(segment_analysis->IndexCSF);
        }
        vector<int> WMDGM_idx;
        WMDGM_idx.push_back(segment_analysis->IndexWM);
        if (segment_analysis->IndexDGM>-1) {
            WMDGM_idx.push_back(segment_analysis->IndexDGM);
        }
        else if(segment_analysis->IndexGM>-1){
            WMDGM_idx.push_back(segment_analysis->IndexGM);
        }
        vector<int> CGM_idx;
        if (segment_analysis->IndexCGM>-1) {
            CGM_idx.push_back(segment_analysis->IndexCGM);
        }
        else if (segment_analysis->IndexGM>-1){
            CGM_idx.push_back(segment_analysis->IndexGM);
        }
        vector<int> DGMVentr_idx;
        if (segment_analysis->IndexDGM>-1) {
            DGMVentr_idx.push_back(segment_analysis->IndexDGM);
        }
        else if(segment_analysis->IndexGM>-1){
            DGMVentr_idx.push_back(segment_analysis->IndexGM);
        }
        if (segment_analysis->IndexICSF >-1) {
            DGMVentr_idx.push_back(segment_analysis->IndexICSF);
        }
        else if (segment_analysis->IndexCSF >-1){
            DGMVentr_idx.push_back(segment_analysis->IndexCSF);
        }
        vector<int> WMDGMVentr_idx;
        WMDGMVentr_idx.push_back(segment_analysis->IndexWM);
        if (segment_analysis->IndexBrainstem>-1) {
            WMDGMVentr_idx.push_back(segment_analysis->IndexBrainstem);
        }
        if (segment_analysis->IndexDGM>-1) {
            WMDGMVentr_idx.push_back(segment_analysis->IndexDGM);
        }
        else if(segment_analysis->IndexGM>-1){
            WMDGMVentr_idx.push_back(segment_analysis->IndexGM);
        }
        if (segment_analysis->IndexICSF >-1) {
            WMDGMVentr_idx.push_back(segment_analysis->IndexICSF);
        }
        else if (segment_analysis->IndexCSF >-1){
            WMDGMVentr_idx.push_back(segment_analysis->IndexCSF);
        }

        vector<int> GMWMVentr_idx;
        GMWMVentr_idx.push_back(segment_analysis->IndexWM);
        if (segment_analysis->IndexBrainstem>-1) {
            GMWMVentr_idx.push_back(segment_analysis->IndexBrainstem);
        }
        if (segment_analysis->IndexDGM>-1 && segment_analysis->IndexCGM >-1) {
            GMWMVentr_idx.push_back(segment_analysis->IndexDGM);
            GMWMVentr_idx.push_back(segment_analysis->IndexCGM);
        }
        else if(segment_analysis->IndexGM>-1){
            GMWMVentr_idx.push_back(segment_analysis->IndexGM);
        }
        if (segment_analysis->IndexICSF >-1) {
            GMWMVentr_idx.push_back(segment_analysis->IndexICSF);
        }
        else if (segment_analysis->IndexCSF >-1){
            GMWMVentr_idx.push_back(segment_analysis->IndexCSF);
        }
        vector<int> GMWM_idx;
        GMWMVentr_idx.push_back(segment_analysis->IndexWM);
        if (segment_analysis->IndexBrainstem>-1) {
            WMDGMVentr_idx.push_back(segment_analysis->IndexBrainstem);
        }
        if (segment_analysis->IndexDGM>-1 && segment_analysis->IndexCGM >-1) {
            GMWMVentr_idx.push_back(segment_analysis->IndexDGM);
            GMWMVentr_idx.push_back(segment_analysis->IndexCGM);
        }
        else if(segment_analysis->IndexGM>-1){
            GMWMVentr_idx.push_back(segment_analysis->IndexGM);
        }
        vector<int> DGM_idx;
        if (segment_analysis->IndexDGM>-1) {
            DGM_idx.push_back(segment_analysis->IndexDGM);
        }
        else if (segment_analysis->IndexGM>-1){
            DGM_idx.push_back(segment_analysis->IndexGM);
        }
        
        vector<float *> Ventr_priors;
        vector<float *> WMDGM_priors;
        vector<float *> CGM_priors;
        vector<float *> DGMVentr_priors;
        vector<float *> GMWMVentr_priors;
        vector<float *> WMDGMVentr_priors;
        vector<float *> GMWM_priors;
        vector<float *> DGM_priors;
        
        nifti_image * DGMPriors=NULL;
        nifti_image * ICSFPriors=NULL;
        nifti_image * CGMPriors=NULL;
        nifti_image * GMPriors=NULL;
        nifti_image * WMPriors=NULL;
        nifti_image * BrainstemPriors=NULL;
        
        if (segment_analysis->flag_inPriorsICSF) {
            ICSFPriors=ReadFromFilename(segment_analysis->filename_inPriorsICSF);
            float * ICSFArray=static_cast<float*>(ICSFPriors->data);
            Ventr_priors.push_back(ICSFArray);
            DGMVentr_priors.push_back(ICSFArray);
            WMDGMVentr_priors.push_back(ICSFArray);
            GMWMVentr_priors.push_back(ICSFArray);
        }
        
        if (segment_analysis->flag_inPriorsDGM) {
            DGMPriors=ReadFromFilename(segment_analysis->filename_inPriorsDGM);
            float * DGMArray=static_cast<float*>(DGMPriors->data);
            DGMVentr_priors.push_back(DGMArray);
            DGM_priors.push_back(DGMArray);
            WMDGM_priors.push_back(DGMArray);
            GMWMVentr_priors.push_back(DGMArray);
            WMDGMVentr_priors.push_back(DGMArray);
            GMWM_priors.push_back(DGMArray);
        }
        
        if (segment_analysis->flag_inPriorsCGM) {
            CGMPriors=ReadFromFilename(segment_analysis->filename_inPriorsCGM);
            float * CGMArray=static_cast<float*>(CGMPriors->data);
            CGM_priors.push_back(CGMArray);
            GMWMVentr_priors.push_back(CGMArray);
            GMWM_priors.push_back(CGMArray);
        }
        else if (segment_analysis->IndexCGM>-1){
            CGMPriors=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexCGM, L2SToUse);
            float * CGMArray=static_cast<float*>(CGMPriors->data);
            GMWMVentr_priors.push_back(CGMArray);
            GMWM_priors.push_back(CGMArray);
            
        }
        else if (segment_analysis->IndexGM>-1){
            GMPriors=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexGM, L2SToUse);
            float * GMArray=static_cast<float*>(GMPriors->data);
            GMWMVentr_priors.push_back(GMArray);
            GMWM_priors.push_back(GMArray);
        }
        if (segment_analysis->flag_inPriorsWM) {
            WMPriors=ReadFromFilename(segment_analysis->filename_inPriorsWM);
            float * WMArray=static_cast<float*>(WMPriors->data);
            WMDGM_priors.push_back(WMArray);
            GMWMVentr_priors.push_back(WMArray);
            WMDGMVentr_priors.push_back(WMArray);
        }
        else{
            WMPriors=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexWM, L2SToUse);
            float * WMArray=static_cast<float*>(WMPriors->data);
            WMDGM_priors.push_back(WMArray);
            GMWMVentr_priors.push_back(WMArray);
            WMDGMVentr_priors.push_back(WMArray);
            if (segment_analysis->IndexBrainstem>-1) {
                BrainstemPriors=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexBrainstem, L2SToUse);
                float * BrainstemArray=static_cast<float*>(BrainstemPriors->data);
                WMDGM_priors.push_back(BrainstemArray);
                GMWMVentr_priors.push_back(BrainstemArray);
                WMDGMVentr_priors.push_back(BrainstemArray);
            }
        }
        
        bool * Ventr_seg=NULL;
        bool * WMDGM_seg=NULL;
        bool * DGM_seg=NULL;
        bool * CGM_seg=NULL;
        bool * DGMVentr_seg=NULL;
        bool * WMDGMVentr_seg=NULL;
        bool * GMWMVentr_seg=NULL;
        bool * GMWM_seg=NULL;
        
        if (Ventr_priors.size()==1) {
            Ventr_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,Ventr_priors,Ventr_idx,L2SToUse,0);
            GMWMVentr_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,GMWMVentr_priors,GMWMVentr_idx,L2SToUse,0);
        }
        if (segment_analysis->IndexDGM>-1) {
            nifti_image * DGMImageSeg=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexDGM, L2SToUse);
            float * DGMData=static_cast<float*>(DGMImageSeg->data);
            DGM_seg=TranscribeArray<float, bool>(DGMData, numel);
            nifti_image_free(DGMImageSeg);
            DGMImageSeg=NULL;
        }
        else if (DGM_priors.size()==1) {
            DGM_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,DGM_priors,DGM_idx,L2SToUse,0);
        }
        if (segment_analysis->IndexCGM>-1) {
            nifti_image * CGMImageSeg=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexCGM, L2SToUse);
            float * CGMData=static_cast<float*>(CGMImageSeg->data);
            CGM_seg=TranscribeArray<float, bool>(CGMData, numel);
            nifti_image_free(CGMImageSeg);
            CGMImageSeg=NULL;
        }
        else if (CGM_priors.size()==1) {
            CGM_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,CGM_priors,CGM_idx,L2SToUse,0);
        }
        else if (DGM_seg!=NULL && segment_analysis->IndexGM >-1){
            nifti_image * GMSeg=HardSegmentationIndex(SummarisedSeg, segment_analysis->IndexGM, L2SToUse);
            float * GMData=static_cast<float *>(GMSeg->data);
            bool * CGM_seg=TranscribeArray<float, bool>(GMData, numel);
            XOROperationBool(CGM_seg, DGM_seg, CGM_seg, numel);
            nifti_image_free(GMSeg);
            GMSeg=NULL;
        }
        if (WMDGMVentr_priors.size()>=3 && DGM_seg!=NULL) {
            WMDGMVentr_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,WMDGMVentr_priors,WMDGMVentr_idx,L2SToUse,0);
        }
        if (WMDGM_priors.size()>=2 && DGM_seg!=NULL) {
            WMDGM_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,WMDGM_priors,WMDGM_idx,L2SToUse,0);
        }
        if (DGMVentr_priors.size()==2) {
            DGMVentr_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,DGMVentr_priors,DGMVentr_idx,L2SToUse,0);
        }
        if (GMWM_priors.size()>=2) {
            GMWM_seg=CreateBoolSegFromSumPriorsIdx(SummarisedSeg,GMWM_priors,GMWM_idx,L2SToUse,0);
        }
        
        //        Correction of WMDGMVentr and WMDGM for CGM inclusion if any
        if (CGM_seg!=NULL && WMDGMVentr_seg!=NULL && WMDGM_seg !=NULL && segment_analysis->IndexDGM<0) {
            for (int i=0; i<numel; i++) {
                if (L2SToUse[i]>=0) {
                    if (CGM_seg[i]) {
                        WMDGM_seg[i]=0;
                        WMDGMVentr_seg[i]=0;
                    }
                }
            }
        }
        
        nifti_image * TestWMDGMVentrNii=CreateNiiFromArray(WMDGMVentr_seg, SummarisedSeg, numel);
        nifti_set_filenames(TestWMDGMVentrNii, FilenameSaveTest.c_str(), 0, 0);
        nifti_image_write(TestWMDGMVentrNii);

        
        
        vector<int> DimVector;
        int  Shift[3];
        float PixDim[3];
        Shift[0]=1;
        for (int d=0; d<3; d++) {
            Dim[d]=SummarisedSeg->dim[d+1];
            DimVector.push_back(Dim[d]);
            PixDim[d]=SummarisedSeg->pixdim[d+1];
            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
        }
        
        
        nifti_image * LaplaceSolImageCGM=NULL;
        nifti_image * LaplaceSolImageWM=NULL;
        nifti_image * LaplaceSolImageWMD=NULL;
        nifti_image * LaplaceSolImageDGM=NULL;
        float * DistanceNormalisedLaplaceCGM=NULL;
        float * DistanceNormalisedLaplaceWMD=NULL;
        float * DistanceNormalisedLaplaceWM=NULL;
        float * DistanceNormalisedLaplaceDGM=NULL;
        
        //        First deal with CGM solution
        float * TmpGMWMVentr=TranscribeArray<bool, float>(GMWMVentr_seg, numel);
        float * TmpWMDGMVentr=TranscribeArray<bool, float>(WMDGMVentr_seg, numel);
        float * TmpDGMVentr=TranscribeArray<bool, float>(DGMVentr_seg, numel);
        float * DilatedGMWMVentr=Erosion_bis(TmpGMWMVentr, 3, DimVector, 0);
        float * DilatedWMDGMVentr=Erosion_bis(TmpWMDGMVentr, 3, DimVector, 0);
        float * DilatedDGMVentr=Erosion_bis(TmpDGMVentr, 3, DimVector, 0);
        
        delete [] TmpGMWMVentr;
        delete [] TmpWMDGMVentr;
        delete [] TmpDGMVentr;
        TmpGMWMVentr=NULL;
        TmpWMDGMVentr=NULL;
        TmpDGMVentr=NULL;

        bool * MaskDilGMWMVentr=TranscribeArray<float, bool>(DilatedGMWMVentr, numel);
        bool * MaskDilWMDGMVentr=TranscribeArray<float, bool>(DilatedWMDGMVentr, numel);
        bool * MaskDilDGMVentr=TranscribeArray<float, bool>(DilatedDGMVentr, numel);
        
        //        SaveTmpResult(ObjectIn, g_scratch_dir + "/ObjectIn.nii.gz", SegToAnalyse);
        //        SaveTmpResult(ObjectTot, g_scratch_dir + "/ObjectOut.nii.gz", SegToAnalyse);
        
        
        
        int * L2S_GMWMVentr=MakeL2S(GMWMVentr_seg, Dim);
        int numelmaskedLapGMWMVentr=0;
        int * S2L_GMWMVentr=MakeS2L(GMWMVentr_seg, Dim, numelmaskedLapGMWMVentr);
        
        int * L2S_WMDGMVentr=MakeL2S(WMDGMVentr_seg, Dim);
        int numelmaskedLapWMDGMVentr=0;
        int * S2L_WMDGMVentr=MakeS2L(WMDGMVentr_seg, Dim, numelmaskedLapWMDGMVentr);
        
        int * L2S_DGMVentr=MakeL2S(DGMVentr_seg, Dim);
        int numelmaskedLapDGMVentr=0;
        int * S2L_DGMVentr=MakeS2L(DGMVentr_seg, Dim, numelmaskedLapDGMVentr);
        
        if (LaplaceSolImageCGM==NULL) { // Meaning that there was nothing saved as Laplace solution
            int * LabelLaplaceGMWMVentr=CreateLaplaceLabelling(WMDGMVentr_seg,GMWMVentr_seg,Dim);
            //                            SaveTmpResult(LabelLaplaceGMWMVentr, g_scratch_dir + "/TestLabel.nii.gz", WMPriors);
            int numelmaskedbGMWMVentr=0;
            int * L2SbGMWMVentr=MakeL2S(MaskDilGMWMVentr, Dim);
            int * S2LbGMWMVentr=MakeS2L(MaskDilGMWMVentr, Dim, numelmaskedbGMWMVentr);
            if (MaskDilGMWMVentr!=NULL) {
                delete [] MaskDilGMWMVentr;
                MaskDilGMWMVentr=NULL;
            }
            int * Label_sGMWMVentr=CreateShort(LabelLaplaceGMWMVentr, S2LbGMWMVentr, numelmaskedbGMWMVentr);
            cout<<GetMax<int>(Label_sGMWMVentr, numelmaskedbGMWMVentr)<<endl;
            //            SaveTmpResult_short(Label_sPar, g_scratch_dir + "/TestLabelB.nii.gz", WMDGMSeg, L2SbPar);
            if (LabelLaplaceGMWMVentr!=NULL) {
                delete [] LabelLaplaceGMWMVentr;
                LabelLaplaceGMWMVentr=NULL;
            }
            
            //            nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/TestLaplace.nii.gz", 1);
            //             nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/FinalLaplace.nii.gz", 1);
            //            float * LaplaceTmpLong=static_cast<float *>(LaplaceTmp->data);
            //            float * LaplaceSol_short=CreateShort(LaplaceTmpLong, S2Lb, numelmaskedb);
            
            
            float * LaplaceSol_shortGMWMVentr=SolvingLaplaceEquation(Label_sGMWMVentr, Dim, Shift, PixDim, 5, L2SbGMWMVentr, S2LbGMWMVentr, numelmaskedbGMWMVentr,WMPriors);
            float * LaplaceSolCGM=CreateLong(LaplaceSol_shortGMWMVentr,L2SbGMWMVentr,numel);
            
            
            //            float * LaplaceSol_shortPar=NULL;
            //            float * LaplaceSolPar=static_cast<float *>(LaplaceTmp->data);
            //                        SaveTmpResult(LaplaceSolPar, g_scratch_dir + "/FinalLaplace.nii.gz", SegLesBasis);
            // Clearing memory for all that won't be needed afterwards
            delete [] L2SbGMWMVentr;
            delete [] S2LbGMWMVentr;
            delete [] Label_sGMWMVentr;
            if (LaplaceSol_shortGMWMVentr!=NULL) {
                delete [] LaplaceSol_shortGMWMVentr;
            }
            
            L2SbGMWMVentr=NULL;
            S2LbGMWMVentr=NULL;
            Label_sGMWMVentr=NULL;
            LaplaceSol_shortGMWMVentr=NULL;
            
            
            
            bool * WMDGMVentr_short=CreateShort(WMDGMVentr_seg, S2L_GMWMVentr, numelmaskedLapGMWMVentr);
            bool * Mask_shortGMWMVentr=CreateShort(GMWMVentr_seg, S2L_GMWMVentr, numelmaskedLapGMWMVentr);
            bool * Border0GMWMVentr=CreateBorderFromBool_L2S(WMDGMVentr_short, Dim, Shift, L2S_GMWMVentr, S2L_GMWMVentr, numelmaskedLapGMWMVentr);
            //            SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", SegToAnalyse, L2S);
            bool * Border1GMWMVentr=CreateBorderFromBool_L2S(Mask_shortGMWMVentr, Dim, Shift, L2S_GMWMVentr, S2L_GMWMVentr, numelmaskedLapGMWMVentr);
            float * LaplaceSol_sGMWMVentr=CreateShort(LaplaceSolCGM, S2L_GMWMVentr, numelmaskedLapGMWMVentr);
            float * DistanceNormalisedLaplace_sGMWMVentr=NormalisedLaplaceLength(LaplaceSol_sGMWMVentr, WMDGMVentr_short, Mask_shortGMWMVentr, L2S_GMWMVentr, S2L_GMWMVentr, numelmaskedLapGMWMVentr, PixDim, Dim, Shift,WMPriors);
            for (int nl=0; nl<nlsize; nl++) {
                stringstream nlName;
                nlName << segment_analysis->numbLaminae[nl];
                int * LayersNormDistCGM=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_sGMWMVentr,segment_analysis->numbLaminae[nl],numelmaskedLapGMWMVentr);
                string FilenamePA_LaplaceLayersCGM=FilenamePA_b+"LaplaceLayersCGM_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                SaveTmpResult_short(LayersNormDistCGM, FilenamePA_LaplaceLayersCGM, WMPriors, L2S_GMWMVentr);
                if (LayersNormDistCGM!=NULL) {
                    delete [] LayersNormDistCGM;
                    LayersNormDistCGM=NULL;
                }
            }

            DistanceNormalisedLaplaceCGM  =CreateLong(DistanceNormalisedLaplace_sGMWMVentr, L2S_GMWMVentr, numel);
            string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplaceCGM_"+FilenamePA_e+".nii.gz";
            nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplaceCGM, WMPriors,numel);
            nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
            nifti_image_write(NormLaplaceDistImage);
            nifti_image_free(NormLaplaceDistImage);
            NormLaplaceDistImage=NULL;
            delete [] DistanceNormalisedLaplace_sGMWMVentr;
            delete [] Mask_shortGMWMVentr;
            delete [] Border0GMWMVentr;
            delete [] Border1GMWMVentr;
            delete [] WMDGMVentr_short;
            Mask_shortGMWMVentr=NULL;
            Border0GMWMVentr=NULL;
            Border1GMWMVentr=NULL;
            WMDGMVentr_short=NULL;
            DistanceNormalisedLaplace_sGMWMVentr=NULL;
        }
        else{
            DistanceNormalisedLaplaceCGM=CopyArray(static_cast<float *>(LaplaceSolImageCGM->data), WMPriors->nvox);
        }
        
        
        if (LaplaceSolImageDGM==NULL) { // Meaning that there was nothing saved as Laplace solution
            int * LabelLaplaceDGMVentr=CreateLaplaceLabelling(Ventr_seg,DGMVentr_seg,Dim);
            //            SaveTmpResult(LabelLaplaceDGMVentr, g_scratch_dir + "/TestLabel.nii.gz", WMPriors);
            int numelmaskedbDGMVentr=0;
            int * L2SbDGMVentr=MakeL2S(MaskDilDGMVentr, Dim);
            int * S2LbDGMVentr=MakeS2L(MaskDilDGMVentr, Dim, numelmaskedbDGMVentr);
            if (MaskDilGMWMVentr!=NULL) {
                delete [] MaskDilGMWMVentr;
                MaskDilGMWMVentr=NULL;
            }
            int * Label_sDGMVentr=CreateShort(LabelLaplaceDGMVentr, S2LbDGMVentr, numelmaskedbDGMVentr);
            cout<<GetMax<int>(Label_sDGMVentr, numelmaskedbDGMVentr)<<endl;
            //            SaveTmpResult_short(Label_sPar, g_scratch_dir + "/TestLabelB.nii.gz", WMDGMSeg, L2SbPar);
            if (LabelLaplaceDGMVentr!=NULL) {
                delete [] LabelLaplaceDGMVentr;
                LabelLaplaceDGMVentr=NULL;
            }
            
            //            nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/TestLaplace.nii.gz", 1);
            //             nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/FinalLaplace.nii.gz", 1);
            //            float * LaplaceTmpLong=static_cast<float *>(LaplaceTmp->data);
            //            float * LaplaceSol_short=CreateShort(LaplaceTmpLong, S2Lb, numelmaskedb);
            
            
            float * LaplaceSol_shortDGMVentr=SolvingLaplaceEquation(Label_sDGMVentr, Dim, Shift, PixDim, 5, L2SbDGMVentr, S2LbDGMVentr, numelmaskedbDGMVentr,WMPriors);
            float * LaplaceSolDGM=CreateLong(LaplaceSol_shortDGMVentr,L2SbDGMVentr,numel);
            
            
            //            float * LaplaceSol_shortPar=NULL;
            //            float * LaplaceSolPar=static_cast<float *>(LaplaceTmp->data);
            //                        SaveTmpResult(LaplaceSolPar, g_scratch_dir + "/FinalLaplace.nii.gz", SegLesBasis);
            // Clearing memory for all that won't be needed afterwards
            delete [] L2SbDGMVentr;
            delete [] S2LbDGMVentr;
            delete [] Label_sDGMVentr;
            if (LaplaceSol_shortDGMVentr!=NULL) {
                delete [] LaplaceSol_shortDGMVentr;
            }
            
            L2SbDGMVentr=NULL;
            S2LbDGMVentr=NULL;
            Label_sDGMVentr=NULL;
            LaplaceSol_shortDGMVentr=NULL;
            
            
            
            bool * Ventr_short=CreateShort(Ventr_seg, S2L_DGMVentr, numelmaskedLapDGMVentr);
            bool * Mask_shortDGMVentr=CreateShort(DGMVentr_seg, S2L_DGMVentr, numelmaskedLapDGMVentr);
            bool * Border0DGMVentr=CreateBorderFromBool_L2S(Ventr_short, Dim, Shift, L2S_DGMVentr, S2L_DGMVentr, numelmaskedLapDGMVentr);
            //            SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", SegToAnalyse, L2S);
            bool * Border1DGMVentr=CreateBorderFromBool_L2S(Mask_shortDGMVentr, Dim, Shift, L2S_DGMVentr, S2L_DGMVentr, numelmaskedLapDGMVentr);
            float * LaplaceSol_sDGMVentr=CreateShort(LaplaceSolDGM, S2L_DGMVentr, numelmaskedLapDGMVentr);
            float * DistanceNormalisedLaplace_sDGMVentr=NormalisedLaplaceLength(LaplaceSol_sDGMVentr, Ventr_short, Mask_shortDGMVentr, L2S_DGMVentr, S2L_DGMVentr, numelmaskedLapDGMVentr, PixDim, Dim, Shift,WMPriors);
            for (int nl=0; nl<nlsize; nl++) {
                stringstream nlName;
                nlName << segment_analysis->numbLaminae[nl];
                int * LayersNormDistDGM=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_sDGMVentr,segment_analysis->numbLaminae[nl],numelmaskedLapDGMVentr);
                string FilenamePA_LaplaceLayersDGM=FilenamePA_b+"LaplaceLayersDGM_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                SaveTmpResult_short(LayersNormDistDGM, FilenamePA_LaplaceLayersDGM, WMPriors, L2S_DGMVentr);
                if (LayersNormDistDGM!=NULL) {
                    delete [] LayersNormDistDGM;
                    LayersNormDistDGM=NULL;
                }
            }
            
            DistanceNormalisedLaplaceDGM  =CreateLong(DistanceNormalisedLaplace_sDGMVentr, L2S_DGMVentr, numel);
            string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplaceDGM_"+FilenamePA_e+".nii.gz";
            nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplaceDGM, WMPriors,numel);
            nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
            nifti_image_write(NormLaplaceDistImage);
            nifti_image_free(NormLaplaceDistImage);
            NormLaplaceDistImage=NULL;
            delete [] DistanceNormalisedLaplace_sDGMVentr;
            delete [] Mask_shortDGMVentr;
            delete [] Border0DGMVentr;
            delete [] Border1DGMVentr;
            delete [] Ventr_short;
            Mask_shortDGMVentr=NULL;
            Border0DGMVentr=NULL;
            Border1DGMVentr=NULL;
            Ventr_short=NULL;
            DistanceNormalisedLaplace_sDGMVentr=NULL;
        }
        else{
            DistanceNormalisedLaplaceDGM=CopyArray(static_cast<float *>(LaplaceSolImageDGM->data), WMPriors->nvox);
        }
        
        
        
        
        if (LaplaceSolImageWM==NULL || LaplaceSolImageWMD==NULL) { // Meaning that there was nothing saved as Laplace solution
            int * LabelLaplaceWMD=CreateLaplaceLabelling(Ventr_seg,WMDGMVentr_seg,Dim);
            int * LabelLaplaceWM=CreateLaplaceLabelling(DGMVentr_seg, WMDGMVentr_seg, Dim);
            //                SaveTmpResult(LabelLaplacePar, g_scratch_dir + "/TestLabel.nii.gz", WMDGMSeg);
            int numelmaskedbWMDGMVentr=0;
            
            int * L2SbWMDGMVentr=MakeL2S(MaskDilWMDGMVentr, Dim);
            int * S2LbWMDGMVentr=MakeS2L(MaskDilWMDGMVentr, Dim, numelmaskedbWMDGMVentr);
            if (MaskDilGMWMVentr!=NULL) {
                delete [] MaskDilWMDGMVentr;
                MaskDilWMDGMVentr=NULL;
            }
            int * Label_sWMD=CreateShort(LabelLaplaceWMD, S2LbWMDGMVentr, numelmaskedbWMDGMVentr);
            int * Label_sWM=CreateShort(LabelLaplaceWM, S2LbWMDGMVentr, numelmaskedbWMDGMVentr);
            cout<<GetMax<int>(Label_sWM, numelmaskedbWMDGMVentr)<<endl;
            //            SaveTmpResult_short(Label_sPar, g_scratch_dir + "/TestLabelB.nii.gz", WMDGMSeg, L2SbPar);
            if (LabelLaplaceWM!=NULL) {
                delete [] LabelLaplaceWM;
                LabelLaplaceWM=NULL;
            }
            if (LabelLaplaceWMD!=NULL) {
                delete [] LabelLaplaceWMD;
                LabelLaplaceWMD=NULL;
            }
            //            nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/TestLaplace.nii.gz", 1);
            //             nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/FinalLaplace.nii.gz", 1);
            //            float * LaplaceTmpLong=static_cast<float *>(LaplaceTmp->data);
            //            float * LaplaceSol_short=CreateShort(LaplaceTmpLong, S2Lb, numelmaskedb);
            
            
            float * LaplaceSol_shortWM=SolvingLaplaceEquation(Label_sWM, Dim, Shift, PixDim, 5, L2SbWMDGMVentr, S2LbWMDGMVentr, numelmaskedbWMDGMVentr,WMPriors);
            float * LaplaceSolWM=CreateLong(LaplaceSol_shortWM,L2SbWMDGMVentr,numel);
            
            float * LaplaceSol_shortWMD=SolvingLaplaceEquation(Label_sWMD, Dim, Shift, PixDim, 5, L2SbWMDGMVentr, S2LbWMDGMVentr, numelmaskedbWMDGMVentr,WMPriors);
            float * LaplaceSolWMD=CreateLong(LaplaceSol_shortWMD,L2SbWMDGMVentr,numel);
            
            
            //            float * LaplaceSol_shortPar=NULL;
            //            float * LaplaceSolPar=static_cast<float *>(LaplaceTmp->data);
            //                        SaveTmpResult(LaplaceSolPar, g_scratch_dir + "/FinalLaplace.nii.gz", SegLesBasis);
            // Clearing memory for all that won't be needed afterwards
            delete [] L2SbWMDGMVentr;
            delete [] S2LbWMDGMVentr;
            delete [] Label_sWM;
            delete [] Label_sWMD;
            if (LaplaceSol_shortWM!=NULL) {
                delete [] LaplaceSol_shortWM;
            }
            if (LaplaceSol_shortWMD!=NULL) {
                delete [] LaplaceSol_shortWMD;
            }
            
            L2SbWMDGMVentr=NULL;
            S2LbWMDGMVentr=NULL;
            Label_sWMD=NULL;
            Label_sWM=NULL;
            LaplaceSol_shortWMD=NULL;
            LaplaceSol_shortWM=NULL;
            
            
            
            bool * Ventr_short=CreateShort(Ventr_seg, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            bool * DGMVentr_short=CreateShort(DGMVentr_seg, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            bool * Mask_shortWMDGMVentr=CreateShort(WMDGMVentr_seg, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            bool * Border0Ventr=CreateBorderFromBool_L2S(Ventr_short, Dim, Shift, L2S_WMDGMVentr, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            bool * Border0DGMVentr=CreateBorderFromBool_L2S(DGMVentr_short, Dim, Shift, L2S_WMDGMVentr, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            //            SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", SegToAnalyse, L2S);
            bool * Border1WMDGMVentr=CreateBorderFromBool_L2S(Mask_shortWMDGMVentr, Dim, Shift, L2S_WMDGMVentr, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            float * LaplaceSol_sWMD=CreateShort(LaplaceSolWMD, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            float * DistanceNormalisedLaplace_sWMD=NormalisedLaplaceLength(LaplaceSol_sWMD, Ventr_short, Mask_shortWMDGMVentr, L2S_WMDGMVentr, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr, PixDim, Dim, Shift,WMPriors);
            for (int nl=0; nl<nlsize; nl++) {
                int * LayersNormDistWMD=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_sWMD,segment_analysis->numbLaminae[nl],numelmaskedLapWMDGMVentr);
                stringstream nlName;
                nlName << segment_analysis->numbLaminae[nl];
                string FilenamePA_LaplaceLayersWMD=FilenamePA_b+"LaplaceLayersWMD_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                SaveTmpResult_short(LayersNormDistWMD, FilenamePA_LaplaceLayersWMD, WMPriors, L2S_WMDGMVentr);
                if (LayersNormDistWMD!=NULL) {
                    delete [] LayersNormDistWMD;
                    LayersNormDistWMD=NULL;
                }
            }
            

            DistanceNormalisedLaplaceWMD  =CreateLong(DistanceNormalisedLaplace_sWMD, L2S_WMDGMVentr, numel);
            string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplaceWMD_"+FilenamePA_e+".nii.gz";
            nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplaceWMD, WMPriors,numel);
            nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
            nifti_image_write(NormLaplaceDistImage);
            nifti_image_free(NormLaplaceDistImage);
            NormLaplaceDistImage=NULL;
            
            
            float * LaplaceSol_sWM=CreateShort(LaplaceSolWM, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr);
            float * DistanceNormalisedLaplace_sWM=NormalisedLaplaceLength(LaplaceSol_sWM, DGMVentr_short, Mask_shortWMDGMVentr, L2S_WMDGMVentr, S2L_WMDGMVentr, numelmaskedLapWMDGMVentr, PixDim, Dim, Shift,WMPriors);
            for (int nl=0; nl<nlsize; nl++) {
                stringstream nlName;
                nlName << segment_analysis->numbLaminae[nl];
                int * LayersNormDistWM=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_sWM,segment_analysis->numbLaminae[nl],numelmaskedLapWMDGMVentr);
                string FilenamePA_LaplaceLayersWM=FilenamePA_b+"LaplaceLayersWM_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                SaveTmpResult_short(LayersNormDistWM, FilenamePA_LaplaceLayersWM, WMPriors, L2S_WMDGMVentr);
                if (LayersNormDistWM!=NULL) {
                    delete [] LayersNormDistWM;
                    LayersNormDistWM=NULL;
                }
            }
            
            DistanceNormalisedLaplaceWM  =CreateLong(DistanceNormalisedLaplace_sWM, L2S_WMDGMVentr, numel);
            string FilenamePA_LaplaceWM=FilenamePA_b+"NormDistLaplaceWM_"+FilenamePA_e+".nii.gz";
            nifti_image * NormLaplaceDistImageWM=CreateNiiFromArray(DistanceNormalisedLaplaceWM, WMPriors,numel);
            nifti_set_filenames(NormLaplaceDistImageWM, FilenamePA_LaplaceWM.c_str(), 0, 0);
            nifti_image_write(NormLaplaceDistImageWM);
            nifti_image_free(NormLaplaceDistImageWM);
            NormLaplaceDistImage=NULL;
            
            
            delete [] DistanceNormalisedLaplace_sWMD;
            delete [] DistanceNormalisedLaplace_sWM;
            delete [] Mask_shortWMDGMVentr;
            delete [] Border0DGMVentr;
            delete [] Border0Ventr;
            delete [] Border1WMDGMVentr;
            delete [] Ventr_short;
            delete [] DGMVentr_short;
            Mask_shortWMDGMVentr=NULL;
            Border0Ventr=NULL;
            Border0DGMVentr=NULL;
            Border1WMDGMVentr=NULL;
            Ventr_short=NULL;
            DGMVentr_short=NULL;
            DistanceNormalisedLaplace_sWM=NULL;
            DistanceNormalisedLaplace_sWMD=NULL;
        }
        else{
            DistanceNormalisedLaplaceWM=CopyArray(static_cast<float *>(LaplaceSolImageWMD->data), WMPriors->nvox);
        }

        
        

        
        // Clearing memory for not needed elements
        delete [] L2S_GMWMVentr;
        delete [] S2L_GMWMVentr;
        delete [] DistanceNormalisedLaplaceCGM;
        
        L2S_GMWMVentr=NULL;
        S2L_GMWMVentr=NULL;
        DistanceNormalisedLaplaceCGM=NULL;
        
        delete [] L2S_WMDGMVentr;
        delete [] S2L_WMDGMVentr;
        delete [] DistanceNormalisedLaplaceWM;
        delete [] DistanceNormalisedLaplaceWMD;
        
        L2S_WMDGMVentr=NULL;
        S2L_WMDGMVentr=NULL;
        DistanceNormalisedLaplaceWMD=NULL;
        DistanceNormalisedLaplaceWM=NULL;
        
        if (SummarisedSeg!=NULL) {
            nifti_image_free(SummarisedSeg);
            SummarisedSeg=NULL;
        }
        if (CGMPriors!=NULL) {
            nifti_image_free(CGMPriors);
            CGMPriors=NULL;
        }
        if (WMPriors!=NULL) {
            nifti_image_free(WMPriors);
            WMPriors=NULL;
        }
        if (BrainstemPriors!=NULL) {
            nifti_image_free(BrainstemPriors);
            BrainstemPriors=NULL;
        }
        if (ICSFPriors!=NULL) {
            nifti_image_free(ICSFPriors);
            ICSFPriors=NULL;
        }
        if (DGMPriors!=NULL) {
            nifti_image_free(DGMPriors);
            DGMPriors=NULL;
        }
        if (WMDGMVentr_seg!=NULL) {
            delete [] WMDGMVentr_seg;
            WMDGMVentr_seg=NULL;
        }
        if (WMDGM_seg!=NULL) {
            delete [] WMDGM_seg;
            WMDGM_seg=NULL;
        }
        if (GMWM_seg!=NULL) {
            delete [] GMWM_seg;
            GMWM_seg=NULL;
        }
        if (GMWMVentr_seg!=NULL) {
            delete [] GMWMVentr_seg;
            GMWMVentr_seg=NULL;
        }
        if (CGM_seg!=NULL) {
            delete [] CGM_seg;
            CGM_seg=NULL;
        }
        if (DGMVentr_seg!=NULL) {
            delete [] DGMVentr_seg;
            DGMVentr_seg=NULL;
        }
        if (Ventr_seg!=NULL) {
            delete [] Ventr_seg;
            Ventr_seg=NULL;
        }
        return EXIT_SUCCESS;
        
        
    }
    
    
    
    
    
    
    // Testing with laplace distance or something else (when do not want to perform anything else)
    if (segment_analysis->flag_test) {
        nifti_image * TestImage=ReadFromFilename(segment_analysis->filename_test);
        // Create binary seg for specific component (currently test filename contains the name for a parcellation file)
        //        float * DataSeg=static_cast<float *>(TestImage->data);
        int numel=TestImage->nvox;
        //        bool * SegBinarisedToUse=new bool[numel];
        //        int * SignFunction=new int[numel];
        int * Dim=new int[3];
        int * Shift=new int[3];
        Shift[0]=1;
        float * PixDim=new float[3];
        for (int d=0; d<3; d++) {
            Dim[d]=TestImage->dim[d+1];
            if (d>0) {
                Shift[d]=Dim[d-1]*Shift[d-1];
            }
            PixDim[d]=TestImage->pixdim[d+1];
        }
        vector<int> ElementsToAssociateOut;
        vector<int> ElementsToAssociateIn;
        ElementsToAssociateOut.push_back(45);
        ElementsToAssociateOut.push_back(46);
        ElementsToAssociateOut.push_back(48);
        ElementsToAssociateOut.push_back(49);
        ElementsToAssociateOut.push_back(37);
        ElementsToAssociateOut.push_back(38);
        ElementsToAssociateOut.push_back(58);
        ElementsToAssociateOut.push_back(59);
        ElementsToAssociateOut.push_back(60);
        ElementsToAssociateOut.push_back(61);
        ElementsToAssociateIn.push_back(50);
        ElementsToAssociateIn.push_back(51);
        ElementsToAssociateIn.push_back(52);
        ElementsToAssociateIn.push_back(53);
        bool * ObjectOut=CreateBoolObjectFromParcellation(TestImage, ElementsToAssociateOut);
        bool * ObjectIn=CreateBoolObjectFromParcellation(TestImage,ElementsToAssociateIn);
        //        SaveTmpResult<bool>(ObjectIn, g_scratch_dir + "/TestOut.nii.gz", TestImage);
        vector<bool *> ObjectVector;
        ObjectVector.push_back(ObjectOut);
        ObjectVector.push_back(ObjectIn);
        bool * Mask=AddArray(ObjectVector, numel);
        bool * MaskDil=Dilation(Mask,3,Dim,Shift);
        //        SaveTmpResult(MaskDil, g_scratch_dir + "/TestMask.nii.gz", TestImage);
        //        SaveTmpResult(Mask, g_scratch_dir + "/TestMaxi.nii.gz", TestImage);
        int * LabelLaplace=CreateLaplaceLabelling(ObjectIn,Mask,Dim);
        //        SaveTmpResult(LabelLaplace, g_scratch_dir + "/TestLabel.nii.gz", TestImage);
        int numelmasked=0;
        int * L2S=MakeL2S(MaskDil, Dim);
        int * S2L=MakeS2L(MaskDil, Dim, numelmasked);
        int * Label_s=CreateShort(LabelLaplace, S2L, numelmasked);
        cout<<GetMax<int>(Label_s, numelmasked)<<endl;
        //        SaveTmpResult_short<int>(Label_s, g_scratch_dir + "/TestLaplaceLabel.nii.gz", TestImage,L2S);
        //      float * LaplaceSol=SolvingLaplaceEquation(Label_s, Dim, Shift, PixDim, 5, L2S, S2L, numelmasked,TestImage);
        //       SaveTmpResult_short(LaplaceSol, g_scratch_dir + "/TestLaplace.nii.gz", TestImage,L2S);
        int * SignIn=InitialisationSignFunction(ObjectIn, Dim);
        //        SaveTmpResult<int>(SignIn, g_scratch_dir + "/TestSignIn.nii.gz", TestImage);
        int * SignOut=InitialisationSignFunction(Mask, Dim);
        //        nifti_image * LaplaceSolutionImage=ReadFromFilename(g_scratch_dir + "/TestLaplace.nii.gz");
        //        int * L2S_b=MakeL2S(Mask, Dim);
        //        int numelmasked_b=0;
        //        int * S2L_b=MakeS2L(Mask, Dim, numelmasked_b);
        //        bool * ObjectIn_short=CreateShort(ObjectIn, S2L_b, numelmasked_b);
        //        bool * Mask_short=CreateShort(Mask, S2L_b, numelmasked_b);
        //        bool * Border0=CreateBorderFromBool_L2S(ObjectIn_short, Dim, Shift, L2S_b, S2L_b, numelmasked_b);
        //        SaveTmpResult_short(Border0, g_scratch_dir + "/TestNormDistLap.nii.gz", TestImage, L2S_b);
        //        bool * Border1=CreateBorderFromBool_L2S(Mask_short, Dim, Shift, L2S_b, S2L_b, numelmasked_b);
        //        SaveTmpResult_short(Border1, g_scratch_dir + "/TestNormDistLap.nii.gz", TestImage, L2S_b);
        //        float * LaplaceSolTmp=static_cast<float *>(LaplaceSolutionImage->data);
        //        float * LaplaceSol_b=CreateShort(LaplaceSolTmp, S2L_b, numelmasked_b);
        //        float * DistanceNormalisedLaplace=NormalisedLaplaceLength(LaplaceSol_b, Border0, Border1, L2S_b, S2L_b, numelmasked_b, PixDim, Dim, Shift,TestImage);
        //        SaveTmpResult_short(DistanceNormalisedLaplace, g_scratch_dir + "/TestNormDistLap.nii.gz", TestImage, L2S_b);
        //        SaveTmpResult<int>(SignOut, g_scratch_dir + "/TestSignOut.nii.gz", TestImage);
        ////        for (int i=0; i<numel; i++) {
        ////            if (DataSeg[i]==60 || DataSeg[i]==61) {
        ////                SegBinarisedToUse[i]=1;
        ////                SignFunction[i]=-1;
        ////            }
        ////            else{
        ////                SignFunction[i]=1;
        ////                SegBinarisedToUse[i]=0;
        ////            }
        ////        }
        //        float * DistanceFinalIn=InitialiseDistanceFromSignFunction(SignIn, Dim, Shift);
        ////        SaveTmpResult<float>(DistanceFinalIn, g_scratch_dir + "/TestOut.nii.gz", TestImage);
        //        float * DistanceFinalOut=InitialiseDistanceFromSignFunction(SignOut, Dim, Shift);
        //        SaveTmpResult<float>(DistanceFinalOut, g_scratch_dir + "/TestDistOut.nii.gz", TestImage);
        //        int * StatusIn=InitialiseStatusFromDistance(DistanceFinalIn, Dim,MaskDil);
        //        int * StatusOut=InitialiseStatusFromDistance(DistanceFinalOut, Dim);
        //        SaveTmpResult<int>(StatusOut, g_scratch_dir + "/TestOut.nii.gz", TestImage);
        //        InitialiseState(StatusIn, DistanceFinalIn, Dim, Shift, PixDim,MaskDil);
        //        InitialiseState(StatusOut, DistanceFinalOut, Dim, Shift, PixDim);
        //        SaveTmpResult<int>(StatusOut, g_scratch_dir + "/TestStatusOut.nii.gz", TestImage);
        //        SaveTmpResult<int>(StatusIn, g_scratch_dir + "/TestStatusIn.nii.gz", TestImage);
        //        SaveTmpResult(DistanceFinalIn, g_scratch_dir + "/TestForDistanceIn.nii.gz",TestImage);
        //        SaveTmpResult(DistanceFinalOut, g_scratch_dir + "/TestForDistanceOut.nii.gz",TestImage);
        //
        ////        vector<DLS *> DLSVector=InitialiseDLSVectorPtr(Status, DistanceFinal, Dim);
        ////        int sizeDLSInit=DLSVector.size();
        ////        HeapSort(DLSVector, sizeDLSInit-1);
        //        int Iteration=1;
        //        IterFastMarchingMap(StatusIn, DistanceFinalIn, Dim, Shift, PixDim,Iteration);
        //        cout<<"Number it for FMM "<<Iteration<<endl;
        //        Iteration=1;
        //        IterFastMarchingMap(StatusOut, DistanceFinalOut, Dim, Shift, PixDim, Iteration);
        //        cout<<"Number it for FMM "<<Iteration<<endl;
        //        SaveTmpResult(DistanceFinalIn, g_scratch_dir + "/TestDistanceFinIn.nii.gz", TestImage);
        //        SaveTmpResult(DistanceFinalOut, g_scratch_dir + "/TestDistanceFinOut.nii.gz", TestImage);
        //        SaveTmpResult<int>(StatusOut, g_scratch_dir + "/TestStatusOut.nii.gz", TestImage);
        nifti_image * DistanceImageIn=ReadFromFilename((g_scratch_dir + "/TestDistanceFinIn.nii.gz").c_str());
        nifti_image * DistanceImageOut=ReadFromFilename((g_scratch_dir + "/TestDistanceFinOut.nii.gz").c_str());
        float * PhiIn=static_cast<float *>(DistanceImageIn->data);
        float * PhiOut=static_cast<float *>(DistanceImageOut->data);
        PARAM_LS * Param_LS=new PARAM_LS();
        Param_LS->numbLaminae=10;
        Param_LS->flag_EquiVolume=0;
        Param_LS->WidthNB=6;
        Param_LS->MaxNC=4;
        Param_LS->EpsCurv=0.0001;
        Param_LS->PropNBSusp=0.8;
        float * SignArrayIn= TranscribeArray<int, float>(SignIn, numel);
        //        SaveTmpResult(SignArrayIn,g_scratch_dir + "/TestSign.nii.gz", TestImage );
        float * SignArrayOut=TranscribeArray<int, float>(SignOut, numel);
        float * SPhiIn=MultiplyElementwise(SignArrayIn, PhiIn, numel);
        float * SPhiOut=MultiplyElementwise(SignArrayOut, PhiOut,numel );
        delete [] SignArrayOut;
        delete [] SignArrayIn;
        delete [] PhiIn;
        delete [] PhiOut;
        SignArrayIn=NULL;
        SignArrayOut=NULL;
        PhiOut=NULL;
        PhiIn=NULL;
        L2S=NULL;
        S2L=NULL;
        numelmasked=0;
        
        vector<float *> PhiFinal=CompleteEquiVolumeLSSolving(ObjectIn, Mask, Dim, PixDim, Param_LS,L2S,S2L,numelmasked,SPhiIn,SPhiOut,TestImage);
        //        SaveTmpResult_short(PhiFinal, g_scratch_dir + "/TestLS.nii.gz", TestImage,L2S);
        //        SaveTmpResult(AMap, g_scratch_dir + "/TestAMap.nii.gz", DistanceImageOut);
        //        float * AMap= AMapCalculation(PhiIn, PhiOut,  Dim, Shift, PixDim,  4);
        //        SaveTmpResult(AMap, g_scratch_dir + "/TestAMap.nii.gz", DistanceImageOut);
        cout<<"distance saved"<<endl;
        
    }
    if (segment_analysis->flag_WMCard && segment_analysis->flag_SegTot) {
        string FilenameWM=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenameWM.find_last_of('/');
        string FilenameWM_b=FilenameWM.substr(0,Index+1);
        string FilenameWM_e=FilenameWM.substr(Index+1,FilenameWM.length());
        FilenameWM=FilenameWM_b+"WMCard_"+FilenameWM_e+".txt";
        segment_analysis->filename_WMCard=const_cast<char*>(FilenameWM.c_str());
        cout<<"Creating name for WMCard..."<<endl;
    }
    if(!segment_analysis->flag_outTxt && segment_analysis->flag_refLes){
        cout<<"Creating name for out result..."<<endl;
        //        if(segment_analysis->flag_TextFile){
        //            string TextFileName=segment_analysis->filename_TextFile.c_str();
        //            int Index=TextFileName.find_last_of('.');

        //            FilenamePA=nifti_makebasename(segment_analysis->filename)
        //        }
        FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        FilenamePA=FilenamePA_b+"SegAnalysis_"+FilenamePA_e+".txt";
        segment_analysis->filename_OutTxt=const_cast<char*>(FilenamePA.c_str());
        cout<<"done";
    }
    
    
    
    
    
    ofstream TxtFile(segment_analysis->filename_OutTxt);
    //    ofstream TxtFileCorr;
    //    if (segment_analysis->flag_SegTot && segment_analysis->flag_TextFile) {
    //        string FilenameSegCorr=nifti_makebasename(segment_analysis->filename_SegTot);
    //        int Index=FilenameSegCorr.find_last_of('/');
    //        string FilenameSegCorr_b=FilenameSegCorr.substr(0,Index+1);
    //        string FilenameSegCorr_e=FilenameSegCorr.substr(Index+1,FilenameSegCorr.length());
    //        FilenameSegCorr=FilenameSegCorr_b+"SegAnalysisCorr_"+FilenameSegCorr_e+".txt";
    //        ofstream TxtFileCorr(FilenameSegCorr);
    //    }

    cout<<segment_analysis->flag_PV<<" "<<segment_analysis->flag_GIFPriors<<" "<<segment_analysis->filename_GIFPriors<<endl;
    if (TreeToAnalyse!=NULL && segment_analysis->flag_PV && segment_analysis->flag_GIFPriors){
        vector<TreeEM*> SelectedPV=SelectionPVClasses(TreeToAnalyse,0,2);
        cout<<SelectedPV.size()<<endl;
        int * L2S=TreeToAnalyse->GetL2S();
        nifti_image * GIFPriors = ReadFromFilename(segment_analysis->filename_GIFPriors);
        float * DataPriorGIF=static_cast<float*>(GIFPriors->data);
        int sizeSel=SelectedPV.size();
        nifti_image * DataImage=TreeToAnalyse->GetDataImage();
        nifti_image * WMINii = HardSegmentationThresholdFromNormResp(TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp(),TreeToAnalyse,0.5);
        nifti_image * MaskNii=HardSegmentationThresholdFromNormResp(TreeToAnalyse->GetNormResp(),TreeToAnalyse,0.5);
        nifti_image * MahalImage = MahalDistMaps(WMINii,MaskNii,DataImage);
        float * MahalData = static_cast<float *>(MahalImage->data);
        string FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string NamePVNii=FilenamePA_b+"PVLes3_"+FilenamePA_e+".nii.gz";

        numel=TreeToAnalyse->GetNumberElements();
        float * GIFCSF=&DataPriorGIF[numel];
        float * GIFCGM=&DataPriorGIF[2*numel];
        float * GIFWM=&DataPriorGIF[3*numel];
        float * PotentialPVLesion=new float[numel];
        for (int i=0;i<numel;i++){
            PotentialPVLesion[i]=0;
            if (GIFCGM[i]>0){
                if (GIFWM[i]>GIFCGM[i]){
                    for(int s=0;s<sizeSel;s++){
                        float weight=MahalData[i+2*numel]/3;
                        weight=weight>1?1:weight;
                        PotentialPVLesion[i]+=weight*SelectedPV[s]->GetNormResp()[L2S[i]];
                        cout<<i;
                    }
                }
                if(GIFWM[i]>0.1 && GIFCGM[i]>GIFCSF[i]){
                    for(int s=0;s<sizeSel;s++){
                        float weight=MahalData[i+2*numel]/10;
                        weight=weight>1?1:weight;
                        PotentialPVLesion[i]+=weight*SelectedPV[s]->GetNormResp()[L2S[i]];
                        cout<<i;
                    }
                }

            }
        }
        nifti_image * PotentialPVNii=CreateNiiFromArray(PotentialPVLesion,GIFPriors,numel);
        nifti_set_filenames(PotentialPVNii,NamePVNii.c_str(),0,0);
        nifti_image_write(PotentialPVNii);
        nifti_image_free(PotentialPVNii);
        nifti_image_free(GIFPriors);
        delete TreeToAnalyse;
        return EXIT_SUCCESS;
    }

    if(TreeToAnalyse!=NULL && segment_analysis->flag_GIFPriors && segment_analysis->flag_ITCheck){
        nifti_image * GIFPriors=ReadFromFilename(segment_analysis->filename_GIFPriors);
        float * GIFData=static_cast<float *>(GIFPriors->data);
        float * GIFBrainstem=&GIFData[5*numel];
        float * NormGM = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();
        int * L2S=TreeToAnalyse->GetL2S();
        float * BrainstemCorr=new float[numel];
        for(int i=0;i<numel;i++){
            BrainstemCorr[i]=0;
            if (GIFBrainstem[i]>0.9){
                BrainstemCorr[i]=NormGM[L2S[i]];
            }
        }
        string FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        nifti_image * ITToCheck=CreateNiiFromArray(BrainstemCorr,GIFPriors,numel);
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameITToCheck=FilenamePA_b+"ITToCheck_"+FilenamePA_e+".nii.gz";
        nifti_set_filenames(ITToCheck, FilenameITToCheck.c_str(),0,0);
        nifti_image_write(ITToCheck);
        nifti_image_free(ITToCheck);
        nifti_image_free(GIFPriors);
        return EXIT_SUCCESS;
    }

    if (TreeToAnalyse!=NULL && segment_analysis->flag_WMCheck && segment_analysis->flag_GIFPriors){
        numel = TreeToAnalyse->GetNumberElements();
        nifti_image * GIFPriors=ReadFromFilename(segment_analysis->filename_GIFPriors);
        float * GIFData=static_cast<float *>(GIFPriors->data);
        float * GIFWM=&GIFData[3*numel];
        float * GIFBrainstem=&GIFData[5*numel];
        float * NormWMI = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexWM)->GetNormResp();
        nifti_image * DataImage = TreeToAnalyse->GetDataImage();
        nifti_image * WMINii = HardSegmentationThresholdFromNormResp(NormWMI,TreeToAnalyse,0.5);
        nifti_image * MahalImage = MahalDistMaps(WMINii,WMINii,DataImage);
        float * ToReconsiderWM = new float [numel];
//        float * WMIData = static_cast<float *>(WMINii->data);
        float * MahalData = static_cast<float *>(MahalImage->data);
        string FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameWMIToCheck=FilenamePA_b+"WMToCheck_"+FilenamePA_e+".nii.gz";
        int * L2S=TreeToAnalyse->GetL2S();
        for (int i=0;i<numel;i++){
            ToReconsiderWM[i]=0;
            if (MahalData[i+2*numel]>1.5 && GIFWM[i]+GIFBrainstem[i]>0.4){
                float weight=MahalData[i+2*numel]/4;
                weight=weight>1?1:weight;
                ToReconsiderWM[i]=weight*NormWMI[L2S[i]];
            }
        }
        nifti_image * WMIToCheck=CreateNiiFromArray(ToReconsiderWM,DataImage,numel);
        nifti_set_filenames(WMIToCheck, FilenameWMIToCheck.c_str(),0,0);
        nifti_image_write(WMIToCheck);
        nifti_image_free(WMIToCheck);
        nifti_image_free(MahalImage);
        nifti_image_free(WMINii);
        nifti_image_free(GIFPriors);
        return EXIT_SUCCESS;
    }
    

    if (TreeToAnalyse!=NULL && segment_analysis->flag_CSFCheck && segment_analysis->flag_GIFPriors){
        numel = TreeToAnalyse->GetNumberElements();
        nifti_image * GIFPriors=ReadFromFilename(segment_analysis->filename_GIFPriors);
        float * GIFData=static_cast<float *>(GIFPriors->data);
        float * GIFWM=&GIFData[3*numel];
        float * GIFBrainstem=&GIFData[5*numel];
        float * NormCSFI = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexCSF)->GetNormResp();
        float * ToReconsiderWM = new float [numel];
        string FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        string FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameWMIToCheck=FilenamePA_b+"CSFToCheck_"+FilenamePA_e+".nii.gz";
        int * L2S=TreeToAnalyse->GetL2S();
        for (int i=0;i<numel;i++){
            ToReconsiderWM[i]=0;
            if (GIFWM[i]+GIFBrainstem[i]>0.9){
                ToReconsiderWM[i]=NormCSFI[L2S[i]];
            }
        }
        nifti_image * WMIToCheck=CreateNiiFromArray(ToReconsiderWM,GIFPriors,numel);
        nifti_set_filenames(WMIToCheck, FilenameWMIToCheck.c_str(),0,0);
        nifti_image_write(WMIToCheck);
        nifti_image_free(WMIToCheck);
        nifti_image_free(GIFPriors);
        return EXIT_SUCCESS;
    }



    if (TreeToAnalyse!=NULL && !(segment_analysis->flag_inLes||segment_analysis->flag_inLesCorr||segment_analysis->flag_inConnect) ) {
        


        //        TreeToAnalyse->SaveTmpResultMasked(TreeToAnalyse->GetChild(0)->GetChild(1)->GetNormResp(), g_scratch_dir + "/TestReconstruct.nii.gz");
        //        TestNormResp=TreeToAnalyse->AreNormRespValid();
        
        //        TreeToAnalyse->SaveTmpResultMasked(TreeToAnalyse->GetDataBFCorrected(), g_scratch_dir + "/TestReconstruct.nii.gz");
        //        TreeToAnalyse->SaveTmpResult(static_cast<float*>(TreeToAnalyse->GetDataImage()->data), g_scratch_dir + "/TestReconstruct2.nii.gz");
        Rule * LesionRule;
        if (segment_analysis->flag_RuleTextFile) {
            LesionRule=BuildRuleFromTextFile(TreeToAnalyse, segment_analysis);
        }
        else{
            LesionRule=BuildingRule(TreeToAnalyse, segment_analysis->LesionRuleType, segment_analysis->LesionUniformRuleType, Modalities);
        }
        // Printing Lesion Rule
        //        int numbOutlierClasses=TreeToAnalyse->GetNodeOutlier()->GetNumberChildren();
        //        cout << "Acceptance Rule"<<endl;
        //        for (int i=0; i<numbOutlierClasses; i++) {
        //            cout << LesionRule->Acceptance[i]<< " ";
        //        }
        //        cout<<endl;
        //
        //        cout << "Acceptance Uniform rule"<<endl;
        //        for (int i=0; i<numbOutlierClasses; i++) {
        //            cout << LesionRule->AcceptanceUniform[i]<< " ";
        //        }
        //        cout<<endl;
        //
        //        cout << "OperationSign "<<endl;
        //        for (int m=0; m<numbmodal; m++) {
        //            for (int i=0; i<numbOutlierClasses; i++) {
        //                cout << LesionRule->SignComparison[m*numbOutlierClasses+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "OperationSignUniform "<<endl;
        //        for (int m=0; m<numbmodal; m++) {
        //            for (int i=0; i<numbOutlierClasses; i++) {
        //                cout << LesionRule->SignComparisonUniform[m*numbOutlierClasses+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "GC Class Comparison "<<endl;
        //        for (int m=0; m<numbmodal; m++) {
        //            for (int i=0; i<numbOutlierClasses; i++) {
        //                cout << LesionRule->CorrespondingGClassComparison[m*numbOutlierClasses+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "GC Class Comparison Uniform "<<endl;
        //        for (int m=0; m<numbmodal; m++) {
        //            for (int i=0; i<numbOutlierClasses; i++) {
        //                cout << LesionRule->CorrespondingGClassComparisonUniform[m*numbOutlierClasses+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "Questionable "<<endl;
        //        int numbQuestionable = LesionRule->Questionable.size();
        //        int numbQuestions=numbQuestionable/3;
        //        for (int q=0; q<numbQuestions; q++) {
        //            for (int i=0; i<3; i++) {
        //                cout << LesionRule->Questionable[q*3+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "RefinedCheck "<<endl;
        //        int numbRefinedChecks = LesionRule->RefinedCheck.size();
        //        int numbChecks=numbRefinedChecks/3;
        //        for (int q=0; q<numbChecks; q++) {
        //            for (int i=0; i<3; i++) {
        //                cout << LesionRule->RefinedCheck[q*3+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "QuestionableUniform "<<endl;
        //        int numbQuestionableUniform = LesionRule->QuestionableUniform.size();
        //        int numbQuestionsUniform=numbQuestionableUniform/3;
        //        for (int q=0; q<numbQuestionsUniform; q++) {
        //            for (int i=0; i<3; i++) {
        //                cout << LesionRule->QuestionableUniform[q*3+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //
        //        cout << "RefinedCheckUniform "<<endl;
        //        int numbRefinedChecksUniform= LesionRule->RefinedCheckUniform.size();
        //        int numbChecksUniform=numbRefinedChecksUniform/3;
        //        for (int q=0; q<numbChecksUniform; q++) {
        //            for (int i=0; i<3; i++) {
        //                cout << LesionRule->RefinedCheckUniform[q*3+i]<< " ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        //        cout <<"CheckSuspicious"<<endl;
        //        int numbSuspicious=LesionRule->CheckSuspiciousInliers.size();
        //        for (int s=0; s<numbSuspicious; s++) {
        //            cout << "Class " << LesionRule->CheckSuspiciousInliers[s][0]<<" checked against ";
        //            int numbChecks=LesionRule->CheckSuspiciousInliers[s].size()-1;
        //            for (int q=0; q<numbChecks; q++) {
        //                cout<<LesionRule->CheckSuspiciousInliers[s][q+1]<<" ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<endl;
        
        vector<TreeEM *> LesionClasses=FindLesionClasses(TreeToAnalyse, LesionRule, Modalities,segment_analysis);
        cout << "Classes found "<< LesionClasses.size()<<endl;
        int numbLesionClasses=LesionClasses.size();
        if (numbLesionClasses>0) {
            TxtFile << "NumbLesionClasses "<<numbLesionClasses<<endl;
            for (int les=0; les<numbLesionClasses; les++) {
                vector<int> HierarchyLesionClass=LesionClasses[les]->GetHierarchyVector();
                int lengthHierarchy=HierarchyLesionClass.size();
                for (int l=0; l<lengthHierarchy; l++) {
                    TxtFile<<HierarchyLesionClass[l]<<" ";
                }
                TxtFile<<endl;
            }
        }
        else{
            TxtFile<<"NumbLesionClasses 0"<<endl;
        }
        nifti_image * LesionClassesImage=ReconstructLesionImage(TreeToAnalyse, LesionRule, Modalities,segment_analysis);
                nifti_set_filenames(LesionClassesImage, (g_scratch_dir + "/StrangeRebuilt.nii.gz").c_str(), 0, 0);
                nifti_image_write(LesionClassesImage);
        nifti_image * LesionPartsFromUniform=LesionFromUniform(TreeToAnalyse, LesionRule, Modalities,LesionClasses,segment_analysis);
        nifti_image * LesionTotOut=NULL;
        nifti_image * SegSecondary=NULL;

        //        if (<#condition#>) {
        //            <#statements#>
        //        }

        if (segment_analysis->flag_simple){
            LesionTotOut=LesionVoxelwise(TreeToAnalyse,LesionRule,Modalities,segment_analysis);

        }
        else{
            LesionTotOut=LesionFromTotOutliers(TreeToAnalyse, LesionRule, Modalities, segment_analysis);

        }

 // Correction to add Ventricular lining there already

        if(segment_analysis->flag_Parc && LesionTotOut!=NULL){
            float * SegData=static_cast<float*>(LesionTotOut->data);
//            int IndexFLAIR;
            bool flag_FLAIR=0;
            std::vector<int>::iterator it;
            it = find (Modalities.begin(), Modalities.end(), 3);
            if (it!=Modalities.end()){
//                IndexFLAIR=std::distance(Modalities.begin(), it)+1;
                flag_FLAIR=1;
            }
//            else{
//                IndexFLAIR=2;
//            }
            if (! flag_FLAIR){
            float * GM_In = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();

            int * L2S = TreeToAnalyse->GetL2S();
            float * GM_In_Long = CreateLong(GM_In,L2S, numel);
            nifti_image * Parc = ReadFromFilename(segment_analysis->filename_Parc);
            float * ParcData = static_cast<float*>(Parc->data);
            for(int i=0;i<numel;i++){
                if(ParcData[i]>65 && ParcData[i]<68){
                    if(SegData[i]<0.5){
                        SegData[i]+=GM_In_Long[i];
                        SegData[i]=SegData[i]>1?1:SegData[i];
                    }
                }
                if (ParcData[i]>100 && GM_In_Long[i]>0.2) {
                    SegData[i] = 0;
                }
            }
            delete [] GM_In_Long;
            nifti_image_free(Parc);
            Parc=NULL;
            }
        }


        //        nifti_image * SegToAnalyse_1=ReconstructLesionImageTot(TreeToAnalyse,LesionRule,Modalities);
        //        SegToAnalyse=HardSegmentationLesionReconstruct(SegToAnalyse_1, TreeToAnalyse);
        switch (segment_analysis->flag_segWeighted) {
        case 0:
            SegToAnalyse=ReconstructLesionImageTot(TreeToAnalyse, LesionRule, Modalities,segment_analysis);

            break;

        default:
            if(!segment_analysis->flag_TO && !segment_analysis->flag_LesWMI){
                SegToAnalyse=ReconstructLesionImageTotWeighted(TreeToAnalyse, LesionClasses, LesionPartsFromUniform,Modalities,segment_analysis,0,CorrectionJuxta);
            }

            else{
                if (LesionTotOut!=NULL) {
                    vector<TreeEM *> LesionClassesNew;
                    FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
                    int Index=FilenamePA.find_last_of('/');
                    FilenamePA_b=FilenamePA.substr(0,Index+1);
                    if(segment_analysis->flag_outputDir){
                        FilenamePA_b=segment_analysis->name_outputDir;
                    }
                    if(segment_analysis->flag_GIFPriors){
                        nifti_image * GIFPriors=ReadFromFilename(segment_analysis->filename_GIFPriors);
                        float * GIFData=static_cast<float*>(GIFPriors->data);
                        float * GIFWM=&GIFData[3*numel];
                        float * GIFDGM=&GIFData[4*numel];
                        float * LesionTotOutData=static_cast<float*> (LesionTotOut->data);
                        nifti_image * MaskFloat=CreateNiiFromArray(static_cast<bool*>(TreeToAnalyse->GetMask()->data),LesionTotOut,numel);
                        nifti_image * InliersGMNii=HardSegmentationThresholdFromNormResp(TreeToAnalyse->GetChild(0)->GetChild(0)->GetNormResp(),TreeToAnalyse,0.5);
                        numel=TreeToAnalyse->GetNumberElements();
                        numbmodal=TreeToAnalyse->GetNumberModalities();
                        float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
                        cout <<"Data corr obtained"<<endl;
                        nifti_image * DataNii=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
                        nifti_image * MahalDist=MahalDistMaps(InliersGMNii,MaskFloat,DataNii);
                        float * MahalData=static_cast<float*>(MahalDist->data);
                        int IndexFLAIR;
                        std::vector<int>::iterator it;
                        it = find (Modalities.begin(), Modalities.end(), 3);
                        if (it!=Modalities.end()){
                            IndexFLAIR=std::distance(Modalities.begin(), it)+1;
                        }
                        else{
                            IndexFLAIR=2;
                        }
                        float ThresholdMahal=1.5;
                        vector<float *> MeanGClasses=TreeToAnalyse->GetMeanGeneralClassesVector();

                        if(MeanGClasses[0][IndexFLAIR-1]<MeanGClasses[1][IndexFLAIR-1]){
                            ThresholdMahal=2;
                        }
                        else{
                            ThresholdMahal=2;
                        }
                        float * MahalToUse=&MahalData[(IndexFLAIR)*numel];
                        float * GIFWMDGM=AddElementwise(GIFWM,GIFDGM,numel);
                        int * Dim=new int[3];
                        vector<int> DimVector;
                        int * Shift=new int[3];
                        float * PixDim=new float[3];
                        Shift[0]=1;
                        for (int d=0; d<3; d++) {
                            Dim[d]=GIFPriors->dim[d+1];
                            DimVector.push_back(Dim[d]);
                            PixDim[d]=GIFPriors->pixdim[d+1];
                            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
                        }
                        bool * GIFBool=ThresholdArray<float,bool>(GIFWMDGM,0.5,numel);
                        int CountChange=0;
                        bool * GIFBorder=CreateBorderFromBool(GIFBool,Dim,Shift);
                        for(int i=0;i<numel;i++){
                            if(GIFBool[i]){
                                if(GIFBorder[i]){
                                    if(MahalToUse[i]>2.5*ThresholdMahal){
                                        LesionTotOutData[i]=1;
                                        CountChange++;
                                    }
                                }
                                else{
                                    if(MahalToUse[i]>ThresholdMahal){
                                        LesionTotOutData[i]=1;
                                        CountChange++;
                                    }
                                }
                            }

                        }
                        int Csize=MeanGClasses.size();
                        for(int i=0;i<Csize;i++){
                            delete [] MeanGClasses[i];
                            MeanGClasses[i]=NULL;
                        }
                        cout<<" Changing due to GIF priors "<<CountChange<<endl;
                        delete [] GIFBool;
                        delete [] GIFBorder;
                        delete [] GIFWMDGM;
                        delete [] Dim;
                        delete [] Shift;
                        delete [] PixDim;
                        nifti_image_free(DataNii);
                        nifti_image_free(GIFPriors);
                        nifti_image_free(MahalDist);
                        nifti_image_free(MaskFloat);
                        nifti_image_free(InliersGMNii);
                    }
                    FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
                    string LesionTotName=FilenamePA_b+"LesionTotInitTest_"+FilenamePA_e+".nii.gz";
                    nifti_set_filenames(LesionTotOut,LesionTotName.c_str(),0,0);
                    nifti_image_write(LesionTotOut);

                    if(FinalLes==NULL){
                        SegToAnalyse=ReconstructLesionImageTotWeighted(TreeToAnalyse, LesionClassesNew, LesionTotOut,Modalities,segment_analysis,0,CorrectionJuxta);
                    }
                    else{
                        SegToAnalyse=CreateNiiFromArray(static_cast<float*>(FinalLes->data),LesionTotOut,numel);
                    }
                }
                //                nifti_set_filenames(SegToAnalyse, (g_scratch_dir + "/TestTotLesion.nii.gz").c_str(), 0, 0);
                //                nifti_image_write(SegToAnalyse);
            }

            int IndexFLAIR;
            bool flag_FLAIR=0;
            std::vector<int>::iterator it;
            it = find (Modalities.begin(), Modalities.end(), 3);
            if (it!=Modalities.end()){
                IndexFLAIR=std::distance(Modalities.begin(), it)+1;
                flag_FLAIR=1;
            }
            else{
                IndexFLAIR=2;
            }
            int CorrectVentr=0;
            if (!flag_FLAIR && (segment_analysis->flag_Parc || segment_analysis->flag_inAuthorised)){
                float * SegData = static_cast<float*>(SegToAnalyse->data);
                if (! segment_analysis->flag_inAuthorised){
                float * CSF_Out = TreeToAnalyse->GetNodeOutlier()->GetChild(segment_analysis->IndexCSF)->GetNormResp();
                float * GM_In = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();

                int * L2S = TreeToAnalyse->GetL2S();
                float *CSF_Out_Long = CreateLong(CSF_Out, L2S,numel);
                float * GM_In_Long = CreateLong(GM_In,L2S, numel);
                nifti_image * Parc = ReadFromFilename(segment_analysis->filename_Parc);
                float * ParcData = static_cast<float*>(Parc->data);
                for(int i=0;i<numel;i++){
                    if (ParcData[i]>46 && ParcData[i]<54){
                        if (CSF_Out_Long[i]>0.5){
                            SegData[i]=0;
                            CorrectVentr++;
                        }
                    }
                    if(ParcData[i]>65 && ParcData[i]<68){
                        if(SegData[i]<0.5){
                            SegData[i]+=GM_In_Long[i];
                            SegData[i]=SegData[i]>1?1:SegData[i];
                        }
                    }
                }

                delete [] GM_In_Long;
                delete [] CSF_Out_Long;
                nifti_image_free(Parc);
                Parc=NULL;
                }
                else{
                    cout<<"Using authorised region"<<endl;
                    nifti_image * Autho = ReadFromFilename(segment_analysis->filename_inAuthorised);
                    float * AuthoData = static_cast<float*>(Autho->data);
                    cout << "Autho read"<<endl;
                    for(int i=0;i<numel;i++){
                        if (AuthoData[i]<0.01 && SegData[i]>0.5 ){
                                SegData[i]=0;
                                CorrectVentr++;

                        }

                }
                    nifti_image_free(Autho);
                    Autho=NULL;
                if(segment_analysis->flag_Parc){
                    float * GM_In = TreeToAnalyse->GetNodeInlier()->GetChild(segment_analysis->IndexGM)->GetNormResp();

                    int * L2S = TreeToAnalyse->GetL2S();
                    float * GM_In_Long = CreateLong(GM_In,L2S, numel);
                    nifti_image * Parc = ReadFromFilename(segment_analysis->filename_Parc);
                    float * ParcData = static_cast<float*>(Parc->data);
                    for(int i=0;i<numel;i++){
                        if(ParcData[i]>65 && ParcData[i]<68){
                            if(SegData[i]<0.5){
                                SegData[i]+=GM_In_Long[i];
                                SegData[i]=SegData[i]>1?1:SegData[i];
                            }
                        }
                    }
                    delete [] GM_In_Long;
                    nifti_image_free(Parc);
                    Parc=NULL;
                }
                }
            }
            cout << "Correction for ventricle is "<< CorrectVentr<<endl;
            FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
            int Index=FilenamePA.find_last_of('/');
            string FilenamePA_b=FilenamePA.substr(0,Index+1);
            if(segment_analysis->flag_outputDir){
                FilenamePA_b=segment_analysis->name_outputDir;
            }
            string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
            string FilenamePA_CheckInit=FilenamePA_b+"CheckInit_"+FilenamePA_e+".nii.gz";
            SaveTmpResult(static_cast<float*>(SegToAnalyse->data),FilenamePA_CheckInit.c_str(),SegToAnalyse);
            string FilenamePA_CheckInit2=FilenamePA_b+"CheckInit2_"+FilenamePA_e+".nii.gz";
            SaveTmpResult(static_cast<float*>(LesionTotOut->data),FilenamePA_CheckInit2.c_str(),SegToAnalyse);
            vector<TreeEM*> LesionClassesNew;
//            SegSecondary=ReconstructLesionImageTotWeighted(TreeToAnalyse, LesionClasses, LesionPartsFromUniform,Modalities,segment_analysis,1,CorrectionJuxta);
            SegSecondary=ReconstructLesionImageTotWeighted(TreeToAnalyse, LesionClassesNew, LesionTotOut,Modalities,segment_analysis,1,CorrectionJuxta);

            break;
        }
        
        //        SegToAnalyse=ReconstructLesionImageTot(TreeToAnalyse,LesionRule,Modalities);
        
        //        nifti_set_filenames(SegToAnalyse_1, (g_scratch_dir + "/CompSegMax.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(SegToAnalyse_1);
        //        nifti_image_free(SegToAnalyse_1);
        
        //        float * MeanLesionToAnalyse=GetMeanLesion(SegToAnalyse, TreeToAnalyse);
        //        TxtFile << "MeanLesionFin ";
        //        for (int m=0; m<numbmodal; m++) {
        //            TxtFile<<MeanLesionToAnalyse[m];
        //        }
        //        TxtFile<<endl;
        nifti_image * SeverityMa=SeverityMap(TreeToAnalyse, SegToAnalyse, Modalities, segment_analysis);
        nifti_image * LesionSuspiciousImage=SuspiciousImage(LesionRule, TreeToAnalyse, Modalities,  SegToAnalyse);
        FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        int Index=FilenamePA.find_last_of('/');
        stringstream sG;
        sG << segment_analysis->LesionRuleType;
//        string cG=sG.str();
        stringstream sU;

        sU << segment_analysis->LesionUniformRuleType;
        string Options="";
        if (segment_analysis->flag_segWeighted>=1) {
            Options+="WS";
            stringstream ws;
            ws << segment_analysis->weightThreshold;
            Options+=ws.str();
            stringstream wt;
            wt <<segment_analysis->flag_segWeighted;
            Options+="WT"+wt.str();
            stringstream wc;
            wc << segment_analysis->weightCompClass;
            Options+="WC"+wc.str();
        }
        stringstream st;
        stringstream CorrLev;
        CorrLev << segment_analysis->flag_correctionLevel;
        stringstream CorrIV;
        CorrIV << segment_analysis->flag_CorrectIV;
        st << segment_analysis->flag_segType;
        string OptionAdded;
        if(segment_analysis->flag_inOptionText){
            OptionAdded=segment_analysis->inOptionText;
        }
        stringstream sFWMI;
        sFWMI<<segment_analysis->flag_LesWMI;
        if(segment_analysis->flag_LesWMI){
            Options+="WMI";
            Options+=sFWMI.str();
        }
        if(segment_analysis->flag_VentrSeg==1){
            Options+="Ventr1";
        }
        if(segment_analysis->flag_Parc){
            Options+="Parc1";
        }
        if(segment_analysis->flag_juxtaCorrection){
            Options+="JC1";
        }
        if(segment_analysis->flag_inArtefact){
            cout<<"AC to put in the options"<<endl;
            Options+="AC1";
        }
        if(segment_analysis->flag_SPCorrection){
            cout<<"SP to put in the options"<<endl;
            Options+="SP1";
        }
        if(segment_analysis->flag_oldLesion){
            Options+="OL1";
        }
        Options+="ST"+st.str()+"CL"+CorrLev.str()+"CIV"+CorrIV.str()+OptionAdded+"_";
        if (segment_analysis->flag_TO) {
            Options+="TO";
        }
        if (segment_analysis->flag_checkInliers) {
            Options+="CI";
        }
        int sizeLA=segment_analysis->vecLeavesToAdd.size();
        if (sizeLA>0) {
            Options+="LA";
        }
        string cG=sG.str();
        string cU=sU.str();
        cG="";
        cU="";
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenamePA_1=FilenamePA_b+"LesionG_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        cout<<FilenamePA_1<<endl;
        string FilenamePA_2=FilenamePA_b+"LesionU_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_3=FilenamePA_b+"LesionTTestb_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_3b=FilenamePA_b+"LesionTb_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_3c=FilenamePA_b+"LesionTc_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_3d=FilenamePA_b+"LesionTd_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_3e=FilenamePA_b+"LesionTe_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_4=FilenamePA_b+"LesionS_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_5=FilenamePA_b+"Severity_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_6=FilenamePA_b+"Connect_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_7=FilenamePA_b+"ConnectR_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_8=FilenamePA_b+"ConnectO_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_9=FilenamePA_b+"LesionCard_"+Options+cG+cU+FilenamePA_e+".txt";
        string FilenamePA_9b=FilenamePA_b+"LesionCard_"+Options+cG+cU+FilenamePA_e+"-b.txt";
        string FilenamePA_9c=FilenamePA_b+"OutlierWMHCard_"+Options+cG+cU+FilenamePA_e+".txt";
        string FilenamePA_9d=FilenamePA_b+"OutlierWMHCard_"+Options+cG+cU+FilenamePA_e+"-b.txt";
        string FilenamePA_9e=FilenamePA_b+"OutlierWMHCardSec_"+Options+cG+cU+FilenamePA_e+".txt";
        string FilenamePA_10=FilenamePA_b+"BorderCSF_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_11=FilenamePA_b+"BorderCSFLabel_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_12=FilenamePA_b+"BorderGMC_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_13=FilenamePA_b+"BorderGMCLabel_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_14=FilenamePA_b+"BorderCSO_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_15=FilenamePA_b+"BorderCSOLabel_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_16=FilenamePA_b+"BorderWMI_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_17=FilenamePA_b+"BorderWMILabel_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_18=FilenamePA_b+"LesionClassif_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_18b=FilenamePA_b+"OutlierWMHClassif_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_19=FilenamePA_b+"LesSegHard_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_20=FilenamePA_b+"TrueLesSeg_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_21=FilenamePA_b+"SummarisedComplete_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_22=FilenamePA_b+"SummarisedCorrected_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_23=FilenamePA_b+"LesionCorrected_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_23Pre=FilenamePA_b+"LesionCorrectedPre_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_24=FilenamePA_b+"LaplaceLayers_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_25=FilenamePA_b+"LesionLabelCorrected_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_26=FilenamePA_b+"LeafCaract_"+Options+cG+cU+FilenamePA_e+".txt";
        string FilenamePA_27=FilenamePA_b+"BinaryNIV_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_28=FilenamePA_b+"Secondary_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_29=FilenamePA_b+"SecondaryCorrected_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenameInfarct = FilenamePA_b+"InfarctDetection_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_ConnectSec=FilenamePA_b+"ConnectSecondary_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_Check1=FilenamePA_b+"Check1_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_Check2=FilenamePA_b+"Check2_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_Check3=FilenamePA_b+"Check3_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        string FilenamePA_Check4=FilenamePA_b+"Check4_"+Options+cG+cU+FilenamePA_e+".nii.gz";
        if (SegSecondary !=NULL){
            nifti_set_filenames(SegSecondary,FilenamePA_28.c_str(),0,0);
            nifti_image_write(SegSecondary);
        }

        if (segment_analysis->flag_Saving) {
            nifti_set_filenames(SeverityMa, FilenamePA_5.c_str(), 0, 0);
            nifti_image_write(SeverityMa);
            nifti_set_filenames(LesionClassesImage, FilenamePA_1.c_str(), 0, 0);
            nifti_image_write(LesionClassesImage);
            if(LesionPartsFromUniform!=NULL){
                nifti_set_filenames(LesionPartsFromUniform, FilenamePA_2.c_str(), 0, 0);
                nifti_image_write(LesionPartsFromUniform);
                //            nifti_image_free(LesionPartsFromUniform);
                //            LesionPartsFromUniform=NULL;
            }
        }
        if(SeverityMa!=NULL){
            nifti_image_free(SeverityMa);
            SeverityMa=NULL;
        }
        if(LesionClassesImage!=NULL){
        nifti_image_free(LesionClassesImage);
        LesionClassesImage=NULL;
        }

        if(LesionSuspiciousImage==NULL){
            cout<<"no suspicious classes"<<endl;
        }
        if(LesionSuspiciousImage!=NULL){
            cout<<"suspicious classes exist"<<endl;
            if (segment_analysis->flag_Saving) {
                nifti_set_filenames(LesionSuspiciousImage, FilenamePA_4.c_str(), 0, 0);
                nifti_image_write(LesionSuspiciousImage);
            }

            nifti_image_free(LesionSuspiciousImage);
            LesionSuspiciousImage=NULL;
        }
        if (segment_analysis->flag_Saving) {
            nifti_set_filenames(SegToAnalyse, FilenamePA_3.c_str(), 0, 0);
            nifti_image_write(SegToAnalyse);
        }
        cout << "Writing SegToAnalyse"<<endl;
        nifti_set_filenames(SegToAnalyse, FilenamePA_3.c_str(), 0, 0);
        nifti_image_write(SegToAnalyse);
        
        nifti_image * RefinedConnectLabelLesions=NULL;
        nifti_image * ImageOrderedLesion=NULL;
        nifti_image * ImageLabelCorr=NULL;
        
        vector<TreeEM *> CSFWMClasses=FindCSFWMOutliers(TreeToAnalyse, Modalities, segment_analysis);
        nifti_image * CSFWMUniform=FindCSFWMOutliersUniform(TreeToAnalyse, Modalities, segment_analysis,CorrectionJuxta);
        
        cout<<"Performing summarised"<<endl;
        if (!segment_analysis->flag_LesWMI){
            SummarisedSeg1=SummarizedSegmentation(TreeToAnalyse, LesionClasses, LesionPartsFromUniform,CSFWMClasses,CSFWMUniform, segment_analysis,ICSF);
        }
        else{
            SummarisedSeg1=SummarizedSegmentation(TreeToAnalyse, LesionClasses, FinalLes,CSFWMClasses,CSFWMUniform, segment_analysis,ICSF);

        }
        cout << "Summarised image obtained "<< endl;
        if (segment_analysis->flag_inPriorsICSF && segment_analysis->flag_CorrectIV) {
            cout << "Reading ICSF priors" << endl;
            nifti_image * PriorsICSF=ReadFromFilename(segment_analysis->filename_inPriorsICSF);
            ICSF=CopyArray(static_cast<float*>(PriorsICSF->data), numel);
            nifti_image_free(PriorsICSF);
        }
        
        if (segment_analysis->flag_Parc){
            cout << "Reading Parcellation"<< endl;
            nifti_image * ParcImage=ReadFromFilename(segment_analysis->filename_Parc);
            int Dim[3];
            int Shift[3];
            for(int d=0;d<3;d++){
                Dim[d]=ParcImage->dim[d+1];
                Shift[d]=d>0?Dim[d-1]*Shift[d-1]:1;
            }
            float * ParcData=static_cast<float*>(ParcImage->data);
            if (numel!=ParcImage->nvox){
                cout<<"Pb with numel !!!"<<endl;
                numel=ParcImage->nvox;
            }
            bool * ParcECSF=UpperThresholdArray<float, bool>(ParcData, 1.5, numel);
            nifti_image * ParcECSFNii=CreateNiiFromArray(ParcECSF, ParcImage, numel);
            float * ParcBrainstem=UpperThresholdArray<float, float>(ParcData, 35.5, numel);
            bool * ParcFinBrainstem=ThresholdArray<float, bool>(ParcBrainstem, 34.5, numel);
            cout << "Brainstem is "<< CountNonZero(ParcFinBrainstem,numel)<<" and numel is "<< numel<<endl;
            delete [] ParcBrainstem;
            ParcBrainstem=NULL;
            bool * BorderBrainstem=CreateBorderFromBool(ParcFinBrainstem, Dim, Shift);
            int * CompLabelECSF=ComponentLabeling(ParcECSFNii, 26);
            bool * MaskData=static_cast<bool*>(TreeToAnalyse->GetMask()->data);
            bool * OpposeMask=OpposeBoolArray(MaskData, numel);
            bool * FullMaskData=new bool[numel];
            OROperationBool(MaskData, OpposeMask, FullMaskData, numel);
            nifti_image * FullMask=CreateNiiFromArray(FullMaskData, ParcImage, numel);
            
                nifti_image * DistanceMask=EuclideanDistanceImage(ParcImage, OpposeMask, FullMask);
            bool * ECSFToSupress=new bool[numel];
            for(int i=0;i<numel;i++){
                ECSFToSupress[i]=0;
            }
            int numbLabelsECSF=GetMaxLabel(CompLabelECSF, numel);
            float * DistanceMaskData=static_cast<float*>(DistanceMask->data);
            for(int l=0;l<numbLabelsECSF;l++){
                bool * ECSFComp=CreateLesionBool(CompLabelECSF, l+1, numel);
                float minDistance=GetMin(DistanceMaskData, ECSFComp, numel);
                if(minDistance<3){
                    OROperationBool(ECSFComp, ECSFToSupress, ECSFToSupress, numel);
                }
                delete [] ECSFComp;
                ECSFComp=NULL;
            }
            int CountDelete=0;
            int CountDeleteBS=0;
            //            Then for all valid elements in ECSFToSupress, Suppress it effectively in SegToAnalyse.
            float * SegToAnalyseData=static_cast<float*>(SegToAnalyse->data);
            for(int i=0;i<numel;i++){

                if(ECSFToSupress[i]==1){
                    if(SegToAnalyseData[i]>0){
                        SegToAnalyseData[i]=0;
                        CountDelete++;
                    }
                }
                if(BorderBrainstem[i]==1){
                    if (SegToAnalyseData[i]>0){
                        SegToAnalyseData[i]=0;
                        CountDeleteBS++;
                    }
                }
            }
            cout << "SegToAnalyse after ECSF Brainstem corr "<< GetSum(SegToAnalyseData, numel)<<endl;
            delete [] BorderBrainstem;
            BorderBrainstem=NULL;
            cout<<"Count delete brainstem "<<CountDeleteBS<<endl;
            cout << "Corrected ECSF is "<<CountDelete<<endl;
            delete [] ECSFToSupress;
            nifti_image_free(ParcImage);
            nifti_image_free(DistanceMask);
            nifti_image_free(ParcECSFNii);
            nifti_image_free(FullMask);
            delete [] CompLabelECSF;
            delete [] OpposeMask;
        }
        nifti_set_filenames(SegToAnalyse, FilenamePA_3b.c_str(), 0, 0);
        //        nifti_image_write(SegToAnalyse);
        std::vector<int>::iterator it;
        it = find (Modalities.begin(), Modalities.end(), 3);
        if (it == Modalities.end()){ // if FLAIR in modalities
            //        If needed Perform the correction with respect to GM border
            bool * ToSuppressGMBorder=ExternalGMBorderCorrection(SegToAnalyse, TreeToAnalyse, segment_analysis);
            CorrectingSegmentation(SegToAnalyse,ToSuppressGMBorder);
            CorrectingSegmentation(SegSecondary,ToSuppressGMBorder);
            nifti_set_filenames(SegToAnalyse, FilenamePA_3c.c_str(), 0, 0);
            //            nifti_image_write(SegToAnalyse);
            bool * ToSuppressWMGM= ExternalGMBorderWM(SegToAnalyse,TreeToAnalyse, segment_analysis);
            CorrectingSegmentation(SegToAnalyse, ToSuppressWMGM);
            CorrectingSegmentation(SegSecondary, ToSuppressWMGM);
            nifti_set_filenames(SegToAnalyse, FilenamePA_3d.c_str(), 0, 0);
            //            nifti_image_write(SegToAnalyse);

            // If needed Perform the correction with respect to GM and ventricle border
            bool * ToSuppressVentricleBorder=VentricleGMBorderCorrection(SegToAnalyse, TreeToAnalyse, segment_analysis);
            CorrectingSegmentation(SegToAnalyse, ToSuppressVentricleBorder);
            CorrectingSegmentation(SegSecondary, ToSuppressVentricleBorder);

        }
        else{
            cout << "FLAIR in segmentation no correction done"<<endl;
        }
        
        SaveTmpResult(static_cast<float*>(SegToAnalyse->data),FilenamePA_Check1.c_str(),SegToAnalyse);
        // Creation of a remaining copy NIV (Non in ventricle corrected)
        SegToAnalyseNIV=CopyFloatNii(SegToAnalyse);
        //        SegToAnalyseNIV=nifti_copy_nim_info(SegToAnalyse);
        //        SegToAnalyseNIV->data=(void*)calloc(SegToAnalyseNIV->nvox, sizeof(float));
        //        float * SegNIVData=static_cast<float*>(SegToAnalyseNIV->data);
        //        float * SegAnalyseData=static_cast<float*>(SegToAnalyse->data);
        //        for (int i=0; i<numel; i++) {
        //            SegNIVData[i]=SegAnalyseData[i];
        //        }
        
        if (SegToAnalyse!=NULL && ICSF!=NULL && segment_analysis->flag_CorrectIV) {
            float * SegToAnalyseData=static_cast<float *>(SegToAnalyse->data);
            for (int i=0; i<numel; i++) {
                if (ICSF[i]>0.45) {
                    SegToAnalyseData[i]=0;
                }
            }
            cout << "After CIV les sum is "<< GetSum(SegToAnalyseData, numel)<<endl;
        }


        //        if (LesionPartsFromUniform!=NULL) {
        //            nifti_image_free(LesionPartsFromUniform);
        //            LesionPartsFromUniform=NULL;
        //        }
        cout<<"Summarised done"<<endl;
        nifti_image * SummarisedSegDec=SummarisedSegDecoupled(SummarisedSeg1, SegToAnalyse, segment_analysis->IndexWM);
        cout<<"Decoupling done"<<endl;
        nifti_image * BinarySegNIV=NULL;
        if(segment_analysis->flag_LesWMI){
            segment_analysis->flag_segType=1;
        }
        //        nifti_image * SummarisedSegIVCorr=CorrectionIV(SummarisedSegDec,TreeToAnalyse,ICSF,segment_analysis)
        switch (segment_analysis->flag_segType) {
        case 1: {// Case where we use a threshold at 0.5
            if(segment_analysis->flag_GIFPriors){
                BinarySeg=HardSegmentationThreshold(SegToAnalyse, 0.75);
                BinarySegNIV=HardSegmentationThreshold(SegToAnalyseNIV,0.75);
            }
            else{
                BinarySeg=HardSegmentationThreshold(SegToAnalyse, 0.5);
                BinarySegNIV=HardSegmentationThreshold(SegToAnalyseNIV,0.5);
            }
        }
            break;
        case 2: {// Case where we use the Initial reconstruct with only comparison to inliers
            BinarySeg=HardSegmentationLesionReconstruct(SegToAnalyse, TreeToAnalyse);
            BinarySegNIV=HardSegmentationLesionReconstruct(SegToAnalyseNIV, TreeToAnalyse);
        }
            break;
        case 3: {// Case where we use the hard segmentation with the temporary summarised result
            BinarySeg=HardSegmentationLesionReconstruct_bis(SegToAnalyse, SummarisedSeg1, TreeToAnalyse, segment_analysis);
            BinarySegNIV=HardSegmentationLesionReconstruct_bis(SegToAnalyseNIV, SummarisedSeg1, TreeToAnalyse, segment_analysis);
        }
            break;
        default: {// Using the summarised second solution as default
            BinarySeg=HardSegmentationLesionReconstruct_bis(SegToAnalyse, SummarisedSeg1, TreeToAnalyse, segment_analysis);
            BinarySegNIV=HardSegmentationLesionReconstruct_bis(SegToAnalyseNIV, SummarisedSeg1, TreeToAnalyse, segment_analysis);
        }
            break;
        }
        cout <<"Hard segmentation lesion done..."<<endl;
        SaveTmpResult(static_cast<float*>(BinarySeg->data),FilenamePA_Check2.c_str(),SegToAnalyse);

        
        //        nifti_set_filenames(BinarySeg, (g_scratch_dir + "/StrangeBinary.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(BinarySeg);

        if (segment_analysis->flag_Saving) {
            nifti_set_filenames(SummarisedSeg1, FilenamePA_21.c_str(), 0, 0);
            nifti_image_write(SummarisedSeg1);
            
            nifti_set_filenames(BinarySegNIV, FilenamePA_27.c_str(), 0, 0);
            nifti_image_write(BinarySegNIV);
            cout<<"Summarised saved"<<endl;
        }
        
        nifti_set_filenames(BinarySeg, FilenamePA_3e.c_str(), 0, 0);
        //        nifti_image_write(BinarySeg);

        
        // MAYBE SAVE THE NIV VERSION AS WELL
        
        
        //        Begin here to perform the DSC TN TP FN FP
        
        
        
        
        //        nifti_image * DGMExtraction=ExtractWMDGM(SummarisedSeg1, segment_analysis, TreeToAnalyse);
        if (segment_analysis->flag_inVentricleSeg) {
            VentricleSeg=ReadFromFilename(segment_analysis->filename_inVentricleSeg);
        }
        else{
            cout<<"Doing Ventricle segmentation ... ";
            if (segment_analysis->flag_inPriorsICSF) {
                PriorsICSF=ReadFromFilename(segment_analysis->filename_inPriorsICSF);
                if (ICSF!=NULL) {
                    delete [] ICSF;
                    ICSF=NULL;
                }
                ICSF=CopyArray(static_cast<float *>(PriorsICSF->data), numel);
                VentricleSeg=VentricleSegmentationPriors(SummarisedSeg1,segment_analysis,TreeToAnalyse,PriorsICSF);
                nifti_image_free(PriorsICSF);
                PriorsICSF=NULL;
            }
            else{
                VentricleSeg=VentricleSegmentation(SummarisedSeg1, segment_analysis, TreeToAnalyse);
            }
            if (segment_analysis->flag_Saving) {
                cout << "Saving Ventricle seg "<<endl;
                cout << VentricleSeg<<endl;
                cout << FilenameVentricle<< endl;
                nifti_set_filenames(VentricleSeg, FilenameVentricle.c_str(), 0, 0);
                nifti_image_write(VentricleSeg);
            }

            cout<<"... done "<<endl;
        }
        cout << VentricleSeg << endl;
        float * VentricleSegData=static_cast<float *>(VentricleSeg->data);
        cout<<VentricleSegData<<endl;
        VentricleBool=TranscribeArray<float, bool>(VentricleSegData, numel);
        cout<<"Ventricle Transcribed "<<endl;
        nifti_image * HardCSF=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexCSF, TreeToAnalyse->GetL2S());
        cout<<"HardCSF created "<<endl;
        float * HardCSFData=static_cast<float *>(HardCSF->data);
        bool * CSFBool=TranscribeArray<float, bool>(HardCSFData, numel);
        if (HardCSF!=NULL) {
            nifti_image_free(HardCSF);
            HardCSF=NULL;
        }
        ECSFBool=new bool[numel];
        XOROperationBool(CSFBool, VentricleBool, ECSFBool, numel);
        delete [] CSFBool;
        CSFBool=NULL;
        cout<<"ECSFBool created"<<endl;
        vector <int> DimVector;
        int Dim[3];
        int Shift[3];
        float PixDim[3];
        for (int d=0; d<3; d++) {
            DimVector.push_back(VentricleSeg->dim[d+1]);
            Dim[d]=DimVector[d];
            PixDim[d]=VentricleSeg->pixdim[d+1];
        }
        Shift[0]=1;
        Shift[1]=Dim[0];
        Shift[2]=Dim[1]*Shift[1];
        
        
        cout<<"Doing DGM region acceptance ... "<<endl;
        if (segment_analysis->flag_inPriorsDGM) {
            nifti_image * PriorsDGM=ReadFromFilename(segment_analysis->filename_inPriorsDGM);
            float * DGMData=static_cast<float *>(PriorsDGM->data);
            DGM=CopyArray(DGMData, numel);
            DGMBool=TranscribeArray<float, bool>(DGMData, numel);
        }
        else{ // If there is no proper definition of DGM we use the ventricular segmentation and extend it by 2 voxels
            int * VentriclesBoundingIndices=FindBoundingBox(VentricleBool, Dim, Shift);
            bool * BoundingBoxVentricles=CreateBoolBoundingBox(VentriclesBoundingIndices, Dim);
            DGMBool=ErosionTemplate(BoundingBoxVentricles, 5, DimVector, 0);
            
        }
        
        cout<<"Doing CGM region acceptance ... "<<endl;
        if (segment_analysis->flag_inPriorsCGM) {
            nifti_image * PriorsCGM=ReadFromFilename(segment_analysis->filename_inPriorsCGM);
            CGMBool=ThresholdArray<float, bool>(static_cast<float*>(PriorsCGM->data), 0.5, numel);
            nifti_image_free(PriorsCGM);
            PriorsCGM=NULL;
        }
        cout<<"flag ECSF "<<segment_analysis->flag_inPriorsECSF<<endl;
        cout<<"Doing SP Potential zone ..."<<endl;
        if (segment_analysis->flag_MNITemplate & segment_analysis->flag_MNITransform) {
            nifti_image * QuadrantResult=CreateQuadrantResult(VentricleSeg, segment_analysis, TreeToAnalyse);
            SPRegion=PotentialSPRegion(QuadrantResult, VentricleBool, CGMBool, segment_analysis->ThresholdSP, TreeToAnalyse->GetMask());
        }

        bool * AcceptedGMBool=NULL;
        if (segment_analysis->flag_AcceptedGM) {
            //        Study of Components of GM
            nifti_image * HardGM=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexGM, TreeToAnalyse->GetL2S());
            float * GMSummarised=static_cast<float *>(HardGM->data);
            int * GMLabels=ComponentLabeling(GMSummarised, 6, Dim, Shift);
            //        int * WMLabels=ComponentLabeling(HardWM, 26);
            float VolumeVox=HardGM->pixdim[1]*HardGM->pixdim[2]*HardGM->pixdim[3];
            int * OrderedGMLabels=OrderedVolumeLabel(GMLabels, 1, numel,VolumeVox);
            int maxGMLabel=GetMaxLabel(OrderedGMLabels, numel);
            
            nifti_image * ComponentsGMLabels=CreateNiiFromArray(OrderedGMLabels, HardGM,numel);
            //            nifti_set_filenames(ComponentsGMLabels, (g_scratch_dir + "/GMComponents.nii.gz").c_str(), 0, 0);
            //            nifti_image_write(ComponentsGMLabels);
            
            delete [] GMLabels;
            GMLabels=NULL;
            nifti_image_free(ComponentsGMLabels);
            ComponentsGMLabels=NULL;

            float SumPixDim=PixDim[0]+PixDim[1]+PixDim[2];
            int IndexWM=segment_analysis->IndexWM;
            int IndexCSF=segment_analysis->IndexCSF;
            for (int l=2; l<=maxGMLabel; l++) {
                bool * GMLabelBool=CreateLesionBool(OrderedGMLabels, l, numel);
                float DistanceECSF=GetDistanceBetweenSeg(ECSFBool, GMLabelBool, Dim, Shift, PixDim);
                float DistanceVentricle=GetDistanceBetweenSeg(VentricleBool, GMLabelBool, Dim, Shift, PixDim);
                float * ProportionGMLabel=ProportionNeighboursLesion(OrderedGMLabels, l, SummarisedSeg1, segment_analysis);
                
                if (DistanceECSF>SumPixDim) {
                    if (DistanceVentricle<SumPixDim) {
                        OROperationBool(GMLabelBool, AcceptedGMBool, AcceptedGMBool, numel);
                    }
                    
                    else if (ProportionGMLabel[IndexWM]>ProportionGMLabel[IndexCSF]) {
                        OROperationBool(GMLabelBool, AcceptedGMBool, AcceptedGMBool, numel);
                    }
                }
                delete [] OrderedGMLabels;
                OrderedGMLabels=NULL;
                delete [] GMLabelBool;
                GMLabelBool=NULL;
                delete [] ProportionGMLabel;
                ProportionGMLabel=NULL;
                if(HardGM!=NULL){
                    nifti_image_free(HardGM);
                    HardGM=NULL;
                }
            }
            
            nifti_image * AcceptedGMImage=CreateNiiFromArray(AcceptedGMBool, VentricleSeg,numel);
            //            nifti_set_filenames(AcceptedGMImage, g_scratch_dir + "/GMAccepted", 0, 0);
            //            nifti_image_write(AcceptedGMImage);
            nifti_image_free(AcceptedGMImage);
        }
        else{
            AcceptedGMBool=new bool[numel];
            for (int i=0; i<numel; i++) {
                AcceptedGMBool[i]=0;
            }
        }
        

        if (ECSFBool!=NULL) {
            delete [] ECSFBool;
            ECSFBool=NULL;
        }
        
        
        ofstream TxtFileLeaf(FilenamePA_26.c_str());
        for (int l=0; l<numbleavesTree; l++) {
            LeafID * LeafIDToPrint=LeafCharacterisation(TreeToAnalyse, VectorLeaves[l], segment_analysis, SummarisedSeg1, 0.1);
            PrintLeafID(LeafIDToPrint, TxtFileLeaf, numbmodal);
            delete LeafIDToPrint;
            LeafIDToPrint=NULL;
            TxtFileLeaf<<endl;
        }
        cout<<"Leaf characterisation  done"<<endl;

        nifti_image * SegToAnalyseHardNIV=NULL;
        if(segment_analysis->flag_connect){
            cout << "Connect is on "<< endl;
            nifti_image * SegToAnalyseHard=NULL;
            if (BinarySeg!=NULL) {
                SegToAnalyseHard=CopyFloatNii(BinarySeg);
                SegToAnalyseHardNIV=CopyFloatNii(BinarySegNIV);
            }
            else{
                if (SummarisedSeg1!=NULL) {
                    cout << "Doing the hard segmentation "<< endl;
                    SegToAnalyseHard=HardSegmentationLesionReconstruct_bis(SegToAnalyse, SummarisedSeg1, TreeToAnalyse, segment_analysis);
                    SegToAnalyseHardNIV=HardSegmentationLesionReconstruct_bis(SegToAnalyseNIV, SummarisedSeg1, TreeToAnalyse, segment_analysis);
                }
                else {
                    SegToAnalyseHard=CopyFloatNii(BinarySeg);
                    SegToAnalyseHardNIV=CopyFloatNii(BinarySegNIV);
                }
            }
            if (BinarySegNIV!=NULL) {
                nifti_image_free(BinarySegNIV);
                BinarySegNIV=NULL;
            }
            
            //          nifti_image *  SegToAnalyseHard=HardSegmentationLesionReconstruct(SegToAnalyse, TreeToAnalyse);
            cout << "Before saving SegToAnalyseHard "<< endl;
            nifti_set_filenames(SegToAnalyseHard, FilenamePA_19.c_str(), 0, 0);
            nifti_image_write(SegToAnalyseHard);
            cout <<"Preparation before SP"<<endl;
            if (SPRegion!=NULL) {
                nifti_image * SPRegionNii=CreateNiiFromArray(SPRegion, SegToAnalyseNIV, numel);
                //                nifti_set_filenames(SPRegionNii, (g_scratch_dir + "/TestSP.nii.gz").c_str(), 0, 0);
                //                nifti_image_write(SPRegionNii);
                nifti_image_free(SPRegionNii);
            }
            cout<<"flag ECSF "<<segment_analysis->flag_inPriorsECSF<<endl;
            cout<<"SP Region written "<<endl;
            //            Correction for SPRegion if it exists
            if (SPRegion!=NULL && segment_analysis->flag_SPCorrection) {
                cout <<"Performing the SP correction"<<endl;
                float * SegHardData=static_cast<float*>(SegToAnalyseHard->data);
                float * SegHardDataNIV=static_cast<float*>(SegToAnalyseHardNIV->data);
                float * SegData=static_cast<float *>(SegToAnalyse->data);
                float * SegDataNIV=static_cast<float *>(SegToAnalyseNIV->data);
                for (int i=0; i<numel; i++) {
                    if (SPRegion[i]) {
                        SegHardData[i]=0;
                        SegHardDataNIV[i]=0;
                        SegData[i]=0;
                        SegDataNIV[i]=0;
                    }
                }
                cout << "New Seg Data after SP is "<< GetSum(SegData,numel)<< endl;
            }


            
            //        ConnectLabelLesions=ImageComponentLabeling(SegToAnalyseHard,segment_analysis);
            
            //RefinedConnectLabelLesions=ImageRefinedLabel(SegToAnalyse, segment_analysis);
            //ImageOrderedLesion=ImageOrderedVolumeLabel(SegToAnalyse, segment_analysis);
            if(segment_analysis->flag_GIFPriors){
                segment_analysis->Neigh=6;
            }
            int * ComponentLabel=ComponentLabeling(SegToAnalyseHard,segment_analysis->Neigh);
            int * ComponentLabelNIV=ComponentLabeling(SegToAnalyseHardNIV, segment_analysis->Neigh);
            cout<<"Connected components labeled done"<<endl;
            float * ComponentLabelFloat=TranscribeArray<int, float>(ComponentLabel, numel);
            ConnectLabelLesions=CreateNiiFromArray(ComponentLabelFloat, SegToAnalyse,numel);
            delete [] ComponentLabelFloat;
            ComponentLabelFloat=NULL;
            float VolumeVox=SegToAnalyse->pixdim[1]*SegToAnalyse->pixdim[2]*SegToAnalyse->pixdim[3];
            OrderedLabels=OrderedVolumeLabel(ComponentLabel, 3, numel,VolumeVox); //WARNING : THIS CHANGE from 3 to 2 at 14/11/14
            delete [] ComponentLabel;
            ComponentLabel=NULL;
            
            int * OrderedLabelsNIV=OrderedVolumeLabel(ComponentLabelNIV, 3, numel,VolumeVox);
            ImageOrderedLesion=nifti_copy_nim_info(SegToAnalyse);
            ImageOrderedLesion->data=(void *) calloc(SegToAnalyse->nvox, sizeof(float));
            float * ImageOrderedLesion_PTR=static_cast<float *>(ImageOrderedLesion->data);
            for (int i=0; i<numel; i++) {
                ImageOrderedLesion_PTR[i]=(float)OrderedLabels[i];
            }

            int * OrderedLabelsSec=NULL;
            int * ComponentLabelSec=NULL;
            if(SegSecondary!=NULL && segment_analysis->flag_Secondary){
                ComponentLabelSec = ComponentLabeling(SegSecondary,segment_analysis->Neigh);
                OrderedLabelsSec=OrderedVolumeLabel(ComponentLabelSec, 3, numel,VolumeVox,segment_analysis->flag_Secondary);
            }

            
            int  Dim[3];
            for (int d=0; d<3; d++) {
                Dim[d]=SegToAnalyse->dim[d+1];
            }
            
            
              vector<Outlier *> OutlierVectorLesions;
              vector<Outlier *> OutlierVectorLesionsSec;
            vector<Lesion*> LesionVector;
            cout<<"PriorsCGM flag "<<segment_analysis->flag_inPriorsCGM;
            cout<<" PriorsECSF flag"<<segment_analysis->flag_inPriorsECSF;
            cout<<" OldLesion flag"<<segment_analysis->flag_oldLesion<<endl;
            if (segment_analysis->flag_inPriorsCGM && segment_analysis->flag_inPriorsECSF && !segment_analysis->flag_oldLesion) {
                cout<<" doing the outlier way" <<endl;
                OutlierVectorLesions=GetVectorOutliers_NIVCorr(OrderedLabels,OrderedLabelsNIV,SegToAnalyse,  SegToAnalyseNIV,TreeToAnalyse,  segment_analysis);
                //                CorrectionForSymmetry(OutlierVectorLesions, 0, 10, 20, OrderedLabels, Dim);
                if(segment_analysis->flag_Secondary){
                    if (OrderedLabelsSec!=NULL){
                        if (CountNonZero(OrderedLabelsSec, numel) > 0){
                        cout << "Possibility to look into individual elements"<< CountNonZero(OrderedLabelsSec,numel) << endl;
                        SaveTmpResult(OrderedLabelsSec,FilenamePA_ConnectSec.c_str(),SegToAnalyse);
                        OutlierVectorLesionsSec=GetVectorOutliers(OrderedLabelsSec,TreeToAnalyse,segment_analysis,2.0/3);
                        }
                        else{
                        cout << "No individual baby lesions to look at" << endl;
                        }
                    }
                }
                int lo1=OutlierVectorLesions.size();
                ofstream TxtFileOutlier1(FilenamePA_9c.c_str());
                ofstream TxtFileOutlierBis1(FilenamePA_9d.c_str());
                int * GravImage1=NULL;
                if (lo1>0) {
                    GravImage1=new int[3];
                    TxtFileOutlier1 <<"Centre gravity ";
                    for (int d=0;d<3; d++) {
                        GravImage1[d]=OutlierVectorLesions[0]->CentreGravity[d]-OutlierVectorLesions[0]->VectorDiffGrav[d];
                        TxtFileOutlier1 << GravImage1[d]<<" ";
                    }
                    TxtFileOutlier1<<endl;
                }
                ImageOrderedLesion_PTR=static_cast<float*>(ImageOrderedLesion->data);
                for (int i=0; i<numel; i++) {
                    ImageOrderedLesion_PTR[i]=OrderedLabels[i];
                }
                cout<<"Printing the outlier way..."<<endl;
                for (int l=0; l<lo1; l++) {
                    cout <<"Label "<<l+1;
                    TxtFileOutlier1<<endl;
                    TxtFileOutlierBis1<<"Label "<<l+1<<" ";
                    TxtFileOutlier1<<"Label "<<l+1<<endl;
                    PrintOutlierCharacteristics(OutlierVectorLesions[l], TxtFileOutlier1, numbmodal, GravImage1, PixDim, segment_analysis);
                    PrintOutlierCharacteristicsForFile(OutlierVectorLesions[l], TxtFileOutlierBis1, numbmodal, GravImage1, PixDim, segment_analysis);
                    TxtFileOutlierBis1<<endl;
                    cout<<"...printed"<<endl;
                }
                cout << "Passed printing "<<endl;


            }
            
            else {
                cout <<"Doing the old way"<<endl;
                vector<Lesion *> LesionVector1=GetVectorLesion_quat(OrderedLabels, SummarisedSegDec,TreeToAnalyse,segment_analysis,VentricleBool,DGMBool,AcceptedGMBool,SPRegion,ICSF,DGM,OrderedLabelsNIV,SegToAnalyse,SegToAnalyseNIV);
                cout<<"First vector Lesion out"<<endl;
                
                //            Since the lesion segmentation is progressively modified, the label of one might have been completely taken by another one, further down the line
                int numbLes=LesionVector1.size();
                for (int l=0; l<numbLes; l++) {
                    delete LesionVector1[l];
                    LesionVector1[l]=NULL;
                }
                cout<<"Void first lesion vector"<<endl;
                vector<Lesion *> LesionVector2=GetVectorLesion_quat(OrderedLabels, SummarisedSegDec,TreeToAnalyse,segment_analysis,VentricleBool,DGMBool,AcceptedGMBool,SPRegion,ICSF,DGM,OrderedLabelsNIV,SegToAnalyse,SegToAnalyseNIV);


                //            Change OrderedLabels to decrease index of those after now void lesion
                for (int i=0; i<numel; i++) {
                    ImageOrderedLesion_PTR[i]=(float)OrderedLabels[i];
                }


                int * VolumeLesion=GetVolumeLabels(OrderedLabels, numel);
                int maxLabel=GetMaxLabel(OrderedLabels, numel);
                int CountNonZeroLesion=CountNonZero(VolumeLesion, maxLabel);

                if(CountNonZeroLesion !=maxLabel){
                    cout<<"Zeroed lesion is "<<maxLabel-CountNonZeroLesion<<endl;
                    int CountZeroTreated=0;
                    for (int l=0; l<maxLabel; l++) {
                        cout<<"Lesion "<<l+1<<" is "<<VolumeLesion[l]<<endl;
                        if (VolumeLesion[l]==0) {
                            for (int i=0; i<numel; i++) {
                                if (OrderedLabels[i]>=l+1-CountZeroTreated) {
                                    OrderedLabels[i]=OrderedLabels[i]-1;
                                }
                            }
                            CountZeroTreated++;
                        }
                    }
                }


                //            vector<Outlier * > OutlierVectorLesion=GetVectorOutliers(OrderedLabels, TreeToAnalyse,segment_analysis);


                if (segment_analysis->flag_inFPTPReclassif) {

                    ImageLesionCorr=LesionReconstructTrueSeg_ter(OrderedLabels, SegToAnalyse, SummarisedSeg1, LesionVector2, TreeToAnalyse, 1,segment_analysis);
                    //                ImageLesionCorr=LesionCorrectedTrueSeg_bis(SegToAnalyse, SummarisedSeg1, TreeToAnalyse, 1, segment_analysis);
                    nifti_set_filenames(ImageLesionCorr, FilenamePA_23Pre.c_str(), 0, 0);
                    nifti_image_write(ImageLesionCorr);
                    nifti_image_free(ImageLesionCorr);
                    ImageLesionCorr=NULL;


                    map<int,int> LesionToChange=ReadReclassifDecisionFromFile(segment_analysis->filename_inFPTPReclassif);
                    //                First check if there are any lesion for which classification need to be modified.
                    int sizeLC=LesionToChange.size();
                    int sizeLI=LesionVector2.size();
                    if (sizeLC!=sizeLI) {
                        cout<<"Impossible to do reclassification correction due to non agreement in sizes"<<endl;
                    }
                    else{
                        int CountToChange=0;
                        for (int l=0; l<sizeLC; l++) {
                            std::map<int,int>::const_iterator pos = LesionToChange.find(l);
                            int value = pos->second;
                            if (value) {
                                CountToChange++;
                                if(LesionVector2[l]->LesionType<0){
                                    LesionVector2[l]->LesionType=RecTP;
                                }
                                else{
                                    LesionVector2[l]->LesionType=RecFP;
                                }
                            }
                        }
                    }
                }
                cout<<GetMaxLabel(OrderedLabels, numel)<<" and in vector"<<LesionVector2.size()<<endl;
                LesionVector=GetVectorLesion_quat(OrderedLabels, SummarisedSeg1, TreeToAnalyse, segment_analysis, VentricleBool, DGMBool, AcceptedGMBool,SPRegion, ICSF, DGM, OrderedLabelsNIV, SegToAnalyse, SegToAnalyseHardNIV);
                
                int * VolumeLesionFin=GetVolumeLabels(OrderedLabels, numel);
                int maxLabelFin=GetMaxLabel(OrderedLabels, numel);
                int CountNonZeroLesionFin=CountNonZero(VolumeLesionFin, maxLabelFin);
                
                if(CountNonZeroLesionFin !=maxLabelFin){
                    cout<<"Zeroed lesion is "<<maxLabelFin-CountNonZeroLesionFin<<endl;
                    int CountZeroTreated=0;
                    for (int l=0; l<maxLabelFin; l++) {
                        cout<<"Lesion "<<l+1<<" is "<<VolumeLesionFin[l]<<endl;
                        if (VolumeLesionFin[l]==0) {
                            for (int i=0; i<numel; i++) {
                                if (OrderedLabels[i]>=l+1-CountZeroTreated) {
                                    OrderedLabels[i]=OrderedLabels[i]-1;
                                }
                            }
                            CountZeroTreated++;
                        }
                    }
                }
                
                cout<<GetMaxLabel(OrderedLabels, numel)<<" and in vector"<<LesionVector.size()<<endl;
                cout<<"Vector final lesion obtained"<<endl;
            }

            
            if(OrderedLabelsNIV!=NULL){
                delete [] OrderedLabelsNIV;
                OrderedLabelsNIV=NULL;
            }
            if(ComponentLabelNIV!=NULL){
                delete [] ComponentLabelNIV;
                ComponentLabelNIV=NULL;
            }
            if(SegToAnalyseHardNIV!=NULL){
                nifti_image_free(SegToAnalyseHardNIV);
                SegToAnalyseHardNIV=NULL;
            }
            int * OutlierClassif=NULL;
            int * OutlierClassifSec=NULL;
            float * SegmentationInfarct=NULL;
            if (OutlierVectorLesions.size()>0) {
                OutlierClassif=OutlierClassification(OrderedLabels,OutlierVectorLesions,Dim);
                if (segment_analysis->flag_Secondary && OrderedLabelsSec!=NULL){
                    OutlierClassifSec=OutlierClassification(OrderedLabelsSec,OutlierVectorLesionsSec,Dim);
                }
            }
            if (SegToAnalyse==NULL){
                cout << "Careful ! SegToAnalyse is NULL "<<endl;
            }
            cout << "Deleting and preparing done after printing"<<endl;
            //            int * OutlierClassif=OutlierClassification(OrderedLabels,OutlierVectorLesions,Dim);
            nifti_image * OutlierClassifNii=NULL;
            if(OutlierClassif!=NULL){
                OutlierClassifNii=CreateNiiFromArray(OutlierClassif, SegToAnalyse, numel);
            }
            if (LesionVector.size()>0 && OutlierClassif==NULL) {
                OutlierClassif=LesionClassification(OrderedLabels, LesionVector, Dim);
                cout<<"Lesion classification obtained"<<endl;
            }
            if (OutlierClassif==NULL){
                cout << "Zeros everywhere in OutlierClassif"<<endl;
                OutlierClassif=new int[numel];
                for (int i=0; i<numel; i++) {
                    OutlierClassif[i]=0;
                }
                cout << "Empty outliers"<<endl;
            }
            if(SegToAnalyse==NULL){
                cout << "SegToAnalyse is NULL"<<endl;
            }
            else{
                cout << SegToAnalyse->nvox << "and numel is "<< numel << "and non zero"<< CountNonZero(static_cast<float*>(SegToAnalyse->data),numel)<<endl;
            }
            nifti_image * ImageLesionClassif=NULL;
            ImageLesionClassif=CreateNiiFromArray(OutlierClassif, SegToAnalyse, numel);
            //                nifti_image * ImageLesionClassif=nifti_copy_nim_info(SegToAnalyse);
            //                ImageLesionClassif->data=(void *) calloc(SegToAnalyse->nvox, sizeof(float));
            //                float * ImageLesionClassif_PTR=static_cast<float *>(ImageLesionClassif->data);
            //                for (int i=0; i<numel; i++) {
            //                    ImageLesionClassif_PTR[i]=(float)LesionClassif[i];
            //                }
            cout <<"Image LesionClassif obtained"<<endl;
            if (segment_analysis->flag_infarcts){
                if(SegToAnalyse!=NULL){
                    float * SegData=static_cast<float*>(SegToAnalyse->data);
                    SegmentationInfarct = new float[numel];
                    for(int i=0;i<numel;i++){
                        SegmentationInfarct[i]=0;
                        if(OutlierClassif[i]==7){
                            SegmentationInfarct[i]=SegData[i];
                        }
                    }
                    nifti_image * InfarctNii = CreateNiiFromArray(SegmentationInfarct,SegToAnalyse,numel);
                    nifti_set_filenames(InfarctNii, FilenameInfarct.c_str(), 0,0);
                    nifti_image_write(InfarctNii);
                    nifti_image_free(InfarctNii);
                }

            }
            if (segment_analysis->flag_Saving) {
                nifti_set_filenames(ImageLesionClassif, FilenamePA_18.c_str(), 0, 0);
                nifti_image_write(ImageLesionClassif);
            }
            

            if(OutlierClassifNii!=NULL && segment_analysis->flag_Saving){
                nifti_set_filenames(OutlierClassifNii, FilenamePA_18b.c_str(), 0, 0);
                nifti_image_write(OutlierClassifNii);
            }
            
            //                if (LesionClassif!=NULL) {
            //                    delete [] LesionClassif;
            //                    LesionClassif=NULL;
            //                }
            if (OutlierClassif!=NULL) {
                delete [] OutlierClassif;
                OutlierClassif=NULL;
            }
            
            //            ImageLesionCorr=CorrectionLesionFromClassif(ImageLesionClassif,SegToAnalyse);
            //            ImageLesionCorr=LesionCorrectedTrueSeg_bis(SegToAnalyse, SummarisedSeg1, TreeToAnalyse, 1, segment_analysis);
            //            nifti_set_filenames(ImageLesionCorr, FilenamePA_23.c_str(), 0, 0);
            //            nifti_image_write(ImageLesionCorr);
            
            ImageLabelCorr=CorrectionLesionFromClassif(ImageLesionClassif, ImageOrderedLesion);
            if (segment_analysis->flag_Saving && ImageLabelCorr!=NULL) {
                nifti_set_filenames(ImageLabelCorr, FilenamePA_25.c_str(), 0, 0);
                nifti_image_write(ImageLabelCorr);
            }
            cout<<"ImageLabelCorr obtained"<<endl;
            if (ImageLesionClassif!=NULL){
                nifti_image_free(ImageLesionClassif);
                ImageLesionClassif=NULL;
            }
            
            ofstream TxtFileLesion(FilenamePA_9.c_str());
            ofstream TxtFileLesionBis(FilenamePA_9b.c_str());
            ofstream TxtFileOutlier(FilenamePA_9c.c_str());
            ofstream TxtFileOutlierBis(FilenamePA_9d.c_str());
            ofstream TxtFileOutlierSec(FilenamePA_9e.c_str());
            int GravImage[3];
            int numbFinalLesions=LesionVector.size();
            if (LesionVector.size()>0) {
                
                int * CountperType= CountLesionType(LesionVector);
                cout<<"Lesion count type obtained"<<endl;
                float * VolumeperType=VolumeLesionType(LesionVector);
                cout<<"Volume lesion type obtained"<<endl;
                float TotVolume=TotVolumeLesion(LesionVector);

                cout<<"NumbFileLesion is "<<numbFinalLesions<<endl;

                cout<<"TextFile opened"<<endl;
                TxtFileLesion<<"TotLesionVolume "<<TotVolume;
                if (numbFinalLesions>0) {
                    for (int d=0;d<3; d++) {
                        GravImage[d]=LesionVector[0]->CentreGravity[d]-LesionVector[0]->VectorDiffGrav[d];
                    }
                }
                cout<<"Preparation to print done"<<endl;
                PrintCountVolperType(CountperType,VolumeperType,TxtFileLesion);
                for (int l=0; l<numbFinalLesions; l++) {
                    cout <<"Label "<<l+1;
                    TxtFileLesion<<endl;
                    TxtFileLesionBis<<"Label "<<l+1<<" ";
                    TxtFileLesion<<"Label "<<l+1<<endl;
                    PrintLesion(LesionVector[l],TxtFileLesion,numbmodal,GravImage,PixDim);
                    PrintLesionForFile(LesionVector[l],TxtFileLesionBis,numbmodal,GravImage,PixDim);
                    TxtFileLesionBis<<endl;
                    cout<<"...printed"<<endl;
                }

                if (CountperType!=NULL) {
                    delete [] CountperType;
                    CountperType=NULL;
                }
                if (VolumeperType!=NULL) {
                    delete [] VolumeperType;
                    VolumeperType=NULL;
                }
            }
            
            int lo=OutlierVectorLesions.size();
            
            

            if (lo>0) {
                cout<<"Printing outlier vector"<<endl;
                for (int d=0;d<3; d++) {
                    GravImage[d]=OutlierVectorLesions[0]->CentreGravity[d]-OutlierVectorLesions[0]->VectorDiffGrav[d];
                }
            }
            for (int l=0; l<lo; l++) {
                cout <<"Label "<<l+1;
                TxtFileOutlier<<endl;
                TxtFileOutlierBis<<"Label "<<l+1<<" ";
                TxtFileOutlier<<"Label "<<l+1<<endl;
                PrintOutlierCharacteristics(OutlierVectorLesions[l], TxtFileOutlier, numbmodal, GravImage, PixDim, segment_analysis);
                PrintOutlierCharacteristicsForFile(OutlierVectorLesions[l], TxtFileOutlierBis, numbmodal, GravImage, PixDim, segment_analysis);
                TxtFileOutlierBis<<endl;
                cout<<"...printed"<<endl;
            }

            
            nifti_set_filenames(SegToAnalyse, FilenamePA_27.c_str(), 0, 0);
            nifti_image_write(SegToAnalyse);
            nifti_image * ImageLesionCorrSec=NULL;
            cout <<"Saved SegToAnalyse"<<endl;
            if (segment_analysis->flag_correct) {
                if(OutlierVectorLesionsSec.size()>0){
                    cout << "Reconstructing Secondary" << endl;
                    SaveTmpResult(OrderedLabelsSec,FilenamePA_ConnectSec.c_str(),SegToAnalyse);
                    ImageLesionCorrSec=LesionReconstructTrueSegOutlier_ter(OrderedLabelsSec,SegSecondary,SummarisedSeg1,OutlierVectorLesionsSec,TreeToAnalyse,1,segment_analysis);
                    int los = OutlierVectorLesionsSec.size();
                    for (int l=0; l<los; l++) {
                        cout <<"Label "<<l+1;
                        TxtFileOutlierSec<<endl;
//                        TxtFileOutlierBis<<"Label "<<l+1<<" ";
                        TxtFileOutlierSec<<"Label "<<l+1<<endl;
                        PrintOutlierCharacteristics(OutlierVectorLesionsSec[l], TxtFileOutlierSec, numbmodal, GravImage, PixDim, segment_analysis);
//                        PrintOutlierCharacteristicsForFile(OutlierVectorLesionsSec[l], TxtFileOutlierBis, numbmodal, GravImage, PixDim, segment_analysis);
//                        TxtFileOutlierBis<<endl;
                        cout<<"...printed"<<endl;
                    }
                }
                if (OutlierVectorLesions.size()>0) {
                    ImageLesionCorr=LesionReconstructTrueSegOutlier_ter(OrderedLabels, SegToAnalyse, SummarisedSeg1, OutlierVectorLesions, TreeToAnalyse, 1,segment_analysis);
                }
                else if(LesionVector.size()>0){
                    ImageLesionCorr=LesionReconstructTrueSeg_ter(OrderedLabels, SegToAnalyse, SummarisedSeg1, LesionVector, TreeToAnalyse, 1,segment_analysis);
                }
                if(ImageLesionCorr==NULL){
                    cout << "Empty Image Lesion corr" << endl;
                    float * EmptyDataLesionCorr=new float[numel];
                    for(int i=0;i<numel;i++){
                        EmptyDataLesionCorr[i]=0;
                    }
                    ImageLesionCorr=CreateNiiFromArray(EmptyDataLesionCorr,TreeToAnalyse->GetDataImage(),numel);
                    delete [] EmptyDataLesionCorr;
                }
                float * LesionCorrData=static_cast<float*>(ImageLesionCorr->data);
                nifti_set_filenames(ImageLesionCorr, FilenamePA_23.c_str(), 0, 0);
                nifti_image_write(ImageLesionCorr);
                if(segment_analysis->flag_Secondary && ImageLesionCorrSec!=NULL){
                    float * LesionCorrDataSec=static_cast<float*>(ImageLesionCorrSec->data);
                    AddElementwiseInPlace(LesionCorrData,LesionCorrDataSec,numel);
                }
                for(int i=0;i<numel;i++){
                    if(LesionCorrData[i]>1){
                        LesionCorrData[i]=1;
                    }
                }

                if (ImageLesionCorrSec==NULL){
                    cout << "Creating empty image for secondary"<< endl;
                    float * EmptyData = new float[numel];
                    for(int i=0;i<numel;i++){
                        EmptyData[i]=0;
                    }
                    ImageLesionCorrSec=CreateNiiFromArray(EmptyData,TreeToAnalyse->GetDataImage(),numel);
                }
                if (ImageLesionCorrSec!=NULL){
                    nifti_set_filenames(ImageLesionCorrSec, FilenamePA_29.c_str(), 0, 0);
                    nifti_image_write(ImageLesionCorrSec);
                }
                //                ImageLesionCorr=LesionCorrectedTrueSeg_bis(SegToAnalyse, SummarisedSeg1, TreeToAnalyse, 1, segment_analysis);

            }
            
            if (ImageLesionCorr!=NULL || OutlierVectorLesions.size()==0){
                nifti_image * SummarisedSeg2=SummarizedSegmentation(TreeToAnalyse, LesionClasses, LesionPartsFromUniform, CSFWMClasses, CSFWMUniform, segment_analysis, ICSF, ImageLesionCorr);
                nifti_set_filenames(SummarisedSeg2, FilenamePA_22.c_str(), 0, 0);
                nifti_image_write(SummarisedSeg2);
                //                if (ImageLabelCorr!=NULL) {
                //                    nifti_image_free(ImageLabelCorr);
                //                    ImageLabelCorr=NULL;
                //                }
                nifti_image_free(SummarisedSeg2);
                SummarisedSeg2=NULL;
                cout << "Freeing Seg 2" << endl;
            }
            //            if (LesionPartsFromUniform!=NULL) {
            //                nifti_image_free(LesionPartsFromUniform);
            //                LesionPartsFromUniform=NULL;
            //            }
            //            Create Ordered Labels based on Corrected image lesion, in order to remove the parts that do not correspond anymore. As the correction is based on the hard
            if (ImageLesionCorr!=NULL) {
                float * ImageLesionCorrData=static_cast<float *>(ImageLesionCorr->data);
                int * OrderLabTemp=new int[numel];
                for (int i=0; i<numel; i++) {
                    if (ImageLesionCorrData[i]>0.5) {
                        OrderLabTemp[i]=OrderedLabels[i];
                    }
                    else{
                        OrderLabTemp[i]=0;
                    }
                }
                float VolumeVox=ImageLesionCorr->pixdim[1]*ImageLesionCorr->pixdim[2]*ImageLesionCorr->pixdim[3];
                OrderedLabelsCorr=OrderedVolumeLabel(OrderLabTemp, 1, numel,VolumeVox);
                delete [] OrderLabTemp;
                OrderLabTemp=NULL;
            }
            
            for(int l=0;l<numbFinalLesions;l++){
                delete LesionVector[l];
                LesionVector[l]=NULL;
            }
            
            if (SegToAnalyseHard!=NULL) {
                nifti_image_free(SegToAnalyseHard);
                SegToAnalyseHard=NULL;
            }
        }
        if(ConnectLabelLesions!=NULL){
            if (segment_analysis->flag_Saving) {
                nifti_set_filenames(ConnectLabelLesions, FilenamePA_6.c_str(), 0, 0);
                nifti_image_write(ConnectLabelLesions);
            }
            nifti_image_free(ConnectLabelLesions);
            ConnectLabelLesions=NULL;
        }
        if(RefinedConnectLabelLesions!=NULL){
            if (segment_analysis->flag_Saving) {
                nifti_set_filenames(RefinedConnectLabelLesions, FilenamePA_7.c_str(), 0, 0);
                nifti_image_write(RefinedConnectLabelLesions);
            }

            nifti_image_free(RefinedConnectLabelLesions);
            RefinedConnectLabelLesions=NULL;
        }
        if(ImageOrderedLesion!=NULL){
            cout << "Saving Image Order";
            nifti_set_filenames(ImageOrderedLesion, FilenamePA_8.c_str(), 0, 0);
            nifti_image_write(ImageOrderedLesion);
            nifti_image_free(ImageOrderedLesion);
            ImageOrderedLesion=NULL;
        }

        float * MeanLesion=GetMeanLesion(SegToAnalyse,TreeToAnalyse);
        if (MeanLesion!=NULL) {
            TxtFile << "MeanLesion ";
            for (int m=0; m<numbmodal; m++) {
                TxtFile << MeanLesion[m]<< " ";
            }
            TxtFile<<endl;
            delete [] MeanLesion;
            MeanLesion=NULL;
        }
        segment_analysis->flag_Les=1;
        cout<<"Lesion rebuilt gaussian components saved"<<endl;
        
        
        

        nifti_image * SummarisedSeg2=SummarizedSegmentation(TreeToAnalyse, LesionClasses, LesionPartsFromUniform, CSFWMClasses, CSFWMUniform, segment_analysis, ICSF, ImageLesionCorr);
        nifti_set_filenames(SummarisedSeg2, FilenamePA_22.c_str(), 0, 0);
        nifti_image_write(SummarisedSeg2);
        cout << "Summarised written "<< endl;
        if (ImageLabelCorr!=NULL) {
            nifti_image_free(ImageLabelCorr);
            ImageLabelCorr=NULL;
        }
        
        //        VentricleSeg=VentricleSegmentation(SummarisedSeg2, segment_analysis, TreeToAnalyse);
        //        nifti_image * DGMExtraction=ExtractWMDGM(SummarisedSeg2, segment_analysis, TreeToAnalyse);
        if (LesionPartsFromUniform!=NULL) {
            nifti_image_free(LesionPartsFromUniform);
            LesionPartsFromUniform=NULL;
        }
        
        if (SummarisedSeg2!=NULL){
            nifti_image_free(SummarisedSeg2);
            SummarisedSeg2=NULL;
        }
        if(SummarisedSeg1!=NULL){
            nifti_image_free(SummarisedSeg1);
            SummarisedSeg1=NULL;
        }
        if(SummarisedSegDec!=NULL){
            nifti_image_free(SummarisedSegDec);
            SummarisedSegDec=NULL;
        }
        if (CSFWMUniform!=NULL) {
            nifti_image_free(CSFWMUniform);
            CSFWMUniform=NULL;
        }
        
        
        if (segment_analysis->flag_WMCard) {
            cout << "Writing WML card"<<endl;
            string FilenameWM=nifti_makebasename(segment_analysis->filename_SegTot);
            int Index=FilenameWM.find_last_of('/');
            string FilenameWM_b=FilenameWM.substr(0,Index+1);
            string FilenameWM_e=FilenameWM.substr(Index+1,FilenameWM.length());
            FilenameWM=FilenameWM_b+"WMCard_"+FilenameWM_e+".txt";
            segment_analysis->filename_WMCard=const_cast<char*>(FilenameWM.c_str());
            WriteWMLIdentityCard(TreeToAnalyse, LesionRule, Modalities,segment_analysis);
        }
        if (LesionRule !=NULL) {
            delete LesionRule;
            LesionRule=NULL;
        }
        cout << "Everything done for lesion construction"<<endl;
    }

    
    
    
    //        Building the laplace solution using Parenchyma extraction and Ventricle segmentation
    if(segment_analysis->flag_LapAnalysis){
        
        FilenamePA=nifti_makebasename(segment_analysis->filename_inSum);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        
        int * L2SToUse;
        if (TreeToAnalyse!=NULL) {
            L2SToUse=TreeToAnalyse->GetL2S();
        }
        else {
            L2SToUse=new int[numel];
            for (int i=0; i<numel; i++) {
                L2SToUse[i]=i;
            }
        }
        cout<<"Going into Laplace analysis"<<endl;
        
        //        First check if there is the Label for the wanted analysis
        nifti_image * SegLesBasis=NULL;
        int * LesionLabel=NULL;
        if (OrderedLabelsCorr==NULL && OrderedLabels==NULL) {
            if (segment_analysis->flag_inLesCorr) {
                SegLesBasis=ReadFromFilename(segment_analysis->filename_inLesCorr);
            }
            else if (segment_analysis->flag_inLes){
                SegLesBasis=ReadFromFilename(segment_analysis->filename_inLes);
            }
            nifti_image * BinarySegForLaplace=HardSegmentationThreshold(SegLesBasis, 0.1);
            int * LesionLabelPrep=ComponentLabeling(BinarySegForLaplace, segment_analysis->Neigh);
            LesionLabel=OrderedVolumeLabel(LesionLabelPrep, 1, numel);
            nifti_image_free(BinarySegForLaplace);
            delete [] LesionLabelPrep;
            LesionLabelPrep=NULL;
            BinarySegForLaplace=NULL;
        }
        else{
            SegLesBasis=CopyFloatNii(SegToAnalyse);
            if (OrderedLabelsCorr!=NULL) {
                LesionLabel=CopyArray(OrderedLabelsCorr, numel);
            }
            else{
                LesionLabel=CopyArray(OrderedLabels, numel);
            }
        }
        cout<<"Lesion presence checked"<<endl;
        
        
        //        First check that everything needed is available and create appropriate element if needed. Needed DGM priors, ICSF priors for DGM and Ventricles segmentation + Summarised
        if(VentricleSeg==NULL){
            if (segment_analysis->flag_inVentricleSeg) {
                cout<<"Reading directly VentricleSeg from file"<<endl;
                VentricleSeg=ReadFromFilename(segment_analysis->filename_inVentricleSeg);
            }
            else{ //        if Ventricle segmentation not available, create it based on summarised and on ICSF
                if(PriorsICSF==NULL) {
                    if (segment_analysis->flag_inPriorsICSF) {
                        PriorsICSF=ReadFromFilename(segment_analysis->filename_inPriorsICSF);
                    }
                    else{
                        cout<<"Pb no PriorsICSF file to build ventricle segmentation for Laplace solution"<<endl;
                    }
                }
                if (SummarisedSeg1==NULL) {
                    if (segment_analysis->flag_inSum) {
                        SummarisedSeg1=ReadFromFilename(segment_analysis->filename_inSum);
                    }
                }
                if (TreeToAnalyse==NULL) {
                    cout<<"Impossible to get Ventricle Seg for Laplace because no Tree available"<<endl;
                }
                VentricleSeg=VentricleSegmentationPriors(SummarisedSeg1, segment_analysis, TreeToAnalyse, PriorsICSF);
                //                Correction for lesion if exist...


            }
        }
        if (DGMSeg==NULL) { // To build DGMSeg, we need SummarisedSeg and DGMPriors //        if DGM segmentation not available, create it based on summarised and on DGM priors.
            if (PriorsDGM==NULL) {
                if (segment_analysis->flag_inPriorsDGM) {
                    PriorsDGM=ReadFromFilename(segment_analysis->filename_inPriorsDGM);
                }
            }
            DGMSeg=DGMSegmentationPriors(SummarisedSeg1,segment_analysis,PriorsDGM,L2SToUse);
        }
        //        Build WM+DGM from DGMSeg
        nifti_image * WMDGMSeg=NULL;
        if (DGMSeg==NULL) {
            WMDGMSeg=ExtractParenchyma(SummarisedSeg1, VentricleSeg, TreeToAnalyse, segment_analysis);
        }
        else{
            WMDGMSeg=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexWM, TreeToAnalyse->GetL2S());
            float * DGMData=static_cast<float *>(DGMSeg->data);
            float * WMDGMSegData=static_cast<float*>(WMDGMSeg->data);
            if (DGMBool==NULL) {
                DGMBool=TranscribeArray<float, bool>(DGMData, numel);
            }
            for (int i=0; i<numel; i++) {
                WMDGMSegData[i]=((WMDGMSegData[i]+DGMBool[i])>0);
            }
        }
        
        //                Correction for potential hole between WM segmented and Ventricles
        
        
        vector<int> DimVec;
        for (int d=0; d<3; d++) {
            DimVec.push_back(WMDGMSeg->dim[d+1]);
        }
        bool * MaskData=static_cast<bool *>(TreeToAnalyse->GetMask()->data);
        bool * OpposedMask=OpposeBoolArray(MaskData, numel);
        bool * DilatedMask=ErosionTemplate(OpposedMask, 3, DimVec, 0);
        
        float * WMDGMData1=static_cast<float*>(WMDGMSeg->data);
        bool * WMDGMBool=TranscribeArray<float, bool>(WMDGMData1, numel);
        
        float * VentricleData=static_cast<float*>(VentricleSeg->data);
        bool * VentrBool=TranscribeArray<float, bool>(VentricleData, numel);
//        nifti_image * DilatedNii=CreateNiiFromArray(DilatedMask, WMSeg, numel);
        nifti_image * DistanceOutNii=EuclideanDistanceImage(VentricleSeg, DilatedMask,TreeToAnalyse->GetMask());
        //        nifti_set_filenames(DistanceOutNii, (g_scratch_dir + "/TestOutDist.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(DistanceOutNii);
        nifti_image * DistanceWMDGM=EuclideanDistanceImage(WMDGMSeg, WMDGMBool, TreeToAnalyse->GetMask());
        nifti_image * DistanceVentr=EuclideanDistanceImage(WMDGMSeg, VentrBool, TreeToAnalyse->GetMask());
        nifti_image * CGMSeg=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexGM, L2SToUse);
        float * CGMData=static_cast<float*>(CGMSeg->data);
        nifti_image * CGMPriors=ReadFromFilename(segment_analysis->filename_inPriorsCGM);
        float * CGMPriorsData=static_cast<float*>(CGMPriors->data);
        for (int i=0; i<numel; i++) {
            if (CGMPriorsData[i]<0.1 && CGMData[i]>0) {
                CGMData[i]=0;
            }
        }
        bool * CGMBool = TranscribeArray<float, bool>(CGMData, numel);
        nifti_image * DistanceCGM=EuclideanDistanceImage(CGMSeg, CGMBool, TreeToAnalyse->GetMask());
        float * DistanceDataOut=static_cast<float*>(DistanceOutNii->data);
        float * DistanceDataWMDGM=static_cast<float*>(DistanceWMDGM->data);
        float * DistanceDataVentr=static_cast<float*>(DistanceVentr->data);
        float * DistanceDataCGM=static_cast<float*>(DistanceCGM->data);
        for (int i=0; i<numel; i++) {
            if (VentricleData[i]+WMDGMData1[i]==0) {
                if (DistanceDataOut[i]>20) {
                    if (DistanceDataOut[i]>DistanceDataWMDGM[i] && DistanceDataVentr[i]<DistanceDataOut[i]) {
                        if(DistanceDataVentr[i]<DistanceDataCGM[i]){
                            VentricleData[i]=1;
                        }
                    }
                }
            }
        }
        
        
        
        cout<<"WMDGMSeg done"<<endl;
        //        nifti_set_filenames(WMDGMSeg, (g_scratch_dir + "/TestWMDGMSeg.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(WMDGMSeg);
        //        nifti_set_filenames(VentricleSeg, (g_scratch_dir + "/TestVentricleSeg.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(VentricleSeg);

        //        Create also parenchyma to have Laplace done both on white matter + DGM and all parenchyma.
        

        
        //        nifti_image * VentricleSeg=VentricleSegmentation(SummarisedSeg1, segment_analysis,TreeToAnalyse);
        //        nifti_image * VentricleSeg=nifti_image_read(g_scratch_dir + "/VentriclesFin.nii.gz", 1);
        nifti_image * Parenchyma=ExtractParenchyma(SummarisedSeg1, VentricleSeg, TreeToAnalyse, segment_analysis);
        cout<<"Extraction parenchyma done "<<endl;
        //        nifti_image * DGMExtraction=ExtractWMDGM(SummarisedSeg1,VentricleSeg, segment_analysis, TreeToAnalyse);
        
        
        int * Dim=new int[3];
        vector<int> DimVector;
        int * Shift=new int[3];
        float * PixDim=new float[3];
        Shift[0]=1;
        for (int d=0; d<3; d++) {
            Dim[d]=SegLesBasis->dim[d+1];
            DimVector.push_back(Dim[d]);
            PixDim[d]=SegLesBasis->pixdim[d+1];
            Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
        }
        
        
//        nifti_image * LaplaceSolImagePar=NULL;
        nifti_image * LaplaceSolImageWMD=NULL;
//        float * DistanceNormalisedLaplacePar=NULL;
        float * DistanceNormalisedLaplaceWMD=NULL;
        //        Create vector for adding VentriclesSeg and Parenchyma extraction and in parallel for adding VentricleSeg and WMDGMSeg
        float * VentriclesData=static_cast<float *>(VentricleSeg->data);
        float * ParenchymaData=static_cast<float *>(Parenchyma->data);
        float * WMDGMData=static_cast<float *>(WMDGMSeg->data);
        
        if (segment_analysis->flag_inLes) {
            nifti_image * LesionSeg=ReadFromFilename(segment_analysis->filename_InLes);
            float * LesionData=static_cast<float*>(LesionSeg->data);
            int CountCorr=0;
            for (int i=0; i<numel; i++) {
                if (LesionData[i]>0&&VentriclesData[i]>0) {
                    VentriclesData[i]=0;
                    WMDGMData[i]=1;
                    CountCorr++;
                }
                
            }
            nifti_image_free(LesionSeg);
            LesionSeg=NULL;
            cout <<"Correction of "<<CountCorr<<" performed for Laplace Ventricular and WMDGM segmentation"<<endl;
        }
        
        
        vector<float *> ToAddPar;
        vector<float *> ToAddWMD;
        ToAddPar.push_back(VentriclesData);
        ToAddPar.push_back(ParenchymaData);
        ToAddWMD.push_back(VentriclesData);
        ToAddWMD.push_back(WMDGMData);
        
        //
        ////        First deal with Parenchyma solution
        //        float * ObjectTotTmpPar=AddArray(ToAddPar, numel);
        //        float * ObjectTotDPar=Erosion_bis(ObjectTotTmpPar, 3, DimVector, 0);
        //        float * ObjectTotCPar=Erosion_bis(ObjectTotDPar, 3, DimVector, 1);
        //        float * ObjectOutPar=SubtractArray(ObjectTotCPar, VentriclesData, numel);
        //
        //        delete [] ObjectTotTmpPar;
        //        bool * ObjectTotPar=TranscribeArray<float, bool>(ObjectTotCPar, numel);
        //        bool * ObjectInPar=TranscribeArray<float , bool>(VentriclesData, numel);
        //        bool * MaskDilPar=TranscribeArray<float, bool>(ObjectTotDPar, numel);
        //
        //        //        SaveTmpResult(ObjectIn, g_scratch_dir + "/ObjectIn.nii.gz", SegToAnalyse);
        //        //        SaveTmpResult(ObjectTot, g_scratch_dir + "/ObjectOut.nii.gz", SegToAnalyse);
        //
        //
        //
        //        int * L2S_Par=MakeL2S(ObjectTotPar, Dim);
        //        int numelmaskedLapPar=0;
        //        int * S2L_Par=MakeS2L(ObjectTotPar, Dim, numelmaskedLapPar);
        //
        //        if (LaplaceSolImagePar==NULL) { // Meaning that there was nothing saved as Laplace solution
        //            int * LabelLaplacePar=CreateLaplaceLabelling(ObjectInPar,ObjectTotPar,Dim);
        ////                SaveTmpResult(LabelLaplacePar, g_scratch_dir + "/TestLabel.nii.gz", WMDGMSeg);
        //            int numelmaskedbPar=0;
        //            int * L2SbPar=MakeL2S(MaskDilPar, Dim);
        //            int * S2LbPar=MakeS2L(MaskDilPar, Dim, numelmaskedbPar);
        //            if (MaskDilPar!=NULL) {
        //                delete [] MaskDilPar;
        //                MaskDilPar=NULL;
        //            }
        //            int * Label_sPar=CreateShort(LabelLaplacePar, S2LbPar, numelmaskedbPar);
        //            cout<<GetMax<int>(Label_sPar, numelmaskedbPar)<<endl;
        ////            SaveTmpResult_short(Label_sPar, g_scratch_dir + "/TestLabelB.nii.gz", WMDGMSeg, L2SbPar);
        //            if (LabelLaplacePar!=NULL) {
        //                delete [] LabelLaplacePar;
        //                LabelLaplacePar=NULL;
        //            }
        //
        ////            nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/TestLaplace.nii.gz", 1);
        ////             nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/FinalLaplace.nii.gz", 1);
        //            //            float * LaplaceTmpLong=static_cast<float *>(LaplaceTmp->data);
        //            //            float * LaplaceSol_short=CreateShort(LaplaceTmpLong, S2Lb, numelmaskedb);
        //
        //
        //            float * LaplaceSol_shortPar=SolvingLaplaceEquation(Label_sPar, Dim, Shift, PixDim, 5, L2SbPar, S2LbPar, numelmaskedbPar,SegLesBasis);
        //            float * LaplaceSolPar=CreateLong(LaplaceSol_shortPar,L2SbPar,numel);
        //
        //
        ////            float * LaplaceSol_shortPar=NULL;
        ////            float * LaplaceSolPar=static_cast<float *>(LaplaceTmp->data);
        ////                        SaveTmpResult(LaplaceSolPar, g_scratch_dir + "/FinalLaplace.nii.gz", SegLesBasis);
        //            // Clearing memory for all that won't be needed afterwards
        //            delete [] L2SbPar;
        //            delete [] S2LbPar;
        //            delete [] Label_sPar;
        //            if (LaplaceSol_shortPar!=NULL) {
        //                delete [] LaplaceSol_shortPar;
        //            }
        //
        //            L2SbPar=NULL;
        //            S2LbPar=NULL;
        //            Label_sPar=NULL;
        //            LaplaceSol_shortPar=NULL;
        //
        //
        //
        //            bool * ObjectIn_shortPar=CreateShort(ObjectInPar, S2L_Par, numelmaskedLapPar);
        //            bool * Mask_shortPar=CreateShort(ObjectTotPar, S2L_Par, numelmaskedLapPar);
        //            bool * Border0Par=CreateBorderFromBool_L2S(ObjectIn_shortPar, Dim, Shift, L2S_Par, S2L_Par, numelmaskedLapPar);
        //            //            SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", SegToAnalyse, L2S);
        //            bool * Border1Par=CreateBorderFromBool_L2S(Mask_shortPar, Dim, Shift, L2S_Par, S2L_Par, numelmaskedLapPar);
        //            float * LaplaceSol_sPar=CreateShort(LaplaceSolPar, S2L_Par, numelmaskedLapPar);
        //            float * DistanceNormalisedLaplace_sPar=NormalisedLaplaceLength(LaplaceSol_sPar, ObjectIn_shortPar, Mask_shortPar, L2S_Par, S2L_Par, numelmaskedLapPar, PixDim, Dim, Shift,SegLesBasis);
        //            int nlsize=segment_analysis->numbLaminae.size();
        //            for (int nl=0; nl<nlsize; nl++) {
        //                int * LayersNormDistPar=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_sPar,segment_analysis->numbLaminae[nl],numelmaskedLapPar);
        //                stringstream nlName;
        //                nlName << segment_analysis->numbLaminae[nl];
        //                string FilenamePA_LaplaceLayersPar=FilenamePA_b+"LaplaceLayersPar_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
        //                SaveTmpResult_short(LayersNormDistPar, FilenamePA_LaplaceLayersPar, SegLesBasis, L2S_Par);
        //                if (LayersNormDistPar!=NULL) {
        //                    delete [] LayersNormDistPar;
        //                    LayersNormDistPar=NULL;
        //                }
        //            }
        //
        //            DistanceNormalisedLaplacePar=CreateLong(DistanceNormalisedLaplace_sPar, L2S_Par, numel);
        //            string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplacePar_"+FilenamePA_e+".nii.gz";
        //            nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplacePar, SegLesBasis,numel);
        //            nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
        //            nifti_image_write(NormLaplaceDistImage);
        //            nifti_image_free(NormLaplaceDistImage);
        //            NormLaplaceDistImage=NULL;
        //            delete [] DistanceNormalisedLaplace_sPar;
        //            delete [] Mask_shortPar;
        //            delete [] Border0Par;
        //            delete [] Border1Par;
        //            delete [] ObjectIn_shortPar;
        //            Mask_shortPar=NULL;
        //            Border0Par=NULL;
        //            Border1Par=NULL;
        //            ObjectIn_shortPar=NULL;
        //            DistanceNormalisedLaplace_sPar=NULL;
        //        }
        //        else{
        //            DistanceNormalisedLaplacePar=CopyArray(static_cast<float *>(LaplaceSolImagePar->data), SegLesBasis->nvox);
        //        }
        //
        //        // Perform the calculation/study using Distance Normalised Laplace and the
        ////        vector<float *> MeanGeneralClasses=TreeToAnalyse->GetMeanGeneralClassesVector();
        //        float * ProportionLesVolDistPar=ProportionLesionVolumeDistance(SegLesBasis, DistanceNormalisedLaplacePar, segment_analysis->numbLaminae[0],PixDim, L2S_Par, S2L_Par, numelmaskedLapPar);
        //        int MaxLabelPar=0;
        //        MaxLabelPar=GetMaxLabel(LesionLabel, numel);
        ////        MeanGeneralClasses=TreeToAnalyse->GetMeanGeneralClassesVector();
        //        int * ExtentLesionPar=ExtentLesionNormDist(LesionLabel, DistanceNormalisedLaplacePar, segment_analysis->numbLaminae[0], L2S_Par, S2L_Par, numelmaskedLapPar);
        ////        MeanGeneralClasses=TreeToAnalyse->GetMeanGeneralClassesVector();
        //        float * LayerLesPar=LayersLesIntensityRelation(SegLesBasis, DistanceNormalisedLaplacePar, TreeToAnalyse, segment_analysis);
        //        // Save in text file the results of ProportionLesVolDist and ExtentLesion
        //        string FilenamePA_NormDistResultPar=FilenamePA_b+"NormDistResultPar_"+FilenamePA_e+".txt";
        //        ofstream TxtNormDistFilePar(FilenamePA_NormDistResultPar.c_str());
        //        PrintNormDistResult(ProportionLesVolDistPar, ExtentLesionPar,LayerLesPar, MaxLabelPar, numbmodal, TxtNormDistFilePar,segment_analysis);
        //
        //        // Clearing memory for not needed elements
        //        delete [] L2S_Par;
        //        delete [] S2L_Par;
        ////        delete [] Dim;
        ////        delete [] Shift;
        ////        delete [] PixDim;
        //        delete [] ObjectInPar;
        //        delete [] ObjectOutPar;
        //        delete [] ObjectTotPar;
        //        delete [] DistanceNormalisedLaplacePar;
        //        delete [] ExtentLesionPar;
        //         delete [] LayerLesPar;
        //        delete [] ProportionLesVolDistPar;
        //
        //        L2S_Par=NULL;
        //        S2L_Par=NULL;
        ////        Dim=NULL;
        ////        Shift=NULL;
        ////        PixDim=NULL;
        //        ObjectOutPar=NULL;
        //        ObjectInPar=NULL;
        //        ObjectTotPar=NULL;
        //        DistanceNormalisedLaplacePar=NULL;
        //        ExtentLesionPar=NULL;
        //        LayerLesPar=NULL;
        //        ProportionLesVolDistPar=NULL;
        //
        //        cout<<"Parenchymal solution for Laplace done "<<endl;
        ////       Then take care of the solution with WMDGMSeg instead of complete Parenchyma
        
        float * ObjectTotTmpWMD=AddArray(ToAddWMD, numel);
        float * ObjectTotDWMD=Erosion_bis(ObjectTotTmpWMD, 3, DimVector, 0);
        float * ObjectTotCWMD=Erosion_bis(ObjectTotDWMD, 3, DimVector, 1);
        float * ObjectOutWMD=SubtractArray(ObjectTotCWMD, VentriclesData, numel);
        
        delete [] ObjectTotTmpWMD;
        ObjectTotTmpWMD=NULL;
        bool * ObjectTotWMD=TranscribeArray<float, bool>(ObjectTotCWMD, numel);
        bool * ObjectInWMD=TranscribeArray<float , bool>(VentriclesData, numel);
        bool * MaskDilWMD=TranscribeArray<float, bool>(ObjectTotDWMD, numel);
        
        //        SaveTmpResult(ObjectIn, g_scratch_dir + "/ObjectIn.nii.gz", SegToAnalyse);
        //        SaveTmpResult(ObjectTot, g_scratch_dir + "/ObjectOut.nii.gz", SegToAnalyse);
        
        
        
        int * L2S_WMD=MakeL2S(ObjectTotWMD, Dim);
        int numelmaskedLapWMD=0;
        int * S2L_WMD=MakeS2L(ObjectTotWMD, Dim, numelmaskedLapWMD);
        
        if (LaplaceSolImageWMD==NULL) { // Meaning that there was nothing saved as Laplace solution beforehand
            int * LabelLaplaceWMD=CreateLaplaceLabelling(ObjectInWMD,ObjectTotWMD,Dim);
            //            SaveTmpResult(LabelLaplace, g_scratch_dir + "/TestLabel.nii.gz", SegToAnalyse);
            int numelmaskedbWMD=0;
            int * L2SbWMD=MakeL2S(MaskDilWMD, Dim);
            int * S2LbWMD=MakeS2L(MaskDilWMD, Dim, numelmaskedbWMD);
            int * Label_sWMD=CreateShort(LabelLaplaceWMD, S2LbWMD, numelmaskedbWMD);
            cout<<GetMax<int>(Label_sWMD, numelmaskedbWMD)<<endl;
            
            //            nifti_image * LaplaceTmpWMD=nifti_image_read(g_scratch_dir + "/TestLaplace.nii.gz", 1);
            //            float * LaplaceTmpLong=static_cast<float *>(LaplaceTmp->data);
            //            float * LaplaceSol_short=CreateShort(LaplaceTmpLong, S2Lb, numelmaskedb);
            
            
            float * LaplaceSol_shortWMD=SolvingLaplaceEquation(Label_sWMD, Dim, Shift, PixDim, 5, L2SbWMD, S2LbWMD, numelmaskedbWMD,SegLesBasis);
            float * LaplaceSolWMD=CreateLong(LaplaceSol_shortWMD,L2SbWMD,numel);
            //                        SaveTmpResult(LaplaceSolWMD, g_scratch_dir + "/FinalLaplaceWMD.nii.gz", SegLesBasis);
            // Clearing memory for all that won't be needed afterwards
            delete [] MaskDilWMD;
            delete [] L2SbWMD;
            delete [] S2LbWMD;
            delete [] Label_sWMD;
            delete [] LaplaceSol_shortWMD;
            MaskDilWMD=NULL;
            L2SbWMD=NULL;
            S2LbWMD=NULL;
            Label_sWMD=NULL;
            LaplaceSol_shortWMD=NULL;
            
            
            bool * ObjectIn_shortWMD=CreateShort(ObjectInWMD, S2L_WMD, numelmaskedLapWMD);
            bool * Mask_shortWMD=CreateShort(ObjectTotWMD, S2L_WMD, numelmaskedLapWMD);
            bool * Border0WMD=CreateBorderFromBool_L2S(ObjectIn_shortWMD, Dim, Shift, L2S_WMD, S2L_WMD, numelmaskedLapWMD);
            //            SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", SegToAnalyse, L2S);
            bool * Border1WMD=CreateBorderFromBool_L2S(Mask_shortWMD, Dim, Shift, L2S_WMD, S2L_WMD, numelmaskedLapWMD);
            float * LaplaceSol_sWMD=CreateShort(LaplaceSolWMD, S2L_WMD, numelmaskedLapWMD);
            float * DistanceNormalisedLaplace_sWMD=NormalisedLaplaceLength(LaplaceSol_sWMD, ObjectIn_shortWMD, Mask_shortWMD, L2S_WMD, S2L_WMD, numelmaskedLapWMD, PixDim, Dim, Shift,SegLesBasis);
            int nlsize=segment_analysis->numbLaminae.size();
            for (int nl=0; nl<nlsize; nl++) {
                stringstream nlName;
                nlName << segment_analysis->numbLaminae[nl];
                int * LayersNormDistWMD=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_sWMD,segment_analysis->numbLaminae[nl],numelmaskedLapWMD);
                string FilenamePA_LaplaceLayersWMD=FilenamePA_b+"LaplaceLayersWMD_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                SaveTmpResult_short(LayersNormDistWMD, FilenamePA_LaplaceLayersWMD, SegLesBasis, L2S_WMD);
                if (LayersNormDistWMD!=NULL) {
                    delete [] LayersNormDistWMD;
                    LayersNormDistWMD=NULL;
                }
            }

            DistanceNormalisedLaplaceWMD=CreateLong(DistanceNormalisedLaplace_sWMD, L2S_WMD, numel);
            string FilenamePA_LaplaceWMD=FilenamePA_b+"NormDistLaplaceWMD_"+FilenamePA_e+".nii.gz";
            nifti_image * NormLaplaceDistImageWMD=CreateNiiFromArray(DistanceNormalisedLaplaceWMD, SegLesBasis,numel);
            nifti_set_filenames(NormLaplaceDistImageWMD, FilenamePA_LaplaceWMD.c_str(), 0, 0);
            nifti_image_write(NormLaplaceDistImageWMD);
            nifti_image_free(NormLaplaceDistImageWMD);
            NormLaplaceDistImageWMD=NULL;
            delete [] DistanceNormalisedLaplace_sWMD;
            delete [] Mask_shortWMD;
            delete [] Border0WMD;
            delete [] Border1WMD;
            delete [] ObjectIn_shortWMD;
            Mask_shortWMD=NULL;
            Border0WMD=NULL;
            Border1WMD=NULL;
            ObjectIn_shortWMD=NULL;
            DistanceNormalisedLaplace_sWMD=NULL;
        }
        else{
            DistanceNormalisedLaplaceWMD=CopyArray(static_cast<float *>(LaplaceSolImageWMD->data), SegToAnalyse->nvox);
        }
        
        // Perform the calculation/study using Distance Normalised Laplace and the
        float * ProportionLesVolDistWMD=ProportionLesionVolumeDistance(SegLesBasis, DistanceNormalisedLaplaceWMD, segment_analysis->numbLaminae[0],PixDim, L2S_WMD, S2L_WMD, numelmaskedLapWMD);
        int MaxLabelWMD=0;
        int * ExtentLesionWMD=ExtentLesionNormDist(LesionLabel, DistanceNormalisedLaplaceWMD, segment_analysis->numbLaminae[0], L2S_WMD, S2L_WMD, numelmaskedLapWMD);
        float * LayerLesWMD=LayersLesIntensityRelation(SegLesBasis, DistanceNormalisedLaplaceWMD, TreeToAnalyse, segment_analysis);
        // Save in text file the results of ProportionLesVolDist and ExtentLesion
        string FilenamePA_NormDistResultWMD=FilenamePA_b+"NormDistResultWMD_"+FilenamePA_e+".txt";
        ofstream TxtNormDistFileWMD(FilenamePA_NormDistResultWMD.c_str());
        PrintNormDistResult(ProportionLesVolDistWMD, ExtentLesionWMD,LayerLesWMD, MaxLabelWMD, numbmodal, TxtNormDistFileWMD,segment_analysis);
        
        // Clearing memory for not needed elements
        delete [] L2S_WMD;
        delete [] S2L_WMD;
        delete [] Dim;
        delete [] Shift;
        delete [] PixDim;
        delete [] ObjectInWMD;
        delete [] ObjectOutWMD;
        delete [] ObjectTotWMD;
        delete [] DistanceNormalisedLaplaceWMD;
        delete [] ExtentLesionWMD;
        delete [] ProportionLesVolDistWMD;
        
        L2S_WMD=NULL;
        S2L_WMD=NULL;
        Dim=NULL;
        Shift=NULL;
        PixDim=NULL;
        ObjectOutWMD=NULL;
        ObjectInWMD=NULL;
        ObjectTotWMD=NULL;
        DistanceNormalisedLaplaceWMD=NULL;
        ExtentLesionWMD=NULL;
        ProportionLesVolDistWMD=NULL;
        
        if (LesionLabel!=NULL) {
            delete [] LesionLabel;
        }
        if (SegLesBasis!=NULL) {
            nifti_image_free(SegLesBasis);
        }
        if(WMDGMSeg!=NULL){
            nifti_image_free(WMDGMSeg);
            WMDGMSeg=NULL;
        }
        
    }
    
    
    
    if (segment_analysis->flag_Parc && segment_analysis->flag_inSum) { // Meaning that supposedly there is a file to use as result of parcellation to build ObjectIn, ObjectOut and the resulting distance files
        FilenamePA=nifti_makebasename(segment_analysis->filename_inSum);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        
        
        
        //        First check if there is the Label for the wanted analysis
        nifti_image * SegLesBasis=NULL;
        int * LesionLabel=NULL;
        if (OrderedLabelsCorr==NULL && OrderedLabels==NULL) {
            if (segment_analysis->flag_inLesCorr) {
                SegLesBasis=ReadFromFilename(segment_analysis->filename_inLesCorr);
            }
            else if (segment_analysis->flag_inLes){
                SegLesBasis=ReadFromFilename(segment_analysis->filename_inLes);
            }
            nifti_image * BinarySegForLaplace=HardSegmentationThreshold(SegLesBasis, 0.1);
            int * LesionLabelPrep=ComponentLabeling(BinarySegForLaplace, segment_analysis->Neigh);
            LesionLabel=OrderedVolumeLabel(LesionLabelPrep, 1, numel);
            nifti_image_free(BinarySegForLaplace);
            delete [] LesionLabelPrep;
            LesionLabelPrep=NULL;
        }
        else{
            SegLesBasis=CopyFloatNii(SegToAnalyse);
            if (OrderedLabelsCorr!=NULL) {
                LesionLabel=CopyArray(OrderedLabelsCorr, numel);
            }
            else{
                LesionLabel=CopyArray(OrderedLabels, numel);
            }
        }
        
        
        nifti_image * ParcellationImage=NULL;
        ParcellationImage=nifti_image_read(segment_analysis->filename_Parc,true);
        if (ParcellationImage!=NULL) {
            if (segment_analysis->vecElementIn.size()==0) {
                int numbIns=DefaultElIn.size();
                for (int e=0; e<numbIns; e++) {
                    segment_analysis->vecElementIn.push_back(DefaultElIn[e]);
                }
            }
            if (segment_analysis->vecElementOut.size()==0) {
                int numbOuts=DefaultElOut.size();
                for (int e=0; e<numbOuts; e++) {
                    segment_analysis->vecElementOut.push_back(DefaultElOut[e]);
                }
            }
            
            nifti_image * LaplaceSolImage=NULL;
            float * DistanceNormalisedLaplace=NULL;
            if (segment_analysis->flag_LaplaceSolImage) { // Case where the Solution for the Laplace Normalised distance is already available
                LaplaceSolImage=nifti_image_read(segment_analysis->filename_LaplaceSolImage,true);
            }
            
            bool * ObjectIn=CreateObjectsFromParcelLesSeg(ParcellationImage, SegLesBasis, segment_analysis->vecElementIn, 2);
            bool * ObjectOut=CreateObjectsFromParcelLesSeg(ParcellationImage,SegLesBasis, segment_analysis->vecElementOut,1);
            //                SaveTmpResult(ObjectIn, g_scratch_dir + "/ObjectIn.nii.gz", ParcellationImage);
            //                SaveTmpResult(ObjectOut, g_scratch_dir + "/ObjectOut.nii.gz", ParcellationImage);
            int numel=ParcellationImage->nvox;
            int * Dim=new int[3];
            int * Shift=new int[3];
            float * PixDim=new float[3];
            Shift[0]=1;
            for (int d=0; d<3; d++) {
                Dim[d]=ParcellationImage->dim[d+1];
                PixDim[d]=ParcellationImage->pixdim[d+1];
                Shift[d]=d>0?Shift[d-1]*Dim[d-1]:1;
            }
            
            vector<bool *> ObjectVector;
            ObjectVector.push_back(ObjectOut);
            ObjectVector.push_back(ObjectIn);
            bool * MaskLap=AddArray(ObjectVector, numel);
            int * L2S=MakeL2S(MaskLap, Dim);
            int numelmasked=0;
            int * S2L=MakeS2L(MaskLap, Dim, numelmasked);
            
            if (LaplaceSolImage==NULL) { // Meaning that there was nothing saved as Laplace solution beforehand
                bool * MaskDil=Dilation(MaskLap,3,Dim,Shift);
                int * LabelLaplace=CreateLaplaceLabelling(ObjectIn,MaskLap,Dim);
                //                    SaveTmpResult(LabelLaplace, g_scratch_dir + "/TestLabel.nii.gz", ParcellationImage);
                int numelmaskedb=0;
                int * L2Sb=MakeL2S(MaskDil, Dim);
                int * S2Lb=MakeS2L(MaskDil, Dim, numelmaskedb);
                int * Label_s=CreateShort(LabelLaplace, S2Lb, numelmaskedb);
                cout<<GetMax<int>(Label_s, numelmaskedb)<<endl;
                
                //                nifti_image * LaplaceTmp=nifti_image_read(g_scratch_dir + "/TestLaplace.nii.gz", 1);
                //                float * LaplaceTmpLong=static_cast<float *>(LaplaceTmp->data);
                //                float * LaplaceSol_short=CreateShort(LaplaceTmpLong, S2Lb, numelmaskedb);
                
                
                float * LaplaceSol_short=SolvingLaplaceEquation(Label_s, Dim, Shift, PixDim, 5, L2Sb, S2Lb, numelmaskedb,ParcellationImage);
                float * LaplaceSol=CreateLong(LaplaceSol_short,L2Sb,numel);
                //                    SaveTmpResult(LaplaceSol, g_scratch_dir + "/TestLabel.nii.gz", ParcellationImage);
                // Clearing memory for all that won't be needed afterwards
                delete [] MaskDil;
                delete [] L2Sb;
                delete [] S2Lb;
                delete [] Label_s;
                delete [] LaplaceSol_short;
                MaskDil=NULL;
                L2Sb=NULL;
                S2Lb=NULL;
                Label_s=NULL;
                LaplaceSol_short=NULL;
                
                
                
                bool * ObjectIn_short=CreateShort(ObjectIn, S2L, numelmasked);
                bool * Mask_short=CreateShort(MaskLap, S2L, numelmasked);
                bool * Border0=CreateBorderFromBool_L2S(ObjectIn_short, Dim, Shift, L2S, S2L, numelmasked);
                //                    SaveTmpResult_short(Border0, g_scratch_dir + "/TestBorder.nii.gz", ParcellationImage, L2S);
                bool * Border1=CreateBorderFromBool_L2S(Mask_short, Dim, Shift, L2S, S2L, numelmasked);
                float * LaplaceSol_s=CreateShort(LaplaceSol, S2L, numelmasked);
                float * DistanceNormalisedLaplace_s=NormalisedLaplaceLength(LaplaceSol_s, ObjectIn_short, Mask_short, L2S, S2L, numelmasked, PixDim, Dim, Shift,ParcellationImage);
                int nlsize=segment_analysis->numbLaminae.size();
                for (int nl=0; nl<nlsize; nl++) {
                    stringstream nlName;
                    nlName << segment_analysis->numbLaminae[nl];
                    int * LayersNormDist=CreateLayersFromNormDistSol(DistanceNormalisedLaplace_s,segment_analysis->numbLaminae[nl],numelmasked);
                    string FilenamePA_LaplaceLayers=FilenamePA_b+"LaplaceLayersGIF_"+nlName.str()+"_"+FilenamePA_e+".nii.gz";
                    SaveTmpResult_short(LayersNormDist, FilenamePA_LaplaceLayers, ParcellationImage, L2S);
                    if (LayersNormDist!=NULL) {
                        delete [] LayersNormDist;
                        LayersNormDist=NULL;
                    }
                }
                
                DistanceNormalisedLaplace=CreateLong(DistanceNormalisedLaplace_s, L2S, numel);
                string FilenamePA_Laplace=FilenamePA_b+"NormDistLaplaceGIF_"+FilenamePA_e+".nii.gz";
                nifti_image * NormLaplaceDistImage=CreateNiiFromArray(DistanceNormalisedLaplace, ParcellationImage,numel);
                nifti_set_filenames(NormLaplaceDistImage, FilenamePA_Laplace.c_str(), 0, 0);
                nifti_image_write(NormLaplaceDistImage);
                nifti_image_free(NormLaplaceDistImage);
                NormLaplaceDistImage=NULL;
                delete [] DistanceNormalisedLaplace_s;
                delete [] Mask_short;
                delete [] Border0;
                delete [] Border1;
                delete [] ObjectIn_short;
                Mask_short=NULL;
                Border0=NULL;
                Border1=NULL;
                ObjectIn_short=NULL;
                DistanceNormalisedLaplace_s=NULL;
            }
            else{
                DistanceNormalisedLaplace=CopyArray(static_cast<float *>(LaplaceSolImage->data), ParcellationImage->nvox);
            }
            
            //            No need anymore for parcellation image
            if(ParcellationImage!=NULL){
                nifti_image_free(ParcellationImage);
                ParcellationImage=NULL;
            }
            
            
            // Perform the calculation/study using Distance Normalised Laplace and the
            float * ProportionLesVolDist=ProportionLesionVolumeDistance(SegLesBasis, DistanceNormalisedLaplace, segment_analysis->numbLaminae[0],PixDim, L2S, S2L, numelmasked);
            int MaxLabel=0;
            MaxLabel=GetMaxLabel(LesionLabel, numel);
            int * ExtentLesion=ExtentLesionNormDist(LesionLabel, DistanceNormalisedLaplace, segment_analysis->numbLaminae[0], L2S, S2L, numelmasked);
            float * LayerLes=LayersLesIntensityRelation(SegLesBasis, DistanceNormalisedLaplace, TreeToAnalyse, segment_analysis);
            // Save in text file the results of ProportionLesVolDist and ExtentLesion
            string FilenamePA_NormDistResult=FilenamePA_b+"NormDistResultGIF_"+FilenamePA_e+".txt";
            ofstream TxtNormDistFile(FilenamePA_NormDistResult.c_str());
            PrintNormDistResult(ProportionLesVolDist, ExtentLesion,LayerLes, MaxLabel, numbmodal, TxtNormDistFile,segment_analysis);
            
            // Clearing memory for not needed elements
            delete [] L2S;
            delete [] S2L;
            delete [] Dim;
            delete [] Shift;
            delete [] PixDim;
            delete [] ObjectIn;
            delete [] ObjectOut;
            delete [] MaskLap;
            delete [] DistanceNormalisedLaplace;
            delete [] LesionLabel;
            delete [] ExtentLesion;
            delete [] ProportionLesVolDist;
            
            L2S=NULL;
            S2L=NULL;
            Dim=NULL;
            Shift=NULL;
            PixDim=NULL;
            ObjectOut=NULL;
            ObjectIn=NULL;
            MaskLap=NULL;
            DistanceNormalisedLaplace=NULL;
            LesionLabel=NULL;
            ExtentLesion=NULL;
            ProportionLesVolDist=NULL;
        }
        if (SegLesBasis!=NULL) {
            nifti_image_free(SegLesBasis);
            SegLesBasis=NULL;
        }
    }
    
    
    if(segment_analysis->flag_Quadrant && segment_analysis->flag_inLesCorr){
        cout<<"Doing Quadrant "<<endl;
        FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenameQuadrant=FilenamePA_b+"QuadrantCode_"+FilenamePA_e+".nii.gz";
        
        //        First check if there is the Label for the wanted analysis
        nifti_image * SegLesBasis=NULL;
        //        int * LesionLabel=NULL;
        if (OrderedLabelsCorr==NULL && OrderedLabels==NULL) {
            if (segment_analysis->flag_inLesCorr) {
                SegLesBasis=ReadFromFilename(segment_analysis->filename_inLesCorr);
            }
            else if (segment_analysis->flag_inLes){
                SegLesBasis=ReadFromFilename(segment_analysis->filename_inLes);
            }
            nifti_image * BinarySegForLaplace=HardSegmentationThreshold(SegLesBasis, 0.1);
            //            int * LesionLabelPrep=ComponentLabeling(BinarySegForLaplace, segment_analysis->Neigh);
            //            LesionLabel=OrderedVolumeLabel(LesionLabelPrep, 1, numel);
            nifti_image_free(BinarySegForLaplace);
            //            delete [] LesionLabelPrep;
            //            LesionLabelPrep=NULL;
        }
        else{
            SegLesBasis=CopyFloatNii(SegToAnalyse);
            //            if (OrderedLabelsCorr!=NULL) {
            //                LesionLabel=CopyArray(OrderedLabelsCorr, numel);
            //            }
            //            else{
            //                LesionLabel=CopyArray(OrderedLabels, numel);
            //            }
        }
        
        
        if (segment_analysis->flag_MNITransform && segment_analysis->flag_MNITemplate) {
            nifti_image * MNITemplate=ReadFromFilename(segment_analysis->filename_MNITemplate);
            
            mat44 MNIMatCorrToInd=MNITemplate->qto_ijk;
            float CentreMNI_xyz[4];
            for (int d=0; d<3; d++) {
                CentreMNI_xyz[d]=0;
            }
            CentreMNI_xyz[3]=1;
            
            mat44 ImageMatCorrToWorld=SegLesBasis->qto_xyz;
            float * TransfoMatrix=ReadAffineFromTextFile(segment_analysis->filename_MNITransform);
            float * MNIMatToInd=new float[16];
            TranscribeMat44(MNIMatCorrToInd, MNIMatToInd);
            float * ImageMatToWorld=new float[16];
            //            int SizeMatVec[4];
            //            SizeMatVec[0]=4;
            //            SizeMatVec[1]=4;
            //            SizeMatVec[2]=4;
            //            SizeMatVec[3]=1;
            //            float * CentreMNI_ijk=TreeToAnalyse->ProductMatrix(MNIMatToInd, CentreMNI_xyz, SizeMatVec);
            
            int SizeTransfo[4];
            SizeTransfo[0]=4;
            SizeTransfo[1]=4;
            SizeTransfo[2]=4;
            SizeTransfo[3]=4;
            TranscribeMat44(ImageMatCorrToWorld, ImageMatToWorld);
            float * Transfo1=TreeToAnalyse->ProductMatrix(TransfoMatrix, ImageMatToWorld,SizeTransfo);
            //            float * TransfoTot=TreeToAnalyse->ProductMatrix(MNIMatToInd, Transfo1, SizeTransfo);
            //            Normally the direct multiplication of the coordinates vector will give the corresponding coordinate vector in the MNI template thus allowing for the separation in quadrants
            
            //            Create the quadrant image according to code for quadrants
            nifti_image * QuadrantResult=QuadrantTransformation(Transfo1,CentreMNI_xyz, SegLesBasis,TreeToAnalyse->GetL2S());
            nifti_set_filenames(QuadrantResult, FilenameQuadrant.c_str(), 0, 0);
            nifti_image_write(QuadrantResult);
            cout<<"Quadrant done"<<endl;
            if (QuadrantResult!=NULL) {
                nifti_image_free(QuadrantResult);
                QuadrantResult=NULL;
            }
            if (Transfo1!=NULL) {
                delete [] Transfo1;
                Transfo1=NULL;
            }
            if (TransfoMatrix!=NULL) {
                delete [] TransfoMatrix;
                TransfoMatrix=NULL;
            }
            if (MNITemplate!=NULL) {
                nifti_image_free(MNITemplate);
                MNITemplate=NULL;
            }

            
        }
        if (SegLesBasis!=NULL) {
            nifti_image_free(SegLesBasis);
            SegLesBasis=NULL;
        }
        //        if (LesionLabel !=NULL){
        //            delete [] LesionLabel;
        //            LesionLabel=NULL;
        //        }
    }

    
    
    
    

    if (segment_analysis->flag_Les && SegToAnalyse==NULL) {
        if (segment_analysis->flag_in) {
            SegToAnalyse=nifti_image_read(segment_analysis->filename_InLes,true);
            cout<<"In les read "<<endl;
        }
        else if(segment_analysis->flag_inConnect){
            ConnectLabelLesions=nifti_image_read(segment_analysis->filename_inConnect, true);
            SegToAnalyse=nifti_image_read(segment_analysis->filename_inConnect, true);
        }
        else if(segment_analysis->flag_inLesCorr){
            SegToAnalyse=nifti_image_read(segment_analysis->filename_inLesCorr, true);
            ImageLesionCorr=nifti_image_read(segment_analysis->filename_inLesCorr, true);
        }
    }
    
    
    nifti_image * SegResult=NULL;
    if (segment_analysis->flag_in && SegToAnalyse==NULL) {
        vector<nifti_image* > VectorImageToSegment;
        int sizeVector=segment_analysis->filename_In.size();
        for (int c=0;c<sizeVector;c++){
            VectorImageToSegment.push_back(nifti_image_read(segment_analysis->filename_In[c],true));
        };
        SegResult=CreateDataImage(VectorImageToSegment);
    }
    
    if (segment_analysis->flag_inConnect) { // Case where we input the already connected result for the
        ConnectLabelLesions=ReadFromFilename(segment_analysis->filename_inConnect);
        if(SegToAnalyse==NULL && segment_analysis->flag_inLesCorr==0){
            SegToAnalyse=ReadFromFilename(segment_analysis->filename_inConnect);
        }
        else if(SegToAnalyse==NULL){
            SegToAnalyse=ReadFromFilename(segment_analysis->filename_inLesCorr);
        }
        float * ConnectedInData=static_cast<float*>(ConnectLabelLesions->data);
        if (BinarySeg==NULL) {
            BinarySeg=HardSegmentationThreshold(ConnectLabelLesions, 0.5);
        }
        int numel=ConnectLabelLesions->nx*ConnectLabelLesions->ny*ConnectLabelLesions->nz;
        int * ComponentLabel=TranscribeArray<float, int>(ConnectedInData, numel);
        OrderedLabels=OrderedVolumeLabel(ComponentLabel, 1, numel);
        if (ComponentLabel!=NULL) {
            delete [] ComponentLabel;
            ComponentLabel=NULL;
        }
        
        
        if (ImageLesionCorr==NULL && segment_analysis->flag_inLesCorr) {
            ImageLesionCorr=ReadFromFilename(segment_analysis->filename_inLesCorr);
        }
        
        if(ImageLesionCorr!=NULL){// Supposedly with the connected input (if it is the non corrected version we need the corrected true lesion seg to avoid the corresponding connected components)
            float * SegData=static_cast<float *>(ImageLesionCorr->data);
            int * ConnectedLabelCorr=TranscribeArray<float, int>(ConnectedInData, numel);
            for (int i=0; i<numel ; i++) {
                if (ConnectedLabelCorr[i]>0) {
                    if (SegData[i]==0) {
                        ConnectedLabelCorr[i]=0;
                    }
                }
            }
            OrderedLabelsCorr=OrderedVolumeLabel(ConnectedLabelCorr, 1, numel);
            delete [] ConnectedLabelCorr;
            ConnectedLabelCorr=NULL;
        }
    }

    nifti_image * Mask=NULL;
    int numbActive=0;
//    if(SegToAnalyse !=NULL){
//        numbActive=SegToAnalyse->nx*SegToAnalyse->ny*SegToAnalyse->ny;
//    }
//    if (SegResult!=NULL) {
//        numbActive=SegResult->nx*SegResult->ny*SegResult->ny;
//    }
    if(segment_analysis->flag_mask){
        Mask=nifti_image_read(segment_analysis->filename_mask,true);
    }
    else{
        Mask=NULL;
    }
    
    if (segment_analysis->flag_outConnect) {
        cout<<"Doing connect out part"<<endl;
        nifti_image * FinalConnect=NULL;
        bool * LesBool=NULL;
        float * LesFloat=NULL;
        int Dim[3];
        int Shift[3];
        float PixDim[3];
        vector<int> DimVector;
        if (segment_analysis->flag_inLesCorr) {
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
            nifti_image * SegLesBasis=ReadFromFilename(segment_analysis->filename_inLesCorr);
            int * ConnectToSave=ComponentLabeling(SegLesBasis, segment_analysis->Neigh);
            int * OrderedToSave=OrderedVolumeLabel(ConnectToSave, 1, numel);
            
            FinalConnect=CreateNiiFromArray(OrderedToSave, SegLesBasis, numel);
            LesBool=TranscribeArray<int, bool>(OrderedToSave, numel);
            
            delete [] OrderedToSave;
            delete [] ConnectToSave;
            nifti_image_free(SegLesBasis);
            OrderedToSave=NULL;
            SegLesBasis=NULL;
            ConnectToSave=NULL;
            
            
        }
        else if ((OrderedLabels !=NULL || OrderedLabelsCorr !=NULL) && SegToAnalyse !=NULL) {
            if (segment_analysis->flag_TextFile) {
                FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
            }
            if (OrderedLabelsCorr !=NULL) {
                FinalConnect=CreateNiiFromArray(OrderedLabelsCorr, SegToAnalyse, numel);
                LesBool=TranscribeArray<int, bool>(OrderedLabelsCorr, numel);
            }
            else{
                FinalConnect=CreateNiiFromArray(OrderedLabels, SegToAnalyse, numel);
                LesBool=TranscribeArray<int, bool>(OrderedLabels, numel);
            }
            
        }
        else if (segment_analysis->flag_inLes) {
            FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
            nifti_image * SegLesBasis=ReadFromFilename(segment_analysis->filename_inLesCorr);
            nifti_image * SegLesThresh=HardSegmentationThreshold(SegLesBasis, 0.2);
            int * ConnectToSave=ComponentLabeling(SegLesThresh, segment_analysis->Neigh);
            int * OrderedToSave=OrderedVolumeLabel(ConnectToSave, 1, numel);
            LesBool=TranscribeArray<int, bool>(OrderedToSave, numel);
            FinalConnect=CreateNiiFromArray(OrderedToSave, SegLesBasis, numel);
            delete [] OrderedToSave;
            nifti_image_free(SegLesBasis);
            nifti_image_free(SegLesThresh);
            OrderedToSave=NULL;
            SegLesBasis=NULL;
            SegLesThresh=NULL;
        }
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenamePA_ConnectLesion=FilenamePA_b+"ConnectFin_"+FilenamePA_e+".nii.gz";
        string FilenamePA_BorderLesion=FilenamePA_b+"BorderLes_"+FilenamePA_e+".nii.gz";
        string FilenamePA_ExtLesion=FilenamePA_b+"ExtLes_"+FilenamePA_e+".nii.gz";
        string FilenamePA_ConnectExt=FilenamePA_b+"ConnectExt_"+FilenamePA_e+".nii.gz";
        if (FinalConnect!=NULL) {
            for (int d=0; d<3; d++) {
                Dim[d]=FinalConnect->dim[d+1];
                DimVector.push_back(Dim[d]);
                PixDim[d]=FinalConnect->pixdim[d+1];
            }
            
            Shift[0]=1;
            Shift[1]=Dim[0];
            Shift[2]=Dim[1]*Shift[1];
            LesFloat=TranscribeArray<bool, float>(LesBool, numel);
            float * ExtendedLesion=AnisotropicMorphologicalChange(LesFloat, 3, DimVector, 0, PixDim);
            nifti_image * ExtLesion=CreateNiiFromArray(ExtendedLesion, FinalConnect, numel);
            int * ConnectExt=ComponentLabeling(ExtLesion, segment_analysis->Neigh);
            int * OrderedConnectExt=OrderedVolumeLabel(ConnectExt, 1, numel);
            nifti_image * ConnectExtImage=CreateNiiFromArray(OrderedConnectExt, FinalConnect, numel);
            nifti_set_filenames(ConnectExtImage, FilenamePA_ConnectExt.c_str(), 0, 0);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(ConnectExtImage);
            }
            
            nifti_image_free(ConnectExtImage);
            nifti_set_filenames(ExtLesion, FilenamePA_ExtLesion.c_str(), 0, 0);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(ExtLesion);
            }
            
            nifti_image_free(ExtLesion);
            ExtLesion=NULL;
            bool * Border=CreateBorderFromBool(LesBool, Dim, Shift);
            nifti_image * BorderConnect=CreateNiiFromArray(Border, FinalConnect, numel);
            nifti_set_filenames(BorderConnect, FilenamePA_BorderLesion.c_str(), 0, 0);
            nifti_set_filenames(FinalConnect, FilenamePA_ConnectLesion.c_str(), 0, 0);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(FinalConnect);
            }
            
            nifti_image_free(FinalConnect);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(BorderConnect);
            }
            
            nifti_image_free(BorderConnect);
            if(LesBool!=NULL){
                delete [] LesBool;
                LesBool=NULL;
            }
            delete [] Border;
            Border=NULL;
            delete [] LesFloat;
            LesFloat=NULL;
            delete [] ExtendedLesion;
            ExtendedLesion=NULL;
            delete [] ConnectExt;
            ConnectExt=NULL;
            delete [] OrderedConnectExt;
            OrderedConnectExt=NULL;
        }
    }
    
    if (segment_analysis->flag_outWMI) {// Consists in getting the WMI segmentation (binary) and the WM deprived of the lesions. Note that the first one is only possible to get if we have the complete tree and not only Sum and Lesion
        cout<<"Trying WM segmentation"<<endl;
        //        Create FilenamePA according to Tree or Summarised
        if (segment_analysis->flag_inSum) {
            FilenamePA=nifti_makebasename(segment_analysis->filename_inSum);
            if (SummarisedSeg1==NULL) {
                SummarisedSeg1=ReadFromFilename(segment_analysis->filename_inSum);
            }
        }
        else if (segment_analysis->flag_TextFile) {
            FilenamePA=nifti_makebasename(segment_analysis->filename_SegTot);
        }
        int Index=FilenamePA.find_last_of('/');
        FilenamePA_b=FilenamePA.substr(0,Index+1);
        if(segment_analysis->flag_outputDir){
            FilenamePA_b=segment_analysis->name_outputDir;
        }
        FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
        string FilenamePA_WMISeg=FilenamePA_b+"WMISeg_"+FilenamePA_e+".nii.gz";
        string FilenamePA_WMmLSeg=FilenamePA_b+"WMmLSeg_"+FilenamePA_e+".nii.gz";
        string FilenamePA_WMSeg=FilenamePA_b+"WMSeg_"+FilenamePA_e+".nii.gz";
        //        Depends on what data is available
        //        Case 1: Summarised and Lesion are availabe
        //        Case 2: TreeToAnalyse and Lesion are available
        //        Check if Lesion data available
        nifti_image * SegLesBasis=NULL;
        if (SegToAnalyse!=NULL) {
            SegLesBasis=CopyFloatNii(SegToAnalyse);
        }
        else if (segment_analysis->flag_inLesCorr){
            SegLesBasis=ReadFromFilename(segment_analysis->filename_inLesCorr);
        }
        else if(segment_analysis->flag_inLes){
            SegLesBasis=ReadFromFilename(segment_analysis->filename_inLes);
        }
        else{
            cout<<"no lesion data a priori available"<<endl;
        }
        nifti_image * WMIHardSeg=NULL;
        nifti_image * WMHardSeg=NULL;
        nifti_image * WMmLHardSeg=NULL;
        int * L2SToUse=NULL;
        bool flag_createL2S=0;
        if (TreeToAnalyse!=NULL) {
            L2SToUse=TreeToAnalyse->GetL2S();
        }
        else if(Mask!=NULL){
            bool * MaskImage = static_cast<bool *>(Mask->data);
            int Dim[3];
            for (int d=0; d<3; d++) {
                Dim[d]=SegLesBasis->dim[d+1];
            }
            L2SToUse=MakeL2S(MaskImage , Dim);
            flag_createL2S=1;
        }
        else{
            L2SToUse=new int[numel];
            for (int i=0; i<numel; i++) {
                L2SToUse[i]=i;
            }
            flag_createL2S=1;
        }
        
        
        if (SegLesBasis!=NULL) {
            cout<<"Lesion Exists"<<endl;
            int Dim[3];
            int Shift[3];
            float PixDim[3];
            vector<int> DimVector;
            string Options="";
            if (segment_analysis->inOptionText!=NULL){
                Options=segment_analysis->inOptionText;
            }
            if (segment_analysis->flag_segWeighted>=1) {
                Options+="WS";
                stringstream ws;
                ws << segment_analysis->weightThreshold;
                Options+=ws.str();
                stringstream wt;
                wt <<segment_analysis->flag_segWeighted;
                Options+="WT"+wt.str();
                stringstream wc;
                wc << segment_analysis->weightCompClass;
                Options+="WC"+wc.str();
            }

            for (int d=0; d<3; d++) {
                Dim[d]=SegLesBasis->dim[d+1];
                DimVector.push_back(Dim[d]);
                PixDim[d]=SegLesBasis->pixdim[d+1];
            }
            Shift[0]=1;
            Shift[1]=Dim[0];
            Shift[2]=Dim[1]*Shift[1];
            if (TreeToAnalyse!=NULL) {
                //                float * MaskValues=CreateLong(TreeToAnalyse->GetNormResp(),L2S,numel);
                //                nifti_image * MaskFloat=CreateNiiFromArray(MaskValues,SegLesBasis,numel);
                //                float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
                //                nifti_image * DataCompImage=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
                //                delete [] DataCorr;
                //                nifti_image * WMSeg=HardSegmentationThresholdFromNormResp(TreeToAnalyse->GetChild(0)->GetChild(1)->GetNormResp(),TreeToAnalyse,0.5);
                //                nifti_image * MahalMap=MahalDistMaps(WMSeg,MaskFloat,DataCompImage);
                if(segment_analysis->flag_CST){
                    cout << "Doing CST correction" <<endl;
                    CorrectionCST(TreeToAnalyse,SegToAnalyse,segment_analysis);
                }
                cout << "CST correction done"<<endl;
                float * DilatedSeg = Erosion_bis(static_cast<float*>(SegLesBasis->data),3,DimVector,0);
                nifti_image * DilatedInit=CreateNiiFromArray(DilatedSeg,SegLesBasis,numel);
                nifti_image * DilatedNii = HardSegmentationThreshold(DilatedInit,0.01);
                vector<TreeEM*> LesionClasses;
                nifti_image_free(DilatedInit);
                nifti_image * FinalWeight = ReconstructLesionImageTotWeighted(TreeToAnalyse,LesionClasses,SegToAnalyse,segment_analysis->vecModalities,segment_analysis);
                nifti_image_free(DilatedNii);
                nifti_image * ThreshFinalWeight = HardSegmentationThreshold(FinalWeight,0.5);
                int * ComponentLabel=ComponentLabeling(ThreshFinalWeight,segment_analysis->Neigh);
                float VolumeVox=FinalWeight->pixdim[1]*FinalWeight->pixdim[2]*FinalWeight->pixdim[3];
                int *OrderedLabels=OrderedVolumeLabel(ComponentLabel,segment_analysis->MiniSize,numel,VolumeVox);
                nifti_image * OrderedNii=CreateNiiFromArray(OrderedLabels,SegLesBasis,numel);
                vector<Outlier*> OutlierCheck=GetVectorOutliers(OrderedLabels,TreeToAnalyse,segment_analysis);
                FilenamePA=nifti_makebasename(segment_analysis->filename_inLesCorr);
                int Index=FilenamePA.find_last_of('/');
                FilenamePA_b=FilenamePA.substr(0,Index+1);
                if(segment_analysis->flag_outputDir){
                    FilenamePA_b=segment_analysis->name_outputDir;
                }
                FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
                string FilenameResultsTxt=FilenamePA_b+"TxtLesion_"+Options+FilenamePA_e+".txt";
                string FilenameConnectNii=FilenamePA_b+"Connect_"+Options+FilenamePA_e+".nii.gz";
                string FilenameCorrect=FilenamePA_b+"Correct_"+Options+FilenamePA_e+".nii.gz";
                string FilenameInfarct=FilenamePA_b+"InfarctCorrect_"+Options+FilenamePA_e+".nii.gz";

                nifti_set_filenames(OrderedNii,FilenameConnectNii.c_str(),0,0);
                nifti_image_write(OrderedNii);
                ofstream TxtFileLesion(FilenameResultsTxt.c_str());
                int numbLabels=OutlierCheck.size();
                int*GravImage=NULL;
                if (numbLabels>0) {
                    GravImage=new int[3];
                    for (int d=0;d<3; d++) {
                        GravImage[d]=(OutlierCheck[0]->CentreGravity[d]-OutlierCheck[0]->VectorDiffGrav[d]/PixDim[d]);
                    }
                }
//                int Grav=GetCenterGravity(static_cast<bool*>(TreeToAnalyse->GetMask()->data),Dim);
                cout<<" Preparing printing..."<<endl;
                nifti_image * ResultingCorrect=CopyFloatNiiImage(ThreshFinalWeight);
                nifti_set_filenames(ResultingCorrect,FilenameCorrect.c_str(),0,0);
                float *ResultData=static_cast<float*>(ResultingCorrect->data);
                float * InfarctDetected = new float[numel];
                for(int i=0;i<numel;i++){
                    InfarctDetected[i] =0;
                }
                for(int l=0;l<numbLabels;l++){
                    TxtFileLesion<<"Label "<<l+1<<endl;
                    PrintOutlierCharacteristics(OutlierCheck[l],TxtFileLesion,TreeToAnalyse->GetNumberModalities(),GravImage,PixDim,segment_analysis);
                    TxtFileLesion<<endl;
                    if(OutlierCheck[l]->OutlierClass<0 && OutlierCheck[l]->OutlierClass>-19){
                        for(int i=0;i<numel;i++){
                            if(OrderedLabels[i]==l+1 || OrderedLabels[i]==0){
                                    ResultData[i]=0;
                            }
                        }
                    }
                    else if(OutlierCheck[l]->OutlierClass==-20){
                        for(int i=0;i<numel;i++){
                            if(OrderedLabels[i]==l+1){
                                    ResultData[i]=ResultData[i]>0.999?1:0;
                            }
                        }
                    }
                    else if(OutlierCheck[l]->OutlierClass==7){
                        for(int i=0;i<numel;i++){
                            if (OrderedLabels[i]==l+1){
                                InfarctDetected[i]=ResultData[i];
                            }
                        }
                    }
                }
                if(segment_analysis->flag_infarcts){
                nifti_image * InfarctNii=CreateNiiFromArray(InfarctDetected, ResultingCorrect,numel);
                nifti_set_filenames(InfarctNii, FilenameInfarct.c_str(),0,0);
                nifti_image_write(InfarctNii);
                nifti_image_free(InfarctNii);
                delete[]InfarctDetected;
                }
                nifti_image_write(ResultingCorrect);
                if (SummarisedSeg1!=NULL) {
                    WMHardSeg=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexWM, L2SToUse);
                }
                vector<TreeEM *> GeneralClasses=TreeToAnalyse->GetGeneralClassesVector();
                int * L2S=TreeToAnalyse->GetL2S();
                float * NormWMI=GeneralClasses[segment_analysis->IndexWM]->GetNormResp();
                vector<TreeEM *> OutlierNodes=TreeToAnalyse->GetNodeOutlier()->GetChildren();
                float * NormWMO=OutlierNodes[segment_analysis->IndexWM]->GetNormResp(); // We assume a symmetric model here
                
                WMIHardSeg=nifti_copy_nim_info(SegLesBasis);
                WMIHardSeg->data=(void *) calloc(numel, sizeof(float));
                float * WMIHardSegData=static_cast<float *>(WMIHardSeg->data);
                WMmLHardSeg=nifti_copy_nim_info(SegLesBasis);
                WMmLHardSeg->data=(void *) calloc(numel, sizeof(float));
                float * WMmLHardSegData=static_cast<float *>(WMmLHardSeg->data);
                float * LesionData=static_cast<float *>(SegLesBasis->data);
                if (SummarisedSeg1!=NULL) {
                    for (int i=0; i<numel; i++) {
                        if (L2S[i]>=0) {
                            if (LesionData[i]+NormWMI[L2S[i]]+NormWMO[L2S[i]]>0.5) {
                                if (LesionData[i]<=NormWMI[L2S[i]] && NormWMO[L2S[i]]<=NormWMI[L2S[i]]) {
                                    WMIHardSegData[i]=1;
                                    WMmLHardSegData[i]=1;
                                }
                                else if(LesionData[i]<=NormWMI[L2S[i]]+NormWMO[L2S[i]]){
                                    WMmLHardSegData[i]=1;
                                    WMIHardSegData[i]=0;
                                }
                                else{
                                    WMIHardSegData[i]=0;
                                    WMmLHardSegData[i]=0;
                                }
                            }
                            else {
                                WMIHardSegData[i]=0;
                                WMmLHardSegData[i]=0;
                            }
                        }
                        else{
                            WMIHardSegData[i]=0;
                            WMmLHardSegData[i]=0;
                        }
                    }

                }
                else{
                    WMHardSeg=nifti_copy_nim_info(SegLesBasis);
                    WMHardSeg->data=(void *) calloc(numel, sizeof(float));
                    float * WMHardSegData=static_cast<float *>(WMmLHardSeg->data);
                    for (int i=0; i<numel; i++) {
                        if (L2SToUse[i]>=0) {
                            if (LesionData[i]+NormWMI[L2SToUse[i]]+NormWMO[L2SToUse[i]]>0.5) {
                                
                                if (LesionData[i]<=NormWMI[L2SToUse[i]] && NormWMO[L2SToUse[i]]<=NormWMI[L2SToUse[i]]) {
                                    WMIHardSegData[i]=1;
                                    WMmLHardSegData[i]=1;
                                    WMHardSegData[i]=1;
                                }
                                else if(LesionData[i]<=NormWMI[L2SToUse[i]]+NormWMO[L2SToUse[i]]){
                                    WMHardSegData[i]=1;
                                    WMmLHardSegData[i]=1;
                                    WMIHardSegData[i]=0;
                                }
                                else{
                                    WMIHardSegData[i]=0;
                                    WMmLHardSegData[i]=0;
                                    if (LesionData[i]>NormWMI[L2SToUse[i]]) {
                                        WMHardSegData[i]=1;
                                    }
                                    else{
                                        WMHardSegData[i]=0;
                                    }
                                }
                            }
                            else {
                                WMIHardSegData[i]=0;
                                WMmLHardSegData[i]=0;
                                WMHardSegData[i]=0;
                            }
                        }
                        else{
                            WMIHardSegData[i]=0;
                            WMHardSegData[i]=0;
                            WMmLHardSegData[i]=0;
                        }
                    }
                }

            }
            else if (SummarisedSeg1!=NULL) {
                WMmLHardSeg=WMIHardSegmentation(SummarisedSeg1,SegLesBasis,segment_analysis,L2SToUse);
                WMHardSeg=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexWM, L2SToUse);
                if (flag_createL2S) {
                    delete [] L2SToUse;
                    L2SToUse=NULL;
                }
            }
        }

        else { // Case where there is no lesion
            if (TreeToAnalyse!=NULL) {
                if (SummarisedSeg1!=NULL) {
                    WMHardSeg=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexWM, L2SToUse);
                    vector<TreeEM *> GeneralClasses=TreeToAnalyse->GetGeneralClassesVector();
                    float * NormWMI=GeneralClasses[segment_analysis->IndexWM]->GetNormResp();
                    vector<TreeEM *> OutlierNodes=TreeToAnalyse->GetOutliersMainNodesVector();
                    float * NormWMO=OutlierNodes[segment_analysis->IndexWM]->GetNormResp(); // We assume a symmetric model here
                    WMIHardSeg=nifti_copy_nim_info(SegLesBasis);
                    WMIHardSeg->data=(void *) calloc(numel, sizeof(float));
                    float * WMIHardSegData=static_cast<float *>(WMIHardSeg->data);
                    float * WMHardSegData=static_cast<float *>(WMHardSeg->data);
                    for (int i=0; i<numel; i++) {
                        if (L2SToUse[i]>=0 && WMHardSegData[i]>0) {
                            if (NormWMI[L2SToUse[i]]>NormWMO[L2SToUse[i]]) {
                                WMIHardSegData[i]=1;
                            }
                            else{
                                WMIHardSegData[i]=0;
                            }
                        }
                        else{
                            WMIHardSegData[i]=0;
                        }
                    }
                }
                else{ // When SummarisedSeg1 is not available the segmentation is based on a threshold instead of being properly done and we do not correct for other errors
                    vector<TreeEM *> GeneralClasses=TreeToAnalyse->GetGeneralClassesVector();
                    float * NormWMI=GeneralClasses[segment_analysis->IndexWM]->GetNormResp();
                    vector<TreeEM *> OutlierNodes=TreeToAnalyse->GetOutliersMainNodesVector();
                    float * NormWMO=OutlierNodes[segment_analysis->IndexWM]->GetNormResp(); // We assume a symmetric model here
                    WMIHardSeg=nifti_copy_nim_info(SegLesBasis);
                    WMIHardSeg->data=(void *) calloc(numel, sizeof(float));
                    float * WMIHardSegData=static_cast<float *>(WMIHardSeg->data);
                    WMHardSeg=nifti_copy_nim_info(SegLesBasis);
                    WMHardSeg->data=(void *) calloc(numel, sizeof(float));
                    float * WMHardSegData=static_cast<float *>(WMHardSeg->data);
                    for (int i=0; i<numel; i++) {
                        if (L2SToUse[i]>=0 && NormWMO[L2SToUse[i]]+NormWMI[L2SToUse[i]]>0.5) {
                            WMHardSegData[i]=1;
                            if (NormWMI[L2SToUse[i]]>NormWMO[L2SToUse[i]]) {
                                WMIHardSegData[i]=1;
                            }
                            else{
                                WMIHardSegData[i]=0;
                            }
                        }
                        else{
                            WMHardSegData[i]=0;
                            WMIHardSegData[i]=0;
                        }
                    }
                }
            }
            else if (SummarisedSeg1!=NULL) {
                WMHardSeg=HardSegmentationIndex(SummarisedSeg1, segment_analysis->IndexWM, L2SToUse);
            }
        }
        if (flag_createL2S) {
            delete [] L2SToUse;
            L2SToUse=NULL;
        }
        if (WMHardSeg!=NULL) {
            nifti_set_filenames(WMHardSeg, FilenamePA_WMSeg.c_str(), 0, 0);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(WMHardSeg);
            }
            nifti_image_free(WMHardSeg);
            WMHardSeg=NULL;
        }
        if (WMmLHardSeg!=NULL) {
            nifti_set_filenames(WMmLHardSeg, FilenamePA_WMmLSeg.c_str(), 0, 0);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(WMmLHardSeg);
            }
            nifti_image_free(WMmLHardSeg);
            WMmLHardSeg=NULL;
        }
        if (WMIHardSeg!=NULL) {
            nifti_set_filenames(WMIHardSeg, FilenamePA_WMISeg.c_str(), 0, 0);
            if (segment_analysis->flag_Saving) {
                nifti_image_write(WMIHardSeg);
            }
            nifti_image_free(WMIHardSeg);
            WMIHardSeg=NULL;
        }
        if (SegLesBasis!=NULL) {
            nifti_image_free(SegLesBasis);
            SegLesBasis=NULL;
        }
    }
    
    
    
    if(segment_analysis->flag_Les && segment_analysis->flag_refLes && BinarySeg==NULL){ // case where we analyse the lesion result and we have no tree to build the comparison
        if(SegResult!=NULL)
        {
            BinarySeg=HardSegmentationCompare(SegToAnalyse,SegResult,1,Mask);

        }
        else{
            BinarySeg=HardSegmentationThreshold(SegToAnalyse,0.5);
        }
        numel=SegToAnalyse->nx*SegToAnalyse->ny*SegToAnalyse->nz;
        if (segment_analysis->flag_connect && OrderedLabels==NULL) {
            int * ComponentLabel=ComponentLabelingNotRefined(BinarySeg,segment_analysis->Neigh);
            cout<<"Non 0 in Binary is "<<CountNonZero(static_cast<float *>(BinarySeg->data), numel);
            cout<<"Non 0 in SegAnalyse is "<<CountNonZero(static_cast<float *>(SegToAnalyse->data), numel);
            cout<<"Non 0 in tested is "<<CountNonZero(ComponentLabel, numel);
            cout<<"Connected components labeled done"<<endl;
            float * ComponentLabelFloat=TranscribeArray<int, float>(ComponentLabel, numel);
            ConnectLabelLesions=CreateNiiFromArray(ComponentLabelFloat, SegToAnalyse,numel);
            delete [] ComponentLabelFloat;
            
            ComponentLabelFloat=NULL;
            //            OrderedLabels=OrderedVolumeLabel(ComponentLabel, 3, numel);
            //            WARNING OTHER METHODS MIGHT ALLOW 1 VOXEL TO STAND ALONE FOR THE ANALYSIS
            OrderedLabels=OrderedVolumeLabel(ComponentLabel, 1, numel);
            cout<<"Non zero seg Tested is "<<CountNonZero(OrderedLabels, numel);
            delete [] ComponentLabel;
            ComponentLabel=NULL;
            if (!segment_analysis->flag_TextFile && segment_analysis->flag_Les) {
                string Filename=nifti_makebasename(segment_analysis->filename_InLes);
                int Index=Filename.find_last_of('/');
                string Filename_b=Filename.substr(0,Index+1);
                string Filename_e=Filename.substr(Index+1,Filename.length());
                string FilenameConnected=Filename_b+"Connect_"+Filename_e+".nii.gz";
                nifti_set_filenames(ConnectLabelLesions, FilenameConnected.c_str(), 0, 0);
                nifti_image_write(ConnectLabelLesions);
            }
        }
    }
    
    
    if (segment_analysis->flag_ConnectRefAnalysis && OrderedLabels!=NULL && OrderedLabelsGT!=NULL) {
        int numel=ConnectLabelLesions->nx*ConnectLabelLesions->ny*ConnectLabelLesions->nz;
        int maxLabelGT=GetMaxLabel(OrderedLabelsGT, numel);
        int * ResultAnalysis=AnalysisPerRefConnected(OrderedLabels, OrderedLabelsGT, numel);
        ofstream TxtFileRefConnect(FilenameRefConnectAnalysis.c_str());
        PrintAnalysisPerRefConnected(ResultAnalysis, maxLabelGT, TxtFileRefConnect);
        if (OrderedLabelsCorr!=NULL) {
            int * ResultAnalysisCorr=AnalysisPerRefConnected(OrderedLabelsCorr, OrderedLabelsGT, numel);
            ofstream TxtFileRefConnectCorr(FilenameRefConnectCorrAnalysis.c_str());
            PrintAnalysisPerRefConnected(ResultAnalysisCorr, maxLabelGT, TxtFileRefConnectCorr);
        }
        if ((segment_analysis->flag_TextFile || segment_analysis->flag_inDataComp) && (segment_analysis->flag_inWMSeg || segment_analysis->flag_TextFile) ) { // Case when we also perform the intensity study of the lesions
            //            TreeEM * TreeTmp=NULL;
            //            nifti_image * DataCompImage=NULL;
            //            nifti_image * WMSeg=NULL;
            //            if (segment_analysis->flag_TextFile) {
            //                TreeTmp=ReadTreeFromTextFile(segment_analysis->filename_TextFile, NULL, segment_analysis->filename_changePath);
            //            }
            
            //            Get WM True Segmentation
            if (segment_analysis->flag_inWMSeg) {
                WMSeg=ReadFromFilename(segment_analysis->filename_inWMSeg);
            }
            else{ // If we do not have the true segmentation we only consider the WM inlier in the following
                vector<TreeEM *> GeneralClassesVector=TreeToAnalyse->GetGeneralClassesVector();
                float * ShortWM=GeneralClassesVector[segment_analysis->IndexWM]->GetNormResp();
                float * LongWM=CreateLong(ShortWM, TreeToAnalyse->GetL2S(), TreeToAnalyse->GetNumberElements());
                WMSeg=CreateNiiFromArray(LongWM, ConnectedGTImage,numel);
            }
            //            Get the Data we have to compare to (either as part of tree or if filename provided as image under filename. Note that consequently the number of modalities might not be the same in one case and another...)
            if (segment_analysis->flag_inDataComp) {
                DataCompImage=ReadFromFilename(segment_analysis->filename_inDataComp);
            }
            else{ // Extract it from TreeTmp
                numbmodal=TreeToAnalyse->GetNumberModalities();
                float * DataCorr=CreateLongPaddingMulti<float>(TreeToAnalyse->GetDataBFCorrected(), 0, TreeToAnalyse->GetL2S(), numel, numbmodal);
                DataCompImage=CreateNiiFromArray(DataCorr, TreeToAnalyse->GetDataImage(), numel*numbmodal);
                delete [] DataCorr;
                //                DataCompImage=CopyFloatNii(TreeToAnalyse->GetDataImage());
            }
            
            //            float * IntensityStudy=AnalysisPerRefConnectedIntensity(OrderedLabels, OrderedLabelsGT, WMSeg,DataCompImage);
            //            ofstream TxtFileIntensityConnect(FilenameIntensityConnectAnalysis);
            //            PrintAnalysisPerRefConnectedIntensity(IntensityStudy, maxLabelGT, TxtFileIntensityConnect);
            //            if (OrderedLabelsCorr!=NULL) {
            //                float * IntensityStudyCorr=AnalysisPerRefConnectedIntensity(OrderedLabelsCorr, OrderedLabelsGT, WMSeg,DataCompImage);
            //                ofstream TxtFileIntensityConnectCorr(FilenameRefConnectCorrAnalysis);
            //                PrintAnalysisPerRefConnectedIntensity(IntensityStudyCorr, maxLabelGT, TxtFileIntensityConnectCorr);
            //            }
            //
            //        }
            //            if (TreeTmp!=NULL) {
            //                delete TreeTmp;
            //                TreeTmp=NULL;
            //            }

        }
    }
    
    nifti_image * BFImage=NULL;
    if (segment_analysis->flag_BF) {
        vector<nifti_image* > VectorImageToSegment;
        int sizeVectorBF=segment_analysis->filename_BF.size();
        for (int c=0;c<sizeVectorBF;c++){
            VectorImageToSegment.push_back(nifti_image_read(segment_analysis->filename_BF[c],true));
        }
        BFImage=CreateDataImage(VectorImageToSegment);
    }
    nifti_image* DataCorrected=NULL;
    if (segment_analysis->flag_data) {
        vector<nifti_image*> VectorImageToSegment;
        int sizeVectorDC=segment_analysis->filename_Data.size();
        for (int c=0;c<sizeVectorDC;c++){
            VectorImageToSegment.push_back(nifti_image_read(segment_analysis->filename_Data[c],true));
        }
        nifti_image * DataImage=CreateDataImage(VectorImageToSegment);
        DataImage=CreateBFCorrectedImage(DataImage, BFImage);
        DataCorrected=NormaliseImage(DataImage,Mask);
        //        nifti_set_filenames(DataCorrected, (g_scratch_dir + "/TestNormalised.nii.gz").c_str(), 0, 0);
        //        nifti_image_write(DataCorrected);
        //        nifti_image_free(DataImage);
        //        DataImage=NULL;
    }
    
    nifti_image*Outliers=NULL;
    if (segment_analysis->flag_outliers) {
        vector<nifti_image* > VectorImageToSegment;
        int sizeVectorO=segment_analysis->filename_Outliers.size();
        for (int c=0;c<sizeVectorO;c++){
            VectorImageToSegment.push_back(nifti_image_read(segment_analysis->filename_Outliers[c],true));
        }
        Outliers=CreateDataImage(VectorImageToSegment);
    }
    
    //    if(!segment_analysis->flag_outTxt){
    //        string FilenamePA=nifti_makebasename(segment_analysis->filename_InLes);
    //        int Index=FilenamePA.find_last_of('/');
    //        string FilenamePA_b=FilenamePA.substr(0,Index+1);
    //        string FilenamePA_e=FilenamePA.substr(Index+1,FilenamePA.length());
    //        FilenamePA=FilenamePA_b+"SegAnalysis_"+FilenamePA_e+".txt";
    //    }
    //    ofstream TxtFile(segment_analysis->filename_OutTxt);
    
    if (segment_analysis->flag_KLD) {
        float KLD=MakeKLDTot(DataCorrected, Mask, SegResult, Outliers);
        vector<float> KLDVector=MakeKLDVector(DataCorrected, Mask, SegResult, Outliers);
        cout<< "Global KLD "<<KLD<<endl;
        cout<<"Modality specific KLD";
        int sizeKLD=KLDVector.size();
        for (int m=0; m<sizeKLD; m++) {
            cout<<" "<<KLDVector[m];
        }
        cout<<endl;
        TxtFile << "KLD "<<KLD;
        for (int m=0; m<sizeKLD; m++) {
            TxtFile<<" "<<KLDVector[m];
        }
        TxtFile<<endl;
    }
    


    
    MaskImage(SegToAnalyse,Mask);
    int Dim[3];
    Dim[0]=SegToAnalyse->nx;
    Dim[1]=SegToAnalyse->ny;
    Dim[2]=SegToAnalyse->nz;
    bool * MaskData=NULL;
    if (Mask!=NULL) {
        MaskData=static_cast<bool*>(Mask->data);
    }
    //    int * L2SAnalysis=MakeL2S(MaskData, Dim);
    int numelmasked=0;
    int * S2LAnalysis=MakeS2L(MaskData, Dim, numelmasked);
    numbActive=NumbActiveVoxels(SegToAnalyse);
    nifti_image * SegToCompare=NULL;
    if (SegToAnalyse !=NULL) {
        if (segment_analysis->flag_refLes) {
            SegToCompare=nifti_image_read(segment_analysis->filename_Ref, true);
            Floatisation(SegToCompare);
        }
        if(segment_analysis->flag_Les && segment_analysis->flag_refLes){ // case where we analyse the lesion result
            if(BinarySeg==NULL){
                if(SegResult!=NULL)
                {
                    BinarySeg=HardSegmentationCompare(SegToAnalyse,SegResult,1,Mask);

                }
                else{
                    BinarySeg=HardSegmentationThreshold(SegToAnalyse,0.5);
                }
            }
            nifti_image * BinarySegCorr=NULL;
            if(ImageLesionCorr!=NULL){
                //         BinarySegCorr=HardSegmentation_bis(ImageLesionCorr,Mask);
                BinarySegCorr=CopyFloatNii(ImageLesionCorr);
            }


            float * TmpBoolSeg=static_cast<float *>(BinarySeg->data);
            //        SaveTmpResult(TmpBoolSeg, g_scratch_dir + "/TestCHB-1.nii.gz", BinarySeg);
            int numel=BinarySeg->nvox;
            bool * BoolSegLong=TranscribeArray<float, bool>(TmpBoolSeg, numel);
            int numelmasked=CalculateNumberMaskedElements(Mask);
            bool * BoolSeg=CreateShort(BoolSegLong, S2LAnalysis, numelmasked);




            //        bool * BoolSeg=SegBoolArray(BinarySeg, Mask);
            bool * BoolSegCorrLong=NULL;
            bool * BoolSegCorr=NULL;
            if (BinarySegCorr!=NULL) {
                float * TmpBoolSegCorr=static_cast<float *>(BinarySegCorr->data);
                BoolSegCorrLong=TranscribeArray<float, bool>(TmpBoolSegCorr, numel);
                BoolSegCorr=CreateShort(BoolSegCorrLong, S2LAnalysis, numelmasked);
            }

            //        bool * BoolSegCorr=SegBoolArray(BinarySegCorr, Mask);
            float * SegToCompareData=static_cast<float *>(SegToCompare->data);
            //        SaveTmpResult(SegToCompareData, g_scratch_dir + "/TestRef.nii.gz", SegToAnalyse);
            nifti_image * RefSeg=HardSegmentation_bis(SegToCompare,Mask);
            bool * BoolRefLong=new bool[numel];
            //        nifti_image * SegToCompareBin=CopyFloatNii(SegToCompare);
            //        Binarisation(SegToCompareBin);
            //        bool * SegToCompareDataBin=static_cast<bool *>(SegToCompareBin->data);
            for (int i=0; i<numel; i++) {
                if (SegToCompareData[i]>0.5) {
                    BoolRefLong[i]=1;
                }
                else{
                    BoolRefLong[i]=0;
                }
            }
            bool * BoolRef=CreateShort(BoolRefLong, S2LAnalysis, numelmasked);
            //        SaveTmpResult_short(BoolRef, g_scratch_dir + "/TestRef.nii.gz", SegToAnalyse, L2SAnalysis);
            //        Lesion segmentation so only  more than threshold is looked at


            //        bool * BoolRef=SegBoolArray(RefSeg, Mask);
            float * FloatSeg=SegFloatArray(SegToAnalyse, Mask);
            float * FloatSegCorr=SegFloatArray(ImageLesionCorr, Mask);
            float * FloatRef=SegFloatArray(SegToCompare, Mask);

            vector<float> MahalDistTP;
            vector<float> MahalDistFP;
            vector<float> MahalDistFN;
            vector<float> MahalDistTPCorr;
            vector<float> MahalDistFPCorr;
            vector<float> MahalDistFNCorr;
            vector<nifti_image *> VectorSegCombine;
            nifti_image * MahalDistCombine=NULL;
            nifti_image * MahalDistImage=MahalDistMaps(WMSeg, SegToAnalyse, DataCompImage);
            nifti_image * MahalDistImageCorr=MahalDistMaps(WMSeg, ImageLesionCorr,DataCompImage);
            nifti_image * MahalDistRef=NULL;
            if (segment_analysis->flag_refLes) {
                MahalDistRef=MahalDistMaps(WMSeg, SegToCompare, DataCompImage);
                VectorSegCombine.push_back(SegToAnalyse);
                VectorSegCombine.push_back(SegToCompare);
                nifti_image * SegCombine=CombineVectorNifti(VectorSegCombine, 0.5);
                nifti_image * DistVentricleCombine=EuclideanDistanceMap(VentricleSeg,SegCombine);
                if (DistVentricleCombine==NULL) {
                    DistVentricleCombine=HardSegmentationThreshold(SegCombine, 0.5);
                }
                //            nifti_set_filenames(DistVentricleCombine, (g_scratch_dir + "/VentricleDist.nii.gz").c_str(), 0, 0);
                //            nifti_image_write(DistVentricleCombine);
                MahalDistCombine=MahalDistMaps(WMSeg, SegCombine, DataCompImage);
                //            Analysis for Lesion intensity and printing
                ofstream TxtFileLesionIntensity(FilenameIntensityConnectAnalysis.c_str());
                cout<<"Non 0 Label GT is "<<CountNonZero(OrderedLabelsGT, numel)<<endl;
                vector<LesionRefIntensity *> VectorLesionIntensity=AnalysisPerRefConnectedIntensity_bis(OrderedLabels, OrderedLabelsGT, MahalDistCombine, segment_analysis,DistVentricleCombine);
                PrintAnalysisPerRefIntensityConnected(VectorLesionIntensity,TxtFileLesionIntensity);
                for (int i=0; i<4; i++) {
                    int Index=FilenameIntensityConnectAnalysis.find_last_of('.');
                    stringstream thresh;
                    thresh << (i*2)*0.5*10;
                    string ThreshStr=thresh.str();
                    string FilenameSave_b=FilenameIntensityConnectAnalysis.substr(0,Index+1);

                    string FilenameThresh=FilenameSave_b+"_"+ThreshStr+".txt";
                    ofstream TxtFileLesionThreshMahal(FilenameThresh.c_str());
                    PrintAnalysisPerRefIntensityConnectedMahalProp(VectorLesionIntensity, i*2+1, TxtFileLesionThreshMahal);
                }




                int numbLesion=VectorLesionIntensity.size();
                for (int l=0; l<numbLesion; l++) {
                    if (VectorLesionIntensity[l]!=NULL) {
                        delete VectorLesionIntensity[l];
                        VectorLesionIntensity[l]=NULL;
                    }
                }

                if (OrderedLabelsCorr!=NULL) {
                    vector<LesionRefIntensity *> VectorLesionIntensityCorr=AnalysisPerRefConnectedIntensity_bis(OrderedLabelsCorr, OrderedLabelsGT, MahalDistCombine, segment_analysis,DistVentricleCombine);
                    ofstream TxtFileLesionIntensityCorr(FilenameIntensityConnectCorrAnalysis.c_str());
                    PrintAnalysisPerRefIntensityConnected(VectorLesionIntensityCorr,TxtFileLesionIntensityCorr);
                    for (int i=0; i<4; i++) {
                        int Index=FilenameIntensityConnectCorrAnalysis.find_last_of('.');
                        stringstream thresh;
                        thresh << (i*2)*0.5*10;
                        string ThreshStr=thresh.str();
                        string FilenameSave_b=FilenameIntensityConnectCorrAnalysis.substr(0,Index+1);

                        string FilenameThresh=FilenameSave_b+"_"+ThreshStr+".txt";
                        ofstream TxtFileLesionThreshMahal(FilenameThresh.c_str());
                        PrintAnalysisPerRefIntensityConnectedMahalProp(VectorLesionIntensityCorr, i*2+1, TxtFileLesionThreshMahal);
                    }


                    int numbLesionCorr=VectorLesionIntensityCorr.size();
                    for (int l=0; l<numbLesionCorr; l++) {
                        if (VectorLesionIntensityCorr[l]!=NULL) {
                            delete VectorLesionIntensityCorr[l];
                            VectorLesionIntensityCorr[l]=NULL;
                        }
                    }
                }

                nifti_image_free(SegCombine);
                SegCombine=NULL;
                if (MahalDistCombine!=NULL) {
                    nifti_image_free(MahalDistCombine);
                    MahalDistCombine=NULL;
                }
                if (DistVentricleCombine!=NULL) {
                    nifti_image_free(DistVentricleCombine);
                    DistVentricleCombine=NULL;
                }

            }
            nifti_image * TPImage=ResultImageFromSeg(SegToAnalyse, RefSeg,Mask,TP);
            if (WMSeg!=NULL && DataCompImage!=NULL) {
                float * MeanTP=GetParamFromImages(TPImage, DataCompImage, 0);
                float * TPSegData=static_cast<float *>(TPImage->data);
                bool * TPSegDataBool=TranscribeArray<float, bool>(TPSegData, numel);
                nifti_image * CopyMahalDist=CopyFloatNii(MahalDistImage);
                float * CopyMahalDistData=static_cast<float *>(CopyMahalDist->data);
                for (int m=0; m<=numbmodal; m++) {
                    MultiplyElementwiseChange<float,bool>(&CopyMahalDistData[m*numel], TPSegDataBool, numel);
                }
                MahalDistFP.push_back(GetMeanData<float ,bool>(CopyMahalDistData,TPSegDataBool,numel));
                MahalDistFP.push_back(GetMax<float,bool>(CopyMahalDistData,TPSegDataBool,numel));
                MahalDistFP.push_back(GetMin<float,bool>(CopyMahalDistData,TPSegDataBool,numel));
                MahalDistFP.push_back(GetProportionAbove<float,bool>(CopyMahalDistData, TPSegDataBool, 4, numel));
                MahalDistFP.push_back(GetProportionUnder<float,bool>(CopyMahalDistData, TPSegDataBool, 2, numel));
                delete [] MeanTP;
                delete [] TPSegDataBool;
                nifti_image_free(CopyMahalDist);
                CopyMahalDist=NULL;
            }
            nifti_image * FPImage=ResultImageFromSeg(SegToAnalyse, RefSeg,Mask, FP);
            if (WMSeg!=NULL && DataCompImage!=NULL) {
                float * MeanFP=GetParamFromImages(FPImage, DataCompImage, 0);
                float * FPSegData=static_cast<float *>(FPImage->data);
                bool * FPSegDataBool=TranscribeArray<float, bool>(FPSegData, numel);
                nifti_image * CopyMahalDist=CopyFloatNii(MahalDistImage);
                float * CopyMahalDistData=static_cast<float *>(CopyMahalDist->data);
                for (int m=0; m<=numbmodal; m++) {
                    MultiplyElementwiseChange<float,bool>(&CopyMahalDistData[m*numel], FPSegDataBool, numel);
                }
                MahalDistFP.push_back(GetMeanData<float ,bool>(CopyMahalDistData,FPSegDataBool,numel));
                MahalDistFP.push_back(GetMax<float,bool>(CopyMahalDistData,FPSegDataBool,numel));
                MahalDistFP.push_back(GetMin<float,bool>(CopyMahalDistData,FPSegDataBool,numel));
                MahalDistFP.push_back(GetProportionAbove<float,bool>(CopyMahalDistData, FPSegDataBool, 4, numel));
                MahalDistFP.push_back(GetProportionUnder<float,bool>(CopyMahalDistData, FPSegDataBool, 2, numel));
                delete [] MeanFP;
                delete [] FPSegDataBool;
                nifti_image_free(CopyMahalDist);
                CopyMahalDist=NULL;
            }
            nifti_image * FNImage=ResultImageFromSeg(SegToAnalyse, RefSeg,Mask, FN);
            if (WMSeg!=NULL && DataCompImage!=NULL) {
                float * MeanFN=GetParamFromImages(FNImage, DataCompImage, 0);
                float * FNSegData=static_cast<float *>(FNImage->data);
                bool * FNSegDataBool=TranscribeArray<float, bool>(FNSegData, numel);
                cout<<"non zero FN "<<CountNonZero(FNSegDataBool, numel)<<endl;
                nifti_image * CopyMahalDist=CopyFloatNii(MahalDistRef);
                float * CopyMahalDistData=static_cast<float *>(CopyMahalDist->data);
                for (int m=0; m<=numbmodal; m++) {
                    MultiplyElementwiseChange<float,bool>(&CopyMahalDistData[m*numel], FNSegDataBool, numel);
                }
                MahalDistFN.push_back(GetMeanData<float ,bool>(CopyMahalDistData,FNSegDataBool,numel));
                MahalDistFN.push_back(GetMax<float,bool>(CopyMahalDistData,FNSegDataBool,numel));
                MahalDistFN.push_back(GetMin<float,bool>(CopyMahalDistData,FNSegDataBool,numel));
                MahalDistFN.push_back(GetProportionAbove<float,bool>(CopyMahalDistData, FNSegDataBool, 4, numel));
                MahalDistFN.push_back(GetProportionUnder<float,bool>(CopyMahalDistData, FNSegDataBool, 2, numel));
                delete [] MeanFN;
                delete [] FNSegDataBool;
                nifti_image_free(CopyMahalDist);
                CopyMahalDist=NULL;
            }

            nifti_image * TPImageCorr=NULL;
            nifti_image * FPImageCorr=NULL;
            nifti_image * FNImageCorr=NULL;
            if (ImageLesionCorr!=NULL) {
                TPImageCorr=ResultImageFromSeg(ImageLesionCorr, RefSeg,Mask,TP);
                if (WMSeg!=NULL && DataCompImage!=NULL) {
                    float * MeanTPCorr=GetParamFromImages(TPImageCorr, DataCompImage, 0);
                    float * TPCorrSegData=static_cast<float *>(TPImageCorr->data);
                    int * IntTPCorrSegData=TranscribeArray<float , int>(TPCorrSegData, numel);
                    MahalDistTPCorr.push_back(GetMahalDist(MeanTPCorr, ParamWM, InvertedCovarianceWM, numbmodal));
                    MahalDistTPCorr.push_back(GetExtremaMahalDistLabel(IntTPCorrSegData, 1, ParamWM, InvertedCovarianceWM, DataCompImage, -1));
                    MahalDistTPCorr.push_back(GetExtremaMahalDistLabel(IntTPCorrSegData, 1, ParamWM, InvertedCovarianceWM, DataCompImage, 1));
                    delete [] MeanTPCorr;
                    delete [] IntTPCorrSegData;
                }
                FPImageCorr=ResultImageFromSeg(ImageLesionCorr, RefSeg,Mask, FP);
                if (WMSeg!=NULL && DataCompImage!=NULL) {
                    float * MeanFPCorr=GetParamFromImages(FPImageCorr, DataCompImage, 0);
                    float * FPCorrSegData=static_cast<float *>(FPImageCorr->data);
                    int * IntFPCorrSegData=TranscribeArray<float , int>(FPCorrSegData, numel);
                    MahalDistFPCorr.push_back(GetMahalDist(MeanFPCorr, ParamWM, InvertedCovarianceWM, numbmodal));
                    MahalDistFPCorr.push_back(GetExtremaMahalDistLabel(IntFPCorrSegData, 1, ParamWM, InvertedCovarianceWM, DataCompImage, -1));
                    MahalDistFPCorr.push_back(GetExtremaMahalDistLabel(IntFPCorrSegData, 1, ParamWM, InvertedCovarianceWM, DataCompImage, 1));
                    delete [] MeanFPCorr;
                    delete [] IntFPCorrSegData;
                }
                FNImageCorr=ResultImageFromSeg(ImageLesionCorr, RefSeg,Mask, FN);
                if (WMSeg!=NULL && DataCompImage!=NULL) {
                    float * MeanFNCorr=GetParamFromImages(FNImageCorr, DataCompImage, 0);
                    float * FNCorrSegData=static_cast<float *>(FNImageCorr->data);
                    int * IntFNCorrSegData=TranscribeArray<float , int>(FNCorrSegData, numel);
                    MahalDistFNCorr.push_back(GetMahalDist(MeanFNCorr, ParamWM, InvertedCovarianceWM, numbmodal));
                    MahalDistFNCorr.push_back(GetExtremaMahalDistLabel(IntFNCorrSegData, 1, ParamWM, InvertedCovarianceWM, DataCompImage, -1));
                    MahalDistFNCorr.push_back(GetExtremaMahalDistLabel(IntFNCorrSegData, 1, ParamWM, InvertedCovarianceWM, DataCompImage, 1));
                    delete [] MeanFNCorr;
                    delete [] IntFNCorrSegData;
                }
                if (DataCompImage!=NULL) {
                    nifti_image_free(DataCompImage);
                    DataCompImage=NULL;
                }
            }


            //        string FilenameSave;
            //        // Setting Filenames for saving images
            //        if(segment_analysis->flag_SegTot){
            //        FilenameSave=nifti_makebasename(segment_analysis->filename_SegTot);
            //        }
            //        else{
            //            FilenameSave=nifti_makebasename(segment_analysis->filename_InLes);
            //        }
            int Index=FilenameSave.find_last_of('/');
            stringstream rG;
            rG << segment_analysis->LesionRuleType;
            stringstream rU;
            rU << segment_analysis->LesionUniformRuleType;
            string G=rG.str();
            string U=rU.str();
            string FilenameSave_b=FilenameSave.substr(0,Index+1);
            string FilenameSave_e=FilenameSave.substr(Index+1,FilenameSave.length());
            string FilenameSave_TP=FilenameSave_b+"TP_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_FN=FilenameSave_b+"FN_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_FP=FilenameSave_b+"FP_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_TPCorr=FilenameSave_b+"TPCorr_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_FNCorr=FilenameSave_b+"FNCorr_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_FPCorr=FilenameSave_b+"FPCorr_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_Mahal=FilenameSave_b+"Mahal_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_MahalCorr=FilenameSave_b+"MahalCorr_"+G+U+FilenameSave_e+".nii.gz";
            string FilenameSave_MahalRef=FilenameSave_b+"MahalRef_"+G+U+FilenameSave_e+".nii.gz";
            //Saving images
            nifti_set_filenames(TPImage, FilenameSave_TP.c_str(), 0, 0);
            nifti_image_write(TPImage);
            nifti_image_free(TPImage);

            nifti_set_filenames(FNImage, FilenameSave_FN.c_str(), 0, 0);
            nifti_image_write(FNImage);
            nifti_image_free(FNImage);

            nifti_set_filenames(FPImage, FilenameSave_FP.c_str(), 0, 0);
            nifti_image_write(FPImage);
            nifti_image_free(FPImage);

            nifti_set_filenames(MahalDistImage, FilenameSave_Mahal.c_str(), 0, 0);
            nifti_image_write(MahalDistImage);
            nifti_image_free(MahalDistImage);
            MahalDistImage=NULL;

            if (MahalDistImageCorr!=NULL) {
                nifti_set_filenames(MahalDistImageCorr, FilenameSave_MahalCorr.c_str(), 0, 0);
                nifti_image_write(MahalDistImageCorr);
                nifti_image_free(MahalDistImageCorr);
                MahalDistImageCorr=NULL;
            }

            if (MahalDistRef!=NULL) {
                nifti_set_filenames(MahalDistRef, FilenameSave_MahalRef.c_str(), 0, 0);
                nifti_image_write(MahalDistRef);
                nifti_image_free(MahalDistRef);
                MahalDistRef=NULL;
            }


            if (TPImageCorr!=NULL) {
                nifti_set_filenames(TPImageCorr, FilenameSave_TPCorr.c_str(), 0, 0);
                nifti_image_write(TPImageCorr);
                nifti_image_free(TPImageCorr);
            }

            if (FNImageCorr!=NULL) {
                nifti_set_filenames(FNImageCorr, FilenameSave_FNCorr.c_str(), 0, 0);
                nifti_image_write(FNImageCorr);
                nifti_image_free(FNImageCorr);
            }

            if (FPImageCorr!=NULL) {
                nifti_set_filenames(FPImageCorr, FilenameSave_FPCorr.c_str(), 0, 0);
                nifti_image_write(FPImageCorr);
                nifti_image_free(FPImageCorr);
            }


            numelmasked=CalculateNumberMaskedElements(Mask);
            if (segment_analysis->flag_numbActive) {
                numbActive=NumbActiveVoxels_bis(Mask);
                TxtFile<<"NumbActives "<<numbActive<<endl;
                if (ImageLesionCorr!=NULL) {
                    TxtFile<<"NumbActives "<<numbActive<<endl;
                }
            }
            if (segment_analysis->flag_TN) {
                int TN=TrueNegatives_bis(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"TN "<<TN<<endl;
                cout<<"TN "<<TN<<" *** ";
                if (ImageLesionCorr!=NULL) {
                    int TNCorr=TrueNegatives_bis(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"TNCorr "<<TNCorr<<endl;
                    cout<<"TNCorr "<<TNCorr<<endl;
                }
            }
            if (segment_analysis->flag_FN){
                int FN=FalseNegatives_bis(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"FN "<<FN<<" ";
                int sizeVector=MahalDistFN.size();
                for(int s=0;s<sizeVector;s++){
                    TxtFile<<MahalDistFN[s]<<" ";
                }
                TxtFile<<endl;
                cout<<"FN "<<FN<< " *** ";
                if (ImageLesionCorr!=NULL) {
                    int FNCorr=FalseNegatives_bis(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"FNCorr "<<FNCorr<<" ";
                    int sizeVector=MahalDistFNCorr.size();
                    for(int s=0;s<sizeVector;s++){
                        TxtFile<<MahalDistFNCorr[s]<<" ";
                    }
                    TxtFile<<endl;
                    cout<<"FNCorr "<<FNCorr<<endl;
                }
            }
            if (segment_analysis->flag_FP) {
                int FP=FalsePositives_bis(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"FP "<<FP<<" ";
                int sizeVector=MahalDistFP.size();
                for(int s=0;s<sizeVector;s++){
                    TxtFile<<MahalDistFP[s]<<" ";
                }
                TxtFile<<endl;
                cout<<"FP "<<FP<< " *** ";
                if (ImageLesionCorr!=NULL) {
                    int FPCorr=FalsePositives_bis(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"FPCorr "<<FPCorr<<" ";
                    int sizeVector=MahalDistFPCorr.size();
                    for(int s=0;s<sizeVector;s++){
                        TxtFile<<MahalDistFPCorr[s]<<" ";
                    }
                    TxtFile<<endl;
                    cout<<"FPCorr "<<FPCorr<<" *** ";
                }
                if (OrderedLabelsGT!=NULL && OrderedLabels!=NULL) {
                    int FPCard=FalsePositives_card(OrderedLabels, OrderedLabelsGT, numel);
                    int maxLesRef=GetMaxLabel(OrderedLabelsGT, numel);
                    float FPRCard=(float)FPCard/(float)maxLesRef;
                    TxtFile<<"FPRCard "<<FPRCard<<endl;
                    cout<<"FPRCard "<<FPRCard<<" *** ";
                    if (OrderedLabelsCorr!=NULL) {
                        int FPCard_corr=FalsePositives_card(OrderedLabelsCorr, OrderedLabelsGT, numel);
                        float FPRCard_corr=(float)FPCard_corr/(float)maxLesRef;
                        TxtFile<<"FPRCardCorr "<<FPRCard_corr<<endl;
                        cout<<"FPRCardCorr "<<FPRCard_corr<<" *** ";
                    }
                }
                float FPRa=FPR(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"FPR "<<FPRa<<endl;
                if (ImageLesionCorr!=NULL) {
                    float FPRCorr=FPR(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"FPRCorr "<<FPRCorr<<endl;
                }
            }
            if (segment_analysis->flag_TP) {
                int TP=TruePositives_bis(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"TP "<<TP<<" ";
                int sizeVector=MahalDistTP.size();
                for(int s=0;s<sizeVector;s++){
                    TxtFile<<MahalDistTP[s]<<" ";
                }
                TxtFile<<endl;
                cout<<"TP "<<TP<< " *** ";
                if (ImageLesionCorr!=NULL) {
                    int TPCorr=TruePositives_bis(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"TPCorr "<<TPCorr<<" ";
                    int sizeVector=MahalDistTPCorr.size();
                    for(int s=0;s<sizeVector;s++){
                        TxtFile<<MahalDistTPCorr[s]<<" ";
                    }
                    TxtFile<<endl;
                    cout<<"TPCorr "<<TPCorr<<" *** ";
                }
                if (OrderedLabelsGT!=NULL && OrderedLabels!=NULL) {
                    int TPCard=TruePositives_card(OrderedLabels, OrderedLabelsGT, numel);
                    int maxLesRef=GetMaxLabel(OrderedLabelsGT, numel);
                    float TPRCard=(float)TPCard/(float)maxLesRef;
                    TxtFile<<"TPRCard "<<TPRCard<<endl;
                    cout<<"TPRCard "<<TPRCard<<" *** ";
                    if (OrderedLabelsCorr!=NULL) {
                        int TPCard_corr=TruePositives_card(OrderedLabelsCorr, OrderedLabelsGT, numel);
                        float TPRCard_corr=(float)TPCard_corr/(float)maxLesRef;
                        TxtFile<<"TPRCardCorr "<<TPRCard_corr<<endl;
                        cout<<"TPRCardCorr "<<TPRCard_corr<<" *** ";
                    }
                }
                float TPRa=TPR(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"TPR "<<TPRa<<endl;
                if (ImageLesionCorr!=NULL) {
                    float TPRCorr=TPR(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"TPRCorr "<<TPRCorr<<endl;
                }


            }

            if (segment_analysis->flag_DE) {
                if (OrderedLabels!=NULL && OrderedLabelsGT!=NULL) {
                    int DE=DetectionError(OrderedLabels, OrderedLabelsGT, numel);
                    TxtFile<<"DE "<<DE<<endl;
                    cout<< "DE "<<DE<<"  ***  ";
                }
                if (OrderedLabelsCorr!=NULL && OrderedLabelsGT!=NULL) {
                    int DECorr=DetectionError(OrderedLabelsCorr, OrderedLabelsGT, numel);
                    TxtFile<<"DECorr "<<DECorr<<endl;
                    cout<<"DECorr "<<DECorr;
                }
                cout<<endl;
            }

            if (segment_analysis->flag_OE) {
                if (OrderedLabels!=NULL && OrderedLabelsGT!=NULL) {
                    int OE=OutlineError(OrderedLabels, OrderedLabelsGT, numel);
                    //                int VolSeg=CountNonZero(BoolSeg, numelmasked);
                    int VolRef=CountNonZero(BoolRef, numelmasked);
                    TxtFile<<"OE "<<OE<<endl;
                    cout<< "OE "<<OE<<"  ***  ";
                    TxtFile<<"OER "<<(float)OE/((VolRef))<<endl;
                    cout<<"OER "<<(float)OE/((VolRef));
                }

                if (OrderedLabelsCorr!=NULL && OrderedLabelsGT!=NULL) {
                    int OECorr=OutlineError(OrderedLabelsCorr, OrderedLabelsGT, numel);
                    //                int VolSegCorr=CountNonZero(OrderedLabelsCorr, numel );
                    int VolRef=CountNonZero(BoolRef, numelmasked);
                    TxtFile<<"OECorr "<<OECorr<<endl;
                    cout<<" OECorr "<<OECorr;
                    TxtFile<<"OERCorr "<<(float)OECorr/((VolRef))<<endl;
                    cout<<"OERCorr "<<(float)OECorr/((VolRef));
                }
                cout<<endl;
            }

            if (segment_analysis->flag_AvDist) {
                float AvDist=AverageDistanceMetric(BoolSegLong, BoolRefLong, SegToAnalyse);
                float MaxAuthorisedDist=0;
                for (int d=0; d<3; d++) {
                    MaxAuthorisedDist+=SegToAnalyse->dim[d+1]*SegToAnalyse->pixdim[d+1]*SegToAnalyse->dim[d+1]*SegToAnalyse->pixdim[d+1];
                }
                MaxAuthorisedDist=sqrt(MaxAuthorisedDist);
                if(AvDist>MaxAuthorisedDist && Mask!=NULL){
                    AvDist=AverageDistanceMetric(MaskData, BoolRefLong, SegToAnalyse);
                }
                TxtFile<<"AvDist "<<AvDist<<endl;
                cout<<"AvDist "<<AvDist<<" *** ";
                if (ImageLesionCorr!=NULL) {
                    float AvDistCorr=AverageDistanceMetric(BoolSegCorrLong, BoolRefLong, SegToAnalyse);
                    if(AvDistCorr>MaxAuthorisedDist && Mask!=NULL){
                        AvDistCorr=AverageDistanceMetric(MaskData, BoolRefLong, SegToAnalyse);
                    }
                    TxtFile<<"AvDistCorr "<<AvDistCorr<<endl;
                    cout<<"AvDistCorr "<<AvDistCorr<<" *** ";
                }
                cout<<endl;
            }

            if (segment_analysis->flag_Sens){
                float Sens=Sensitivity(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"Sens "<<Sens<<endl;
                if (ImageLesionCorr!=NULL) {
                    float SensCorr=Sensitivity(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"SensCorr "<<SensCorr<<endl;
                }
            }
            if (segment_analysis->flag_Spec) {
                float Spec=Specificity(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"Spec "<<Spec<<endl;
                if (ImageLesionCorr!=NULL) {
                    float SpecCorr=Specificity(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"SpecCorr "<<SpecCorr<<endl;
                }
            }
            if (segment_analysis->flag_Acc) {
                float Acc=Accuracy(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"Acc "<<Acc<<endl;
                if (ImageLesionCorr!=NULL) {
                    float AccCorr=Accuracy(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"AccCorr "<<AccCorr<<endl;
                }
            }
            if (segment_analysis->flag_FPR) {
                float FPRa=FPR(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"FPR "<<FPRa<<endl;
                if (ImageLesionCorr!=NULL) {
                    float FPRCorr=FPR(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"FPRCorr "<<FPRCorr<<endl;
                }
            }
            if (segment_analysis->flag_TPR) {
                float TPRa=TPR(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"TPR "<<TPRa<<endl;
                if (ImageLesionCorr!=NULL) {
                    float TPRCorr=TPR(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"TPRCorr "<<TPRCorr<<endl;
                }
            }
            if (segment_analysis->flag_VD) {
                float VDa=VD(BoolSeg, BoolRef, numelmasked);
                TxtFile<<"VD "<<VDa<<endl;
                if (ImageLesionCorr!=NULL) {
                    float VDCorr=VD(BoolSegCorr, BoolRef, numelmasked);
                    TxtFile<<"VDCorr "<<VDCorr<<endl;
                }
            }
            if (segment_analysis->flag_PSI) {
                float PSI=ProbabilitySimilarityIndex(FloatSeg, BoolRef, numelmasked);
                TxtFile<<"PSI "<<PSI<<endl;
                if (ImageLesionCorr!=NULL) {
                    float PSICorr=ProbabilitySimilarityIndex(FloatSegCorr, BoolRef, numelmasked);
                    TxtFile<<"PSICorr "<<PSICorr<<endl;
                }
            }
            if (segment_analysis->flag_DCh) {
                float TP=TruePositives_bis(BoolSeg, BoolRef, numelmasked);
                float FN=FalseNegatives_bis(BoolSeg, BoolRef, numelmasked);
                float FP=FalsePositives_bis(BoolSeg, BoolRef, numelmasked);
                float DCS=2.0*TP/(FP+FN+2.0*TP);
                cout<<"DCS "<<DCS<<" *** ";
                TxtFile<<"DCS "<<DCS<<endl;

                if (ImageLesionCorr!=NULL) {
                    float TPCorr=TruePositives_bis(BoolSegCorr, BoolRef, numelmasked);
                    float FNCorr=FalseNegatives_bis(BoolSegCorr, BoolRef, numelmasked);
                    float FPCorr=FalsePositives_bis(BoolSegCorr, BoolRef, numelmasked);
                    float DCSCorr=2.0*TPCorr/(FPCorr+FNCorr+2.0*TPCorr);
                    TxtFile<<"DCSCorr "<<DCSCorr<<endl;
                    cout<<"DCSCorr "<<DCSCorr<<endl;;
                }
            }
            if (segment_analysis->flag_Corr){
                float Corr=Correlation(SegToAnalyse, SegToCompare, Mask);
                TxtFile<<"Corr "<<Corr<<endl;
                if (ImageLesionCorr!=NULL) {
                    float CorrCorr=Correlation(ImageLesionCorr, SegToCompare, Mask);
                    TxtFile<<"CorrCorr "<<CorrCorr<<endl;
                }
            }
            if (segment_analysis->flag_TLLSoft) {
                float TLLSoft=TLL(FloatSeg, numelmasked);
                TxtFile<<"TLLSoft "<<TLLSoft<<endl;
                if (ImageLesionCorr!=NULL) {
                    float TLLSoftCorr=TLL(FloatSegCorr, numelmasked);
                    TxtFile<<"TLLSoftCorr "<<TLLSoftCorr<<endl;
                }
            }
            if (segment_analysis->flag_TLLHard) {
                float TLLHard=TLL(BoolSeg, numelmasked);
                TxtFile<<"TLLHard "<<TLLHard<<endl;
                if (ImageLesionCorr!=NULL) {
                    float TLLHardCorr=TLL(BoolSegCorr, numelmasked);
                    TxtFile<<"TLLHardCorr "<<TLLHardCorr<<endl;
                }

            }

            delete [] BoolRef;
            if (BoolRefLong!=NULL) {
                delete [] BoolRefLong;
                BoolRefLong=NULL;
            }
            delete [] BoolSeg;
            if (BoolSegLong!=NULL) {
                delete [] BoolSegLong;
                BoolSegLong=NULL;
            }
            if (BoolSegCorrLong !=NULL) {
                delete [] BoolSegCorrLong;
                BoolSegCorrLong=NULL;
            }
            delete [] FloatRef;
            delete [] FloatSeg;

            if(S2LAnalysis!=NULL){
                delete [] S2LAnalysis;
                S2LAnalysis=NULL;
            }
            if (ConnectedGTImage!=NULL) {
                nifti_image_free(ConnectedGTImage);
                ConnectedGTImage=NULL;
            }
            if (WMSeg !=NULL) {
                nifti_image_free(WMSeg);
                WMSeg=NULL;
            }
            if (OrderedLabels!=NULL) {
                delete [] OrderedLabels;
                OrderedLabels=NULL;
            }
            if (OrderedLabelsCorr!=NULL) {
                delete [] OrderedLabelsCorr;
                OrderedLabelsCorr=NULL;
            }
            if (OrderedLabelsGT!=NULL) {
                delete [] OrderedLabelsGT;
                OrderedLabelsGT=NULL;
            }
            if (BoolSegCorr!=NULL) {
                delete [] BoolSegCorr;
                BoolSegCorr=NULL;
            }
            if (FloatSegCorr!=NULL) {
                delete [] FloatSegCorr;
                FloatSegCorr=NULL;
            }
            if (Mask!=NULL) {
                nifti_image_free(Mask);
                Mask=NULL;
            }
            if (RefSeg!=NULL) {
                nifti_image_free(RefSeg);
            }
            if (BinarySeg!=NULL) {
                nifti_image_free(BinarySeg);
                BinarySeg=NULL;
            }

            if(BinarySegCorr!=NULL){
                nifti_image_free(BinarySegCorr);
                BinarySegCorr=NULL;
            }
            if (SegToAnalyse!=NULL) {
                nifti_image_free(SegToAnalyse);
                SegToAnalyse=NULL;
            }
            if(SegToAnalyseNIV!=NULL){
                nifti_image_free(SegToAnalyseNIV);
                SegToAnalyseNIV=NULL;
            }
            if (ImageLesionCorr!=NULL) {
                nifti_image_free(ImageLesionCorr);
                ImageLesionCorr=NULL;
            }
            if (SegToCompare!=NULL) {
                nifti_image_free(SegToCompare);
                SegToCompare=NULL;
            }
            BoolRef=NULL;
            BoolSeg=NULL;
            FloatRef=NULL;
            FloatSeg=NULL;
            cout<< "Analysis done for Lesion ref"<<endl;
            return EXIT_SUCCESS;
        }
        else if(segment_analysis->flag_Analysis){
            if(segment_analysis->flag_compMult){
                if (segment_analysis->numbMult != SegToAnalyse->nu*SegToAnalyse->nt){
                    if(SegToAnalyse->nu*SegToAnalyse->nt<segment_analysis->numbMult){
                        fprintf(stderr,"Err:\t The number of classes to analyse cannot be lower than the number of classes in the reference.\n");
                        Usage(argv[0]);
                        return EXIT_SUCCESS;
                    }
                    else{
                        cout<<"Outlier model used"<<endl;
                        vector<nifti_image*> PartialSeg;
                        for (int c=0;c<segment_analysis->numbMult;c++){
                            PartialSeg.push_back(nifti_image_read(segment_analysis->filename_compMult[c],true));
                        }
                        SegToCompare=PrepareCompareSoftSegmentationOutlierMult(PartialSeg,Mask);
                        MaskImage(SegToAnalyse,Mask);
                        numbActive=NumbActiveVoxels(SegToAnalyse);
                    }

                }
                else{
                    vector<nifti_image*> PartialSeg;
                    for (int c=0;c<segment_analysis->numbMult;c++){
                        PartialSeg.push_back(nifti_image_read(segment_analysis->filename_compMult[c],true));
                    }
                    SegToCompare=PrepareCompareSoftSegmentations(PartialSeg,Mask);
                    //          for(int c=0;c<segment_analysis->numbMult;c++){
                    //              nifti_image_free(PartialSeg[c]);
                    //          }
                    //          PartialSeg.clear();
                }
            }
            if(segment_analysis->flag_compJoint){
                SegToCompare=nifti_image_read(segment_analysis->filename_compJoint,true);
                int numbMultComp=SegToCompare->nu*SegToCompare->nt;
                int numbMultAna=SegToAnalyse->nu*SegToAnalyse->nt;
                if(numbMultAna==numbMultComp){
                    MaskImage(SegToCompare,Mask);
                }
                else if(numbMultAna<numbMultComp){
                    fprintf(stderr,"Err:\t The number of classes to analyse cannot be lower than the number of classes in the reference.\n");
                    Usage(argv[0]);
                    return EXIT_SUCCESS;
                }
                else{
                    cout << "Outlier model to handle"<<endl;
                    PrepareCompareSoftSegmentationOutlierJoint(SegToCompare,Mask);
                    MaskImage(SegToAnalyse,Mask);
                    numbActive=NumbActiveVoxels(SegToAnalyse);
                }

                //      nifti_image_free(TempSegToCompare);
            }
            int numbclasses = SegToAnalyse->nu*SegToAnalyse->nt;
            if(segment_analysis->flag_HS){
                nifti_image * SegHSResult=HardSegmentation(SegToAnalyse);
                nifti_set_filenames(SegHSResult, segment_analysis->filename_Out, 0, 0);
                nifti_image_write(SegHSResult);
                nifti_image_free(SegHSResult);
            }
            //  nifti_set_filenames(SegToCompare, (g_scratch_dir + "/SegToCompare.nii.gz").c_str(),0,0);
            //  nifti_image_write(SegToCompare);
            if(segment_analysis->flag_Distance==1){
                nifti_image * SoftSeg2=nifti_image_read(segment_analysis->filename_Dist,1);
                vector<nifti_image*> VectDistanceImages=DistanceFromClassesHS(SegToAnalyse,SoftSeg2,segment_analysis->ExploRadius);
                int lenDistanceImages=VectDistanceImages.size();
                string PrimeName=nifti_makebasename(segment_analysis->filename_Dist);
                for(int i=0;i<lenDistanceImages;i++){
                    ostringstream Number;
                    Number << i;
                    string FinalName=PrimeName+"_Dist"+Number.str()+".nii.gz";
                    cout<<FinalName<<endl;
                    nifti_set_filenames(VectDistanceImages[i],FinalName.c_str(),0,0);
                    nifti_image_write(VectDistanceImages[i]);
                    nifti_image_free(VectDistanceImages[i]);
                }
                VectDistanceImages.clear();
            }
            if(segment_analysis->flag_NW==1){
                nifti_image * NWResult=NeighbourhingWeight(SegToAnalyse,segment_analysis->ExploRadius);
                string PrimeName=nifti_makebasename(segment_analysis->filename_In[0]);
                string FinalName=PrimeName+"_NW.nii.gz";
                nifti_set_filenames(NWResult,FinalName.c_str(),0,0);
                nifti_image_write(NWResult);
                nifti_image_free(NWResult);
            }
            if(segment_analysis->flag_outTxt && segment_analysis->flag_Les==0){
                ofstream TxtFile(segment_analysis->filename_OutTxt);
                if(segment_analysis->flag_numbActive){
                    TxtFile<< "ActiveNumber "<<numbActive<<endl;
                }
                if(segment_analysis->flag_DCh && SegToCompare!=NULL){
                    float * DCh=HardDiceScores(SegToAnalyse,SegToCompare);
                    TxtFile << "HDC ";
                    for(int c=0;c<numbclasses+1;c++){
                        TxtFile<<DCh[c]<<" ";
                    }
                    TxtFile<< endl;
                    delete [] DCh;
                }
                if(segment_analysis->flag_DCs && SegToCompare!=NULL){
                    float * DCs=SoftDiceScores(SegToAnalyse,SegToCompare);
                    TxtFile << "SDC ";
                    for(int c=0;c<numbclasses+1;c++){
                        TxtFile<<DCs[c]<<" ";
                    }
                    TxtFile<< endl;
                    delete[]DCs;
                    DCs=NULL;

                }
                if(segment_analysis->flag_TP && SegToCompare!=NULL){
                    float * TP=TruePositives(SegToAnalyse,SegToCompare);
                    TxtFile << "TP ";
                    for(int c=0;c<numbclasses+1;c++){
                        TxtFile<<TP[c]<<" ";

                    }
                    TxtFile<< endl;
                    delete [] TP;
                    TP=NULL;
                }
                if(segment_analysis->flag_FP && SegToCompare!=NULL){
                    float * FP=FalsePositives(SegToAnalyse,SegToCompare);
                    TxtFile << "FP ";
                    for(int c=0;c<numbclasses+1;c++){
                        TxtFile<<FP[c]<<" ";
                    }
                    TxtFile<< endl;
                    delete[] FP;
                    FP=NULL;
                }
                if(segment_analysis->flag_TN && SegToCompare!=NULL){
                    float * TN=TrueNegatives(SegToAnalyse,SegToCompare);
                    TxtFile << "TN ";
                    for(int c=0;c<numbclasses+1;c++){
                        TxtFile<<TN[c]<<" ";

                    }
                    TxtFile<< endl;
                    delete [] TN;
                    TN=NULL;
                }
                if(segment_analysis->flag_FN && SegToCompare!=NULL){
                    float * FN=FalseNegatives(SegToAnalyse,SegToCompare);
                    TxtFile << "FN ";
                    for(int c=0;c<numbclasses+1;c++){
                        TxtFile<<FN[c]<<" ";

                    }
                    TxtFile<< endl;
                    delete [] FN;
                    FN = NULL;
                }

            }


            else if(segment_analysis->flag_Les==0){
                if(segment_analysis->flag_numbActive){
                    cout<< "ActiveNumber "<<numbActive<<endl;
                }
                if(segment_analysis->flag_DCh && SegToCompare!=NULL){
                    float * DCh=HardDiceScores(SegToAnalyse,SegToCompare);
                    cout << "HDC ";
                    for(int c=0;c<numbclasses+1;c++){
                        cout<<DCh[c]<<" ";
                    }
                    cout<< endl;
                    delete [] DCh;
                    DCh=NULL;

                }
                if(segment_analysis->flag_DCs && SegToCompare!=NULL){
                    float * DCs=SoftDiceScores(SegToAnalyse,SegToCompare);
                    cout << "SDC ";
                    for(int c=0;c<numbclasses+1;c++){
                        cout<<DCs[c]<<" ";
                    }
                    cout<< endl;
                    delete[]DCs;
                    DCs=NULL;
                }
                if(segment_analysis->flag_TP && SegToCompare!=NULL){
                    float * TP=TruePositives(SegToAnalyse,SegToCompare);

                    cout << "TP ";
                    for(int c=0;c<numbclasses+1;c++){
                        cout<<TP[c]<<" ";

                    }
                    cout<< endl;
                    delete [] TP;
                    TP=NULL;
                }
                if(segment_analysis->flag_FP && SegToCompare!=NULL){
                    float * FP=FalsePositives(SegToAnalyse,SegToCompare);
                    cout << "FP ";
                    for(int c=0;c<numbclasses+1;c++){
                        cout<<FP[c]<<" ";

                    }
                    cout<< endl;
                    delete [] FP;
                    FP=NULL;
                }
                if(segment_analysis->flag_TN && SegToCompare!=NULL){
                    float * TN=TrueNegatives(SegToAnalyse,SegToCompare);
                    cout << "TN ";
                    for(int c=0;c<numbclasses+1;c++){
                        cout<<TN[c]<<" ";

                    }
                    cout<< endl;
                    delete [] TN;
                    TN = NULL;
                }
                if(segment_analysis->flag_FN && SegToCompare!=NULL){
                    float * FN=FalseNegatives(SegToAnalyse,SegToCompare);
                    cout << "FN ";
                    for(int c=0;c<numbclasses+1;c++){
                        cout<<FN[c]<<" ";

                    }
                    cout<< endl;
                    delete [] FN;
                    FN=NULL;
                }

            }
        }
        nifti_image_free(SegToAnalyse);
        SegToAnalyse=NULL;
        nifti_image_free(SegToCompare);
    }

    if (TreeToAnalyse!=NULL) {
        delete TreeToAnalyse;
        TreeToAnalyse=NULL;
    }
    if(ParamWM!=NULL){
        delete [] ParamWM;
        ParamWM=NULL;
    }
    if(InvertedCovarianceWM!=NULL){
        delete [] InvertedCovarianceWM;
        InvertedCovarianceWM=NULL;
    }
    if(ConnectedGTImage!=NULL){
        nifti_image_free(ConnectedGTImage);
        ConnectedGTImage=NULL;
    }
    if (VentricleSeg!=NULL) {
        nifti_image_free(VentricleSeg);
        VentricleSeg=NULL;
    }
    if (ICSF!=NULL) {
        delete [] ICSF;
        ICSF=NULL;
    }
    if (DGM!=NULL) {
        delete [] DGM;
        DGM=NULL;
    }
    if(DGMBool!=NULL){
        delete [] DGMBool;
        DGMBool=NULL;
    }
    if (SummarisedSeg1!=NULL) {
        nifti_image_free(SummarisedSeg1);
        SummarisedSeg1=NULL;
    }
    if (PriorsDGM!=NULL) {
        nifti_image_free(PriorsDGM);
        PriorsDGM=NULL;
    }
    if (PriorsICSF!=NULL) {
        nifti_image_free(PriorsICSF);
        PriorsICSF=NULL;
    }
    if(VentricleBool!=NULL){
        delete [] VentricleBool;
        VentricleBool=NULL;
    }
    if (SPRegion!=NULL) {
        delete [] SPRegion;
        SPRegion=NULL;
    }
    if(ImageLesionCorr!=NULL){
        nifti_image_free(ImageLesionCorr);
        ImageLesionCorr=NULL;
    }
    if (Mask!=NULL) {
        nifti_image_free(Mask);
        Mask=NULL;
    }
    if(SegToAnalyse!=NULL){
        nifti_image_free(SegToAnalyse);
        SegToAnalyse=NULL;
    }
    if(BinarySeg!=NULL){
        nifti_image_free(BinarySeg);
        BinarySeg=NULL;
    }
    if(SegToAnalyseNIV!=NULL){
        nifti_image_free(SegToAnalyseNIV);
        SegToAnalyseNIV=NULL;
    }
    if(SegToCompare!=NULL){
        nifti_image_free(SegToCompare);
        SegToCompare=NULL;
    }
    if (S2LAnalysis!=NULL) {
        delete [] S2LAnalysis;
        S2LAnalysis=NULL;
    }
    if(OrderedLabels!=NULL){
        delete [] OrderedLabels;
        OrderedLabels=NULL;
    }
    if (OrderedLabelsCorr!=NULL) {
        delete [] OrderedLabelsCorr;
        OrderedLabelsCorr=NULL;
    }
    if (OrderedLabelsGT!=NULL) {
        delete [] OrderedLabelsGT;
        OrderedLabelsGT=NULL;
    }



    std::cout << "Analysis done!\n";
    delete [] segment_analysis;

    return EXIT_SUCCESS;
}

