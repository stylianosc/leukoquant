#include "_seg_GIF.h"




void Usage(char *exec)
{
#ifdef _OPENMP
    int omp_true_max_threads=omp_get_num_procs();
#endif
    printf("\n  GIF ");
#ifdef _OPENMP
    printf("(OpenMP x%d)",omp_get_num_procs());
#endif
    printf(":\n  Usage -> %s <mandatory> <options>\n\n",exec);
    printf("\n  * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    printf("    -in <filename>\t| Input target image filename\n");
    printf("    -db <XML>   \t| Path to database <XML> file\n");
    printf("\n  * * * * * * * * * * * General Options * * * * * * * * * * * * * * * * * * * * * *\n\n");
    printf("    -mask <filename>\t| Mask over the input image\n");
    printf("    -out <path> \t| Output folder [./]\n");
    printf("    -cpp <cpp_path>\t| Read/store the cpps in <cpp_path>, rather than in memory \n");
    printf("    -geo \t\t| Save Geo to output folder\n");
    printf("    -upd \t\t| Update label given the previous db result.\n");
    printf("    -dbt <fname> <tar>\t| Sets a database (<fname>) specific target (<tar>) [in]\n");
    printf("    -v <int>    \t| Verbose level (0 = off, 1 = on, 2 = debug) [0]\n");
#if defined (_OPENMP)
    int omp_default_num_threads=omp_get_num_procs();
    if(getenv("OMP_NUM_THREADS")!=NULL)
       omp_default_num_threads=atoi(getenv("OMP_NUM_THREADS"));
    printf("    -omp <int>    \t| Use a maximum of <int> openMP threads [default=%d, max=%d]\n",omp_default_num_threads,omp_true_max_threads);
    printf("    -ompj <int>    \t| Run <int> regs in parallel [default=%d, max=%d]\n",omp_default_num_threads,omp_true_max_threads);
#endif
#ifdef _GIT_HASH
   printf("    -version    \t| Print current source code git hash key and exit");
   printf("                  |%s" ,_GIT_HASH);
#endif
    printf("\n  * * * * * * * * * * * Fusion Options * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    printf("    -lssd_ker <float>\t| SSD kernel stdev in voxels (mm if negative) [-2.5]\n");
    printf("    -ldef_ker <float>\t| DEF kernel stdev in voxels (mm if negative) [-2.5]\n");
    printf("    -lncc_ker <float>\t| NCC kernel stdev in voxels (mm if negative) [-2.5]\n");
    printf("    -t1dti_ker <float>\t| T1DTI kernel stdev in voxels (mm if negative) [-2.5]\n");
    printf("    -lssd_weig <float>\t| SSD distance weight <float> [0.0]\n");
    printf("    -ldef_weig <float>\t| DEF distance weight <float> [0.0]\n");
    printf("    -lncc_weig <float>\t| NCC distance weight <float> [1.0]\n");
    printf("    -t1dti_weig <float>\t| T1DTI distance weight <float> [0.0]\n");
    printf("    -temper <float>\t| Kernel emperature <float> [0.15]\n");
    printf("    -sort_beta <float>\t| The beta scaling factor [0.5]\n");
    printf("    -sort_numb <char>\t| The number of elements in the sort [7]\n");
    printf("\n  * * * * * * * * * * * Reg Options * * * * * * * * * * * * * * * * * * * * * *\n\n");
    printf("    -regAff <aff.txt>\t| Input affine file from database to target\n");
    printf("    -regNMI \t\t| Ust NMI as a registration similarity, instead of LNCC\n");
    printf("    -regBE <float>\t| Bending energy value for the registration [0.005]\n");
    printf("    -regJL <float>\t| Jacobian log value for the registration [0.0001]\n");
    printf("    -regSL\t\t| Skip the second Level non-rigid registration\n");
    printf("    -regSpeed\t\t| Faster (but less accurate) affine registration\n");
    printf("\n  * * * * * * * * * * * Seg Options * * * * * * * * * * * * * * * * * * * * * *\n\n");
    printf("    -segRF <fl1> <fl1>\t| Relax Priors (fl1 = relax factor, fl2 = gauss std) [0,0]\n");
    printf("    -segMRF <float> \t| The segmentation MRF beta value [0.15]\n");
    printf("    -segPT <float> \t| The segmentation prior threhsold [0.2]\n");
    printf("    -segIter <int> \t| The minimum number of iterationss [3]\n");
    printf("\n  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");


    return;
}

int main(int argc, char **argv)
{

#ifndef NDEBUG
    printf("[NiftySeg DEBUG] NiftySeg has been compiled in DEBUG mode\n");
#endif
    if (argc <= 1)
    {
        Usage(argv[0]);
        return 0;
    }

#ifdef _OPENMP
    int omp_curr_num_threads=omp_get_num_procs();
    if(getenv("OMP_NUM_THREADS")!=NULL)
       omp_curr_num_threads=atoi(getenv("OMP_NUM_THREADS"));
    omp_set_num_threads(omp_curr_num_threads);
#endif
    string cpp_folder;
    string out_folder="./";
    string target_file;
    string target_mask;
    string db_xml;
    string input_affine;
    std::map <string, string> DatabaseTargetMap;

    float ssd_ker=-2.5f;
    float def_ker=-2.5f;
    float ncc_ker=-2.5f;
    float t1dti_ker=-2.5f;
    float ssd_weig=0.0f;
    float def_weig=0.0f;
    float ncc_weig=1.0f;
    float t1dti_weig=0.0f;
    float temperature=0.15;

    char nrrNMI=0;
    bool nrrRunAllLevels=true;
    float nrrBE=0.005;
    float nrrJL=0.0001;
    int aladinSpeed=0;
    int platf=0;

    float segRFsmo=0;
    float segRFalpha=0;
    float segMRFbeta=0.15;
    float probthresh=0.2;
    int segMinNumbIter=3;
    int ompj=1;

    bool saveGeo=false;
    bool updateLab=false;

    unsigned char sort_number=7;
    float sort_beta=0.5f;

    bool mandatory_target_defined=false;
    bool mandatory_database_defined=false;
    bool nipypeDoOnlyAffine=false;
    int verbose=0;


    for(int i=1; i<argc; i++)
    {
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }
        else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc)
        {
            target_file=argv[++i];
            mandatory_target_defined=true;
        }
        else if(strcmp(argv[i], "-db") == 0 && (i+1)<argc)
        {
            db_xml = argv[++i];
            mandatory_database_defined=true;
        }
        else if(strcmp(argv[i], "-mask") == 0 && (i+1)<argc)
        {
            target_mask=argv[++i];
        }
        else if(strcmp(argv[i], "-cpp") == 0 && (i+1)<argc)
        {
            cpp_folder = argv[++i];
        }
        else if(strcmp(argv[i], "-out") == 0 && (i+1)<argc)
        {
            out_folder = argv[++i];
        }
        else if(strcmp(argv[i], "-dbt") == 0 && (i+2)<argc)
        {
            string stringfirst=string(argv[++i]);
            string stringsecond=string(argv[++i]);
            DatabaseTargetMap.insert(std::pair<string, string>(stringfirst,stringsecond));
        }
        else if(strcmp(argv[i], "-geo") == 0)
        {
            saveGeo = true;
        }
        else if(strcmp(argv[i], "-upd") == 0)
        {
            updateLab = true;
        }
        else if(strcmp(argv[i], "-v") == 0 && (i+1)<argc)
        {
            verbose = atoi(argv[++i]);
        }
        // Parameters
        else if(strcmp(argv[i], "-lssd_ker") == 0 && (i+1)<argc)
        {
            ssd_ker = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-ldef_ker") == 0 && (i+1)<argc)
        {
            def_ker = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-lncc_ker") == 0 && (i+1)<argc)
        {
            ncc_ker = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-t1dti_ker") == 0 && (i+1)<argc)
        {
            t1dti_ker = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-lssd_weig") == 0 && (i+1)<argc)
        {
            ssd_weig = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-ldef_weig") == 0 && (i+1)<argc)
        {
            def_weig = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-lncc_weig") == 0 && (i+1)<argc)
        {
            ncc_weig = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-t1dti_weig") == 0 && (i+1)<argc)
        {
            t1dti_weig = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-temper") == 0 && (i+1)<argc)
        {
            temperature = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-regAff") == 0 && (i+1)<argc)
        {
            input_affine=argv[++i];
        }
        else if(strcmp(argv[i], "-regNMI") == 0 && (i)<argc)
        {
            nrrNMI=1;
        }
        else if(strcmp(argv[i], "-regBE") == 0 && (i+1)<argc)
        {
            nrrBE=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-regJL") == 0 && (i+1)<argc)
        {
            nrrJL=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-regSL") == 0 && (i)<argc)
        {
            nrrRunAllLevels=false;
        }
        else if(strcmp(argv[i], "-segRF") == 0 && (i+2)<argc)
        {
            segRFalpha=atof(argv[++i]);
            segRFsmo=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-segMRF") == 0 && (i+1)<argc)
        {
            segMRFbeta=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-segPT") == 0 && (i+1)<argc)
        {
            probthresh=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-segIter") == 0 && (i+1)<argc)
        {
            segMinNumbIter=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-sort_beta") == 0 && (i+1)<argc)
        {
            sort_beta = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-sort_numb") == 0 && (i+1)<argc)
        {
            sort_number = (unsigned char)(atoi(argv[++i]));
        }
        else if( (strcmp(argv[i], "-affOnly") == 0 || strcmp(argv[i], "--aff_only") == 0 ) && (i)<argc) // the --aff_only option is for Nico (nipype)
        {
          nipypeDoOnlyAffine=true;
        }
        else if(strcmp(argv[i], "-regSpeed") == 0 && (i)<argc)
        {
            aladinSpeed=1;
        }
        else if(strcmp(argv[i], "-platf") == 0 && (i+1)<argc)
        {
            // 0=CPU (default), 1=CUDA -- requires this seg_GIF build to have
            // been linked against a CUDA-enabled NiftyReg (bin_gpu build).
            platf=atoi(argv[++i]);
        }
#if defined (_OPENMP)
        else if(strcmp(argv[i], "-omp")==0 || strcmp(argv[i], "--omp")==0)
        {
            omp_set_num_threads(atoi(argv[++i]));
        }
        else if(strcmp(argv[i], "-ompj")==0 || strcmp(argv[i], "--omp")==0)
        {
            ompj=atoi(argv[++i]);
        }
#endif
#ifdef _GIT_HASH
      else if( strcmp(argv[i], "-version")==0 ||
          strcmp(argv[i], "-Version")==0 ||
          strcmp(argv[i], "-V")==0 ||
          strcmp(argv[i], "-v")==0 ||
          strcmp(argv[i], "--v")==0 ||
          strcmp(argv[i], "--version")==0)
      {
         printf("%s\n",_GIT_HASH);
         return EXIT_SUCCESS;
      }
#endif
        else
        {
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }

    if(mandatory_database_defined==false)
    {
        fprintf(stderr,"Err:\tMandatory -db parameter is missing\n\n");
        Usage(argv[0]);
        return 1;
    }
    if(mandatory_target_defined==false )
    {
        fprintf(stderr,"Err:\tMandatory -in parameter is missing\n\n");
        Usage(argv[0]);
        return 1;
    }

    //    if(ssd_weig==0 && def_weig==0 && ncc_weig==0 && t1dti_weig==0) // if no weight is set up... then assume lncc=1
    //    {
    //        ncc_weig=1;
    //    }

    GIF GIF_object(db_xml.c_str(), target_file.c_str());
    GIF_object.setResultsPath(out_folder.c_str());
    GIF_object.setVerbose(verbose);

    if(cpp_folder.length()>0)
        GIF_object.setTransformationsFolder(cpp_folder.c_str());

    if(target_mask.length()>0)
        GIF_object.setRootTargetMask(target_mask.c_str());

    if(input_affine.length()>0)
        GIF_object.setInputAffineFromDatabase(input_affine);


    GIF_object.set_kernelstd_def(def_ker);
    GIF_object.set_kernelstd_ssd(ssd_ker);
    GIF_object.set_kernelstd_lncc(ncc_ker);
    GIF_object.set_kernelstd_t1dti(t1dti_ker);

    GIF_object.set_weight_def(def_weig);
    GIF_object.set_weight_ssd(ssd_weig);
    GIF_object.set_weight_lncc(ncc_weig);
    GIF_object.set_weight_t1dti(t1dti_weig);

    GIF_object.set_sort_number(sort_number);
    GIF_object.set_sort_beta(sort_beta);

    GIF_object.set_temperature(temperature);

    GIF_object.set_nrr_metric(nrrNMI);
    GIF_object.set_regJacobianLogWeight(nrrJL);
    GIF_object.set_regBendingEnergyWeight(nrrBE);
    GIF_object.set_regRunSecondLevel(nrrRunAllLevels);

    GIF_object.set_segRFsmo(segRFsmo);
    GIF_object.set_segRFalpha(segRFalpha);
    GIF_object.set_segMRFbeta(segMRFbeta);
    GIF_object.set_probthresh(probthresh);
    GIF_object.set_segMinNumIter(segMinNumbIter);
    GIF_object.set_regSpeed(aladinSpeed);
    GIF_object.set_regPlatf(platf);

    GIF_object.set_saveGeo(saveGeo);
    GIF_object.set_updateLab(updateLab);
    GIF_object.set_ompj(ompj);




    // For each defined database, set its target
    std::map<string,string>::iterator DatabaseIterator = DatabaseTargetMap.begin();
    while(DatabaseIterator != DatabaseTargetMap.end())
    {
        GIF_object.setDatabaseTarget(DatabaseIterator->first,DatabaseIterator->second);
        DatabaseIterator++;
    }


    GIF_object.Run();



    return 0;
}



