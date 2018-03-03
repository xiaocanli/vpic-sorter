#include "stdlib.h"
#include "hdf5.h"
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include "configuration.h"

void print_help();

/******************************************************************************
 * Get the analysis configuration from command line arguments
 ******************************************************************************/
int get_configuration(int argc, char **argv, int mpi_rank, config_t *config)
{
    int c;
    static const char *options="a:b:c:ef:g:hi:k:l:m:o:pqrst:u:vw";
    static struct option long_options[] =
    {
        {"tmax", required_argument, 0, 'b'},
        {"tinterval", required_argument, 0, 'i'},
        {"tstep", required_argument, 0, 11},
        {"filepath", required_argument, 0, 1},
        {"species", required_argument, 0, 2},
        {"filename_traj", required_argument, 0, 3},
        {"nptl_traj", required_argument, 0, 4},
        {"ratio_emax", required_argument, 0, 5},
        {"tmin", required_argument, 0, 6},
        {"is_recreate", required_argument, 0, 7},
        {"nsteps", required_argument, 0, 8},
        {"reduced_tracer", required_argument, 0, 9},
        {"load_tracer_meta", required_argument, 0, 'r'},
        {0, 0, 0, 0},
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    extern char *optarg;

    /* Default values */
    config->key_index = 1;
    config->sort_key_only = 0;
    config->skew_data = 0;
    config->verbose = 0;
    config->write_result = 1;
    config->collect_data = 1;
    config->weak_scale_test = 0;
    config->weak_scale_test_length = 1000000;
    config->local_sort_threaded = 0;
    config->local_sort_threads_num = 16;
    config->meta_data = 0;
    config->multi_tsteps = 0;
    config->ux_kindex = 0;
    config->nptl_traj = 10;
    config->ratio_emax = 1;
    config->tracking_traj = 0;
    config->load_tracer_meta = 0;
    config->is_recreate = 0;    // Do not recreate a HDF5 file when it exists
    config->nsteps = 1;         // # steps are saved in each time interval
    config->reduced_tracer = 0; // original tracer

    while ((c = getopt_long (argc, argv, options, long_options, &option_index)) != -1){
        switch (c){
            case 'a':
                strcpy(config->filename_attribute, optarg);
                break;
            case 'b':
                config->tmax = atoi(optarg);
                break;
            case 'c':
                config->collect_data = 0;
                break;
            case 'e':
                config->sort_key_only = 1;
                break;
            case 'f':
                strcpy(config->filename, optarg);
                break;
            case 'g':
                strcpy(config->group_name, optarg);
                break;
            case 'i':
                config->tinterval = atoi(optarg);
                break;
            case 'k':
                config->key_index = atoi(optarg);
                break;
            case 'l':
                config->weak_scale_test = 1;
                config->weak_scale_test_length = atoi(optarg);
                break;
            case 'm':
                config->meta_data = 1;
                strcpy(config->filename_meta, optarg);
                break;
            case 'o':
                strcpy(config->filename_sorted, optarg);
                break;
            case 'p':
                config->multi_tsteps = 1;
                break;
            case 'q':
                config->tracking_traj = 1;
                break;
            case 'r':
                config->load_tracer_meta = 1;
                break;
            case 's':
                config->skew_data = 1;
                break;
            case 't':
                config->local_sort_threaded = 1;
                config->local_sort_threads_num = atoi(optarg);
                break;
            case 'u':
                config->ux_kindex = atoi(optarg);
                break;
            case 'v':
                config->verbose = 1;
                break;
            case 'w':
                config->write_result = 0;
                break;
            case 1:
                strcpy(config->filepath, optarg);
                break;
            case 2:
                strcpy(config->species, optarg);
                break;
            case 3:
                strcpy(config->filename_traj, optarg);
                break;
            case 4:
                config->nptl_traj = atoi(optarg);
                break;
            case 5:
                config->ratio_emax = atof(optarg);
                break;
            case 6:
                config->tmin = atoi(optarg);
                break;
            case 7:
                config->is_recreate = atoi(optarg);
                break;
            case 8:
                config->nsteps = atoi(optarg);
                break;
            case 9:
                config->reduced_tracer = atoi(optarg);
                break;
            case 11:
                config->tstep = atoi(optarg);
                break;
            case 'h':
                if (mpi_rank == 0) {
                    print_help();
                }
                return 1;
            default:
                printf("Error option [%s]\n", optarg);
                exit(-1);
        }
    }
    return 0;
}


/******************************************************************************
 * Print help information.
 ******************************************************************************/
void print_help(){
    char *msg="Usage: %s [OPTION] \n\
               -a name of the attribute file to store sort table  \n\
               -b the particle output maximum time step \n\
               -c whether to collect sorted data \n\
               -e only sort the key  \n\
               -f name of the file to sort \n\
               -g group path within HDF5 file to data set \n\
               -h help (--help)\n\
               -i the particle output time interval \n\
               -k the index key of the file \n\
               -l whether to do weak scale test \n\
               -m whether to use meta data to determine particle position \n\
               -o name of the file to store sorted results \n\
               -p run sorting for multiple time steps \n\
               -q tracking the trajectories of particles\n\
               -r whether to load tracer meta data \n\
               -s the data is in skew shape \n\
               -t whether to use threads in local sorting \n\
               -u the key index of ux \n\
               -v verbose  \n\
               -w won't write the sorted data \n\
               --tmin the particle output minimum time step \n\
               --tstep current time step \n\
               --filepath file path saving the particle tracing data \n\
               --species particle species for sorting \n\
               --filename_traj output file for particle trajectories \n\
               --nptl_traj number of particles for trajectory tracking \n\
               --ratio_emax ratio of Emax of all particles to that of tracked ones \n\
               --is_recreate whether to recreate a HDF5 file \n\
               --nsteps # of steps are save in each time interval \n\
               --reduced_tracer whether to use reduced tracer \n\
               example: ./h5group-sorter -f testf.h5p  -g /testg  -o testg-sorted.h5p -a testf.attribute -k 0 \n";
    fprintf(stdout, msg, "h5group-sorter");
}
