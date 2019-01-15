/*
 * Generally, this file sort the dataset inside a HDF5 group
 * The original code is developed by Bin Dong at LBNL.
 * Modified by Xiaocan Li.
 *
 */

#include "stdlib.h"
#include "hdf5.h"
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "qsort-parallel.h"
#include "mpi_io.h"
#include "vpic_data.h"
#include "configuration.h"
#include "get_data.h"
#include "meta_data.h"
#include "tracked_particle.h"
#include "particle_tags.h"

void set_filenames(int tstep, config_t *config_init, config_t *config);
void set_filenames_reduced(int tstep, char *filepath, char *species,
        char *filename, char *group_name, char *filename_sorted,
        char *filename_attribute, char *filename_meta);

/******************************************************************************
 * Main of the parallel sampling sort
 ******************************************************************************/
int main(int argc, char **argv){
    int mpi_size, mpi_rank;
    double t0, t1;
    int is_help;
    char *final_buff;
    unsigned long long rsize;

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    t0 = MPI_Wtime();

    config_t *config = (config_t *)malloc(sizeof(config_t));
    config_t *config_init = (config_t *)malloc(sizeof(config_t));
    init_configuration(config);
    init_configuration(config_init);

    is_help = get_configuration(argc, argv, mpi_rank, config_init);

    // when -h flag is set to seek help of how to use this program
    if (is_help) {
        free_configuration(config);
        free_configuration(config_init);
        MPI_Finalize();
        return 1;
    }

    copy_configuration(config, config_init);

    int ntf, mtf, tstep, nsteps_tot;
    mtf = config->tmin / config->tinterval;
    ntf = config->tmax / config->tinterval + 1;
    if (config->nsteps == 1) {
        nsteps_tot = ntf;
    } else {
        nsteps_tot = config->tmax;
    }

    // Get the particle tags from sorted-by-energy data of the last time frame,
    // and then sort the tags
    int *tags, qindex;
    int row_size, dataset_num, max_type_size, key_value_type;
    char *tracked_particles, *tracked_particles_sum, *package_data;
    dset_name_item *dname_array;
    hsize_t my_data_size, rest_size;
    if (config->tracking_traj) {
        tags = (int *)malloc(config->nptl_traj * sizeof(int));
        char filename_ene[MAX_FILENAME_LEN];
        tstep = (ntf - 1) * config->tinterval;
        snprintf(filename_ene, MAX_FILENAME_LEN, "%s%s%d%s%s%s",
                config->filepath, "/T.", tstep, "/", config->species,
                "_tracer_energy_sorted.h5p");
        if (config->nsteps != 1) tstep--;
        // get particle tags w.r.t some energy
        get_particle_tags(filename_ene, tstep, config->ratio_emax,
                config->nptl_traj, tags);

        // sort the tags
        qsort(tags, config->nptl_traj, sizeof(int), CompareInt32Value);

        // get the indices for Ux and q in the HDF5 file
        snprintf(config->group_name, MAX_FILENAME_LEN, "%s%d", "/Step#", tstep);
        dname_array = (dset_name_item *)malloc(MAX_DATASET_NUM * sizeof(dset_name_item));
        package_data = get_vpic_pure_data_h5(mpi_rank, mpi_size, filename_ene,
            config->group_name, &row_size, &my_data_size, &rest_size, &dataset_num,
            &max_type_size, &key_value_type, dname_array);
        free(package_data);
        qindex = get_dataset_index("q", dname_array, dataset_num);
        config->ux_kindex = get_dataset_index("Ux", dname_array, dataset_num);

        tracked_particles = (char *)malloc(nsteps_tot * config->nptl_traj * row_size);
        for (int j = 0; j < nsteps_tot*config->nptl_traj*row_size; j++) {
            tracked_particles[j] = 0;
        }
        if (mpi_rank == 0) {
            tracked_particles_sum = (char *)malloc(nsteps_tot * config->nptl_traj * row_size);
            for (int j = 0; j < nsteps_tot*config->nptl_traj*row_size; j++) {
                tracked_particles_sum[j] = 0;
            }
        }
    } // if (config->tracking_traj)

    if (config->multi_tsteps) {
        for (int i = mtf; i < ntf; i++) {
            tstep = i * config->tinterval;
            if (mpi_rank == 0) printf("Time Step: %d\n", tstep);
            if (config->reduced_tracer) {
                set_filenames_reduced(tstep, config->filepath, config->species,
                        config->filename, config->group_name, config->filename_sorted,
                        config->filename_attribute, config->filename_meta);
            } else {
                set_filenames(tstep, config_init, config);
            }
            if (config->nsteps == 1) {
                final_buff = sorting_single_tstep(mpi_size, mpi_rank, config, &rsize);
                if (config->tracking_traj) {
                    get_tracked_particle_info(final_buff, qindex, row_size,
                            rsize, i, ntf, tags, config->nptl_traj, tracked_particles);
                }
                if(config->collect_data == 1) {
                    free(final_buff);
                }
            } else {
                for (tstep = (i-1) * config->tinterval; tstep < i*config->tinterval; tstep++) {
                    snprintf(config->group_name, MAX_FILENAME_LEN, "%s%d", "/Step#",
                            tstep);
                    final_buff = sorting_single_tstep(mpi_size, mpi_rank, config, &rsize);
                    if (config->tracking_traj) {
                        get_tracked_particle_info(final_buff, qindex, row_size,
                                rsize, tstep, nsteps_tot, tags, config->nptl_traj,
                                tracked_particles);
                    }
                    if(config->collect_data == 1) {
                        free(final_buff);
                    }
                }
            }
        }
    } else {
        if (mpi_rank == 0) {
            printf("Input filename: %s\n", config->filename);
            printf("Group name: %s\n", config->group_name);
            printf("Output filename: %s\n", config->filename_sorted);
            printf("Meta data filename: %s\n", config->filename_meta);
        }
        set_filenames(config->tstep, config_init, config);

        final_buff = sorting_single_tstep(mpi_size, mpi_rank, config, &rsize);
        if(config->collect_data == 1) {
            free(final_buff);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (config->tracking_traj) {
        MPI_Reduce(tracked_particles, tracked_particles_sum,
                nsteps_tot*config->nptl_traj*row_size, MPI_CHAR, MPI_SUM, 0,
                MPI_COMM_WORLD);

        /* Save the particle data. */
        if (mpi_rank == 0) {
            save_tracked_particles(config->filename_traj, tracked_particles_sum,
                    nsteps_tot, config->nptl_traj, row_size, dataset_num, max_type_size,
                    dname_array, tags);
        }
        free(tracked_particles);
        if (mpi_rank == 0) {
            free(tracked_particles_sum);
        }
        free(tags);
        free(dname_array);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    if(mpi_rank == 0) {
        printf("Overall time is [%f]s \n", (t1 - t0));
    }

    free_configuration(config);
    free_configuration(config_init);

    MPI_Finalize();
    return 0;
}

/******************************************************************************
 * Set filenames to include file directories
 ******************************************************************************/
void set_filenames(int tstep, config_t *config_init, config_t *config)
{
    char *tempname = (char *)malloc(MAX_FILENAME_LEN * sizeof(char));
    if (config->single_h5) {
        snprintf(config->group_name, MAX_FILENAME_LEN, "%s%d", "/Step#", tstep);
        snprintf(config->subgroup_name, MAX_FILENAME_LEN, "%s%d%s%s%s",
                "/Step#", tstep, "/", config->species, "_tracer");
        snprintf(config->meta_group_name, MAX_FILENAME_LEN, "%s%d%s",
                "/Step#", tstep, "/grid_metadata");
        printf("%s\n", config->subgroup_name);
        printf("%s\n", config->meta_group_name);
    } else {
        snprintf(config->group_name, MAX_FILENAME_LEN, "%s%d", "/Step#", tstep);
    }

    strcpy(tempname, config_init->filename);
    snprintf(config->filename, MAX_FILENAME_LEN, "%s%s%d%s%s",
            config->filepath, "/T.", tstep, "/", tempname);

    strcpy(tempname, config_init->filename_sorted);
    snprintf(config->filename_sorted, MAX_FILENAME_LEN, "%s%s%d%s%s",
            config->filepath, "/T.", tstep, "/", tempname);

    snprintf(config->filename_attribute, MAX_FILENAME_LEN, "%s", "attribute");

    strcpy(tempname, config_init->filename_meta);
    snprintf(config->filename_meta, MAX_FILENAME_LEN, "%s%s%d%s%s",
            config->filepath, "/T.", tstep, "/", tempname);
    free(tempname);
}

/******************************************************************************
 * Set filenames for reduced tracer.
 ******************************************************************************/
void set_filenames_reduced(int tstep, char *filepath, char *species,
        char *filename, char *group_name, char *filename_sorted,
        char *filename_attribute, char *filename_meta)
{
    snprintf(group_name, MAX_FILENAME_LEN, "%s%d", "/Step#", tstep);
    snprintf(filename, MAX_FILENAME_LEN, "%s%s%d%s%s%s", filepath,
            "/T.", tstep, "/", species, "_tracer_reduced_sorted.h5p");
    snprintf(filename_sorted, MAX_FILENAME_LEN, "%s%s%d%s%s%s",
            filepath, "/T.", tstep, "/", species, "_tracer_energy_sorted.h5p");
    snprintf(filename_attribute, MAX_FILENAME_LEN, "%s", "attribute");
    snprintf(filename_meta, MAX_FILENAME_LEN, "%s%s%d%s%s%s", filepath,
            "/T.", tstep, "/grid_metadata_", species, "_tracer_reduced.h5p");
}
