#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "basic_functions.h"
#include "percolation_functions.h"
#include "file_io.h"

int main(int argc, char **argv) {

    clock_t start, end;
    double cpu_time_used, mean_cpu_time=0.;

    unsigned long **pairs, **copy_pairs;
    unsigned long num_pairs;
    unsigned long long *update_time;
    int **bond, **edges, **copy_bond;
    int i, j, x, num_instance, N, E, C, num_iter, pid_id, file_id, avg_k, event_index, max_file_length, sampling_factor;
    double **largest_cluster;
    double removed_edges;
    char dir_path[200], file_name[200], EE_file_name[200], complexity_fname[200];
	char *f_ext = ".csv";
	FILE *f, *g, *h;

    pid_id = time(NULL) * getpid();
    init_genrand64((unsigned)pid_id);

    N = atoi(argv[1]); /* size of network */
    avg_k = atoi(argv[2]); /* mean degree */
    C = atoi(argv[3]); /* maximum cost available */ 
    num_iter = atoi(argv[4]); /* number of iterations */
	num_instance = atoi(argv[5]); /* number of graph instances */

    ////////////////////////////
    E = (int)(0.5 * (double) N * avg_k);

    bond = (int **)malloc(sizeof(int *) * (N + 1));
    bond[0] = (int *)malloc(sizeof(int) * (N + 1));
    for (i = 1; i <= N; i++) {bond[i] = (int *)malloc(sizeof(int) * 1);}

	copy_bond = (int **)malloc(sizeof(int *) * (N + 1));
	copy_bond[0] = (int *)malloc(sizeof(int) * (N + 1));

    edges = (int **)malloc(sizeof(int *) * 2);
    edges[0] = (int *)malloc(sizeof(int) * (E + 1));
    edges[1] = (int *)malloc(sizeof(int) * (E + 1));

    update_time = (unsigned long long *)calloc(E + 1, sizeof(unsigned long long)); 

    largest_cluster = (double **)malloc(sizeof(double *) * 6);
	for (i = 0; i <= 5; i++) {largest_cluster[i] = (double *)calloc(E + 1, sizeof(double));}

    start = clock();
    for (x = 0; x < num_instance; x++) {
        generate_ERgraph(N, bond, E);

        for (i = 0; i < num_iter; i++) {
           
            copy_network (N, bond, copy_bond);
            clean_time_stamp(update_time, E);
            efficient_pair_removal(C, N, copy_bond, edges, E, update_time, (int) (7 * pow(N, 0.44))); 
            modified_NZ_algorithm(N, edges, largest_cluster);

            for (j = 1; j <= N; j++) {free(copy_bond[j]);}

            /* some tedious jobs to get file names*/
            if (C > 0) {
                sprintf(dir_path, "./data/N%d/k%d/C%d/", N, avg_k, C);
                file_id = get_max_num_from_filenames(dir_path, f_ext) + 1;
                sprintf(file_name, "./data/N%d/k%d/C%d/full_data_%d.csv", N, avg_k, C, file_id);
            } else {
                sprintf(dir_path, "./data/N%d/k%d/C%d/", N, avg_k, N);
                file_id = get_max_num_from_filenames(dir_path, f_ext) + 1;
                sprintf(file_name, "./data/N%d/k%d/C%d/full_data_%d.csv", N, avg_k, N, file_id);
            }
            

			f = fopen(file_name, "w");
			for (j = 1; j <= edges[0][0]; j+=1) {
				removed_edges = (double) j / edges[0][0];
				fprintf(f, "%.10f %.10f %.10f %.10f %.10f %.10f %d %d %lld\n", 
						removed_edges, 
						largest_cluster[1][j], /* 1st moment of percolation strength */
						largest_cluster[2][j], /* 2st moment of percolation strength */
						largest_cluster[3][j], /* 3st moment of percolation strength */
						largest_cluster[4][j], /* 4st moment of percolation strength */
						largest_cluster[5][j], /* average cluster size */
						edges[0][j], /* edge1*/
						edges[1][j], /* edge2 */
						update_time[j] /* time */
						);
			}
			fclose(f);
        }    
    }
    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("# mean simulation time for N=%d C=%d calculated from %d instances and %d iterations : %g\n", N, C, num_instance, num_iter, cpu_time_used/(num_instance * num_iter)); fflush(stdout);
    

     /* Free memories */
    for (i = 0; i <= N; i++) {free(bond[i]);}
    free(bond);
    free(copy_bond[0]);
	free(copy_bond);
    free(edges[0]);
    free(edges[1]);
    free(edges);
    for (i = 0; i <= 5; i++) {free(largest_cluster[i]);}
    free(largest_cluster);    
    free(update_time);
    
}
