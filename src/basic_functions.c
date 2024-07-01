#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "basic_functions.h"

int find_node (int i, int j, int **bond) {
    int k;
    for (k = 1; k <= bond[0][i]; k++) {
        if (bond[i][k] == j){
            return k;
        }
    }
    return -1;
}

void generate_ERgraph (int N, int **bond, int E) {
    int i, n, m, control;
    for (i = 1; i <= N; i++) {bond[0][i] = 0;}

    for (i = 1; i <= E; i++) {
      
        n = (int)(genrand64_real3() * (double) N) + 1;
        m = (int)(genrand64_real3() * (double) N);
        control = 0;

        while (control == 0) {
            m += 1;
            if (m > N) {m = 1;}
            if (n!=m  && find_node (n, m, bond) < 0) {control = 1;}
        }

        bond[0][n] += 1;
        bond[n] = (int *)realloc(bond[n], (bond[0][n] + 1) * sizeof(int));
        bond[n][bond[0][n]] = m;
        bond[0][m] += 1;
        bond[m] = (int *)realloc(bond[m], (bond[0][m] + 1) * sizeof(int));
        bond[m][bond[0][m]] = n;
    }
}

void generate_2DSQL_w_PB(int L, int **bond) {
    int i, j, idx; // i: x-coordinate, y: y-coordinate

    for (i = 1; i <= L*L; i++) {
        bond[i] = (int *)realloc(bond[i], sizeof(int) * 5);
        bond[0][i] = 4;
    }

    for (i = 1; i <= L; i++) {
        for (j = 1; j <= L; j++) {
            idx = (i - 1) + (j - 1) * L + 1;
            bond[idx][1] = L*(j-1)+((i-1)+1+L)%L + 1; // right
            bond[idx][2] = L*(j-1)+((i-1)-1+L)%L + 1; // left
            bond[idx][3] = L*(((j-1)+1+L)%L)+(i-1) + 1; // up
            bond[idx][4] = L*(((j-1)-1+L)%L)+(i-1) + 1; // down
        }
    }
}


void copy_network (int N, int **origin_bond, int **copy_bond) {
	int i, j;
	for (i = 1; i <= N; i++) {
		copy_bond[i] = (int *)malloc(sizeof(int) * (origin_bond[0][i] + 1));
        copy_bond[0][i] = origin_bond[0][i]; // degree
		for (j = 1; j <= origin_bond[0][i]; j++) {
			copy_bond[i][j] = origin_bond[i][j]; // adj nodes
		}
	}
    /* check validity*/
    // for (i = 1; i <= N; i++) {
    //     for (j = 1; j <= copy_bond[0][i]; j++) {
    //         printf("copy_bond[%d][%d] = %d\n", i, j, copy_bond[i][j]);
    //     }
    // }
}

void copy_list_of_pair (unsigned long num_pairs, unsigned long **origin_pairs, unsigned long **copy_pairs) {
    unsigned long i;
    copy_pairs[0] = (unsigned long *)malloc(sizeof(unsigned long) * (num_pairs + 1));
    copy_pairs[1] = (unsigned long *)malloc(sizeof(unsigned long) * (num_pairs + 1));
    for (i = 1; i <= num_pairs; i++) {
        copy_pairs[0][i] = origin_pairs[0][i];
        copy_pairs[1][i] = origin_pairs[1][i];
    }
    copy_pairs[0][0] = num_pairs;


    // chunk to check whether the copy is successful
    // for (i = 1; i <= num_pairs; i++) {
    //     if (copy_pairs[0][i] != origin_pairs[0][i]) {
    //         printf("copy is not done properly!\n");
    //         break;
    //     }
    //     if (copy_pairs[1][i] != origin_pairs[1][i]) {
    //         printf("copy is not done properly!\n");
    //         break;
    //     }
    // }
}

unsigned long long geometric_distribution(double P) {
    
    double u;
    unsigned long long interval;

    if (P <=0 ) {return -1;}
    if (P == 1 ) {return 1;}

    u = genrand64_real3();
    interval = (unsigned long long) (1 + floor(log(u) / log(1. - P)));

    return interval;
}

int warning_message(int find_pairs, int C) {
    
    if (find_pairs != 0 && find_pairs != 1) {
        printf("Error!\nCheck your pair-search parameter agian! We received %d\n!", find_pairs);
        return -1;
    }
    if (C == -1 && find_pairs == 1) {
        printf("Warning!\nFinding all pairs with parameter C = N might consume all the available memory you have!\n");
        printf("We will run the simulation without finding the pairs in advance!\n");
        return -1;
    }

    return 0;
}

void clean_time_stamp(unsigned long long *time_stamp, int E) {
    int i;

    for (i = 1; i <= E; i++) {
        time_stamp[i] = 0;
    }
}

int get_sampling_factor(int N, int criteria) {
    /* utility function to sub-sample the full phase trajectory */
    int file_length = 2 * N;
    int sampling_factor = 1;

    while ((int) (file_length / sampling_factor) >= criteria) {
        sampling_factor *= 2;
    }

    return sampling_factor;
}


