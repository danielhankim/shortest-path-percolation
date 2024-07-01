#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "percolation_functions.h"
#include "basic_functions.h"

int simple_bfs(int C, int N, int **bond, int source, int target, int *vec, int *tmp_vec, int *visited, int *reset, int **dag, int *nr_sp) {
    
    /*
     * This function implements the BFS algorithm especially designed for the 
     * shortest-path percolation model. To be specific, we focus on whether
     * we can reach the target node from the source node with the traveling cost C.
     *
     * This function returns one of the followings:
     *     control = 0: we should still run the bfs 
     *     control = 1: we reach the maximal traveling cost 
     *                  before finding the target node
     *     control = 2: we find the target node
     */


    int i, j, m, n, distance, control;

    vec[0] = 1;
    vec[1] = source;
    visited[source] = 0;
    dag[0][source] = 0;
    distance = 0;
    nr_sp[source] = 1;

    reset[0] = 1;
    reset[1] = source;
    control = 0;

    while (vec[0] > 0 && control == 0) {
        
        tmp_vec[0] = 0;
        distance += 1;

        if (distance == C) {control = 1;} /* We're out of traveling cost */
        
        for (i = 1; i <= vec[0]; i++) {
            n = vec[i];
            for (j = 1; j <= bond[0][n]; j++) {
                m = bond[n][j];
                
                if (m == target) {control = 2;} /* We found our target node */
                
                if (visited[m] < 0) { /* If node m is not visited yet */
                    tmp_vec[0] += 1;
                    tmp_vec[tmp_vec[0]] = m;
                    visited[m] = visited[n] + 1; 
                    reset[0] += 1;
                    reset[reset[0]] = m;
                }

                if (visited[m] == distance) {
                    dag[0][m] += 1;
                    dag[m][dag[0][m]] = n;
                    nr_sp[m] += nr_sp[n];
                }

            }
        }
        for (i = 0; i <= tmp_vec[0]; i++) {vec[i] = tmp_vec[i];}
    }

    return control;
}



int find_all_possible_pairs(int C, int N, int **bond, unsigned long **pairs) {

    /*
     * We use the breadth-first search (BFS) to find all the possible pairs between 
     * the source and target node. While performing BFS, we construct a 
     * directed acyclic graph (DAG) that contains the entire information 
     * about all eventual shortest pathes between the source and target.
     * In addition, this function returns the number of all possible pairs.
     */

    unsigned long num_pairs;
    int **dag, **list_pairs;
    int *vec, *tmp_vec, *visited, *reset, *nr_sp;
    int i, j, n, control, source, target;

    vec = (int *)malloc(sizeof(int) * (N + 1));
    tmp_vec = (int *)malloc(sizeof(int) * (N + 1));
    visited = (int *)malloc(sizeof(int) * (N + 1)); /* check whether we visited the node */
    reset = (int *)malloc(sizeof(int) * (N + 1));
    dag = (int **)malloc(sizeof(int *) * (N + 1));
    dag[0] = (int *)malloc(sizeof(int) * (N + 1));
    nr_sp = (int *)malloc(sizeof(int) * (N + 1)); /* list to save the number of shortest path */
    list_pairs = (int **)malloc(sizeof(int *) * (N + 1));
    list_pairs[0] = (int *)malloc(sizeof(int) * (N + 1));    

    for (i = 1; i <= N; i++) {
        dag[i] = (int *)malloc(sizeof(int) * (bond[0][i] + 1));
        visited[i] = -1;
        dag[0][i] = 0;
        nr_sp[i] = 0;
        list_pairs[i] = (int *)malloc(sizeof(int) * 1);
        list_pairs[0][i] = 0;
    }

    num_pairs = 0;
    target = -1; /* To find all the possible pairs, we don't designate a target */

    for (source = 1; source <= N; source++) {
        control = simple_bfs (C, N, bond, source, target, vec, tmp_vec, visited, reset, dag, nr_sp);
        for (i = 1; i <= reset[0]; i++) {
            n = reset[i];
            if (source < n && find_node (source, n, list_pairs) < 0) {
                list_pairs[0][source] += 1;
                list_pairs[source] = (int *)realloc(list_pairs[source], sizeof(int) * (list_pairs[0][source] + 1));
                list_pairs[source][list_pairs[0][source]] = n; /* save pair (source, n) */
                num_pairs += 1;
            }
        }
        /* reset some data */
        for (i = 1; i <= reset[0]; i++) {
            n = reset[i];
            visited[n] = -1;
            dag[0][n] = 0;
            nr_sp[n] = 0;
        }
    }


    pairs[0] = (unsigned long *)realloc(pairs[0], sizeof(unsigned long) * (num_pairs + 1));
    pairs[1] = (unsigned long *)realloc(pairs[1], sizeof(unsigned long) * (num_pairs + 1));

    num_pairs = 0;

    for (i = 1; i <= N; i++) {
        for (j = 1; j <= list_pairs[0][i]; j++) {
            num_pairs += 1;
            pairs[0][num_pairs] = (unsigned long) i;
            pairs[1][num_pairs] = (unsigned long) list_pairs[i][j];
        }
    }
    pairs[0][0] = num_pairs;

    free(vec);
    free(tmp_vec);
    free(visited);
    free(reset);
    for (i = 0; i <= N; i++) {free(dag[i]);}
    free(dag);
    free(nr_sp);
    for (i = 0; i <= N; i++) {free(list_pairs[i]);}
    free(list_pairs);   

    return num_pairs;
}

void remove_pair_wo_pair_info(int C, int N, int **bond, int **edges, int E, unsigned long long *update_time) {

    /*
     * 1. We randomly choose two nodes from the network and 
     * find all shortest paths with length shorter or equal to C.
     * 2. If no shortest paths are found at point 1, we proceed to point 4.
     * 3. If at least one shortest path is found at point 1, we select one at random.
     * 4. We increase t <- t+1 and go back to point 1.
     */ 

    unsigned long long tau;
    int p, n, m, i, j, s, t, v1, v2, count_E, tmp;
    int **dag;
    int *path, *vec, *tmp_vec, *visited, *reset, *nr_sp;


    path = (int *)malloc((N + 1) * sizeof(int)); 
    vec = (int *)malloc(sizeof(int) * (N + 1));
    tmp_vec = (int *)malloc(sizeof(int) * (N + 1));
    visited = (int *)malloc(sizeof(int) * (N + 1)); /* check whether we visited the node */
    reset = (int *)malloc(sizeof(int) * (N + 1));
    dag = (int **)malloc(sizeof(int *) * (N + 1));
    dag[0] = (int *)malloc(sizeof(int) * (N + 1));
    nr_sp = (int *)malloc(sizeof(int) * (N + 1)); /* list to save the number of shortest path */
    
    
    for (i = 1; i <= N; i++) {
        dag[i] = (int *)malloc(sizeof(int) * (bond[0][i] + 1));
        visited[i] = -1;
        dag[0][i] = 0;
        nr_sp[i] = 0;
    }

    count_E = 0;
    tau = 1;
    while(1) { 
        /* randomly choose two nodes */
        while (1) { 
            n = (int)(genrand64_real3() * (double)N) + 1;
            m = (int)(genrand64_real3() * (double)N) + 1;
            if (n != m) {
                if (n > m) {tmp = n; n = m; m = tmp;} /* swap variables: Do I need this?*/
                break;
            }
            
        }
        if (C == -1) {
            p = get_path_of_pair (n, m, N, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        } else {
            p = get_path_of_pair (n, m, C, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        }
        if (p > 0) {
            // printf("p = %d\n", p);
            // printf("removing pair: (%d, %d)\n", n, m);
            for (j = 1; j < path[0]; j++) {
                s = path[j];
                t = path[j + 1];
                // printf("removing edge (%d, %d)\n", s, t);
                // printf("# of bonds of %d = %d\n", s, bond[0][s]);
                // printf("# of bonds of %d = %d\n", t, bond[0][t]);

                // for (i = 1; i <= bond[0][s]; i++) {
                //     printf("bond[%d][%d]=%d\n", s, i, bond[s][i]);
                // } 
                v1 = find_node (s, t, bond);
                bond[s][v1] = bond[s][bond[0][s]];
                bond[0][s] -= 1;
                // for (i = 1; i <= bond[0][s]; i++) {
                //     printf("bond[%d][%d]=%d\n", s, i, bond[s][i]);
                // }

                v2 = find_node (t, s, bond);
                bond[t][v2] = bond[t][bond[0][t]];
                bond[0][t] -= 1;
                
                count_E += 1;
                if (s < t) {
                    edges[0][count_E] = s;
                    edges[1][count_E] = t;
                } else {
                    edges[0][count_E] = t;
                    edges[1][count_E] = s;
                }
                update_time[count_E] += tau;
                // printf("count_E=%d\n", count_E);
                // printf("E = %d\n", E);
            }
        }


        // printf("Line Number %s->%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
        tau += 1;

        if (count_E == E) {break;}

    }
    // printf("Line Number %s->%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
    edges[0][0] = count_E;

    free(path);
    free(vec);
    free(tmp_vec);
    free(visited);
    free(reset);
    for(i = 0; i <= N; i++) {free(dag[i]);}
    free(dag);
    free(nr_sp);
}

void remove_pair_wo_pair_info_v2(int C, int N, int **bond, int **edges, int E, unsigned long long *update_time, int threshold) {

    /*
     * 1. We randomly choose two nodes from the network and 
     * find all shortest paths with length shorter or equal to C.
     * 2. If no shortest paths are found at point 1, we proceed to point 4.
     * 3. If at least one shortest path is found at point 1, we select one at random.
     * 4. We increase t <- t+1 and go back to point 1.
     * 5. If we fail to find a path between for a given pair of nodes more than `threshold`, we terminate the program
     */ 

    unsigned long long tau;
    int p, n, m, i, j, s, t, v1, v2, count_E, tmp, num_failures, source, target;
    int **dag;
    int *path, *vec, *tmp_vec, *visited, *reset, *nr_sp;


    path = (int *)malloc((N + 1) * sizeof(int)); 
    vec = (int *)malloc(sizeof(int) * (N + 1));
    tmp_vec = (int *)malloc(sizeof(int) * (N + 1));
    visited = (int *)malloc(sizeof(int) * (N + 1)); /* check whether we visited the node */
    reset = (int *)malloc(sizeof(int) * (N + 1));
    dag = (int **)malloc(sizeof(int *) * (N + 1));
    dag[0] = (int *)malloc(sizeof(int) * (N + 1));
    nr_sp = (int *)malloc(sizeof(int) * (N + 1)); /* list to save the number of shortest path */
    
    
    for (i = 1; i <= N; i++) {
        dag[i] = (int *)malloc(sizeof(int) * (bond[0][i] + 1));
        visited[i] = -1;
        dag[0][i] = 0;
        nr_sp[i] = 0;
    }

    count_E = 0;
    tau = 1;
    num_failures = 0;
    while(num_failures < threshold) { 
        /* randomly choose two nodes that are not the same */
        while (1) { 
            n = (int)(genrand64_real3() * (double)N) + 1;
            m = (int)(genrand64_real3() * (double)N) + 1;
            if (n != m) {break;} 
            // if (n > m) {tmp = n; n = m; m = tmp;} /* swap variables: Do I need this?*/
        }
        /*to expediate the path searching process*/
        source = n;
        target = m;
        if (bond[0][source] > bond[0][target]) {
            n = target;
            m = source;
        }
        if (C == -1) {
            p = get_path_of_pair (n, m, N, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        } else {
            p = get_path_of_pair (n, m, C, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        }

        if (p > 0) {
            num_failures = 0;
            for (j = 1; j < path[0]; j++) {
                s = path[j];
                t = path[j + 1];
                
                v1 = find_node (s, t, bond);
                bond[s][v1] = bond[s][bond[0][s]];
                bond[0][s] -= 1;

                v2 = find_node (t, s, bond);
                bond[t][v2] = bond[t][bond[0][t]];
                bond[0][t] -= 1;
                
                count_E += 1;
                if (s < t) {
                    edges[0][count_E] = s;
                    edges[1][count_E] = t;
                } else {
                    edges[0][count_E] = t;
                    edges[1][count_E] = s;
                }
                update_time[count_E] += tau;
            }
        } else { 
            /* if we fail to find a path */
            num_failures += 1;
        }

        tau += 1;
    }

    for (i = count_E+1; i <= E; i++) {
        update_time[i] += tau;
    }
    edges[0][0] = count_E;

    free(path);
    free(vec);
    free(tmp_vec);
    free(visited);
    free(reset);
    for(i = 0; i <= N; i++) {free(dag[i]);}
    free(dag);
    free(nr_sp);
}

void remove_pair_w_pair_info(int C, int N, int **bond, unsigned long **pairs, int **edges) {
    
    /*
     * This is the modified version of the function `remove_pair_wo_pair_info`.
     * The main purpose of the modification is to accelerate the speed of the
     * simulation. For this, we find all pairs of nodes which have at least
     * one shortest path equal to or shorter than C at the beginning. This 
     * modification makes the program to consume a lot of memory. However, 
     * this allows us to save time when finding the available pairs of nodes. 
     * Regarding memory consumption, I highly recommend not to use this function
     * with C = -1 on large networks.
     */

    unsigned long num_pairs, q;
    unsigned long *order;

    int p, n, m, i, j, E, s, t, v;
    int **dag;
    int *path, *vec, *tmp_vec, *visited, *reset, *nr_sp;

    num_pairs = pairs[0][0];

    path = (int *)malloc((N + 1) * sizeof(int)); 
    vec = (int *)calloc(N + 1, sizeof(int));
    tmp_vec = (int *)malloc(sizeof(int) * (N + 1));
    visited = (int *)malloc(sizeof(int) * (N + 1)); /* check whether we visited the node */
    reset = (int *)malloc(sizeof(int) * (N + 1));
    nr_sp = (int *)malloc(sizeof(int) * (N + 1)); /* list to save the number of shortest path */
    order = (unsigned long *)malloc(sizeof(unsigned long) * (num_pairs + 1));

    dag = (int **)malloc((N + 1) * sizeof(int *));
	dag[0] = (int *)malloc((N + 1) * sizeof(int));
	for (i = 1; i <= N; i++) {dag[i] = (int *)malloc((bond[0][i] + 1) * sizeof(int));}
    
    for (i = 1; i <= N; i++) {
        visited[i] = -1;
        dag[0][i] = 0;
        nr_sp[i] = 0;
    }

    for (q = 1; q <= num_pairs; q++) {order[q] = q;}
    randomize_vector(order, num_pairs);

    for (E = 0, q = 1; q <= num_pairs; q++) {
        n = (int)pairs[0][order[q]];
        m = (int)pairs[1][order[q]];

        if (C == -1) {
            p = get_path_of_pair (n, m, N, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        } else {
            p = get_path_of_pair (n, m, C, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        }

        if (p > 0) {
            for (j = 1; j < path[0]; j++) {
                s = path[j];
                t = path[j + 1];
                    
                v = find_node (s, t, bond);
                bond[s][v] = bond[s][bond[0][s]];
                bond[0][s] -= 1;
                
                v = find_node (t, s, bond);
                bond[t][v] = bond[t][bond[0][t]];
                bond[0][t] -= 1;

                E += 1;
                if (s < t) {
                    edges[0][E] = s;
                    edges[1][E] = t;
                } else {
                    edges[0][E] = t;
                    edges[1][E] = s;
                }
            }
        }

    }

    edges[0][0] = E;

    free(order);
    free(path);
    free(vec);
    free(tmp_vec);
    free(visited);
    free(reset);
    for(i = 0; i <= N; i++) {free(dag[i]);}
    free(dag);
    free(nr_sp);
}

void remove_pair_w_pair_info_v2(int C, int N, int **bond, unsigned long **pairs, int **edges, unsigned long long *update_time) {
    
    /*
     * This is the modified version of the function `remove_pair_w_pair_info`.
     * The main purpose of this modification is to save the "update time" 
     * information properly while using the available pair information. 
     * Instead of choosing a random pair of nodes with shortest path length  
     * shorter than or euqal to C and saving the time we actually find such pair, 
     * we just sample the update time from the geometric distribution. 
     * Still, this approach might consume huge memory when C = -1 on large networks
     * and highly recommend not to use this option with large networks.
     */

    unsigned long long interval;
    unsigned long num_pairs, q;
    unsigned long *order;

    int p, n, m, i, j, E, s, t, v, control, source, target;
    int **dag;
    int *path, *vec, *tmp_vec, *visited, *reset, *nr_sp;
    double link_prob;

    num_pairs = pairs[0][0];
    path = (int *)malloc((N + 1) * sizeof(int)); 
    vec = (int *)calloc(N + 1, sizeof(int));
    tmp_vec = (int *)malloc(sizeof(int) * (N + 1));
    visited = (int *)malloc(sizeof(int) * (N + 1)); /* check whether we visited the node */
    reset = (int *)malloc(sizeof(int) * (N + 1));
    nr_sp = (int *)malloc(sizeof(int) * (N + 1)); /* list to save the number of shortest path */

    dag = (int **)malloc((N + 1) * sizeof(int *));
	dag[0] = (int *)malloc((N + 1) * sizeof(int));
	for (i = 1; i <= N; i++) {dag[i] = (int *)malloc((bond[0][i] + 1) * sizeof(int));}
    
    for (i = 1; i <= N; i++) {
        visited[i] = -1;
        dag[0][i] = 0;
        nr_sp[i] = 0;
    }



    interval = 0;
    E = edges[0][0];
    while (pairs[0][0] > 0) { 
        
        /* 1. We extract a random number from the geometric distribution */
        num_pairs = pairs[0][0];
        link_prob = 2.0 * (double)num_pairs / (double)N / (double)(N - 1);
        // printf("link_prob = %.20f\n", link_prob);
        interval += geometric_distribution(link_prob);

        /* 2. We select a random pair from the pair list 
         * and apply the BFS to find and remove the edges
         * from one of the shortest-paths that has length
         * shorter than or equal to C.
         */

        q = (unsigned long) (genrand64_real3() * num_pairs + 1.); // I can also start by just randomizing the pair list?

        n = (int)pairs[0][q];
        m = (int)pairs[1][q];

        /*to expediate the path searching process*/
        source = n;
        target = m;
        if (bond[0][source] > bond[0][target]) {
            n = target;
            m = source;
        }

        if (C == -1) {
            p = get_path_of_pair (n, m, N, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        } else {
            p = get_path_of_pair (n, m, C, N, bond, path, vec, tmp_vec, visited, reset, dag, nr_sp);
        }
        if (p > 0) { // if there exists a path!
            for (j = 1; j < path[0]; j++) {
                s = path[j];
                t = path[j + 1];
                    
                v = find_node (s, t, bond);
                bond[s][v] = bond[s][bond[0][s]];
                bond[0][s] -= 1;
                
                v = find_node (t, s, bond);
                bond[t][v] = bond[t][bond[0][t]];
                bond[0][t] -= 1;

                E += 1;
                if (s < t) {
                    edges[0][E] = s;
                    edges[1][E] = t;
                } else {
                    edges[0][E] = t;
                    edges[1][E] = s;
                }
                update_time[E] += interval;
            }

        	/* 3. Verify that the selected pair in step 2. can still belong to the
         	* pair list by measuring the shortest-path distance. If such a distance
         	* is larger than C, we remove the pair from the pair list.
         	*/

        	control = simple_bfs(C, N, bond, n, m, vec, tmp_vec, visited, reset, dag, nr_sp);
        	// if (control == 1) {remove_and_resize(pairs, q);}
            if (control == 1) {adjust_pair(pairs, q);}

        } else { 
            // we could not find a path! thus we remove the pair
			// remove_and_resize(pairs, q);
            adjust_pair(pairs, q);
		}

        /* reset some values */
        for (i = 1; i <= reset[0]; i++) {
            n = reset[i];
            visited[n] = -1;
            dag[0][n] = 0;
            nr_sp[n] = 0;
        }

        // if (E == num_edges) {break;} /* If we use all the edges, we end the program */
    }

    edges[0][0] = E;

    free(path);
    free(vec);
    free(tmp_vec);
    free(visited);
    free(reset);
    for(i = 0; i <= N; i++) {free(dag[i]);}
    free(dag);
    free(nr_sp);
}

int get_path_of_pair(int source, int target, int C, int N, int **bond, int *path, int *vec, int *tmp_vec, int *visited, int *reset, int **dag, int *nr_sp) {
    int i, n, control;

    control = simple_bfs(C, N, bond, source, target, vec, tmp_vec, visited, reset, dag, nr_sp);
    path[0] = 0;
    if (control == 2) {sample_random_path(C, N, dag, source, target, path, nr_sp);}

    /* reset some values */
    for (i = 1; i <= reset[0]; i++) {
        n = reset[i];
        visited[n] = -1;
        dag[0][n] = 0;
        nr_sp[n] = 0;
    }

    return path[0];
}

void sample_random_path(int C, int N, int **dag, int source, int target, int *path, int *nr_sp) {
    
    int i, n, found, control;
    unsigned long q, T;
    
    control = 0;
    path[0] = 1;
    path[1] = target;

    while (control == 0) {
        n = path[path[0]];
        if (dag[0][n] > 0) {
            T = 0;
            for (i = 1; i <= dag[0][n]; i++) {
                T += nr_sp[dag[n][i]];
            }

            q = (unsigned long)(genrand64_real3() * (double) T) + 1;
            if (q > T) {q = 1;}
            
            i = 0;
            T = 0;
            found = 0;

            while (found == 0) {
                i += 1;
                T += nr_sp[dag[n][i]];
                if (q <= T) {
                    found = dag[n][i];
                }
            }

            path[0] += 1;
            path[path[0]] = found;

            if (path[0] == C + 1) {control = 1;} // if path length is longer than cost
            if (path[path[0]] == source) {control = 1;} // if we reach the source
        } else {
            path[0] = 0;
            return;
        }
    }
    return;
}

void randomize_vector(unsigned long *vector, unsigned long size) { // suspected bottleneck
    
    unsigned long i, q, k, tmp;
    
    for (i = 1; i <= 2 * size; i++) {
        q = (unsigned long)(genrand64_real3() * (double)size) + 1;
        if (q > size) {q =1;}
        k = (unsigned long)(genrand64_real3() * (double)size) + 1;
        if (k > size) {k =1;}
        tmp = vector[q];
        vector[q] = vector[k];
        vector[k] = tmp;
    }
}

void remove_and_resize(unsigned long **pairs, unsigned long q) {
    unsigned long i, j, num_pairs;
    unsigned long **cpy_pairs;

    num_pairs = pairs[0][0];

    cpy_pairs = (unsigned long **)calloc(2, sizeof(unsigned long *));
    cpy_pairs[0] = (unsigned long *)calloc(num_pairs, sizeof(unsigned long));
    cpy_pairs[1] = (unsigned long *)calloc(num_pairs, sizeof(unsigned long));

    for (i = 1, j = 1; i <= num_pairs; i++) {
        if (i != q) {
            cpy_pairs[0][j] = pairs[0][i];
            cpy_pairs[1][j] = pairs[1][i];
            j += 1;
        }
    }

    for (i = 1; i <= num_pairs - 1; i++) {
        pairs[0][i] = cpy_pairs[0][i];
        pairs[1][i] = cpy_pairs[1][i];
    }

    pairs[0] = (unsigned long *)realloc(pairs[0], num_pairs * sizeof(unsigned long));
    pairs[1] = (unsigned long *)realloc(pairs[1], num_pairs * sizeof(unsigned long));
    pairs[0][0] -= 1;

    free(cpy_pairs[0]);
    free(cpy_pairs[1]);
    free(cpy_pairs);
}

void adjust_pair(unsigned long **pairs, unsigned long q) {
    unsigned long node1, node2;

    node1 = pairs[0][q];
    node2 = pairs[0][q];

    // bring the last pairs to the position of used pair
    pairs[0][q] = pairs[0][pairs[0][0]]; 
    pairs [1][q] = pairs[1][pairs[0][0]];

    // send the used pair to the end of the list
    pairs[0][pairs[0][0]] = node1;
    pairs[1][pairs[0][0]] = node2;

    pairs[0][0] -= 1;
}


void tree_find_root (int n, int **root, int *res) {

    if (root[0][n] == n) {
        res[0] = root[0][n];
        res[1] = root[1][n];
        return;
    }

    tree_find_root (root[0][n], root, res);
  
}

void modified_NZ_algorithm (int N, int **edges, double **largest_cluster) {

    int **root;
    int *res;
    int i, j, n, m, e, rn, rm, sn, sm, gcc, count, control, larg;
    double av, av2, tmp;
    
    av2 = N - 1;

    /* allocate memory */
    root = (int **)malloc(sizeof(int *) * 2);
    root[0] = (int *)malloc(sizeof(int) * (N + 1));
    root[1] = (int *)malloc(sizeof(int) * (N + 1));
    res = (int *)malloc(sizeof(int) * 2); /* helper array to save root node & cluster size info */

    for (i = 1; i <= N; i++) {
        root[0][i] = i; /* root node id */
        root[1][i] = 1; /* size of the cluster which shares the same root node */
    }

    larg = 1;
    /* it seems there is an indexing issue when applying NZ algorithm */
    for (e = edges[0][0]; e >= 1; e--) { 
        /* choose two nodes at the end of an edge */
        n = edges[0][e]; 
        m = edges[1][e];

        /* find the root node and the size of the clusters with the same root node */
        tree_find_root(n, root, res);
        rn = res[0]; 
        sn = res[1];
        tree_find_root(m, root, res); 
        rm = res[0];
        sm = res[1];

        av2 += (double)(larg) * (double) (larg);

        /* 
         * Merge two clusters related to node n, and node m
         * only if the tow nodes do not share the same root node
         */

        if (rn != rm) {
            
            av2 -= (double) (sn) * (double) (sn);
            av2 -= (double) (sm) * (double) (sm);

            if (sn > sm) {
                root[0][rm] = rn;
                root[1][rm] = 1; /* Is this proper? YES*/ 
                root[1][rn] = sn + sm;
            } else {
                root[0][rn] = rm;
                root[1][rn] = 1; 
                root[1][rm] = sn + sm;
            }
            if (sn + sm > larg) {larg = sn + sm;}

        }

        tmp = (double)larg / (double) N;
        
        largest_cluster[1][e] = tmp; /* 1st moment */
        largest_cluster[2][e] = tmp * tmp; /* 2nd moment */
        largest_cluster[3][e] = tmp * tmp * tmp; /* 3rd moment*/
        largest_cluster[4][e] = tmp * tmp * tmp * tmp; /* 4th moment */

        av = (double) (N) - (double) (larg); /* number of nodes belonging to finite clusters */
        
        if (rm != rn) {av2 += (double) (sn + sm) * (double) (sn + sm);}
        av2 -= (double) (larg) * (double) (larg);
        if (av > 0) {largest_cluster[5][e] = av2 / av;} /* average cluster size */
    }

    free(root[0]);
    free(root[1]);
    free(root);
    free(res);
}

void efficient_pair_removal(int C, int N, int **bond, int **edges, int E, unsigned long long *update_time, int threshold) {

    unsigned long **pairs;
    unsigned long num_pairs;

    /* phase 1 */
    remove_pair_wo_pair_info_v2 (C, N, bond, edges, E, update_time, threshold);

    /* phase 2 - (a)*/
    pairs = (unsigned long **)malloc(sizeof(unsigned long *) * 2);
    pairs[0] = (unsigned long *)malloc(sizeof(unsigned long) * 1);
    pairs[1] = (unsigned long *)malloc(sizeof(unsigned long) * 1);
    num_pairs = find_all_possible_pairs (C, N, bond, pairs);
    remove_pair_w_pair_info_v2 (C, N, bond, pairs, edges, update_time);

    free(pairs[0]); 
    free(pairs[1]); 
    free(pairs);
}


int find_event_index(double **largest_cluster, int num_edges) {

    int i, event_index;
    double max_gap, gap;

    max_gap = -1.;
    for (i = 1; i <= num_edges - 1; i++) {
        gap = largest_cluster[1][i] - largest_cluster[1][i+1];
        if (gap > max_gap) {
            max_gap = gap;
            event_index = i;
        }
    }
    return event_index;
}

int find_transition_index(double **largest_cluster, double criteria, int start) {
    int index = start;
    double value;
    
    while(1) {
        value = largest_cluster[1][index];
        if (value < criteria) {break;}
        index += 1;
    }
    index -= 1;
    return index;
}

unsigned long long find_longest_shortest_path_length(int N, int avg_k, unsigned long long *update_time) {
    int max_index, i;
    unsigned long long time_old, time_new, path_length, max_length;

    max_index = (int) (N * avg_k * 0.5);
    time_old = 0;
    path_length = 1;
    max_length = 1;

    for (i = 1; i <= max_index; i++) {
        time_new = update_time[i];
        if (time_old == time_new) {
            path_length += 1;
        } else {
            if (path_length > max_length) {
                max_length = path_length;
            }
            path_length = 1;
        }
        time_old = time_new;
    } 

    /* check one more time */
    if (path_length > max_length) {
        max_length = path_length;
    }

    return max_length;
}
