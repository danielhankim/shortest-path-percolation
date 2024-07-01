int simple_bfs (int C, int N, int **bond, int source, int target, int *vec, int *tmp_vec, int *visited, int *reset, int **dag, int *nr_sp);
int find_all_possible_pairs (int C, int N, int **bond, unsigned long **pairs);
void remove_pair_w_pair_info (int C, int N, int **bond, unsigned long **pairs, int **edges);
void remove_pair_w_pair_info_v2 (int C, int N, int **bond, unsigned long **pairs, int **edges, unsigned long long *update_time);
void remove_pair_wo_pair_info (int C, int N, int **bond, int **edges, int E, unsigned long long *update_time);
void remove_pair_wo_pair_info_v2 (int C, int N, int **bond, int **edges, int E, unsigned long long *update_time, int threshold);
int get_path_of_pair (int source, int target, int C, int N, int **bond, int *path, int *vec, int *tmp_vec, int *visited, int *reset, int **dat, int *nr_sp);
void sample_random_path (int C, int N, int **dag, int source, int target, int *path, int *nr_sp);
void randomize_vector (unsigned long *vector, unsigned long size);
void remove_and_resize(unsigned long **pairs, unsigned long q);
void tree_find_root (int n, int **root, int *res);
void modified_NZ_algorithm (int N, int **edges, double **largest_cluster);
void new_modified_NZ_algorithm (int N, int **edges, double **largest_cluster);
void adjust_pair(unsigned long **pairs, unsigned long q);
void efficient_pair_removal(int C, int N, int **bond, int **edges, int E, unsigned long long *update_time, int threshold);
void efficient_pair_removal_v2(int C, int N, int **bond, int **edges, int E, unsigned long long *update_time, int threshold);
void efficient_pair_removal_v3(int C, int N, int **bond, int **edges, int E, unsigned long long *update_time, int threshold, char *cost_name);
int find_event_index(double **largest_cluster, int num_edges);
int find_transition_index(double **largest_cluster, double criteria, int start);
unsigned long long find_longest_shortest_path_length(int N, int avg_k, unsigned long long *update_time);