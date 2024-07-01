int find_node (int i, int j, int **bond);
void generate_ERgraph (int N, int **bond, int E);
void generate_2DSQL_w_PB(int L, int **bond);
void copy_network (int N, int **origin_bond, int **copy_bond);
void copy_list_of_pair (unsigned long num_pairs, unsigned long **origin_pairs, unsigned long **copy_pairs);
unsigned long long geometric_distribution(double P);
int warning_message(int find_pairs, int C);
void clean_time_stamp(unsigned long long *time_stamp, int E);
int get_sampling_factor(int N, int criteria);
