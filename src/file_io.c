#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "file_io.h"

int get_max_num_from_filenames(const char* dir_path, char* f_ext) {
    DIR *d;
    struct dirent *dir;
    int max_num = -1;

    d = opendir(dir_path);
    if (d) {
        while ((dir = readdir(d)) != NULL) {
            char *extension = strstr(dir->d_name, f_ext);
            if (extension != NULL) {
                char *num_str = extension;
                while (num_str > dir->d_name && *(num_str - 1) != '_') {
                    num_str--;  // Go back until '_' or start of string
                }

                *extension = '\0';  // Null-terminate before extension
                int num = atoi(num_str);
                if (num > max_num) {
                    max_num = num;
                }
            }
        }
        closedir(d);
    }

    return max_num;
}

int read_network_file(char *fn, int **bond) {
	FILE *fp;
	int i, j, N;
	fp = fopen(fn, "r");

	if (fp == NULL) {
        printf("!!! FAIL TO OPEN YOUR FILE !!! \n");
        printf("!!! CHECK YOUR FILE AGAIN  !!! \n");
        return -1;
    } else {
        do{
            fscanf(fp, "%d;%d\n", &i, &j);
            bond[0][i]++;
            bond[i] = (int *)realloc(bond[i], sizeof(int) * (bond[0][i] + 1));
            bond[i][bond[0][i]] = j;

            bond[0][j]++;
            bond[j] = (int *)realloc(bond[j], sizeof(int) * (bond[0][j] + 1));
            bond[j][bond[0][j]] = i;


        }while(!feof(fp));
        fclose(fp);
    }
	return 0;
}