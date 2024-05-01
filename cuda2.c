%%cuda
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATCH_SCORE 3
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2

__device__ int max(int a, int b, int c) {
    int max_val = 0;
    if (a > max_val) max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    return max_val;
}

__global__ void smithWatermanKernel(char *t, char *q, int l1, int l2, int *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < l1 * l2) {
        int i = idx / l1;
        int j = idx % l1;
        int score = 0;
        if (q[i] == t[j]) 
            score = max(result[(i - 1) * l1 + j - 1] + MATCH_SCORE,
                        result[(i - 1) * l1 + j] + GAP_PENALTY,
                        result[i * l1 + j - 1] + GAP_PENALTY);
        else
            score = max(result[(i - 1) * l1 + j - 1] + MISMATCH_SCORE,
                        result[(i - 1) * l1 + j] + GAP_PENALTY,
                        result[i * l1 + j - 1] + GAP_PENALTY);
        result[i * l1 + j] = score;
    }
}

int main(int argc, char const *argv[]) {
    double time;
    clock_t start_time, end_time;
    int l2 = 1000; // Length of sequence
    start_time = clock();
    char *q = (char *)malloc(l2 * sizeof(char));
    FILE *fp;
    fp = fopen("sequences.txt", "r");
    if (fp == NULL) {
        printf("Error opening file.\n");
        return 1;
    }
    fgets(q, l2 + 1, fp); // Read one sequence
    fclose(fp);

    int n = 3; // Number of target sequences
    int *maxScores = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        char *t = q; // Use the same sequence for all targets
        int *d_result;
        cudaMalloc(&d_result, l2 * l2 * sizeof(int));
        cudaMemcpy(&d_result, maxScores, l2 * l2 * sizeof(int), cudaMemcpyHostToDevice);
        smithWatermanKernel<<<(l2 * l2 + 255) / 256, 256>>>(t, q, l2, l2, d_result);
        cudaMemcpy(&maxScores, d_result, l2 * l2 * sizeof(int), cudaMemcpyDeviceToHost);

        int max_score = 0;
        for (int k = 0; k < l2 * l2; k++) {
            if (maxScores[k] > max_score) {
                max_score = maxScores[k];
            }
        }
        maxScores[i] = max_score;
    }

    for (int i = 0; i < n; i++) {
        printf("the max alignment at %d is: %d\n", i, maxScores[i]);
    }

    end_time = clock();
    time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Max alignment score: %d\n", maxScores[0]); // All scores should be the same
    printf("The average execution time is: %f s\n", time);

    free(q);
    free(maxScores);
    return 0;
}

