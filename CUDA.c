%%cuda
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATCH_SCORE 3
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2

__device__ int max_score_device = 0;

__device__ int max(int a, int b, int c) {
    int max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    return max_val;
}

__global__ void smithWatermanKernel(char *t, char *q, int l1, int l2, int *matrix) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (i <= l2 && j <= l1) {
        if (q[i - 1] == t[j - 1]) 
            matrix[i*(l1+1) + j] = max(matrix[(i - 1)*(l1+1) + (j - 1)] + MATCH_SCORE, matrix[(i - 1)*(l1+1) + j] + GAP_PENALTY, matrix[i*(l1+1) + (j - 1)] + GAP_PENALTY);
        else
            matrix[i*(l1+1) + j] = max(matrix[(i - 1)*(l1+1) + (j - 1)] + MISMATCH_SCORE, matrix[(i - 1)*(l1+1) + j] + GAP_PENALTY, matrix[i*(l1+1) + (j - 1)] + GAP_PENALTY);
        atomicMax(&max_score_device, matrix[i*(l1+1) + j]);
    }
}

int smithWaterman(char t[], char q[], int l1, int l2) {
    char *d_t, *d_q;
    int *d_matrix;
    int max_score = 0;

    cudaMalloc((void**)&d_t, l1 * sizeof(char));
    cudaMalloc((void**)&d_q, l2 * sizeof(char));
    cudaMalloc((void**)&d_matrix, (l1+1)*(l2+1) * sizeof(int));

    cudaMemcpy(d_t, t, l1 * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_q, q, l2 * sizeof(char), cudaMemcpyHostToDevice);

    dim3 blockDim(16, 16);
    dim3 gridDim((l2 + blockDim.x - 1) / blockDim.x, (l1 + blockDim.y - 1) / blockDim.y);

    smithWatermanKernel<<<gridDim, blockDim>>>(d_t, d_q, l1, l2, d_matrix);

    cudaMemcpyFromSymbol(&max_score, max_score_device, sizeof(int), 0, cudaMemcpyDeviceToHost);

    cudaFree(d_t);
    cudaFree(d_q);
    cudaFree(d_matrix);

    return max_score;
}

int main() {
    double time;
    clock_t start_time, end_time;

    int l2, max_score = 0, longestSimilarSeq = 0;

    start_time = clock();

    // Open input file
    FILE *file = fopen("sequences.txt", "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    // Read query sequence
    char q[1000];
    fscanf(file, "%*s %s", q);
    l2 = strlen(q);

    // Read number of target sequences
    int n;
    fscanf(file, "%*s %d", &n);

    // Read target sequences
    char t[n][1000];
    for (int i = 0; i < n; i++) {
        fscanf(file, "%*s %s", t[i]);
    }

    fclose(file);

    // Compute Smith-Waterman score for each target sequence
    for (int i = 0; i < n; i++) {
        int score = smithWaterman(t[i], q, strlen(t[i]), l2);
        if (score > max_score) {
            max_score = score;
            longestSimilarSeq = i;
        }
    }

    end_time = clock();
    time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Max alignment score: %d\n", max_score);
    printf("The average execution time is: %f s\n", time);

    return 0;
}