#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h> // Include stdlib.h for malloc and free

#define MATCH_SCORE 3
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2

int max(int a, int b, int c) {
    int max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    return max_val;
}

int smithWaterman(char t[], char q[], int l1, int l2) {
    int max_score = 0;
    int* matrix = (int*)malloc((l1+1)*(l2+1) * sizeof(int)); // allocate memory for matrix

    // Initialize first row and first column
    for (int i = 0; i <= l2; i++) matrix[i*(l1+1)] = 0;
    for (int i = 0; i <= l1; i++) matrix[i] = 0;

    #pragma omp parallel for reduction(max:max_score)
    for (int i = 1; i <= l2; i++) {
        for (int j = 1; j <= l1; j++) {
            if (q[i - 1] == t[j - 1]) 
                matrix[i*(l1+1) + j] = max(matrix[(i - 1)*(l1+1) + (j - 1)] + MATCH_SCORE, matrix[(i - 1)*(l1+1) + j] + GAP_PENALTY, matrix[i*(l1+1) + (j - 1)] + GAP_PENALTY);
            else
                matrix[i*(l1+1) + j] = max(matrix[(i - 1)*(l1+1) + (j - 1)] + MISMATCH_SCORE, matrix[(i - 1)*(l1+1) + j] + GAP_PENALTY, matrix[i*(l1+1) + (j - 1)] + GAP_PENALTY);
            if (matrix[i*(l1+1) + j] > max_score)
                max_score = matrix[i*(l1+1) + j];
        }
    }

    free(matrix); // free allocated memory
    return max_score;
}

int main(int argc, char const *argv[]) {
    double time;
    clock_t start_time, end_time;

    int n, l2, max_score=0,longestSimilarSeq=0;
    
    start_time=clock();

    printf("Enter the length query sequence\n");
    scanf("%d",&l2); // length of sequence q
    char q[l2];
    printf("Enter query sequence\n");
    for(int i=0;i<l2;i++){
        scanf(" %c",&q[i]);
    }
    printf("Enter number of target sequences\n");
    scanf("%d",&n);
    char t[n][l2];
    int maxScores[n], maxScoresFromSlave[n];

    printf("Enter %d target sequence(s)\n", n);
    for(int i=0;i<n;i++){
        printf("Sequence %d: ", i+1);
        for(int j=0;j<l2;j++){
            scanf(" %c",&t[i][j]);
        }
    }

    #pragma omp parallel for shared(maxScoresFromSlave) num_threads(n)
    for(int i = 0; i < n; i++) {
        maxScoresFromSlave[i] = smithWaterman(t[i], q, l2, l2);
    }

    for(int i = 0; i < n; i++) {
        if(maxScoresFromSlave[i] > max_score) {
            max_score = maxScoresFromSlave[i];
            longestSimilarSeq = i;
        }
    }

    end_time = clock();
    time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Max alignment score: %d\n", max_score);
    printf("The average execution time is: %f s\n", time);

    return 0;
}
