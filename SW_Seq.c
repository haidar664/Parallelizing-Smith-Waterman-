#include <stdio.h>
#include <time.h>

#define MATCH_SCORE 3
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2

int max(int a, int b, int c) {
    int max_val = 0; 
    if (a > max_val) max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    return max_val;
}

int smithWaterman(char t[], char q[], int l1, int l2) {
    int matrix[l2 + 1][l1 + 1];
    int max_score = 0;
    for (int i = 0; i < l2 + 1; i++) {
        matrix[i][0] = 0;
    }
    for (int i = 0; i < l1 + 1; i++) {
        matrix[0][i] = 0;
    }
    for (int i = 1; i < l2 + 1; i++) {
        for (int j = 1; j < l1 + 1; j++) {
            if (q[i - 1] == t[j - 1]) 
                matrix[i][j] = max(matrix[i - 1][j - 1] + MATCH_SCORE, matrix[i - 1][j] + GAP_PENALTY, matrix[i][j - 1] + GAP_PENALTY);
            else
                matrix[i][j] = max(matrix[i - 1][j - 1] + MISMATCH_SCORE, matrix[i - 1][j] + GAP_PENALTY, matrix[i][j - 1] + GAP_PENALTY);
            if (matrix[i][j] > max_score)
                max_score = matrix[i][j];
        }
    }
    // for (int i = 0; i < l2 + 1; i++) {
    //     for (int j = 0; j < l1 + 1; j++) {
    //         printf("%d ", matrix[i][j]);
    //     }
    //     printf("\n");
    // }
    return max_score;
}

int main(int argc, char const *argv[]) {
    double time;
    clock_t start_time, end_time;
    int n,l2,max_score;
    start_time=clock();
    printf("Enter the length query sequence\n");
    scanf("%d",&l2); // length of sequence q
    char q[l2];
    printf("Enter query sequence\n");
    for(int i=0;i<l2;i++){
        scanf(" %c",&q[i]);
    }
    while ((getchar()) != '\n');
    printf("Enter number of target sequences\n");
    scanf("%d",&n);
    int maxScores[n];
    for(int i=0;i<n;i++){
        // printf("Enter the length target sequence\n");
        // scanf("%d",&l1); // length of sequence t
        char t[l2];
        printf("Enter the %d target sequence\n",i);
        for(int i=0;i<l2;i++){
            scanf(" %c",&t[i]);
        }
        while ((getchar()) != '\n');
        max_score = smithWaterman(t, q, l2, l2);
        maxScores[i]=max_score;
    }
    for(int i=0;i<n;i++){
    	printf("the max alignment at %d is: %d\n",i,maxScores[i]);
        if(maxScores[i]>max_score) max_score=maxScores[i];
    }
    end_time=clock();
    time=((double) (end_time-start_time))/CLOCKS_PER_SEC;
    printf("Max alignment score: %d\n", max_score);
    printf("The average execution time is: %f s\n", time);
    return 0;
}
