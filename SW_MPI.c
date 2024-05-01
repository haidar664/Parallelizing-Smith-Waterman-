#include <stdio.h>
#include <mpi.h>
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
    int matrix[l2 + 1][l2 + 1];
    int max_score = 0;
    for (int i = 0; i < l2 + 1; i++) {
        matrix[i][0] = 0;
    }
    for (int i = 0; i < l2 + 1; i++) {
        matrix[0][i] = 0;
    }
    for (int i = 1; i < l2 + 1; i++) {
        for (int j = 1; j < l2 + 1; j++) {
            if (q[i - 1] == t[(l1-l2)+(j-1)]) 
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

    int n, l2,longestSimilarSeq=0;
    int nbSubsets,package[2];
    start_time=clock();

    MPI_Init(NULL,NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    if(world_rank==0){
        printf("Enter the length query sequence\n");
        scanf("%d",&l2); // length of sequence q
        char q[l2];
        printf("Enter query sequence\n");
        for(int i=0;i<l2;i++){
            scanf(" %c",&q[i]);
        }
        printf("Enter number of target sequences\n");
        scanf("%d",&n);
        nbSubsets=n/(world_size-1);
        char t[l2*nbSubsets];
        int maxScores[n],maxScoresFromSlave[1+nbSubsets];
        package[0]=l2;
        package[1]=nbSubsets;
        for(int slave=1;slave<world_size;slave++){
            printf("%d- Enter %d target sequence(s)\n",slave,nbSubsets);
            for(int i=0;i<l2*nbSubsets;i++){
                scanf(" %c",&t[i]);
            }
            MPI_Send(&package,2*sizeof(int),MPI_INT, slave, 0,MPI_COMM_WORLD);            
            printf("package is sented ");
            MPI_Send(&q,l2,MPI_CHAR, slave, 0,MPI_COMM_WORLD);
            printf("q is sented ");
            MPI_Send(&t,l2*nbSubsets,MPI_CHAR, slave, 0,MPI_COMM_WORLD);
            printf("t is sented\n");
        }
        int max_score = 0;
        for(int slave=1;slave<world_size/*+1*/;slave++){
            MPI_Recv(&maxScoresFromSlave,1+nbSubsets*sizeof(int),MPI_INT, MPI_ANY_SOURCE, 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   
            
            int index=maxScoresFromSlave[0];
            printf("index= %d\n", index );
            for(int i=0;i<nbSubsets;i++){
                maxScores[i+index]=maxScoresFromSlave[i+1];
                printf("maxScore of %d is: %d\n",i+index, maxScores[i+index]);
                if(max_score<maxScores[i+index]){
                    max_score=maxScores[i+index];
                }
            }
            // offset += nbSubsets;
            printf("I am the master, it is %d\n", slave );
        }
        
        printf("max score= %d\n", max_score );
        
        end_time=clock();
        time=((double) (end_time-start_time))/CLOCKS_PER_SEC;
        printf("Max alignment score: %d\n", max_score);
        printf("The average execution time is: %f   s\n", time);
    }

    else{
        printf("I am slave %d\n", world_rank);
        
        MPI_Recv(&package,2*sizeof(int),MPI_INT, 0, 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        l2=package[0];nbSubsets=package[1];
        char q[l2];
        char t[l2*nbSubsets];
        int maxScores[nbSubsets+1];
        maxScores[0]=world_rank;
        printf("I am slave %d\n", maxScores[0]);
        MPI_Recv(&q,l2,MPI_CHAR, 0, 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&t,l2*nbSubsets,MPI_CHAR, 0, 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int i=0;i<nbSubsets;i++){
            maxScores[i+1]=smithWaterman(t,q,l2*nbSubsets,l2);
            printf("similiarity is calculated - %d",world_rank);
        }
        MPI_Send(&maxScores,nbSubsets*sizeof(int)+1,MPI_INT, 0, 0,MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
