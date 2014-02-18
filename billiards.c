#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "billiardslib.h"

#define XSCAN_UNITQUAD {2.0*i/N-1, 2.0*j/N-1, height}
#define XSCAN_UNITLINE {1.0*(i+j*N)/(N*N), 1.0*(i+j*N)/(N*N), height}
#define XSCAN_HEIGHT {1.0*i/N, 1.0*i/N, height+ 0.5*(double)j/N}
#define XSCAN_HEXAGON {6.0*i/N-3.0, 6.0*j/N-2.0, 0.75 + 5*(double)j/N}

int main(int argc, char **argv){
    if (argc != 5){
        printf("<N> <feltspeed> <wallrestore> <filename>\n");
        return 1;
    }

    int N;
    double eta = 0.0;
    double xi = 1.0;
    char filename[1024];

    N = atoi(argv[1]);
    eta = atof(argv[2]);
    xi = atof(argv[3]);
    strncpy(filename, argv[4], 1024);

    double ttotal = 0.0;
    int *bounces = malloc(sizeof(int)*N*N);

    double rate = 0.0;
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    int ci = 545;
    int cj = 276;
    int steps = 0;
    int userinformed = 0;
    #pragma omp parallel shared(userinformed)
    for (int i=0; i<N; i++){

        userinformed = 0;
        #pragma omp for nowait schedule(dynamic,N/16) reduction(+:steps)
        for (int j=0; j<N; j++){
            steps++;

            double xin[2] = {ci, cj};
            double vin[2] = {i-N/2, j-N/2};

            int b = trackBall(xin, vin, eta, xi, &ttotal);
            bounces[i+j*N] = b;
        }

        if (!userinformed)
        {
            userinformed = 1;
            clock_gettime(CLOCK_REALTIME, &end);
            rate = steps/((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9);
            printf("done: %0.4f \t rate: %0.2f\r", (float)steps/(N*N), rate);
            fflush(stdout);
        }
    }
    printf("done: %0.4f \t rate: %0.2f\n", (float)steps/(N*N), rate);

    FILE *f = fopen(filename, "wb");
    fwrite(bounces, sizeof(int), N*N, f);
    fclose(f);

    free(bounces);
    return 0;
}
