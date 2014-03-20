#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "billiardslib.h"

int main(int argc, char **argv){
    if (argc != 6){
        printf("<N> <feltspeed> <wallrestore> <n-samles> <filename>\n");
        return 1;
    }
    ran_seed(19238);
    int N, samples;
    double eta = 0.0;
    double xi = 1.0;
    char filename[1024];

    N = atoi(argv[1]);
    eta = atof(argv[2]);
    xi = atof(argv[3]);
    samples = atoi(argv[4]);
    strncpy(filename, argv[5], 1024);

    double ttotal = 0.0;
    int *bounces = malloc(sizeof(int)*N*N/2);

    double rate = 0.0;
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    int tbounces, lastact;
    int imin = 65, imax = 1023;
    int jmin = 61, jmax = 532;
    int steps = 0;
    int userinformed = 0;
    #pragma omp parallel shared(userinformed)
    for (int i=0; i<N; i++){

        userinformed = 0;
        #pragma omp for nowait schedule(dynamic,N/16) reduction(+:steps)
        for (int j=0; j<N/2; j++){
            steps++;

            double xin[2] = {(double)((imax-imin)*i)/N + imin,
                             (double)((jmax-jmin)*j)/(N/2) + jmin};

            for (int kk=0; kk<samples; kk++){
                double vin[2] = {ran_ran2()*N-N/2, ran_ran2()*N-N/2};

                int rr = trackBall(xin, vin, eta, xi, &ttotal, &tbounces, &lastact);
                if (rr == RESULT_INPOCKET)
                    bounces[i+j*N] += 1;//tbounces;
            }
        }

        if (!userinformed)
        {
            userinformed = 1;
            clock_gettime(CLOCK_REALTIME, &end);
            rate = steps/((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9);
            printf("done: %0.4f \t rate: %0.2f\r", (float)steps/(N*N/2), rate);
            fflush(stdout);
        }
    }
    printf("done: %0.4f \t rate: %0.2f\n", (float)steps/(N*N), rate);

    FILE *f = fopen(filename, "wb");
    fwrite(bounces, sizeof(int), N*N/2, f);
    fclose(f);

    free(bounces);
    return 0;
}
