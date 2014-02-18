#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cupgamelib.h"

int main(int argc, char **argv){
    if (argc != 6){
        printf("<restoring> <x> <y> <z> <filename>\n");
        return 1;
    }

    //double h = 0.95; double r = 0.3; restore = 0.75; parameters for a real game
    double h = 0.7;
    double r = 0.3;
    double x, y, z;
    double restore;
    char filename[1024];

    restore = atof(argv[1]);
    x = atof(argv[2]);
    y = atof(argv[3]);
    z = atof(argv[4]);
    strncpy(filename, argv[5], 1024);

    double xin[3] = {x, y, z};
    double vin[3] = {0.0, 0.0, -1e-1};

    const int MAXLEN = 10000;
    double *traj = (double*)malloc(sizeof(double)*MAXLEN);
    int tracklen = trackTrajectory(3, xin, 3, vin, h, r, restore, 
                                   MAXLEN/TSAMPLES, MAXLEN, traj);

    FILE *f = fopen(filename, "wb");
    fwrite(traj, sizeof(double), tracklen, f);
    fclose(f);

    return 0;
}
