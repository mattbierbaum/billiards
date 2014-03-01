#ifndef __BILLIARDS_H__
#define __BILLIARDS_H__

//========================================================
// global constants for the calculation
//========================================================
#define XTOL 1e-15

#define MAX_COLLISIONS 10000
#define TSAMPLES 25

#define RESULT_NOTHING   0
#define RESULT_COLLISION 1
#define RESULT_INPOCKET  2

#define RADIUS_POCKET 60
#define RADIUS_CORNER 25

#define UNUSED  -1
#define T_RAIL   0
#define T_CORNER 1
#define T_POCKET 2
#define T_BALL   3

#define M_PI 3.14159265358979323846

#define SIGN(x) ((x) <= 0 ? -1 : 1)

//========================================================
/* These are functions that should be called externally */
int trackBall(double *p, double *v, double eta, double xi, double *t);
int trackTrajectory(double *p, double *v, double eta, double xi, double *traj,
        int *tlen, int maxlen);

//========================================================
/* internal use functions only */
double mymod(double a, double b);
double dot(double *r1, double *r2);
double cross(double *a, double *b);
void position(double *x0, double *v0, double t, double *out);
void velocity(double *v0, double t, double eta, double *out);

//========================================================
/* actual dynamics functions */
void normal(double *pos, int obj, double *out);
void reflect(double *r, double *normal, double *out, double restore);
double collide_rail(double *p, double *v, double *r0, double *r1);
double collide_circle(double *p, double *v, double *c, double rad);
int docollision(double *p, double *v, int last, double *tcoll, int *next);
#endif
