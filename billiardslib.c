#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "billiardslib.h"

#define NOBJS 24

int TOBJS[NOBJS] = {
    T_RAIL, T_RAIL, T_RAIL, T_RAIL, T_RAIL, T_RAIL,
    T_CORNER, T_CORNER, T_CORNER, T_CORNER, T_CORNER, T_CORNER,
    T_CORNER, T_CORNER, T_CORNER, T_CORNER, T_CORNER, T_CORNER,
    T_POCKET, T_POCKET, T_POCKET, T_POCKET, T_POCKET, T_POCKET,
    /*T_BALL*/
};


double OBJS[NOBJS][6] = {
    /* rail,    x0   y0    x1   y1   normal */
    {   64,  87,   64, 508,  1,  0}, 
    {   88, 533,  518, 533,  0, -1},
    {  573, 533, 1002, 533,  0, -1},
    { 1026, 508, 1026,  87, -1,  0},
    { 1003,  60,  518,  60,  0,  1},
    {  574,  60,   88,  60,  0,  1},

    /* corner     cx   cy   radius */
    {   88,  48, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {   52,  87, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {   52, 508, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {   88, 545, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {  518, 545, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {  574, 545, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    { 1003, 545, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    { 1038, 508, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    { 1038,  87, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    { 1003,  48, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {  518,  48, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},
    {  574,  48, RADIUS_CORNER, UNUSED, UNUSED, UNUSED},

    /* pocket     cx   cy   radius */
    {   45,  44, RADIUS_POCKET, UNUSED, UNUSED, UNUSED},
    {   45, 551, RADIUS_POCKET, UNUSED, UNUSED, UNUSED},
    {  545, 573, RADIUS_POCKET, UNUSED, UNUSED, UNUSED},
    { 1044, 551, RADIUS_POCKET, UNUSED, UNUSED, UNUSED},
    { 1044,  45, RADIUS_POCKET, UNUSED, UNUSED, UNUSED},
    {  545,  21, RADIUS_POCKET, UNUSED, UNUSED, UNUSED},

    /*{  -100,  276,  180, UNUSED, UNUSED, UNUSED},*/
};

//============================================================================
// These are helper functions for the basics
//============================================================================
inline double mymod(double a, double b){
      return a - b*(int)(a/b) + b*(a<0);
}

inline double dot(double *r1, double *r2){
    return r1[0]*r2[0] + r1[1]*r2[1];
}

inline double cross(double *a, double *b){
    return a[0]*b[1]-a[1]*b[0];
}

inline void position(double *x0, double *v0, double t, double *out){
    out[0] = x0[0] + v0[0]*t;
    out[1] = x0[1] + v0[1]*t;
}

inline void velocity(double *v0, double t, double eta, double *out){
    out[0] = v0[0] - eta*t;
    out[1] = v0[1] - eta*t;
}

//============================================================================
// a little bit of geometry
//============================================================================
void normal(double *pos, int obj, double *out){
    if (TOBJS[obj] == T_RAIL){
        out[0] = OBJS[obj][4];
        out[1] = OBJS[obj][5];
    } else {
        out[0] = OBJS[obj][0] - pos[0];
        out[1] = OBJS[obj][1] - pos[1];

        double invlen = 1./sqrt(dot(out, out));
        out[0] *= invlen;
        out[1] *= invlen;
    } 
}

void reflect(double *r, double *normal, double *out, double restore){
    double dt = dot(r, normal);
    out[0] = 2*dt*normal[0] - r[0];
    out[1] = 2*dt*normal[1] - r[1];
    out[0] *= -restore; out[1] *= -restore;
}

double collide_rail(double *p, double *v, double *r0, double *r1){
    double vr[2] = {r1[0] - r0[0], r1[1] - r0[1]};
    double vcross = cross(vr, v);

    if (fabs(vcross) < XTOL)
        return -1;

    // based on the notations above, t0 is the time for 
    // the rail velocity and t1 is the time in terms of the ball 
    double t0 = (vr[0]*(r0[1] - p[1]) + vr[1]*(p[0] - r0[0])) / vcross;
    double t1 = ( v[1]*(p[0] - r0[0]) +  v[0]*(r0[1] - p[1])) / vcross;

    if (t1 > 0 && t1 < 1)
        return t0;
    return -1;
}

double collide_circle(double *p, double *v, double *cr, double rad){
    double d[2] = {p[0] - cr[0], p[1] - cr[1]};
    double a = dot(v, v);
    double b = 2*dot(v, d);
    double c = dot(d, d) - rad*rad;
    double desc = b*b - 4*a*c;

    if (desc > 0){
        double t0 = (-b + sqrt(desc)) / (2*a);
        double t1 = (-b - sqrt(desc)) / (2*a);

        if (t0 > 0 && t0 < t1) return t0;
        if (t1 > 0 && t1 < t0) return t1;
        return t0;
    }
    return -1;
}

//============================================================================
// finds the collision time for an initial condition
// by finding the roots of a poly and finding the nearest collision
//============================================================================
int docollision(double *pos, double *vel, double eta, double xi, int last,
        double *tcoll, int *next){
    /* 
     * This functions determines whether a particular trajectory collides
     * with different parts of the table and returns:
     *  0 : There was no collision (how odd!)
     *  1 : There was a collision and it occured at tcoll
     *  2 : The ball landed in a pocket
     * It does not modify the values of pos, vel; modified tcoll
    */
    int outcome = RESULT_NOTHING;
    double ttmp = 1e10;
    double tpos[2], tvel[2];
    memcpy(tpos, pos, sizeof(double)*2);
    memcpy(tvel, vel, sizeof(double)*2);

    *tcoll = NAN;
    // loop over everything that it could possibly interact with
    // and keep track of the nearest interaction
    for (int i=0; i<NOBJS; i++){
        if (i == last) continue;

        double *o = (double*)&OBJS[i];
        int type = TOBJS[i];

        if (type == T_RAIL)
            ttmp = collide_rail(tpos, tvel, &o[0], &o[2]);
        //if (type == T_CORNER)
        //    ttmp = collide_circle(tpos, tvel, &o[0], o[2]);
        if (type == T_POCKET)
            ttmp = collide_circle(tpos, tvel, &o[0], o[2]);
        if (type == T_BALL)
            ttmp = collide_circle(tpos, tvel, &o[0], o[2]);

        if (ttmp > 0 && (isnan(*tcoll) || ttmp < *tcoll)){
            *tcoll = ttmp;// + eta + xi;
            *next = i;
            if (type == T_RAIL || type == T_CORNER)
                outcome = RESULT_COLLISION;
            else
                outcome = RESULT_INPOCKET;
        }
    }
    return outcome;
}

int trackBall(double *pos, double *vel, double eta, double xi, double *t){
    int result, actlast, actnext;
    double tcoll, vlen, ttotal;

    int tbounces = 0;
    double tpos[2], tvel[2], norm[2];
    memcpy(tpos, pos, sizeof(double)*2);
    memcpy(tvel, vel, sizeof(double)*2);

    ttotal = 0;
    actlast = -1;
    actnext = -1;
    while (tbounces < MAX_COLLISIONS){
        // get the next collision
        actlast = actnext;
        result = docollision(tpos, tvel, eta, xi, actlast, &tcoll, &actnext);

        if (result == RESULT_NOTHING)
            return 0;
        if (result == RESULT_INPOCKET)
            break;

        ttotal += tcoll;
        // figure out where it hit and what speed
        position(tpos, tvel, tcoll, tpos);
        velocity(tvel, tcoll, eta, tvel);

        vlen = dot(tvel, tvel);
        if (vlen < XTOL) break;

        normal(tpos, actnext, norm);
        reflect(tvel, norm, tvel, xi);
        tbounces++;
    }

    *t = ttotal;
    return tbounces;
}

