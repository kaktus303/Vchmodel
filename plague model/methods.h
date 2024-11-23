#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
struct vector {
    double x;
    double y;
    double s;
    double t;
}typedef vector;
void runhe_kutte_4(double step, double *points_x, double *points_y, int n);
double runhe_kutte_4_1(double step, double x, double y);
void ABM(double step, double *points_x, double *points_y, int n);