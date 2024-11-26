#include<stdio.h>
#include<stdlid.h>
#define p1 -11
#define p2 9
#define u10 5
#define u20 -3
#define start 0
#define finish 10
#define n 10000
#define lamda1 -20
#define lamda2 -2
double f1_test(double u1, double u2)
{
    return p1 * u1 + p2 * u2;
}
double f2_test(double u1, double u2)
{
    return p2 * u1 + p1 * u2;
}
double Eiler(int n, double *x_points, double *u1_points, double *u2_points)
{
    double step = (finish - start)/(n-1);
    double x = 0, u1 = u10, u2 = u20;
    u1_points[0] = u10;
    u2_points[0] = u20;
    for(int i = 1; i < n; i++)

}
void x_filling(double *x_points)
{
    double step = (finish - start)/(n-1), x = start, x_points[0] = start;
    for(int i = 1;i < n;++i)
    {
        x+=step;
        x_points[i] = x;
    }
}
double runhe_kutte_4(double x, double u1, double u2)
{

}

int main()
{
    double *x_points,*u1_points, *u2_points;   
}