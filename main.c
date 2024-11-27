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
void Eiler(int n, double *x_points, double *u1_points, double *u2_points)
{
    double step = (finish - start)/(n-1);
    double x = 0, u1 = u10, u2 = u20;
    u1_points[0] = u10;
    u2_points[0] = u20;
    for(int i = 1; i < n; i++)
    {
        u1_points[i] = u1_points[i-1] + step * f1_test(u1_points[i-1], u2_points[i-1]);
        u2_points[i] = u2_points[i-1] + step * f2_test(u1_points[i-1], u2_points[i-1]);
    }

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
void file_write(FILE *file_graph, FILE *file_error, double *points_x, double *points_y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        fprintf(file_graph, "%.15lf %.15lf\n", points_x[i], points_y[i]);
        fprintf(file_error, "%.15lf %.15lf\n", points_x[i], fabs(points_y[i] - test_function_fi(points_x[i])));
    }
}
double runhe_kutte_4(double x, double u1, double u2)
{

}
void answer(double *points_answer, double points_x)
{
    for(int i = 0; i < n; ++i)
    {

    }
}

int main()
{
    double *x_points = malloc(sizeof(double) * n),*u1_points = malloc(sizeof(double) * n), *u2_points= malloc(sizeof(double) * n);  
    FILE *f = fopen("test.txt", w);
    
}