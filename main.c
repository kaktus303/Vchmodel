#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define p1 -11
#define p2 9
#define u10 -3
#define u20 5
#define start 0.0
#define finish 10.0
#define n 1000
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
void Eiler(double *x_points, double *u1_points, double *u2_points)
{
    double step = (finish - start)/(n-1);
    double x = 0, u1 = u10, u2 = u20;
    u1_points[0] = u10;
    u2_points[0] = u20;
    for(int i = 1; i < n; i++)
    {
        u1_points[i] = u1_points[i-1] + step * f1_test(u1_points[i-1], u2_points[i-1]);
        // printf("%lf\n",f1_test(u1_points[i-1], u2_points[i-1]));
        u2_points[i] = u2_points[i-1] + step * f2_test(u1_points[i-1], u2_points[i-1]);
    }

}
void x_filling(double *x_points)
{
    double step = (finish - start)/(n-1), x = start;
    x_points[0] = start;
    printf("%.15lf", step);
    for(int i = 1;i < n;++i)
    {
        x+=step;
        x_points[i] = x;
    }
}
void file_write(FILE *file_graph, FILE *file_error, double *points_x, double *points_y)
{
    for (int i = 0; i < n; ++i)
    {
        fprintf(file_graph, "%.15lf %.15lf\n", points_x[i], points_y[i]);
        // fprintf(file_error, "%.15lf %.15lf\n", points_x[i], fabs(points_y[i] - test_function_fi(points_x[i])));
    }
}
double runhe_kutte_4(double x, double u1, double u2)
{

}
void answer(double *points_answer, double  *points_x)
{
    for(int i = 0; i < n; ++i)
    {
        points_answer[i] = 0.5 * (u10 + u20) * exp(lamda1 * points_x[i])+  0.5 * (u10 - u20) * exp(lamda2 * points_x[i]);
    }
}

int main()
{
    double *x_points = malloc(sizeof(double) * n),*u1_points = malloc(sizeof(double) * n), *u2_points= malloc(sizeof(double) * n), *answer_points = malloc(sizeof(double) * n);  
    FILE *f = fopen("test.txt", "w"), *f1 = fopen("answer.txt", "w");
    x_filling(x_points);
    Eiler(x_points,u1_points,u2_points);
    for(int i = 0;i<n;++i)
    {
       //  printf("%lf", x_points[i]);
    }
    answer(answer_points, x_points);
    // x_points[0] = 0;
    file_write(f,f,x_points,u1_points);
    file_write(f1,f1,x_points, answer_points);
    
}