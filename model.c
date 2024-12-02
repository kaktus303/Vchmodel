#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define tau0 0
#define x0 0
#define y0 0
#define s0 10
#define Alfas 2.7
#define Alfay 0.8
#define Betax 0.9
#define Betay 1.2
#define Betas 0.9
#define Alfaxy 1.5
#define n 1000
#define start 0.0
#define finish 10.0
#define number 4
double ksi(double t)
{
    return 1.0;
}

double f1(double t, double x, double y, double s, double tau)
{
    return Betax * (pow(s, number) / (1.0 + pow(s, number))) * ksi(t) - Alfaxy * x * y;
}
double f2(double t, double x, double y, double s, double tau)
{
    return Betay * x * (t - tau) * ksi(t - tau) - Alfay * y;
}
double f3(double t, double x, double y, double s, double tau)
{
    return Betas - Alfas * y * s;
}
void t_filling(double *t_points)
{
    double step = (finish - start) / (n - 1), x = start;
    t_points[0] = start;
    printf("%.15lf\n", step);
    for (int i = 1; i < n; ++i)
    {
        x += step;
        t_points[i] = x;
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
double runhe_kutte_4(double *t_points, double *x_points, double *y_points, double *s_points, int stoper)
{

    double step = (finish - start) / (n - 1);
    double k1x, k1y, k1s, k2x, k2y, k2s, k3x, k3y, k3s, k4x, k4y, k4s;
    x_points[0] = x0;
    y_points[0] = y0;
    s_points[0] = s0;
    for (int i = 1; i < n && i < stoper; ++i)
    {
        k1x = step * f1(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0);
        k1y = step * f2(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0);
        k1s = step * f3(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0);

        k2x = step * f1(t_points[i - 1] + (step / 2.0), x_points[i - 1] + k1x / 2, y_points[i - 1] + k1y / 2, s_points[i - 1] + k1s / 2, tau0);
        k2y = step * f2(t_points[i - 1] + (step / 2.0), x_points[i - 1] + k1x / 2, y_points[i - 1] + k1y / 2, s_points[i - 1] + k1s / 2, tau0);
        k2s = step * f3(t_points[i - 1] + (step / 2.0), x_points[i - 1] + k1x / 2, y_points[i - 1] + k1y / 2, s_points[i - 1] + k1s / 2, tau0);

        k3x = step * f1(t_points[i - 1] + (step / 2.0), x_points[i - 1] + k2x / 2, y_points[i - 1] + k2y / 2, s_points[i - 1] + k2s / 2, tau0);
        k3y = step * f2(t_points[i - 1] + (step / 2.0), x_points[i - 1] + k2x / 2, y_points[i - 1] + k2y / 2, s_points[i - 1] + k2s / 2, tau0);
        k3s = step * f3(t_points[i - 1] + (step / 2.0), x_points[i - 1] + k2x / 2, y_points[i - 1] + k2y / 2, s_points[i - 1] + k2s / 2, tau0);

        k4x = step * f1(t_points[i - 1] + (step), x_points[i - 1] + k3x, y_points[i - 1] + k3y, s_points[i - 1] + k3s, tau0);
        k4y = step * f2(t_points[i - 1] + (step), x_points[i - 1] + k3x, y_points[i - 1] + k3y, s_points[i - 1] + k3s, tau0);
        k4s = step * f3(t_points[i - 1] + (step), x_points[i - 1] + k3x, y_points[i - 1] + k3y, s_points[i - 1] + k3s, tau0);

        x_points[i] = x_points[i - 1] + (1.0 / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        y_points[i] = y_points[i - 1] + (1.0 / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
        s_points[i] = s_points[i - 1] + (1.0 / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
    }
}
double ABM(double *t_points, double *x_points, double *y_points, double *s_points, int stoper)
{
    double step = (finish - start) / (n - 1);
    x_points[0] = x0;
    y_points[0] = y0;
    s_points[0] = s0;
    double x_pre = 0, y_pre = 0, s_pre= 0;
    int addams_counter = 0;
    runhe_kutte_4(t_points, x_points, y_points, s_points, 4);
    for (int i = 4; i < n; i++)
    {
        x_pre = x_points[i - 1] + (step / 24.0) * (55 * f1(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0) - 59 * f1(t_points[i - 2], x_points[i - 2], y_points[i - 2], s_points[i - 2], tau0) + 37 * f1(t_points[i - 3], x_points[i - 3], y_points[i - 3], s_points[i - 3], tau0) - 9 * f1(t_points[i - 4], x_points[i - 4], y_points[i - 4], s_points[i - 4], tau0));
        y_pre = y_points[i - 1] + (step / 24.0) * (55 * f2(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0) - 59 * f2(t_points[i - 2], x_points[i - 2], y_points[i - 2], s_points[i - 2], tau0) + 37 * f2(t_points[i - 3], x_points[i - 3], y_points[i - 3], s_points[i - 3], tau0) - 9 * f2(t_points[i - 4], x_points[i - 4], y_points[i - 4], s_points[i - 4], tau0));
        s_pre = s_points[i - 1] + (step / 24.0) * (55 * f3(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0) - 59 * f3(t_points[i - 2], x_points[i - 2], y_points[i - 2], s_points[i - 2], tau0) + 37 * f3(t_points[i - 3], x_points[i - 3], y_points[i - 3], s_points[i - 3], tau0) - 9 * f3(t_points[i - 4], x_points[i - 4], y_points[i - 4], s_points[i - 4], tau0));


        // x = x + step;
        x_points[i] = x_points[i - 1] + (step / 24.0) * (9 * f1(t_points[i], x_pre, y_pre, s_pre, tau0) + 19 * f1(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0) - 5 * f1(t_points[i - 2], x_points[i - 2], y_points[i - 2], s_points[i - 2], tau0) + f1(t_points[i - 3], x_points[i - 3], y_points[i - 3], s_points[i - 3], tau0));
        y_points[i] = y_points[i - 1] + (step / 24.0) * (9 * f2(t_points[i], x_pre, y_pre, s_pre, tau0) + 19 * f2(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0) - 5 * f2(t_points[i - 2], x_points[i - 2], y_points[i - 2], s_points[i - 2], tau0) + f2(t_points[i - 3], x_points[i - 3], y_points[i - 3], s_points[i - 3], tau0));
        s_points[i] = s_points[i - 1] + (step / 24.0) * (9 * f3(t_points[i], x_pre, y_pre, s_pre, tau0) + 19 * f3(t_points[i - 1], x_points[i - 1], y_points[i - 1], s_points[i - 1], tau0) - 5 * f3(t_points[i - 2], x_points[i - 2], y_points[i - 2], s_points[i - 2], tau0) + f3(t_points[i - 3], x_points[i - 3], y_points[i - 3], s_points[i - 3], tau0));
    }
}
int main()
{
    double *t_points = malloc(sizeof(double) * n), *x_points = malloc(sizeof(double) * n), *y_points = malloc(sizeof(double) * n), *s_points = malloc(sizeof(double) * n);
    FILE *f = fopen("test.txt", "w"), *f1 = fopen("answer.txt", "w");
    t_filling(t_points);
    // Eiler(x_points,u1_points,u2_points);
    // ABM(x_points, u1_points, u2_points);
    runhe_kutte_4(t_points, x_points, y_points, s_points, n);
    // ABM(t_points, x_points, y_points, s_points, n);
    for (int i = 0; i < n; ++i)
    {
        //  printf("%lf", x_points[i]);
    }
    // answer(answer_points, x_points);
    // x_points[0] = 0;
    file_write(f, f, t_points, s_points);
    // file_write(f1, f1, x_points, answer_points);
}