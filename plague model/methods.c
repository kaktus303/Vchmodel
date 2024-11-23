#include"methods.h"

void runhe_kutte_4(double step, double *points_x, double *points_y, int n)
{
    double x = start, y = test_function_fi(start), k1 = 0, k2 = 0, k3 = 0, k4 = 0;
    for (int i = 0; i < n; i++)
    {
        // file_write(file_graph, file_error, x, y);
        points_x[i] = x;
        points_y[i] = y;
        k1 = step * f(x, y);
        k2 = step * f(x + step / 2.0, y + k1 / 2.0);
        k3 = step * f(x + step / 2.0, y + k2 / 2.0);
        k4 = step * f(x + step, y + k3);
        y = y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        x += step;
    }
}
double runhe_kutte_4_1(double step, double x, double y)
{
    double k1 = 0, k2 = 0, k3 = 0, k4 = 0;
    k1 = step * f(x, y);
    k2 = step * f(x + step / 2.0, y + k1 / 2.0);
    k3 = step * f(x + step / 2.0, y + k2 / 2.0);
    k4 = step * f(x + step, y + k3);
    y = y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    return y;
}


void ABM(double step, double *points_x, double *points_y, int n)
{
    double x = start, y = test_function_fi(start), y_pre = 0;
    int addams_counter = 0;
    for (int i = 0; i < 4; ++i)
    {
        points_x[i] = x;
        points_y[i] = y;
        y = runhe_kutte_4_1(step, x, y);
        x += step;
    }
    // file_write(file_graph, file_error, x4, )
    for (int i = 4; i < n; i++)
    {
        points_x[i] = x;
        points_y[i] = y;
        y_pre = y + (step / 24.0) * (55 * f(x, y) - 59 * f(points_x[i - 1], points_y[i - 1]) + 37 * f(points_x[i - 2], points_y[i - 2]) - 9 * f(points_x[i - 3], points_y[i - 3]));
        x = x + step;
        y = y + (step / 24.0) * (9 * f(x, y_pre) + 19 * f(points_x[i], points_y[i]) - 5 * f(points_x[i - 1], points_y[i - 1]) + f(points_x[i - 2], points_y[i - 2]));
    }
}