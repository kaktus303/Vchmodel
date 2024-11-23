#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define start 0.0
#define finish 10.0
#define Koshi 1.00
#define Const 1.0
#define counter 20
#define lamda 1000
#define start_n 8
void arrayequal(int n, double *a, double *b)
{
    for (int i = 0; i < n; ++i)
    {
        a[i] = b[i];
    }
}
char *give_file_way(int i)
{
    char *file_way = malloc(sizeof(char) * 10);
    strcpy(file_way, "graph");
    file_way[5] = i + '0';
    strcat(file_way, ".txt");
    return file_way;
}
char *give_error_way(int i)
{
    char *file_way = malloc(sizeof(char) * 10);
    strcpy(file_way, "error");
    file_way[5] = i + '0';
    strcat(file_way, ".txt");
    return file_way;
}
double test_function_fi(double x)
{
    return (sin(2 * x) + cos(x));
}
double test_function_fi_diff(double x)
{
    return (2 * cos(2 * x) - sin(x));
}
void file_write(FILE *file_graph, FILE *file_error, double *points_x, double *points_y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        fprintf(file_graph, "%.15lf %.15lf\n", points_x[i], points_y[i]);
        fprintf(file_error, "%.15lf %.15lf\n", points_x[i], fabs(points_y[i] - test_function_fi(points_x[i])));
    }
}
double f(double x, double y)
{
    return (test_function_fi_diff(x) + Const * (y - test_function_fi(x)));
}
void eiler(double step, double *points_x, double *points_y, int n)
{
    double x = start, y = test_function_fi(start);
    for (int i = 0; i < n; i++)
    {
        points_x[i] = x;
        points_y[i] = y;
        // file_write(file_graph, file_error, x, y);
        y += step * f(x, y);
        x += step;
    }
}
double eiler_1(double step, double x, double y)
{
    // file_write(file_graph, file_error, x, y);
    y += step * f(x, y);
    return y;
}
void runhe_kutte(double step, double *points_x, double *points_y, int n)
{
    double x = start, y = test_function_fi(start), k1 = 0, k2 = 0, k3 = 0;
    for (int i = 0; i < n; i++)
    {
        // file_write(file_graph, file_error, x, y);
        points_x[i] = x;
        points_y[i] = y;
        k1 = step * f(x, y);
        k2 = step * f(x + step / 2.0, y + k1 / 2.0);
        k3 = step * f(x + step, y - k1 + 2.0 * k2);
        y = y + (k1 + 4.0 * k2 + k3) / 6.0;
        x += step;
    }
}
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
void iterat(double step, double *points_x, double *points_y, int n, double epsilon, int iterations)
{
    double x = start, y = test_function_fi(start), x_old= 0 , y_old= 0, y_iterational;
    int addams_counter = 0;
    for (int i = 0; i < n; i++)
    {
        // file_write(file_graph, file_error, x, y);
        points_x[i] = x;
        points_y[i] = y;
        y_old = y;
        y += step * f(x, y);
        x_old = x;
        x += step;
        y_iterational = y;
        addams_counter = 0;
        // printf("Значение i = %d y_old = %e y_iterational = %e, y = %e\n", i, y_old, y_iterational, y);
        do
        {
            addams_counter++;
            y_iterational = y;
            y = y_old + (f(x_old, y_old) + f(x, y)) * (step / 2.0);
            // printf("Iteration counter = %d y_iterational = %e y = %e,       %e \n", addams_counter, y_iterational, y, (fabs(y_iterational - y) / (fabs(y_iterational))));
        } while (addams_counter < iterations && (fabs(y_iterational - y) / (fabs(y_iterational))) > epsilon);
        // printf("Value i = %d y_old = %e y_iterational = %e, y = %e x = %e x_old %e\n", i, y_old, y_iterational, y, x, x_old);
        // printf("%lf\n", y_iterational);
        // printf("%d\n", addams_counter);
    }
}

void addams_4(double step, double *points_x, double *points_y, int n)
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
void addams_3(double step, double *points_x, double *points_y, int n)
{
    double x = start, y = test_function_fi(start);
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
        y = y + (step / 24.0) * (55 * f(x, y) - 59 * f(points_x[i - 1], points_y[i - 1]) + 37 * f(points_x[i - 2], points_y[i - 2]) - 9 * f(points_x[i - 3], points_y[i - 3]));
        x = x + step;
    }
}
void Gir(double step, double *points_x, double *points_y, int n)
{
    double x = start, y = test_function_fi(start), func;
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
        x += step;
        y = (0.48 * step * (test_function_fi_diff(x) - Const * test_function_fi(x)) + 1.92 * y - 1.44 * points_y[i - 1] + 0.64 * points_y[i - 2] - points_y[i - 3] * 0.12) / (1.0 - 0.48 * Const * step);
    }
}
double max_error(double *points_x, double *points_y, int n)
{
    double max = 0;
    for (int i = 0; i < n; ++i)
    {
        if (fabs(points_y[i] - test_function_fi(points_x[i])) > max)
            max = fabs(points_y[i] - test_function_fi(points_x[i]));
    }
    return max;
}
void execuse_methods(int n, double step, double *points_x, double *points_y, char *number, int i)
{
    if (*number == '1')
    {
        if (i == 1)
            printf("Eiler method\n");
        eiler(step, points_x, points_y, n);
    }
    if (*number == '2')
    {
        if (i == 1)
            printf("Runhe-Kutte 3 method\n");
        runhe_kutte(step, points_x, points_y, n);
    }
    if (*number == '3')
    {
        if (i == 1)
            printf("Runhe-Kutte 4 method\n");
        runhe_kutte_4(step, points_x, points_y, n);
    }
    if (*number == '4')
    {
        if (i == 1)
            printf("Iterat method\n");
        iterat(step, points_x, points_y, n, 0.00001, 100);
    }
    if (*number == '5')
    {
        if (i == 1)
            printf("Addams 4 method\n");
        addams_4(step, points_x, points_y, n);
    }
    if (*number == '6')
    {
        if (i == 1)
            printf("Addams 3 method\n");
        addams_3(step, points_x, points_y, n);
    }
    if (*number == '7')
    {
        if (i == 1)
            printf("Gir 3 method\n");
        Gir(step, points_x, points_y, n);
    }
}
void allmethods(int n)
{
    FILE *graph, *error;
    graph = fopen("graphics/eiler_graph.txt", "w");
    error = fopen("graphics/eiler_error.txt", "w");
    double step = (finish - start) / (n - 1);
    double *points_x = malloc(sizeof(double) * n);
    double *points_y = malloc(sizeof(double) * n);
    eiler(step, points_x, points_y, n);
    file_write(graph, error, points_x, points_y, n);
    graph = fopen("graphics/runhe3_graph.txt", "w");
    error = fopen("graphics/runhe3_error.txt", "w");
    runhe_kutte(step, points_x, points_y, n);
    file_write(graph, error, points_x, points_y, n);
    graph = fopen("graphics/runhe4_graph.txt", "w");
    error = fopen("graphics/runhe4_error.txt", "w");
    runhe_kutte_4(step, points_x, points_y, n);
    file_write(graph, error, points_x, points_y, n);
    graph = fopen("graphics/iterat_graph.txt", "w");
    error = fopen("graphics/iterat_error.txt", "w");
    iterat(step, points_x, points_y, n, 0.000001, 100);
    file_write(graph, error, points_x, points_y, n);
    graph = fopen("graphics/addams3_graph.txt", "w");
    error = fopen("graphics/addams3_error.txt", "w");
    addams_3(step, points_x, points_y, n);
    file_write(graph, error, points_x, points_y, n);
    graph = fopen("graphics/addams4_graph.txt", "w");
    error = fopen("graphics/addams4_error.txt", "w");
    addams_4(step, points_x, points_y, n);
    file_write(graph, error, points_x, points_y, n);
    graph = fopen("graphics/gir_graph.txt", "w");
    error = fopen("graphics/gir_error.txt", "w");
    Gir(step, points_x, points_y, n);
    file_write(graph, error, points_x, points_y, n);
    free(points_x);
    free(points_y);
    fclose(graph);
    fclose(error);
}
int main(int argc, char *argv[])
{
    printf("Hello, World!\n");
    FILE *graph, *error, *information, *pogr;
    double step, *points_x, *points_y, epsilon_old = 1, *points_y_old, max_pogr = 0;
    int n = start_n, main_counter = 0;
    char y[10] = "graph";
    char copy[10] = "graph";
    char c;
    if (*argv[1] == '1')
    {
        for (int i = 1; i < counter + 1; ++i)
        {

            c = '0' + i;
            // graph = fopen(give_file_way(i), "w");
            // error = fopen(give_error_way(i), "w");
            graph = fopen("addams_with_eiler_graph.txt", "w");
            error = fopen("addams_with_eiler_error.txt", "w");
            step = (finish - start) / (n - 1);
            points_x = malloc(sizeof(double) * n);
            points_y = malloc(sizeof(double) * n);
            execuse_methods(n, step, points_x, points_y, argv[2], i);
            file_write(graph, error, points_x, points_y, n);
            printf("%d    %e\n", n, max_error(points_x, points_y, n));
            // printf("%d\n", n);
            // printf("%.6lf\n", log2(epsilon_old/max_error(points_x, points_y, n)));
            // fclose(graph);
            // fclose(error);
            // epsilon_old = max_error(points_x, points_y, n);
            // printf("222");
            // if (i > 2)
            // {
            //    /// printf("666");
            //     max_pogr = 0;
            //     for (int i = 0; i < (n / 2); ++i)
            //     {
            //         if ((fabs(points_y[i * 2] - points_y_old[i])) > max_pogr)
            //         {
            //             max_pogr = fabs(points_y[i * 2] - points_y_old[i]);
            //         }
            //     }
            //     printf("%e      %e\n", max_pogr*(pow(2.0, 4)/(pow(2.0,4) - 1.0)), max_error(points_x,points_y, n));
            // }

            // free(points_y_old);
            // points_y_old = malloc(sizeof(double) * n);
            // arrayequal(n, points_y_old, points_y);
            free(points_x);
            free(points_y);
            n *= 2;
        }
    }
    if (*argv[1] == '2')
    {
        allmethods(n);
    }
    if (*argv[1] == '3')
    {
        for (int i = 1; i < counter + 1; ++i)
        {
            if (i != 1)
            {
                max_pogr = 0;
                // printf("%d      %e      %e    ", n, step, max_error(points_x, points_y, (n - 1) / 2));
                // printf("%d\n", n);
                // printf("%e \n",  max_error(points_x, points_y, (n - 1) / 2));
                step = (finish - start) / (n - 1);
                points_x = malloc(sizeof(double) * n);
                points_y = malloc(sizeof(double) * n);
                execuse_methods(n, step, points_x, points_y, argv[2], i);
                for (int i = 0; i <= ((n - 1) / 2); ++i)
                {
                    if ((fabs(points_y[i * 2] - points_y_old[i])) > max_pogr)
                    {
                        max_pogr = fabs(points_y[i * 2] - points_y_old[i]);
                    }
                }
                printf("%e\n", max_pogr * (pow(2.0, 1) / (pow(2.0, 1) - 1.0)));

                // printf("%.6lf\n", log2(epsilon_old/max_error(points_x, points_y, n)));
                points_y_old = malloc(sizeof(double) * n);
                arrayequal(n, points_y_old, points_y);
                epsilon_old = max_error(points_x, points_y, n);
                n = n * 2 - 1;
            }
            else
            {
                step = (finish - start) / (n - 1);
                points_x = malloc(sizeof(double) * n);
                points_y = malloc(sizeof(double) * n);
                execuse_methods(n, step, points_x, points_y, argv[2], i);
                points_y_old = malloc(sizeof(double) * n);
                arrayequal(n, points_y_old, points_y);
                epsilon_old = max_error(points_x, points_y, n);
                n = n * 2 - 1;
            }
        }
    }
    if (*argv[1] == '4')
    {
        n = 8;
        FILE *iterat_file = fopen("iterat_file7.txt", "w"), *h = fopen("h.txt", "w"), *h2 = fopen("h^2.txt", "w");
        for (int i = 0; i < counter; ++i)
        {
            step = (finish - start) / (n - 1);
            points_x = malloc(sizeof(double) * n);
            points_y = malloc(sizeof(double) * n);
            iterat(step, points_x, points_y, n, 0.0000001, 1000);
            fprintf(iterat_file, "%d %.15lf\n", n, max_error(points_x, points_y, n));
            // sfprintf(h, "%d %.15lf\n",n, 10000*step);
            // fprintf(h2, "%d %.15lf\n",n, 3500*step*step);
            iterat(step, points_x, points_y, n, 0.00001, 0);
            fprintf(h2, "%d %.15lf\n", n, max_error(points_x, points_y, n));
            free(points_x);
            free(points_y);
            n *= 2;
        }
    }
    return 0;
}