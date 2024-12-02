#define Alfas 2.7
#define Alfay 0.8
#define Betax 0.9
#define Betay 1.2
#define Betas 0.9
#define Alfaxy 1.5
#define n 4
double ksi(double t)
{
    return 1.0;
}

double f1(double t, double x, double y, double s, double tau)
{
    return Betax * (pow(s, n) / (1.0 + pow(s, n))) * ksi(t) - Alfaxy * x * y;
}
double f2(double t, double x, double y, double s, double tau)
{
    return Betay * x * (t - tau) * ksi(t - tau) - Alfay * y;
}
double f3(double t, double x, double y, double s, double tau)
{
    return Betas - Alfas * y * s;
}
int main()
{
    double *x_points = malloc(sizeof(double) * n), *y_points = malloc(sizeof(double) * n), *s_points = malloc(sizeof(double) * n), *t_points = malloc(sizeof(double) * n);
    FILE *f = fopen("test.txt", "w"), *f1 = fopen("answer.txt", "w");
    x_filling(x_points);
    // Eiler(x_points,u1_points,u2_points);
    ABM(x_points, u1_points, u2_points);
    for (int i = 0; i < n; ++i)
    {
        //  printf("%lf", x_points[i]);
    }
    answer(answer_points, x_points);
    // x_points[0] = 0;
    file_write(f, f, x_points, u1_points);
    file_write(f1, f1, x_points, answer_points);
}
