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