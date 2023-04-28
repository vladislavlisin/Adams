#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

class differ
{
public:
    differ(double tau) : tau(tau), result_x(1), result_Vx(1), result_Vy(1), result_y(1)
    {
        result_x[0] = 7 * pow(10, -2);
        result_Vx[0] = 7 * pow(10, -2);
        result_y[0] = 1;
        result_Vy[0] = -1;
        Run();
    }

private:
    double tau;

    double ny = -0.5f;

    std::vector<double> result_Vx;
    std::vector<double> result_x;
    std::vector<double> result_Vy;
    std::vector<double> result_y;

// двухшаговый метод Адамса (для нахождения y_(n+1) использую метод Эйлера)

    void calc_x(int num)
    {
        double t = tau * num;
        double V1,V2,x1,x2;

        V1 = tau * f_Vx(t, result_x[num - 1], result_Vx[num - 1]) + result_Vx[num - 1];
        x1 = tau * f_x(t, result_x[num - 1], result_Vx[num - 1]) + result_x[num - 1];
        V2 = tau * ((3.0/2.0)*f_Vx(t + tau, x1, V1) - (1.0/2.0)*f_Vx(t, result_x[num - 1], result_Vx[num - 1])) + V1;
        x2 = tau * ((3.0/2.0)*f_x(t + tau, x1, V1) - (1.0/2.0)*f_x(t, result_x[num - 1], result_Vx[num - 1])) + x1;

        result_Vx.push_back(V2);
        result_x.push_back(x2);
    }

    void calc_y(int num)
    {
        double t = tau * num;
        double V1, V2, y1, y2;

        V1 = tau * f_Vy(t, result_y[num - 1], result_Vy[num - 1]) + result_Vy[num - 1];
        y1 = tau * f_y(t, result_y[num - 1], result_Vy[num - 1]) + result_y[num - 1];
        V2 = tau * ((3.0/2.0)*f_Vy(t + tau, y1, V1) - (1.0/2.0)*f_Vy(t, result_y[num - 1], result_Vy[num - 1])) + V1;
        y2 = tau * ((3.0/2.0)*f_y(t + tau, y1, V1) - (1.0/2.0)*f_y(t, result_y[num - 1], result_Vy[num - 1])) + y1;

        result_Vy.push_back(V2);
        result_y.push_back(y2);
    }

    double f_y(double t, double y, double v)
    {
        return v;
    }

    double f_Vy(double t, double y, double v)
    {
        return ny * (-2 * y);
    }

    double f_x(double t, double x, double v)
    {
        return v;
    }

    double f_Vx(double t, double x, double v)
    {
        return ny * (-2 * x);
    }

    void Run()
    {
        int i = 0;
        while (result_x[i] <= 1 && result_y[i] <= 1)
        {
            i++;
            calc_x(i);
            calc_y(i);
        }

        double cur_point = 0;
        double maxX = 0;
        double maxY = 0;
        double maxC = 0;

        for (int j = 0; j < i; j++)
        {
            maxX = std::max(maxX, pogr(result_x[j], analitX(cur_point)));
            maxY = std::max(maxY, pogr(result_y[j], analitY(cur_point)));
            maxC = std::max(maxC, pogr(result_x[j] * result_y[j], analitY(cur_point) * analitX(cur_point)));
            std::cout << std::fixed << std::setprecision(9)
            << pogr(result_x[j] * result_y[j], analitY(cur_point) * analitX(cur_point)) << std::endl;
            cur_point += tau;
        }

        cur_point = 0;
        std:: cout << "Approx points" << "   ---   " << "Exact points" << std:: endl;
        for(int k = 0; k < result_x.size(); k++)
        {
            std:: cout << "(" << result_x[k] << "," << result_y[k] << ")" << " --- "
            << "(" << analitX(cur_point) << "," << analitY(cur_point) << ")" <<  std:: endl;
            cur_point += tau;
        }

        std:: cout<< maxX << " " << maxY << " " << maxC << std::endl;
        std:: cout << "it" << i << std:: endl;
    }

    double analitX(double t)
    {
        return (7.0/100.0)*exp(t);
    }

    double analitY(double t)
    {
        return exp(-t);
    }

    double pogr(double cn, double ca)
    {
        return std::abs(((cn - ca) / ca)) * 100;
    }

};

// sol

#define ETA (-1/(double)2)
#define E_X (-2*x)
#define E_Y (-2*y)

// x^2 + y^2 = R^2
double Ex(double x, double y) {
    return ETA*E_X;
}

double Ey(double x, double y) {
    return ETA*E_Y;
}

double analitX(double t)
{
    return (2.0/100.0)*exp(t);
}

double analitY(double t)
{
    return exp(-t);
}

void adams (const double& n, const double& tau, const double& y_0, const double& x_0, const double& v_y_0,
                 const double& v_x_0)
{
    double y_half, y_n = 0;
    double y_old = y_0;
    double x_half, x_n = 0;
    double x_old = x_0;
    double v_y_half, v_y_n, v_x_half, v_x_n;
    double v_y_old = v_y_0;
    double v_x_old = v_x_0;
    double c_n = 0;
    int i = 0;
    double c_a = x_0*y_0;

    while(x_n <= 1 && y_n <= 1)
    {
        // первый этап метод эйлера первого порядка точности для координаты и скорости
        v_x_half = v_x_old + tau*Ex(x_old, y_old);
        v_y_half = v_y_old + tau*Ey(x_old, y_old);
        x_half = x_old + tau*v_x_old;
        y_half = y_old + tau*v_y_old;

        // сам Адомс
        v_x_n = (tau/2) * (3.0*Ex(x_half, y_half) - Ex(x_old, y_old)) + v_x_half;
        v_y_n = (tau/2) * (3.0*Ex(x_half, y_half) - Ex(x_old, y_old)) + v_y_half;
        x_n = (tau/2) * (3.0*v_x_half - v_x_old) + x_half;
        y_n = (tau/2) * (3.0*v_y_half - v_y_old) + y_half;

        y_old = y_half;
        x_old = x_half;
        v_x_old = v_x_half;
        v_y_old = v_y_half;

        c_n = x_n * y_n;
        std::cout << "(x,y): " << "(" << x_n << "," << y_n << ")" << " " << "(v_x,v_y): "
                  << "(" << v_x_n << "," << v_y_n << ")" << std::endl;
        std::cout << "error: " << std::abs((c_n - c_a)/c_a) << " %" << std::endl;
        i++;
    }

    std::cout << "it number: " << i << std::endl;
    x_n = 0;
    y_n = 0;
    i = 0;
    v_y_old = v_y_0;
    v_x_old = v_x_0;
    y_old = y_0;
    x_old = x_0;
    std::cout << std::endl << std::endl;

    while (x_n <= 1 && y_n <= 1)
    {
        v_y_half = v_y_old + tau*Ey(x_old, y_old);
        y_half = y_old + tau*v_y_old;
        v_x_half = v_x_old + tau*Ex(x_old, y_old);
        x_half = x_old + tau*v_x_old;

        v_y_n = v_y_old + (tau/2) * (Ey(x_old, y_old) + Ey(x_half, y_half));
        y_n = y_old + (tau/2) * (v_y_old + v_y_half);
        v_x_n = v_x_old + (tau/2) * (Ex(x_old, y_old) + Ex(x_half, y_half));
        x_n = x_old + (tau/2)*(v_x_old + v_x_half);

        y_old = y_n;
        x_old = x_n;
        v_x_old = v_x_n;
        v_y_old = v_y_n;

        c_n = x_n * y_n;
        std::cout << "(x,y): " << "(" << x_n << "," << y_n << ")" << " " << "(v_x,v_y): "
                  << "(" << v_x_n << "," << v_y_n << ")" << std::endl;
        std::cout << "error: " << 100 * std::abs((c_n - c_a)/c_a) << "%" << std::endl;
        i++;
    }

    std::cout << "it number: " << i << std::endl;
}

int main()
{
    adams(20, 0.2, 1, 0.012, -1, 0.012);
    return 0;
}






//    differ dif(0.2);