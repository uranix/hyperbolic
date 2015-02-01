#include <iostream>

#include "BaseProblem.h"
#include "Limiters.h"
#include "Methods.h"
#include "Solver.h"

struct Euler : public BaseProblem<double, 3, Euler> {
    double gamma;
    Euler(double gamma) : gamma(gamma) { }
    double pressure(double rho, double eps) const {
        return (gamma - 1) * rho * eps;
    }
    double sound(double rho, double eps) const {
        (void)rho;
        return sqrt(gamma * (gamma - 1) * eps);
    }
    vec F(const vec &U) const {
        double rho = U[0];
        double rhomin = static_cast<double>(1e-6);
        double epsmin = static_cast<double>(1e-6);
        if (rho < rhomin)
            rho = rhomin;
        double u = U[1] / rho;
        double eps = U[2] / rho - u * u / 2;
        if (eps < epsmin)
            epsmin = eps;
        double p = pressure(rho, eps);
        return vec(rho * u, rho * u * u + p, (rho * (eps + u * u / 2) + p) * u);
    }
    vec lambda(const vec &Umid) const {
        double rho = Umid[0];
        double rhomin = static_cast<double>(1e-6);
        double epsmin = static_cast<double>(1e-6);
        if (rho < rhomin)
            rho = rhomin;
        double u = Umid[1] / rho;
        double eps = Umid[2] / rho - u * u / 2;
        if (eps < epsmin)
            epsmin = eps;

        double c = sound(rho, eps);

        return vec(u - c, u, u + c);
    }
    mat omega(const vec &Umid) const {
        mat ret;
        double rho = Umid[0];
        double rhomin = static_cast<double>(1e-6);
        double epsmin = static_cast<double>(1e-6);
        if (rho < rhomin)
            rho = rhomin;
        double u = Umid[1] / rho;
        double eps = Umid[2] / rho - u * u / 2;
        if (eps < epsmin)
            epsmin = eps;

        double c = sound(rho, eps);
        double s = c / (gamma - 1);
        ret <<
            u * (u / 2 + s),   -s - u, static_cast<double>(1),
            u * u / 2 - s * c,     -u, static_cast<double>(1),
            u * (u / 2 - s),    s - u, static_cast<double>(1);
        return ret;
    }
};

int main() {
    const int N = 50;
    const double h = 1.0f / N;
    const double C = 0.5f;

    Solver<Euler> solver(Euler(1.4));

    std::vector<Cell<Euler>> cells(N);
    std::vector<CellExtra<Euler>> cex(N);
    std::vector<Face<Euler>> faces(N + 1);

    for (size_t i = 0; i < N; i++) {
        double x = (i + 0.5) * h;
        cells[i].U[0] = x < 0.5 ? 1 : 0.1;
        cells[i].U[1] = x < 0.5 ? 0 : 0;
        cells[i].U[2] = x < 0.5 ?
            (  1 / (1.4 - 1)) :
            (0.1 / (1.4 - 1));
    }

    for (int steps = 0; steps < N; steps ++) {
        double amax = 0;
        for (size_t i = 1; i < N; i++) {
            double a = solver.prepareAndComputeAmax(
                    cells[i - 1], cells[i],
                    cex[i - 1], cex[i],
                    faces[i]);
            if (a > amax)
                amax = a;
        }
        cex.front().L.fill(0);
        cex.back().R.fill(0);
        const double dt_h = C / amax;

        for (size_t i = 1; i < N; i++)
            solver.computeLimitedFluxes<VanLeer, CourantIsaaksonRees, LaxWendroff>(
                    cells[i - 1], cells[i],
                    cex[i - 1], cex[i],
                    faces[i], dt_h);

        for (size_t i = 1; i < N - 1; i++)
            solver.integrateFluxes(faces[i], faces[i+1], cells[i], cells[i], dt_h);
    }

    for (size_t i = 0; i < N; i++) {
        double x = (i + 0.5) * h;
        std::cout << x << " " << cells[i].U[0] << std::endl;
    }

    return 0;
}
