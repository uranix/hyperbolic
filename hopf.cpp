#include <iostream>

#include "BaseProblem.h"
#include "Limiters.h"
#include "Methods.h"
#include "Solver.h"

struct Hopf : public BaseProblem<float, 1, Hopf> {
    vec F(const vec &U) const {
        vec ret;
        ret << 0.5 * U[0] * U[0];
        return ret;
    }
    vec lambda() const {
        vec ret;
        ret << U()[0];
        return ret;
    }
    mat omega() const {
        mat ret;
        ret << 1;
        return ret;
    }
};

int main() {
    const int N = 50;
    const float h = 1.0f / N;
    const float C = 0.5f;

    Solver<Hopf> solver{Hopf()};

    std::vector<Solver<Hopf>::Cell> cells(N);
    std::vector<Solver<Hopf>::Face> faces(N + 1);

    for (size_t i = 0; i < N; i++) {
        float x = (i + 0.5) * h;
        cells[i].U[0] = sin(8 * atan(1) * x);
    }

    for (int steps = 0; steps < 20; steps ++) {
        float amax = 0;
        for (size_t i = 1; i < N; i++) {
            float a = solver.prepareAndComputeAmax(cells[i - 1], cells[i], faces[i]);
            if (a > amax)
                amax = a;
        }
        cells.front().Q.fill(0);
        cells.back().P.fill(0);
        const float dt_h = C / amax;

        for (size_t i = 1; i < N; i++)
            solver.computeLimitedFluxes<VanAlbada, LaxFriedrichs, LaxWendroff>(cells[i - 1], cells[i], faces[i], dt_h);

        for (size_t i = 1; i < N - 1; i++)
            solver.integrateFluxes(faces[i], faces[i+1], cells[i], dt_h);
    }

    for (size_t i = 0; i < N; i++) {
        float x = (i + 0.5) * h;
        std::cout << x << " " << cells[i].U[0] << std::endl;
    }

    return 0;
}
