#ifndef __SOLVER_H__
#define __SOLVER_H__

template<class Problem>
struct Solver {
    struct Cell {
        typename Problem::vec U, P, Q;
    };

    struct Face {
        typename Problem::vec wU, F;
    };

    const Problem prob;

    Solver(const Problem &prob) : prob(prob) { }

    typename Problem::real prepareAndComputeAmax(
            Cell &l, Cell &r, Face &f
        )
    {
        const auto &wU = static_cast<typename Problem::real>(0.5) * (l.U + r.U);
        prob.assign(wU);

        f.wU = wU;
        const auto &Omega = prob.omega();
        const auto &dW = Omega * (r.U - l.U);
        l.P = dW;
        r.Q = dW;

        return prob.lambda().cwiseAbs().maxCoeff();
    }

    template<
        template<typename T> class Limiter,
        template<typename Prob> class MethodLow,
        template<typename Prob> class MethodHigh
    >
    void computeLimitedFluxes(
            const Cell &l, const Cell &r, Face &f,
            typename Problem::real dt_h,
            const Limiter<typename Problem::real> &limiter = Limiter<typename Problem::real>(),
            const MethodLow <Problem> &methodLow  = MethodLow <Problem>(),
            const MethodHigh<Problem> &methodHigh = MethodHigh<Problem>()
        )
    {
        prob.assign(f.wU);

        const auto delta = static_cast<typename Problem::real>(1e-6);
        const auto d2 = delta * delta;

        const auto &lam = prob.lambda();
        typename Problem::vec phi;
        for (int i = 0; i < Problem::n; i++) {
            if (lam[i] > 0)
                phi[i] = l.P[i] / (l.P[i] * l.P[i] + d2) * limiter.psi(l.P[i], l.Q[i]);
            else
                phi[i] = r.Q[i] / (r.Q[i] * r.Q[i] + d2) * limiter.psi(r.P[i], r.Q[i]);
        }
        const auto &Flo = methodLow .F(prob, l.U, r.U, dt_h);
        const auto &Fhi = methodHigh.F(prob, l.U, r.U, dt_h);

        const auto &Phi = prob.fromDiag(phi);

        f.F = Flo + Phi * (Fhi - Flo);
    }

    void integrateFluxes(
            const Face &l, const Face &r, Cell &c,
            typename Problem::real dt_h
        )
    {
        c.U += dt_h * (l.F - r.F);
    }
};

#endif
