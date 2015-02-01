#ifndef __SOLVER_H__
#define __SOLVER_H__

template<class Problem>
struct Cell {
    typename Problem::vec U;
};

template<class Problem>
struct CellExtra {
    typename Problem::vec R, L;
};

template<class Problem>
struct Face {
    typename Problem::vec wU, F;
};


template<class Problem, class Cell = Cell<Problem>, class CellExtra = CellExtra<Problem>, class Face = Face<Problem> >
struct Solver {
    const Problem prob;
    typedef Cell cell_type;
    typedef Face face_type;

    Solver(const Problem &prob) : prob(prob) { }

    typename Problem::real prepareAndComputeAmax(
            const Cell &l, const Cell &r,
            CellExtra &le, CellExtra &re, Face &f
        )
    {
        const auto &wU = static_cast<typename Problem::real>(0.5) * (l.U + r.U);

        f.wU = wU;
        const auto &Omega = prob.omega(wU);
        const auto &dW = Omega * (r.U - l.U);
        le.R = dW;
        re.L = dW;

        return prob.lambda(wU).cwiseAbs().maxCoeff();
    }

    template<
        template<typename T> class Limiter,
        template<typename Prob> class MethodLow,
        template<typename Prob> class MethodHigh
    >
    void computeLimitedFluxes(
            const Cell &l, const Cell &r,
            const CellExtra &le, const CellExtra &re,
            Face &f,
            typename Problem::real dt_h,
            const Limiter<typename Problem::real> &limiter = Limiter<typename Problem::real>(),
            const MethodLow <Problem> &methodLow  = MethodLow <Problem>(),
            const MethodHigh<Problem> &methodHigh = MethodHigh<Problem>()
        )
    {
        const auto delta = static_cast<typename Problem::real>(1e-6);
        const auto d2 = delta * delta;

        const auto &lam = prob.lambda(f.wU);
        typename Problem::vec phi;
        for (int i = 0; i < Problem::n; i++) {
            if (lam[i] > 0)
                phi[i] = le.R[i] / (le.R[i] * le.R[i] + d2) * limiter.psi(le.R[i], le.L[i]);
            else
                phi[i] = re.L[i] / (re.L[i] * re.L[i] + d2) * limiter.psi(re.R[i], re.L[i]);
        }
        const auto &Flo = methodLow .F(prob, l.U, r.U, dt_h);
        const auto &Fhi = methodHigh.F(prob, l.U, r.U, dt_h);

        const auto &Phi = prob.fromDiag(f.wU, phi);

        f.F = Flo + Phi * (Fhi - Flo);
    }

    void integrateFluxes(
            const Face &l, const Face &r, const Cell &from, Cell &to,
            typename Problem::real dt_h
        )
    {
        to.U = from.U + dt_h * (l.F - r.F);
    }
};

#endif
