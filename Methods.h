#ifndef __METHODS_H__
#define __METHODS_H__

template<class Problem>
struct BaseMethod {
    typedef typename Problem::vec vec;
    typedef typename Problem::real real;
    vec F(const Problem &prob, const vec &U1, const vec &U2, const real dt_h) const;
};

template<class Problem>
struct LaxWendroff : public BaseMethod<Problem> {
    typedef typename Problem::vec vec;
    typedef typename Problem::real real;
    vec F(const Problem &prob, const vec &U1, const vec &U2, const real dt_h) const {
        vec U = static_cast<real>(0.5) * (U1 + U2 + dt_h * (prob.F(U1) - prob.F(U2)));
        return prob.F(U);
    }
};

template<class Problem>
struct CourantIsaaksonRees : public BaseMethod<Problem> {
    typedef typename Problem::vec vec;
    typedef typename Problem::real real;
    vec F(const Problem &prob, const vec &U1, const vec &U2, const real) const {
        const auto half = static_cast<real>(0.5);
        const auto &wU = half * (U1 + U2);
        const auto &lam = prob.lambda(wU);
        const auto &Aabs = prob.fromDiag(wU, lam.cwiseAbs());
        return half * (prob.F(U1) + prob.F(U2) + Aabs * (U1 - U2));
    }
};

template<class Problem>
struct LaxFriedrichs : public BaseMethod<Problem> {
    typedef typename Problem::vec vec;
    typedef typename Problem::real real;
    vec F(const Problem &prob, const vec &U1, const vec &U2, const real) const {
        const auto half = static_cast<real>(0.5);
        const auto &wU = half * (U1 + U2);
        const auto &lam = prob.lambda(wU);
        const auto Amax = lam.cwiseAbs().maxCoeff();
        return half * (prob.F(U1) + prob.F(U2) + Amax * (U1 - U2));
    }
};

#endif
