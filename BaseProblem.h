#ifndef __BASE_PROBLEM_H__
#define __BASE_PROBLEM_H__

#include <Eigen/Core>
#include <Eigen/LU>

template<typename T, int _n, class Problem>
struct BaseProblem {
    static constexpr int n = _n;
    typedef Eigen::Matrix<T, n, n> mat;
    typedef Eigen::Matrix<T, n, 1> vec;
    typedef T real;
public:
    BaseProblem() { }
    mat omega(const vec &Umid) const;
    vec lambda(const vec &Umid) const;
    mat fromDiag(const vec &Umid, const vec &diag) const {
        const Problem &prob = static_cast<const Problem &>(*this);
        const auto &Omega = prob.omega(Umid);
        const auto &Diag = diag.asDiagonal();
        const auto &iOmega = Omega.inverse();
        return iOmega * Diag * Omega;
    }
    vec F(const vec &U) const;
};

#endif
