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

private:
    mutable vec wU;
protected:
    vec U() const {
        return wU;
    }
public:
    BaseProblem() { }
    void assign(const vec &v) const { wU = v; }
    mat omega() const;
    vec lambda() const;
    mat fromDiag(const vec &diag) const {
        const Problem &prob = static_cast<const Problem &>(*this);
        const auto &Omega = prob.omega();
        const auto &Diag = diag.asDiagonal();
        const auto &iOmega = Omega.inverse();
        return iOmega * Diag * Omega;
    }
    vec F(const vec &U) const;
};

#endif
