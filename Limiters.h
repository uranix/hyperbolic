#ifndef __LIMITERS_H__
#define __LIMITERS_H__

#include <algorithm>

template<typename T, template <typename Q> class Limiter>
struct BaseLimiter {
    T phi(T r) const {
        return static_cast<const Limiter<T> *>(this)->phi(r);
    }
    T psi(T a, T b) const {
        if (a * b <= 0)
            return 0;
        if (fabs(a) < fabs(b))
            return b * phi(a / b);
        else
            return a * phi(b / a);
    }
};

template<typename T>
struct Minmod : public BaseLimiter<T, Minmod> {
    T phi(T r) const { return r; }
};

template<typename T>
struct MC : public BaseLimiter<T, MC> {
    T phi(T r) const { return std::min(2 * r, (1 + r) / 2); }
};

template<typename T>
struct Sweeby : public BaseLimiter<T, Sweeby> {
    const T beta;
    Sweeby(T beta) : beta(beta) { }
    T phi(T r) const { return std::min(beta * r, static_cast<T>(1)); }
};

template<typename T>
struct SuperBee : public Sweeby<T> {
    SuperBee() : Sweeby<T>(2) { }
};

template<typename T>
struct Ospre : public BaseLimiter<T, Ospre> {
    T phi(T r) const { return 3 * (r * r + r) / 2 / (r * r + r + 1); }
};

template<typename T>
struct VanAlbada : public BaseLimiter<T, VanAlbada> {
    T phi(T r) const { return (r * r + r) / (r * r + 1); }
};

template<typename T>
struct VanLeer : public BaseLimiter<T, VanLeer> {
    T phi(T r) const { return 2 * r / (1 + r); }
};

#endif
