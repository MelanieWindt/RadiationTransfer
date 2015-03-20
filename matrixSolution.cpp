#include "matrixSolution.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

void solve(const int n, double a[], double b[], double x[]) {
    for (int k = 0; k < n; k++) {
        int ks = 0;
        double piv = 0;
        for (int i = k, ik = k * n + k; i < n; i++, ik += n)
            if (fabs(a[ik]) > piv) {
                piv = fabs(a[ik]);
                ks = i;
            }
        if (piv == 0)
            throw std::invalid_argument("Matrix is singular");
        std::swap(b[k], b[ks]);
        for (int ki = k * n, ksi = ks * n, i = 0; i < n; i++, ki++, ksi++)
            std::swap(a[ki], a[ksi]);
        double akk = a[k * n + k];
        for (int i = k + 1, ik = (k + 1) * n + k; i < n; i++, ik += n) {
            a[ik] /= akk;
            for (int j = k + 1, ij = i * n + k + 1, kj = k * n + k + 1; j < n; j++, ij++, kj++)
                a[ij] -= a[ik] * a[kj];
        }
    }
 
    for (int i = 0; i < n; i++)
        for (int j = 0, ij = i * n; j < i; j++, ij++)
            b[i] -= a[ij] * b[j];
 
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1, ij = i * n + i + 1; j < n; j++, ij++)
            b[i] -= a[ij] * b[j];
        b[i] /= a[i * n + i];
    }
 
 
    if (x != b)
        std::copy(b, b + n, x);
}