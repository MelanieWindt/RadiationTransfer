#ifndef __FACESOLUT_H__
#define __FACESOLUT_H__

#include "meshProcessor/mesh.h"

#include <algorithm>

#define ORDER 2
#define NFREQ 4

double correctValue(double solutP, double I1, double I2) {
	double min, max;
	if (I1 > I2) {
		max = 0.25*(3*I1 + I2);
		min = 0.25*(3*I2 + I1);
	}
	else {
		max = 0.25*(3*I2 + I1);
		min = 0.25*(3*I1 + I2);
	}

	if (solutP < min) 
		return min;
	else if (solutP > max)
		return max;
	else 
		return solutP;

}

namespace basis {

template<int n>
double l(double x, double y);

#if ORDER == 15
template<> double l<0>(double x, double y) {
	if (x+y > 0.5) 
		return 0;
	else 
		return 1 - 2*(x + y);
}

template<> double l<1>(double x, double y) {
	if ( x < 0.5)
		return 0;
	else 
		return 2*x - 1;
}

template<> double l<2>(double x, double y) { return l<1>(y, x); }

template<> double l<3>(double x, double y) {
	if (x + y < 0.5) 
		return 0;
	if (x > 0.5)
		return 2*y;
	if (y > 0.5)
		return 2*x;
	return 2*x + 2*y - 1;
}


template<> double l<4>(double x, double y) {
	if (x > 0.5) 
		return 0;
	if (x + y < 0.5)
		return 2*y;
	if (y > 0.5)
		return -2*x - 2*y + 2;
	return -2*x + 1;
}

template<> double l<5>(double x, double y) { return l<4>(y, x); } 
#endif

#if ORDER == 2
template<> double l<0>(double x, double y) { return (x + y - 1) * (2 * x + 2 * y - 1); }
template<> double l<1>(double x, double y) { return x * (2 * x - 1); }
template<> double l<2>(double x, double y) { return y * (2 * y - 1); }
template<> double l<3>(double x, double y) { return 4 * x * y; }
template<> double l<4>(double x, double y) { return -4 * y * (x + y - 1); }
template<> double l<5>(double x, double y) { return -4 * x * (x + y - 1); }
#endif

#if ORDER == 1
template<> double l<0>(double x, double y) { return 1 - x - y; }
template<> double l<1>(double x, double y) { return x; }
template<> double l<2>(double x, double y) { return y; }
template<> double l<3>(double x, double y) { return 0; }
template<> double l<4>(double x, double y) { return 0; }
template<> double l<5>(double x, double y) { return 0; }
#endif

template<typename T>
void put(std::ostream &o, T v) {
	union {
		T x;
		char buf[sizeof(T)];
	} w;
	w.x = v;
	std::reverse(w.buf, w.buf + sizeof(T));
	o.write(w.buf, sizeof(T));
}

void save() {
	const int M = 101;
	const double h = 1. / (M - 1);
	std::fstream vtk("basis.vtk", std::ios::out | std::ios::binary);
	vtk << "# vtk DataFile Version 3.0\n";
	vtk << "Finite element basis\n";
	vtk << "BINARY\n";
	vtk << "DATASET UNSTRUCTURED_GRID\n";
	vtk << "POINTS " << M * (M + 1) / 2 << " float\n";
	for (int i = 0; i < M; i++)
		for (int j = 0; j < M - i; j++) {
			put<float>(vtk, i * h);
			put<float>(vtk, j * h);
			put<float>(vtk, 0);
		}
	vtk << "\nCELLS " << (M - 1) * (M - 1) << " " << 4 * (M - 1) * (M - 1) << "\n";
	int k = 0;
	for (int i = 0; i < M - 1; i++) {
		for (int j = 0; j < M - i - 1; j++) {
			put<int>(vtk, 3);
			put<int>(vtk, k);
			put<int>(vtk, k + 1);
			put<int>(vtk, k + M - i);
			k++;
			if (j < M - i - 2) {
				put<int>(vtk, 3);
				put<int>(vtk, k);
				put<int>(vtk, k + M - i);
				put<int>(vtk, k + M - i - 1);
			}
		}
		k++;
	}
	vtk << "\nCELL_TYPES " << (M - 1) * (M - 1) << "\n";
	for (int i = 0; i < (M - 1) * (M - 1); i++)
		put<int>(vtk, 5);
	vtk << "\nPOINT_DATA " << M * (M + 1) / 2 << "\n";
#define PUT(z) \
	vtk << "SCALARS l" << z << " float 1\nLOOKUP_TABLE default\n"; \
	for (int i = 0; i < M; i++) \
		for (int j = 0; j < M - i; j++) { \
			double x = i * h; \
			double y = j * h; \
			put<float>(vtk, l<z>(x,y)); \
		}

	PUT(0);
	PUT(1);
	PUT(2);
	PUT(3);
	PUT(4);
	PUT(5);
}

}

struct faceSolut {
	double v [6][NFREQ];

	double evaluate (const mesh3d::vector &p, const mesh3d::vector points [], int ifreq) const {
		double ATr [2], ATAinv [2][2], detInv, x, y;
		const mesh3d::vector &p0 = points [0];
		const mesh3d::vector &p1 = points [1];
		const mesh3d::vector &p2 = points [2];

		const mesh3d::vector &p20 = p2 - p0;
		const mesh3d::vector &p10 = p1 - p0;
		const mesh3d::vector &r = p - p0;

		ATr [0] = p10.dot(r);
		ATr [1] = p20.dot(r);

		double ad = p10.norm2()*p20.norm2();
		double bc = std::pow(p20.dot(p10), 2);
		detInv = ad - bc;
		double costheta2 = bc / ad;

		detInv = 1 / detInv;

		ATAinv[1][1] = detInv * p10.dot(p10) ;
		ATAinv[0][1] = ATAinv [1][0] = -detInv *p10.dot(p20) ;
		ATAinv[0][0] = detInv * p20.dot(p20) ;

		x = ATAinv [0][0] * ATr [0] + ATAinv [0][1] * ATr [1];
		y = ATAinv [1][0] * ATr [0] + ATAinv [1][1] * ATr [1];

		double eps = 1e-6;
		if (x < -eps || y < -eps || x + y > 1 + eps) {
			std::cout << "Problem p* = " << p
				<< ", pcomp-p* = " << p0 - p + x * p10 + y * p20
				<< ", sin^2 theta = " << 1 - costheta2
				<< ", vol = " << r.dot(p10 % p20)
				<< ", xi = " << x
				<< ", eta = " << y 
				<< std::endl;
			__builtin_trap();
		}

		return 	v[0][ifreq]*basis::l<0>(x, y) + 
				v[1][ifreq]*basis::l<1>(x, y) + 
				v[2][ifreq]*basis::l<2>(x, y) + 
				v[3][ifreq]*basis::l<3>(x, y) + 
				v[4][ifreq]*basis::l<4>(x, y) +
				v[5][ifreq]*basis::l<5>(x, y);
	}

	void setZero () {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < NFREQ; j++) {
				v[i][j] = 0;
			}
		}
	}

	void set (double w [6][NFREQ]) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < NFREQ; j++) {
				v[i][j] = w[i][j]; 
			}
		}

#if ORDER == 2
		for (int ifreq = 0; ifreq < NFREQ; ifreq ++) {
			v[5][ifreq] = correctValue(v[5][ifreq], v[0][ifreq], v[1][ifreq]);
			v[3][ifreq] = correctValue(v[3][ifreq], v[1][ifreq], v[2][ifreq]);
			v[4][ifreq] = correctValue(v[4][ifreq], v[0][ifreq], v[2][ifreq]);
		}
#endif
	}

	void copy_from (const faceSolut &origSolut, const mesh3d::vector pointsOrig [3], const mesh3d::vector pointsFlip [6]) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < NFREQ; j++) {
				v[i][j] = origSolut.evaluate (pointsFlip[i], pointsOrig, j);
			}
		}	
	}

};

#endif