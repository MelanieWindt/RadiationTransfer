#ifndef __FACESOLUT_H__
#define __FACESOLUT_H__

#include "meshProcessor/mesh.h"

#include <algorithm>

namespace basis {

template<int n>
double l(double x, double y);

#define SECOND_ORD 1

#if SECOND_ORD
template<> double l<0>(double x, double y) { return (x + y - 1) * (2 * x + 2 * y - 1); }
template<> double l<1>(double x, double y) { return x * (2 * x - 1); }
template<> double l<2>(double x, double y) { return y * (2 * y - 1); }
template<> double l<3>(double x, double y) { return 4 * x * y; }
template<> double l<4>(double x, double y) { return -4 * y * (x + y - 1); }
template<> double l<5>(double x, double y) { return -4 * x * (x + y - 1); }
#else
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
	double v [6];

	double evaluate (const mesh3d::vector &p, const mesh3d::vector points []) const {
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

		return 	v[0]*basis::l<0>(x, y) + 
				v[1]*basis::l<1>(x, y) + 
				v[2]*basis::l<2>(x, y) + 
				v[3]*basis::l<3>(x, y) + 
				v[4]*basis::l<4>(x, y) +
				v[5]*basis::l<5>(x, y);
	}

	void setZero () {
		for (int i = 0; i < 6; i++) {
			v[i] = 0;
		}
	}

	void set (double w []) {
		for (int i = 0; i < 6; i++) {
			v[i] = w[i]; 
		}
	}

	void copy_from (const faceSolut &origSolut, const mesh3d::vector pointsOrig [3], const mesh3d::vector pointsFlip [6]) {
		for (int i = 0; i < 6; i++) {
			v[i] = origSolut.evaluate (pointsFlip[i], pointsOrig);
		}	
	}

};

#endif