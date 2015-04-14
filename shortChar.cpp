#include "meshProcessor/vol_mesh.h"
#include "meshProcessor/vtk_stream.h"
#include "meshProcessor/mesh.h"
#include <iostream>
#include <cstdlib>

#include "coloring.h"
#include "OuterCrit.h"
#include "LebedevQuad.h"
#include "faceSolut.h"
#include "traceInfo.h"
#include "matrixSolution.h"
#include "timer.h"

using namespace mesh3d;

double containTest (const face &f, vector point) {
	double S1 = 0.5 * ((f.p(1).r() - point) % (f.p(2).r() - point)).norm();     
	double S2 = 0.5 * ((f.p(0).r() - point) % (f.p(2).r() - point)).norm();   
	double S3 = 0.5 * ((f.p(0).r() - point) % (f.p(1).r() - point)).norm();   
	double S0 = 0.5 * ((f.p(1).r() - f.p(0).r()) % (f.p(2).r() - f.p(0).r())).norm();

	return (S1+S2+S3)/S0;

}

double traceFromPoint(const face &f, const vector &omega, vector initPoint, vector &Q) {
	vector norm = f.normal();
	double num = norm.dot((initPoint - f.p(0).r()));
	double den = norm.dot(omega);
	if (fabs(den) < 1e-8) 
		return -1;
	double l = num/ den ;
	Q = initPoint - l*omega;
	return l;
}

traceInfo traceTet (const tetrahedron &tet, const vector &omega, const vector &initPoint, const OuterCrit &cr) {
	double sigma[4];
	vector Qs[4];

	for (int i = 0; i < 4; i++) {
		if (cr.isOuter(tet.f(i)))  {
			sigma[i] = 2;
			continue;
		}
		traceFromPoint(tet.f(i), omega, initPoint, Qs[i]);
		sigma[i] = containTest(tet.f(i), Qs[i]);
	}

	double x = sigma[0];
	int num = 0;
	for (int i = 1; i < 4; i++) {
		if (sigma[i] < x) {
			x = sigma[i];
			num = i;
		}
	}

	return traceInfo(num, Qs[num]);
}

double kappa_by_color (index color, int ifreq) {
	if (color == 1)
		return 10;
	return 0.01/(1 + ifreq);
}

double Ieq_by_color (index color, int ifreq) {
	(void)ifreq;
	if (color == 1)
		return 1;
	return 0;
}

void solveEq(const vector points[7], double v[10]) {
	double A [10][10];
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			A[i][j] = 0;
		}
	}
	for (int i = 0; i < 7; i++) {
		double x1 = points[i].x;
		double x2 = points[i].y;
		double x3 = points[i].z;
		A[i][0] = 0.5*x1*x1;
		A[i][1] = 0.5*x2*x2;
		A[i][2] = 0.5*x3*x3;
		A[i][3] = x2*x3;
		A[i][4] = x1*x3;
		A[i][5] = x2*x1;
		A[i][6] = x1;
		A[i][7] = x2;
		A[i][8] = x3;
		A[i][9] = 1;
	}

	A[7][0] = A[8][5] = A[9][4] = points[6].x;
	A[7][5] = A[8][1] = A[9][3] = points[6].y;
	A[7][4] = A[8][3] = A[9][2] = points[6].z;

	A[7][6] = A[8][7] = A[9][8] = 1;
	solve (10, (double *)A, v, v);
}

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

void one_dir(const mesh &m, const vector &omega, std::vector<double> &one_dir_sol) {

	Timer t;

	std::cout << "Processing direction " << omega << std::endl;

	OuterCrit cr (omega);

	std::vector<int> colors;

	Timer d;

	paintGraph(m, cr, colors);

	std::cout << "one_dir function, painting graph " << d.stopAndGetElapsedTime () << std::endl;

	int maxColor = 0;
	for (size_t i = 0; i < colors.size() ; i++) {
		if (colors[i] > maxColor ) { 
			maxColor = colors[i];
		}
	}

	std::vector<std::vector<int> > order(maxColor+1);

	for (size_t i = 0; i < colors.size(); i++) {
		order[colors[i]].push_back(i);
	}

	std::vector<std::vector<faceSolut> > solution(m.tets().size(), std::vector<faceSolut>(4));

	Timer b;

	for (int color = 0; color < maxColor+1; color++){
		for(size_t j = 0; j < order[color].size(); j++){
			int tetNum = order[color][j];
			const tetrahedron &tet = m.tets(tetNum);
			one_dir_sol[tetNum] = 0;
			for (int k = 0; k < 4; k++) {
				if (cr.isOuter(tet.f(k)))
					continue;
				const face &flipped = tet.f(k).flip();
				if (flipped.is_border()) {
					solution[tetNum][k].setZero();
				}
				else {
					int flipTetNum = flipped.tet().idx();
					int faceNum = flipped.face_local_index();

					vector pointsOrig [3];
					vector pointsFlip [6];

					for (int bar = 0; bar < 3; bar ++) {
						pointsOrig [bar] = tet.f(k).p(bar).r();
						pointsFlip [bar] = tet.f(k).flip().p(bar).r();
					}

					pointsFlip [5] = 0.5*(pointsFlip[0]+pointsFlip[1]);
					pointsFlip [3] = 0.5*(pointsFlip[1]+pointsFlip[2]);
					pointsFlip [4] = 0.5*(pointsFlip[0]+pointsFlip[2]);

					solution[tetNum][k].copy_from(solution[flipTetNum][faceNum], pointsOrig, pointsFlip);
				}
			}
			for (int k = 0; k < 4; k++) {
				if (!cr.isOuter(tet.f(k)))
					continue;

				vector points [6];
				double v[6][NFREQ];

				for (int i = 0; i < 3; i++) {
					points[i] = tet.f(k).p(i).r();
				}

				points [5] = 0.5*(points[0]+points[1]);
				points [3] = 0.5*(points[1]+points[2]);
				points [4] = 0.5*(points[0]+points[2]);
				 
				for(int i = 0; i < 6; i++){
					const vector &p = points[i];

					traceInfo qInfo = traceTet(tet, omega, p, cr);

					vector basePoints[3];
					for (int bk = 0; bk < 3; bk++) {
						basePoints[bk] = tet.f(qInfo.num).p(bk).r();
					}
					
					for (int ifreq = 0; ifreq < NFREQ; ifreq++) {
						double solutQ = solution[tetNum][qInfo.num].evaluate(qInfo.point, basePoints, ifreq);
						double delta = norm(p - qInfo.point);
						double eKappeDelta = exp(-kappa_by_color(tet.color(), ifreq)*delta);
						double Ieq = Ieq_by_color(tet.color(), ifreq);
						double solutP = solutQ*eKappeDelta + Ieq*(1 - eKappeDelta);

						v[i][ifreq] = solutP;
					}

				}

#if ORDER == 2
				for (int ifreq = 0; ifreq < NFREQ; ifreq ++) {
					v[5][ifreq] = correctValue(v[5][ifreq], v[0][ifreq], v[1][ifreq]);
					v[3][ifreq] = correctValue(v[3][ifreq], v[1][ifreq], v[2][ifreq]);
					v[4][ifreq] = correctValue(v[4][ifreq], v[0][ifreq], v[2][ifreq]);
				}
#endif

					solution[tetNum][k].set(v);
	
			}

			for (int k = 0; k < 4; k++) {
				vector points [3];

				for (int i = 0; i < 3; i++) {
					points[i] = tet.f(k).p(i).r();
				}

				for (int ifreq = 0; ifreq < NFREQ; ifreq ++) {
					one_dir_sol[tetNum] += 0.25*solution[tetNum][k].evaluate(tet.f(k).center(), points, ifreq);
				}
			}
		}
	}
	std::cout << "one_dir function, cycle by colors " << b.stopAndGetElapsedTime () << std::endl;
	std::cout << "steps = " << maxColor + 1 << std::endl;
	std::cout << "one_dir function " << t.stopAndGetElapsedTime () << std::endl;
}

int main() {
	try {
		basis::save();
		vol_mesh vm("../newMesh300k.vol");
		std::cout << "opening mesh" << std::endl;

		mesh m(vm);
		bool res = m.check(&std::cout);
		std::cout << "Mesh check: " << (res ? "OK" : "failed") << std::endl;
		
		LebedevQuad quad(5);
		std::cout << "Using " << quad.order << " directions" << std::endl;
		std::vector<double> U(m.tets().size(), 0);
		std::vector<double> I(m.tets().size());

		for (int s = 0; s < quad.order; s++) {
			one_dir(m, vector(quad.x[s], quad.y[s], quad.z[s]), I);
			for (size_t i = 0; i < m.tets().size(); i++)
				U[i] += quad.w[s] * I[i];
			if (s == 0) {
				vtk_stream vtk("onedir.vtk");
				vtk.write_header(m, "I_1");
				vtk.append_cell_data(I.data(), "I");
				vtk.close(); 
			}
		}
	
		vtk_stream vtk("mesh.vtk");
		vtk.write_header(m, "U");
		vtk.append_cell_data(U.data(), "U");
		vtk.close(); 
	}
	catch (std::exception &e) {
		std::cerr << "Exception occured: " << e.what() << std::endl;
	} 
	return 0;
}
