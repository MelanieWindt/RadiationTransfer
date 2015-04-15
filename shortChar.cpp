#include "meshProcessor/vol_mesh.h"
#include "meshProcessor/vtk_stream.h"
#include "meshProcessor/mesh.h"
#include <iostream>
#include <fstream>
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

double kappa_by_coord (const mesh3d::vector &p, int ifreq) {
	if (p.norm() < 0.35)
		return 10;
	return 0.5/(1 + ifreq) * (.5 + .25 * std::atan(p.x) / std::atan(1));
}

double Ieq_by_coord (const mesh3d::vector &p, int ifreq) {

	if (p.norm() < 0.35)
		return 1 + 0.1 * ifreq;
	return 0;
}

std::vector<double> one_dir(const mesh &m, const vector &omega, std::vector< std::vector<double> > &one_dir_sol) {

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
						double eKappeDelta = exp(-kappa_by_coord(tet.center(), ifreq)*delta);
						double Ieq = Ieq_by_coord(tet.center(), ifreq);
						double solutP = solutQ*eKappeDelta + Ieq*(1 - eKappeDelta);

						v[i][ifreq] = solutP;
					}

				}

				solution[tetNum][k].set(v);
	
			}

			for (int k = 0; k < 4; k++) {

				for (int i = 0; i < 3; i++) {
					index v = tet.f(k).p(i).idx();
					for (int ifreq = 0; ifreq < NFREQ; ifreq++) {
						one_dir_sol[v][ifreq] = solution[tetNum][k].v[i][ifreq];
					}
				}
			}
		}
	}

	std::vector<double> integral(NFREQ);

	for (size_t i = 0; i < m.tets().size(); i++) {
		const tetrahedron &tet = m.tets(i);
		for (int j = 0; j < 4; j++) {
			const face &f = tet.f(j);
			if (!f.flip().is_border())
				continue;
			vector norm = f.normal();
			double normOmega = norm.dot(omega);
			if (normOmega > 0)
				continue;
			vector points [3];
 
			for (int k = 0; k < 3; k++) {
				points[k] = f.p(k).r();
			}
			for (int ifreq = 0; ifreq < NFREQ; ifreq++) {
				integral[ifreq] += -f.surface()*normOmega*solution[i][j].evaluate(f.center(), points, ifreq);
			}
		}
	}

	std::cout << "one_dir function, cycle by colors " << b.stopAndGetElapsedTime () << std::endl;
	std::cout << "steps = " << maxColor + 1 << std::endl;
	std::cout << "one_dir function " << t.stopAndGetElapsedTime () << std::endl;

	return integral;
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
		std::vector< std::vector <double> > U(m.vertices().size(), std::vector<double>(NFREQ, 0));
		std::vector< std::vector <double> > I(m.vertices().size(), std::vector<double>(NFREQ));

		std::vector <std::vector<double> > integral (quad.order);

		for (int s = 0; s < quad.order; s++) {
			integral [s] = one_dir(m, vector(quad.x[s], quad.y[s], quad.z[s]), I);
			for (size_t i = 0; i < I.size(); i++){
				for (int ifreq = 0; ifreq < NFREQ; ifreq ++) {
					U[i][ifreq] += quad.w[s] * I[i][ifreq];
				}
			}
			if (s == 0) {
				vtk_stream vtk("onedir.vtk");
				vtk.write_header(m, "I_1");
				for (int ifreq = 0; ifreq < NFREQ; ifreq ++) {
					std::vector<double> If(I.size());
					for (size_t i = 0; i < If.size(); i++)
						If[i] = I[i][ifreq];
					char name[10];
					sprintf(name, "I%d", ifreq);
					vtk.append_point_data(If.data(), name);
				}
				vtk.close(); 
			}
		}
	
		vtk_stream vtk("mesh.vtk");
		vtk.write_header(m, "U");
		for (int ifreq = 0; ifreq < NFREQ; ifreq ++) {
			std::vector<double> Uf(U.size());
				for (size_t i = 0; i < Uf.size(); i++)
					Uf[i] = U[i][ifreq];
				char name[10];
				sprintf(name, "U%d", ifreq);
				vtk.append_point_data(Uf.data(), name);
		}

		vtk.close();

		std::ofstream myfile;
	  	myfile.open ("1.txt");
	  	for (int i = 0; i < quad.order; i++) {
	  		for (int j = 0; j < NFREQ; j++) {
	  			myfile << integral[i][j] << " ";
	  		}
	  		myfile << std::endl;
	  	}
	  	myfile.close();
	}


	catch (std::exception &e) {
		std::cerr << "Exception occured: " << e.what() << std::endl;
	} 
	return 0;
}
