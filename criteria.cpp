#include "meshProcessor/vol_mesh.h"
#include "meshProcessor/vtk_stream.h"
#include "meshProcessor/mesh.h"
#include <iostream>
#include <memory>
#include <algorithm>

using namespace mesh3d;

double det (double const A [3][3]) {
	return A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[1][0]*A[2][1]*A[0][2] 
			- A[0][2]*A[1][1]*A[2][0] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] ;

}

double findDelta (double const A [3][3], double const b [], int numb) {
	int num = numb -1;
	double matr [3][3];
	for (int i = 0; i<3; i++) {
		for (int j = 0; j < 3; j++) {
			matr[i][j] = A[i][j];
		}
		matr[i][num] = b[i];
	}
	return det (matr) ;

}


void solveSLE (double const A [3][3], double b [3]) {
	double det0 = det (A);
	double det1 = findDelta(A, b, 1);
	double det2 = findDelta(A, b, 2);
	double det3 = findDelta(A, b, 3);

	b[0] = det1/det0;
	b[1] = det2/det0;
	b[2] = det3/det0;

}

struct sortPrior {
	typedef std::pair <double, int> pair;
	bool operator () (const pair &a, const pair &b){ 
		return a.first < b.first;
	}

};

int main() {
	try {
		vol_mesh vm("../mesg.vol");

		mesh m(vm);
		bool res = m.check(&std::cout);
		std::cout << "Mesh check: " << (res ? "OK" : "failed") << std::endl;
		
		std::vector<double> projection (m.tets().size());
		std::vector<double> numbers (m.tets().size());
		std::vector<int>  swapNums (m.tets().size());
		std::vector<std::pair<double, int> > projOmg (m.tets().size());

		vector omega (1,1,1);

		for (index i = 0; i < m.tets().size(); i++) {
			vector p1 = m.tets()[i].p(0).r() ;
			vector p2 = m.tets()[i].p(1).r() ;
			vector p3 = m.tets()[i].p(2).r() ;
			vector p4 = m.tets()[i].p(3).r() ;

			double matr [3][3]; 
			vector c;
			double b [3];
			double radius;

			matr[0][0] = (p1 - p4).x;
			matr[0][1] = (p1 - p4).y;
			matr[0][2] = (p1 - p4).z;

			matr[1][0] = (p2 - p4).x;
			matr[1][1] = (p2 - p4).y;
			matr[1][2] = (p2 - p4).z;

			matr[2][0] = (p3 - p4).x;
			matr[2][1] = (p3 - p4).y;
			matr[2][2] = (p3 - p4).z;

			b[0] = 0.5*(p1.norm2() - p4.norm2());
			b[1] = 0.5*(p2.norm2() - p4.norm2());
			b[2] = 0.5*(p3.norm2() - p4.norm2());

			solveSLE(matr, b);
			c.x = b[0];
			c.y = b[1];
			c.z = b[2];

			projOmg [i].first = c.dot(omega);
			projection [i] = projOmg [i].first;
			projOmg [i].second = i;
		}

		std::sort(projOmg.begin(), projOmg.end(), sortPrior ());

		for (int i = 0; i < m.tets().size(); i++) {
			swapNums[projOmg[i].second] = i;

		}

		for (int i = 0; i < m.tets().size(); i++) {
			numbers[i] = swapNums[i];
		}

		for (index i = 0; i < m.tets().size(); i++) {
			const tetrahedron &tet = m.tets()[i];
			int s0 = swapNums[i];
			for (int j = 0; j < 4; j++) {
				const face &f = tet.f(j);
				double crit = f.normal().dot(omega);
				const face &f1 = f.flip();
				if(!f1.is_border()){
					int k = f1.tet().idx();
					int s1 = swapNums[k];
					if (crit < 0)
						std::cerr << s0 << " "<< s1 <<" 1"<< std::endl;
					if (((crit< 0 ) && (s1 < s0)) || ((crit > 0) && (s1 > s0))) {
					}					
				}
			}
		}

		vtk_stream vtk("mesh.vtk");
		vtk.write_header(m, "Test");
		vtk.append_cell_data(numbers.data(), "numbers");
		vtk.append_cell_data(projection.data(), "projection");
		vtk.close();
	}
	catch (std::exception &e) {
		std::cerr << "Exception occured: " << e.what() << std::endl;
	} 
	return 0;
}
