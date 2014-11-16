#include "meshProcessor/vol_mesh.h"
#include "meshProcessor/vtk_stream.h"
#include "meshProcessor/mesh.h"
#include <iostream>
#include <memory>

using namespace mesh3d;

int main() {
	try {
		vol_mesh vm("../mesg.vol");

		mesh m(vm);
		bool res = m.check(&std::cout);
		std::cout << "Mesh check: " << (res ? "OK" : "failed") << std::endl;
		
		std::vector<float> u(m.vertices().size());
		for (index i = 0; i < m.vertices().size(); i++) {
			u[i] = m.vertices()[i].r().norm();
		}

		std::vector<vec<double> > w(m.tets().size());
		for (index i = 0; i < m.tets().size(); i++) {
			const vector &r = m.tets()[i].center();
			w[i].x = r.z;
			w[i].y = -r.y;
			w[i].z = r.x;
		}

		vtk_stream vtk("mesh.vtk");
		vtk.write_header(m, "Test");
		vtk.append_cell_data(&w[0], "w");
		vtk.append_point_data(&u[0], "u");
		vtk.close();
	} catch (std::exception &e) {
		std::cerr << "Exception occured: " << e.what() << std::endl;
	}
	return 0;
}
