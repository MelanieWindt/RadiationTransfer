#ifndef __FACESOLUT_H__
#define __FACESOLUT_H__

#include "meshProcessor/mesh.h"

struct faceSolut {
	double a;

	double operator () (const mesh3d::vector &point) {
		return a;
	}

};

#endif