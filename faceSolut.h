#ifndef __FACESOLUT_H__
#define __FACESOLUT_H__

#include "meshProcessor/mesh.h"

struct faceSolut {
	mesh3d::vector a;
	double b;


	double evaluate (const mesh3d::vector &point) {
		return a.dot(point) + b;
	}

};

#endif