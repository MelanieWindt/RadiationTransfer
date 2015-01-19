#ifndef __FACESOLUT_H__
#define __FACESOLUT_H__

#include "meshProcessor/mesh.h"

struct faceSolut {
	double d, c11, c22, c33, c13, c23, c12, b1, b2, b3;
 
	double evaluate (const mesh3d::vector &p) {
		double bx = b1*p.x + b2*p.y + b3*p.z;
		double xCx = c11*p.x*p.x + c22*p.y*p.y + c33*p.z*p.z + 2*(c12*p.x*p.y + c13*p.x*p.z + c23*p.z*p.y);
		return  0.5*xCx + bx + d;
	}

	void setZero () {
		d = c11 = c22 = c33 = c13 = c23 = c12 = b1 = b2 = b3 = 0;
	}

	void set(const double v[10]) {
		c11 = v[0];
		c22 = v[1];
		c33 = v[2];
		c23 = v[3];
		c13 = v[4];
		c12 = v[5];
		b1 = v[6];
		b2 = v[7];
		b3 = v[8];
		d =  v[9];
	}

};


#endif