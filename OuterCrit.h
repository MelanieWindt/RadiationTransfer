#ifndef __OUTERCRIT_H__
#define __OUTERCRIT_H__

#include "meshProcessor/mesh.h"

class OuterCrit {
	mesh3d::vector omega;
public:
	OuterCrit (mesh3d::vector omega) : omega (omega) {
	}

	bool isOuter (const mesh3d::face & f) const {
		int num1 = f.idx () ;
		mesh3d::vector norm1 = f.normal();
		int num2 = f.flip () .idx ();
		mesh3d::vector norm2 = f.flip().normal();
		
		if (num1 > num2) {
			if (omega.dot(norm1)<0)
				return true;
			else return false; 
		} else {
			if (omega.dot(norm2)<0)
 				return false;
			else return true; 
		}
	}
};

#endif