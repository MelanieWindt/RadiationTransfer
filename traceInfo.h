#ifndef __TRACEINFO_H__
#define __TRACEINFO_H__

#include "meshProcessor/mesh.h"

struct traceInfo {
	int num;
	mesh3d::vector point;
	traceInfo (const int num,const mesh3d::vector &point): num (num), point(point) {
	}

};

#endif