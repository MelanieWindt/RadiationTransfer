#ifndef __LEBEDEVQUAD_H__
#define __LEBEDEVQUAD_H__

struct LebedevQuad {
	int order;
	double *x;
	double *y;
	double *z;
	double *w;
	LebedevQuad(int mindegree);
	~LebedevQuad();
};

#endif
