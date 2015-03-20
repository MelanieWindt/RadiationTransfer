#ifndef __COLORING_H__
#define __COLORING_H__

#include <vector>
#include "meshProcessor/mesh.h"
#include "OuterCrit.h"

void paintGraph (const mesh3d::mesh & m, const OuterCrit &cr, std::vector<int> &colors);

#endif