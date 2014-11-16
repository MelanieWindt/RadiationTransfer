#ifndef __COLORING_H__
#define __COLORING_H__

#include <vector>
#include "meshProcessor/mesh.h"
#include "OuterCrit.h"

bool dfs (int v, const std::vector<std::vector<int> > & g, std::vector<int> & cl);
void paintGraph (const mesh3d::mesh & m, const OuterCrit &cr, std::vector<int> &colors);

#endif