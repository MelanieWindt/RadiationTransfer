#include "coloring.h"
#include <queue>
#include <cassert>
#include <iostream>

using namespace mesh3d;
void print_queue(const std::queue<int> &qorig);

int maximum (const std::vector <int> &vect) {
	int max = vect[0];
	for (size_t i = 1; i < vect.size (); i++) {
		if (max < vect[i]) 
			max = vect [i];
	}
	return max;
}

bool isColorable (const tetrahedron &tet, const std::vector<int> &colors, const OuterCrit &cr) {

	for (int i = 0; i < 4; i++) {
		const face & f = tet.f(i);
		if (cr.isOuter(f)) {
			continue;
		}
		
		const face & flipped = f.flip ();
		if (flipped.is_border()) {
			continue;
		}

		if (colors [flipped.tet ().idx ()] == -2)
			return false;
	}
	return true;
}

void paint (const tetrahedron &tet, std::vector<int> &colors, const OuterCrit &cr, std::queue <int> &q, int &cnt) {
	if (colors[tet.idx()] != -2) {
		return;
	}

	if (!isColorable(tet, colors, cr)) {
		q.push (tet.idx ());
		cnt--;
		//std::cout << "cnt = " << cnt << std::endl;
		//print_queue(q);
		return;
	}


	std::vector<int> nums (4);
	for (int i = 0; i < 4; i++) {
		const face & f = tet.f(i);

		if (cr.isOuter(f)) {
			nums[i] = -1;
			continue;
		}

		const face & flipped = f.flip ();
		if (flipped.is_border()) {
			nums[i] = -1;
			continue;
		}

		nums[i] = colors[flipped.tet ().idx ()];
	}


	for (int i = 0; i < 4; i++) {
		const face & f = tet.f(i);
		const face & flipped = f.flip();
		if (!cr.isOuter (f)) 
			continue; 

		if (flipped.is_border())
			continue;

		if (colors[flipped.tet ().idx () ] == -2 ) {
			q.push (flipped.tet ().idx () );
		}
	}

	cnt = q.size();
	//std::cout << "cnt = " << cnt << std::endl;
	//print_queue(q);

	assert(colors[tet.idx()] == -2);
	colors[tet.idx ()] = maximum(nums) + 1;
} 


void initQueue (const mesh & m, std::queue <int> &q, const OuterCrit &cr) {
	 for (size_t i = 0; i < m.faces () .size () ; i++)  {
	 	const face &f = m.faces()[i];
	 	if (!f.is_border()) 
	 		continue;
	 	if (!cr.isOuter(f))
	 		continue;
	 	q.push (f.flip ().tet ().idx()) ;
	 }
}

void print_queue(const std::queue<int> &qorig) {
	std::queue<int> q(qorig);

	std::cout << '[' << std::endl;
	while (!q.empty()) {
		std::cout << q.front() << ", ";
		q.pop();
	}
	std::cout << ']' << std::endl;
}

void paintGraph (const mesh & m, const OuterCrit &cr, std::vector<int> &colors) {
	colors.assign (m.tets().size(), -2) ;
	std::queue <int> q; 
	initQueue (m, q, cr);
	int cnt = q.size();
	while (!q.empty ()) {
		int current = q.front() ;
		q.pop ();
		paint (m.tets()[current], colors, cr, q, cnt) ;

		if (cnt == -1) {
			std::cout << "Cycles found" << std::endl;
			//std::cout << "queue size = " << q.size();
			exit (1);
		}

	}

	for (size_t i = 0; i < m.tets().size(); i++)  {
		if (colors [i] == -2) 
			std::cerr << "Uncolored tets found" << std::endl;
	}
}