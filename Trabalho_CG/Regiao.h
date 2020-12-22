#ifndef REGIAO_H_INCLUDED
#define REGIAO_H_INCLUDED

#include "Coord2D.h"
class Regiao
{
public:
	Regiao() {

	}
	Regiao(Coord2D mi, Coord2D mx) {
		cMin = mi;
		cMax = mx;
	}
	Coord2D cMin, cMax;
};

#endif
