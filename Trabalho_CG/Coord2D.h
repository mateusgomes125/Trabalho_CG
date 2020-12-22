#ifndef COORD2D_H_INCLUDED
#define COORD2D_H_INCLUDED
class Coord2D
{
public:
	Coord2D() {
		x = 0;
		y = 0;
	}
	Coord2D(int _x, int _y) {
		x = _x;
		y = _y;
	}
	int x, y;
};
#endif
