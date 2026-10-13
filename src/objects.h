#include <iostream>
#include <vector>
#include <math.h>

struct LATTICE
{
	float v_up = 0;
	float v_down = 0;
	float v_left = 0;
	float v_right = 0;

	int size;
	
	float positions[8];

	bool is_border = false;

	LATTICE(int size, int x, int y, int posinitx, int posinity);

	float Divergence_Vitesse_Cell();
};

class FIELD
{
public:
	int resx;
	int resy;
	int size;

	int InPosx;
	int InPosy;

	std::vector<std::vector<LATTICE>> vector_field;

	FIELD(int pos1x,int pos1y,int Resx,int Resy,int Size);

	void CreateField();

	void EnableBorder();

	void InFlow(float speed_ms);
};