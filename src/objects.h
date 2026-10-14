#include <iostream>
#include <vector>
#include <math.h>

struct LATTICE
{
	int rel_x;
	int rel_y;

	int size;

	float divergence;
	
	float positions[8];

	bool is_border = false;

	LATTICE(int size, int x, int y, int posinitx, int posinity);
};

class FIELD
{
private:
	int resx;
	int resy;
	int size;

	int InPosx;
	int InPosy;

	float* VelocityX;
	float* VelocityY;

	std::vector<std::vector<LATTICE>> vector_field;
public:
	FIELD(int pos1x,int pos1y,int Resx,int Resy,int Size);
	~FIELD();

	void CreateField();

	void EnableBorder();

	void Divergence_Vitesse_Cell(int x,int y);

	void InFlow(float speed_ms);
};