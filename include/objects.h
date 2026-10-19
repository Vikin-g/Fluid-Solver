#include <iostream>
#include <vector>
#include <math.h>

struct CELL
{
	int rel_x;
	int rel_y;

	int size;

	float divergence;
	
	float positions[8];

	bool is_border = false;

	CELL(int size = 1, int x = 1, int y = 1, int posinitx = 1, int posinity = 1);
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

public:
	
	 CELL **vector_field;
	
	FIELD(int pos1x,int pos1y,int Resx,int Resy,int Size);
	~FIELD();

	void CreateField();

	void EnableBorder();

	void Divergence_Vitesse_Cell(int x,int y);

	void InFlow(float speed_ms);

	float GetVelocityAtPoint(int x,int y, bool horizontal);

	void ChangeVelocityAtPoint(int x, int y, bool horizontal, float n_v);	//to the left if x or under if y
};