#include "objects.h"

unsigned int indices[] = {
	 0,1,2,
	 2,3,0
};

float pos[8] = {
	0.0f,0.0f,
	1.0f,0.0f,
	1.0f,1.0f,
	0.0f,1.0f
};

LATTICE::LATTICE(int size, int x, int y, int posinitx, int posinity) 
{
	int i = 0;
	while (i < 8)
	{
		if ((i + 1) % 2 != 0) {
			positions[i] = posinitx + x + (size * pos[i]);	//translate position in field to position in window for quad
		}
		else {
			positions[i] = posinity + y + (size * pos[i]);
		}
		i += 1;
	}
	size = size;
}

LATTICE::Divergence_Vitesse_Cell()
{

	// v par rapport à la cellule----------------- Finite difference approximation
	float gradient_x = (v_right - v_left) / size;
	float gradient_y = (v_down - v_up) / size;

	float gradient = gradient_x + gradient_y;

	return gradient;
}

FIELD::FIELD(int pos1x, int pos1y, int Resx, int Resy, int Size)
{
	InPosx = pos1x; InPosy = pos1y; resx = Resx; resy = Resy; size = Size;
}

void FIELD::CreateField() 
{
	for (int x = 0;x < resx;x++)
	{
		vector_field.push_back(std::vector<LATTICE>());
		for (int y = 0; y < resy; y++)
		{
			vector_field[x].push_back(LATTICE(size, x, y, InPosx, InPosy)); //intialize field
		}
	}
}

void FIELD::EnableBorder() 
{
	for (int i = 0; i<resx;i++) 
	{
		vector_field[i][0].is_border = true;
		vector_field[i][resy - 2].is_border = true;
	}
}

void FIELD::InFlow(float speed_ms)
{
	for (int i = 1; i < resy - 1;i++)
	{
		vector_field[0][i].v_right = speed_ms;
	}
}