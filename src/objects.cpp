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

	rel_x = x;
	rel_y = y;
}

FIELD::FIELD(int pos1x, int pos1y, int Resx, int Resy, int Size)
{
	InPosx = pos1x; InPosy = pos1y; resx = Resx; resy = Resy; size = Size;
}

void FIELD::CreateField() 
{

	VelocityX = new float[resx + 1][resy];
	VelocityY = new float[resx][resy + 1];

	for (int x = 0;x < resx;x++)
	{
		vector_field.push_back(std::vector<LATTICE>());
		for (int y = 0; y < resy; y++)
		{
			vector_field[x].push_back(LATTICE(size, x, y, InPosx, InPosy)); //intialize field
			VelocityY[x][y] = 0.0f;					//initialize array of velocities in-between cells
		}
		VelocityY[x][-1] = 0.0f;
	}

	for (int i = 0; i < VelocityX.size();i++)
	{
		for (int j = 0; j < resy; j++)
			VelocityX[i][j] = 0.0f;
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

void FIELD::Divergence_Vitesse_Cell(int x, int y)		//divergence at cell index x and y
{
	LATTICE& cell = vector_field[x][y];
	float gradient_x = (VelocityX[x+1][y] - VelocityX[x][y]) / cell.size;
	float gradient_y = (VelocityY[x][y + 1] - VelocityY[x][y])/cell.size;

	cell.divergence = gradient_x + gradient_y;
}

void FIELD::InFlow(float speed_ms)
{
	for (int i = 1; i < resy - 1;i++)
	{
		VelocityX[0][i] = speed_ms;
	}
}

FIELD::~FIELD()
{
	delete[] VelocityX;
	delete[] VelocityY;
}