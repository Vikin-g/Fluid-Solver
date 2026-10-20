#include "../include/objects.h"

#define HORIZONTAL true
#define VERTICAL false

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

CELL::CELL(int size, int x, int y, int posinitx, int posinity, float ro)
	:rel_x(x), rel_y(y), size(size), ro(ro)
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
}

FIELD::FIELD(int pos1x, int pos1y, int Resx, int Resy, int Size, float rho)
	:InPosx(pos1x),InPosy(pos1y), resx(Resx), resy(Resy), size(Size)
{ 
	CreateField(rho);
}

void FIELD::CreateField(float ro) 
{
	vector_field = new CELL* [resx];
	VelocityX = new float[(resx + 1) * resy];
	VelocityY = new float[resx * (resy + 1)];

	for (int y = 0;y < resy;y++)
	{
		for (int x = 0; x < resx; x++)
		{
			vector_field[x] = new CELL [resy];								//intialize field
			vector_field[x][y] = CELL(size, x, y, InPosx, InPosy, ro);

			VelocityY[x + y * resx] = 0.0f;										//initialize array of velocities in-between cells
			VelocityY[x + resx * resy] = 0.0f;									// add extra top value 

			VelocityX[x + y * (resx+1)] = 0.0f;
		}
		VelocityX[resx + y * resx] = 0.0f;										// add extra right value
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
	CELL& cell = vector_field[x][y];
	float gradient_x = (GetVelocityAtPoint(x + 1, y, HORIZONTAL) - GetVelocityAtPoint(x, y, HORIZONTAL)) / cell.size;
	float gradient_y = (GetVelocityAtPoint(x, y + 1, VERTICAL) - GetVelocityAtPoint(x, y, VERTICAL)) / cell.size;

	cell.divergence = gradient_x + gradient_y;
}

void FIELD::InFlow(float speed_ms)
{
	for (int i = 1; i < resy - 1;i++)
	{
		ChangeVelocityAtPoint(0, i, HORIZONTAL, speed_ms);
	}
	dt = floor(5 * size / speed_ms);						//time interval calculations forward Euler
}

float FIELD::GetVelocityAtPoint(int x, int y, bool horizontal)
{
	if (horizontal)
		return VelocityX[x + y * resx];
	else
		return VelocityY[x + y * resx];
}

void FIELD::ChangeVelocityAtPoint(int x, int y, bool horizontal, float n_v)
{
	if (horizontal)
		VelocityX[x + y * resx] = n_v;
	else
		VelocityY[x + y * resx] = n_v;
}

FIELD::~FIELD()
{
	delete[] VelocityX;
	delete[] VelocityY;

	for (int i = 0; i < resy; i++)
		delete[] vector_field[i];
	delete[] vector_field;
}