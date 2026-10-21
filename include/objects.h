#pragma once

#include <visualize.h>

enum class DIRECTION {
	HORIZONTAL,
	VERTICAL
};

struct CELL
{
	int rel_x;
	int rel_y;

	int size;

	float ro;
	float divergence;
	
	float positions[8];

	bool is_border = false;

	CELL(int size = 1, int x = 1, int y = 1, int posinitx = 1, int posinity = 1, float rho = 1);
};

struct VEC
{
	int pos;

	float prev_x;
	float prev_y;

	float vel;

	VEC(int Pos = 0, float Velocity = 0.0f);
};

class FIELD
{
private:
	int resx;
	int resy;
	int size;

	int InPosx;
	int InPosy;

	VEC* TVelocityX;
	VEC* TVelocityY;

	float* VelocityX;
	float* VelocityY;

	void CreateField(float ro);

	float& GetVelocityAtPoint(int x, int y, DIRECTION direction);	//to the left if x or under if y

	int Index(int x, int y);
public:
	float dt;
	CELL** vector_field;

	FIELD(int pos1x, int pos1y, int Resx, int Resy, int Size, float rho);
	~FIELD();

	void EnableBorder();

	void Divergence_Vitesse_Cell(int x, int y);

	void InFlow(float speed_ms);

	void VirtualLagrangian(int x, int y);					//find value of particle before current position

	void ViewVelocity(int x, int y, DIRECTION direction);
};