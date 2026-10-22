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
	
	float divergence;

	bool is_border = false;

	float pressure;

	CELL(int rel_x = 0, int rel_y = 0);
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

	float rho;

	float u_p;					// uncertainty for pressure

	VEC* TVelocityX;
	VEC* TVelocityY;

	float* VelocityX;
	float* VelocityY;

	float* NewV_x;
	float* NewV_y;

	void CreateField();

	float& GetVelocityAtPoint(int x, int y, DIRECTION direction);	//to the left if x or under if y

	int Index(int x, int y);

	void SolvePressureCell(int x, int y);

	void UpdateVelocities();

	void UpdateAll();

	void Gauss_Seidel();

	
public:
	float dt;
	CELL** vector_field;

	FIELD(int pos1x, int pos1y, int Resx, int Resy, int Size, float rho, float U_p);
	~FIELD();

	void EnableBorder();

	void Divergence_Vitesse_Cell(int x, int y);
	
	float CalcR_p(int x, int y);

	float CalcS_p(int x, int y);
	
	void InFlow(float speed_ms);

	void VirtualLagrangian(int x, int y);					//find value of particle before current position

	void ViewVelocity(int x, int y, DIRECTION direction);
};