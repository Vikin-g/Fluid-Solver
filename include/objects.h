#pragma once

#include <visualize.h>

enum class DIRECTION {
	HORIZONTAL,
	VERTICAL
};

enum class Material {
	FREE_WALL,
	WATER,
	OBJECT
};

struct CELL
{
	int rel_x;
	int rel_y;
	
	float divergence;

	Material material;

	float pressure;

	CELL(int rel_x = 0, int rel_y = 0);
};

template<typename T>
void SwapBuffers(T*& front, T*& back)
{
	T* temp = front;
	front = back;
	back = temp;
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

	float in_v;

	float* TVelocityX;
	float* TVelocityY;

	float* VelocityX;
	float* VelocityY;

	void CreateField();

	float& GetVelocityAtPoint(int x, int y, DIRECTION direction);	//to the left if x or under if y

	float& GetPressureAtPoint(int x, int y);

	int Index(int x, int y);

	void SolvePressureCell(int x, int y);

	void UpdateVelocities();

	void Solve();

	void Gauss_Seidel();

	void VirtualLagrangian();					//find value of particle before current position
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

	Material IsBorder(int x, int y);

	

	void ViewVelocity(int x, int y, DIRECTION direction);
};