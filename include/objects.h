#pragma once

#include <visualize.h>
#include <board.h>
#include <string>

enum class DIRECTION {
	HORIZONTAL,
	VERTICAL
};

enum class Material {
	FREE_WALL,
	WATER,
	OBJECT
};

struct Pressure 
{
	int flow_up = 0;
	int flow_left = 0;
	int flow_down = 0;
	int flow_right = 0;

	int fluidEdgeCount = 0;
};

struct CELL
{
	int rel_x = 0;
	int rel_y = 0;
	
	float divergence = 0;

	Material material;

	float pressure = 0;
	
	CELL() = default;
	CELL(int rel_x, int rel_y);
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
	float size;

	int InPosx;
	int InPosy;

	float rho;

	float CFL = 0.5;

	int time;
	float dt;

	float u_p;					// uncertainty for pressure

	float in_v;

	float mu = 1e-3f;

	float* TVelocityX;
	float* TVelocityY;

	float* VelocityX;
	float* VelocityY;

	void CreateField();

	float& GetVelocityAtPoint(int x, int y, DIRECTION direction);	//to the left if x or under if y

	float& GetPressureAtPoint(int x, int y);

	int IndexU(int x, int y);
	int IndexV(int x, int y);

	void SolvePressureCell(int x, int y);

	void UpdateVelocities();

	void Gauss_Seidel();

	void VirtualLagrangian();					//find value of particle before current position

	Pressure PrepPressure(int x, int y);

	void Diffuse();

	void Solve3D();

	float BilinearPos(float f[], int w, int height, float x, float y);


public:

	CELL** vector_field;

	FIELD(int pos1x, int pos1y, int Resx, int Resy, float Size, float rho, float U_p);
	~FIELD();

	void Solve();

	void EnableBorder();

	void Divergence_Vitesse_Cell(int x, int y);
	
	float CalcR_p(int x, int y);

	float CalcS_p(int x, int y);
	
	void InFlow(float speed_ms);

	Material IsBorder(int x, int y);

	float IntegrateFriction();

	void ViewVelocity(int x, int y, DIRECTION direction);

	void FindBoard(ObjectMask, int, int);

	float SampleX(float x, float y);
	float SampleY(float x, float y);

	bool IsSolidV(int x, int y);
	bool IsSolidU(int x, int y);
};

//------------------------------------FIELD 3D----------------------------------

class FIELD3D
{
public:

	int intersec;
	int pos1x;
	int pos1y;
	int resx;
	int resy;
	int size;

	float v_in;

	std::string path;

	FIELD3D(int pos1x, int pos1y, int Resx, int Resy, int Size, std::string Path, float v_in);

	float SolveAll();
};