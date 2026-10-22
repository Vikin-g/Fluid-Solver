#include <objects.h>

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

CELL::CELL(int rel_x, int rel_y)
	:rel_x(rel_x), rel_y(rel_y)
{}
/*	int i = 0;
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
*/
VEC::VEC(int Pos, float Velocity)
	:pos(Pos), vel(Velocity)
{ }

//----------------------------- FIELD----------------------------------------------

FIELD::FIELD(int pos1x, int pos1y, int Resx, int Resy, int Size, float rho, float U_p)
	:InPosx(pos1x),InPosy(pos1y), resx(Resx), resy(Resy), size(Size), rho(rho), u_p(U_p)
{ 
	CreateField();
}

void FIELD::CreateField() 
{
	vector_field = new CELL* [resx];

	VelocityX = new float[(resx + 1) * resy];
	VelocityY = new float[resx * (resy + 1)];

	TVelocityX = new VEC[(resx + 1) * resy];
	TVelocityY = new VEC[resx * (resy + 1)];

	for (int y = 0;y < resy;y++)
	{
		for (int x = 0; x < resx; x++)
		{
			vector_field[x] = new CELL [resy];								//intialize field
			vector_field[x][y] = CELL(x, y);

			VelocityY[x + y * resx] = 0.0f;										//initialize array of velocities in-between cells
			VelocityY[x + resx * resy] = 0.0f;									// add extra top value 

			VelocityX[x + y * (resx+1)] = 0.0f;
		}
		VelocityX[resx + y * resx] = 0.0f;										// add extra right value
	}
}

void FIELD::SolvePressureCell(int x, int y)
{
	float v_left =	GetVelocityAtPoint(x + 0, y + 0, DIRECTION::HORIZONTAL);
	float v_right = GetVelocityAtPoint(x + 1, y + 0, DIRECTION::HORIZONTAL);
	float v_down =	GetVelocityAtPoint(x + 0, y + 1, DIRECTION::VERTICAL);
	float v_up =	GetVelocityAtPoint(x + 1, y + 1, DIRECTION::VERTICAL);

	float p_left =	vector_field[x - 1][y + 0].pressure;
	float p_right = vector_field[x + 1][y + 0].pressure;
	float p_down =	vector_field[x + 0][y - 1].pressure;
	float p_up =	vector_field[x + 0][y + 1].pressure;

	float DeltaVelocity = v_up - v_down + v_right - v_left;
	float PressureSum = p_up + p_down + p_left + p_right;
	vector_field[x][y].pressure = ((size*rho/dt)*DeltaVelocity + PressureSum)/4;
}

void FIELD::UpdateAll()
{
	
}

void FIELD::Divergence_Vitesse_Cell(int x, int y)		//divergence at cell index x and y
{
	CELL& cell = vector_field[x][y];
	float gradient_x = (GetVelocityAtPoint(x + 1, y, DIRECTION::HORIZONTAL) - GetVelocityAtPoint(x, y, DIRECTION::HORIZONTAL)) / size;
	float gradient_y = (GetVelocityAtPoint(x, y + 1, DIRECTION::VERTICAL) - GetVelocityAtPoint(x, y, DIRECTION::VERTICAL)) / size;

	cell.divergence = gradient_x + gradient_y;
}

void FIELD::VirtualLagrangian(int x, int y)
{
	float v_left = GetVelocityAtPoint(x + 0, y + 0, DIRECTION::HORIZONTAL);
	float v_right = GetVelocityAtPoint(x + 1, y + 0, DIRECTION::HORIZONTAL);
	float v_down = GetVelocityAtPoint(x + 0, y + 1, DIRECTION::VERTICAL);
	float v_up = GetVelocityAtPoint(x + 1, y + 1, DIRECTION::VERTICAL);

	int index = Index(x, y);

	TVelocityX[index].pos = x;
	TVelocityX[index].vel = v_left;
	TVelocityX[index].prev_x = x - v_left * dt;

	TVelocityY[index].pos = y;
	TVelocityY[index].vel = v_down;
	TVelocityY[index].prev_y = y - v_down * dt;
}

void FIELD::Gauss_Seidel()
{
	float r_rel = 1.0;			//relative residual

	int i = 0;

	while (r_rel > u_p)  // 
	{

		if (i == 50)
			break;
		else if (i % 5 == 0)
		{
			float r_2 = 0.0f;			// residual norm (sqrt sum of r^2)
			float s_2 = 0.0f;			//source term

			for (int x = 0; x < resx; x++)
			{
				for (int y = 0; y < resy; y++)
				{
					r_2 += std::pow(CalcR_p(x, y), 2);
					s_2 += std::pow(CalcS_p(x, y), 2);

					SolvePressureCell(x, y);
				}
			}

			r_rel = std::pow(r_2, 0.5f) / std::pow(s_2, 0.5f);
		}
		else
		{
			for (int x = 0; x < resx;x++)
			{
				for (int y = 0;y < resy;y++)
				{
					SolvePressureCell(x, y);
				}
			}
		}
		i++;
	}
}

float FIELD::CalcR_p(int x, int y)
{
	float p_left = vector_field[x - 1][y + 0].pressure;
	float p_right = vector_field[x + 1][y + 0].pressure;
	float p_down = vector_field[x + 0][y - 1].pressure;
	float p_up = vector_field[x + 0][y + 1].pressure;

	float this_divergence = vector_field[x][y].divergence;
	float this_p = vector_field[x][y].pressure;

	float SumPressures = p_left + p_right + p_down + p_up - 4 * this_p;

	return (SumPressures) / (size * size) - CalcR_p(x,y); //sum of pressures divided by area - sourceterm of particle
}

float FIELD::CalcS_p(int x, int y)
{
	float divergence = vector_field[x][y].divergence;

	return rho * divergence / dt;
}

void FIELD::InFlow(float speed_ms)
{
	for (int i = 1; i < resy - 1;i++)
	{
		GetVelocityAtPoint(0, i, DIRECTION::HORIZONTAL) = speed_ms;
	}
	dt = floor(5 * size / speed_ms);						//time interval calculations forward Euler
}

void FIELD::EnableBorder() 
{
	for (int i = 0; i<resx;i++) 
	{
		vector_field[i][0].is_border = true;
		vector_field[i][resy - 2].is_border = true;
	}
}

float& FIELD::GetVelocityAtPoint(int x, int y, DIRECTION direction)
{
	int index = Index(x,y);

	if (direction == DIRECTION::HORIZONTAL)
		return VelocityX[index];
	else
		return VelocityY[index];
}



void FIELD::ViewVelocity(int x, int y, DIRECTION direction)
{
	std::cout << GetVelocityAtPoint(x, y, direction) << std::endl;
}

int FIELD::Index(int x, int y)
{
	return x + y * resx;
}

FIELD::~FIELD()
{
	delete[] VelocityX;
	delete[] VelocityY;

	delete[] TVelocityX;
	delete[] TVelocityY;

	delete[] NewV_x;
	delete[] NewV_y;

	for (int i = 0; i < resy; i++)
		delete[] vector_field[i];
	delete[] vector_field;
}