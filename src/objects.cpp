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

	TVelocityX = new float[(resx + 1) * resy];
	TVelocityY = new float[resx * (resy + 1)];

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

void FIELD::Solve()
{
	for (int i = 0; i < time; i++)
	{
		Gauss_Seidel();
		UpdateVelocities();
		VirtualLagrangian();
	}
}

float FIELD::IntegrateFriction()			//Integrate friction in x-component
{
	float Fx;
	for (int x = 0;x < resx; x++)
	{
		for (int y = floor(resy/2); y < resy; y++)
		{
			Pressure p = PrepPressure(x, y);
			if (IsBorder(x, y) == Material::OBJECT && p.flow_right == 1)
			{
				float gradient_x = (GetVelocityAtPoint(x + 1, y, DIRECTION::HORIZONTAL) - GetVelocityAtPoint(x, y, DIRECTION::HORIZONTAL)) / size;
				float sigma_x = -GetPressureAtPoint(x + 1, y) + 2 * mu * gradient_x;
				float sigma_y = mu * vector_field[x + 1][y].divergence;

				Fx += (sigma_x * n + sigma_y * n) * size;
			}
		}
	}
}

void FIELD::SolvePressureCell(int x, int y)
{
	if (vector_field[x][y].material == Material::WATER)
	{
		Pressure p = PrepPressure(x, y);

		float v_left = GetVelocityAtPoint(x + 0, y + 0, DIRECTION::HORIZONTAL);
		float v_right = GetVelocityAtPoint(x + 1, y + 0, DIRECTION::HORIZONTAL);
		float v_down = GetVelocityAtPoint(x + 0, y + 1, DIRECTION::VERTICAL);
		float v_up = GetVelocityAtPoint(x + 1, y + 1, DIRECTION::VERTICAL);

		float p_left = vector_field[x - 1][y + 0].pressure * p.flow_left;
		float p_right = vector_field[x + 1][y + 0].pressure * p.flow_right;
		float p_down = vector_field[x + 0][y - 1].pressure * p.flow_down;
		float p_up = vector_field[x + 0][y + 1].pressure * p.flow_up;

		float VelSum = (v_up - v_down + v_right - v_left);
		float PressureSum = p_up + p_down + p_left + p_right;
		GetPressureAtPoint(x, y) = (PressureSum - (size * rho / dt) * VelSum) /p.fluidEdgeCount;
	}
	else
		GetPressureAtPoint(x, y) = 0;
}

Pressure FIELD::PrepPressure(int x, int y)
{
	Pressure p;
	if (IsBorder(x - 1, y + 0) == Material::WATER)
	{
		p.fluidEdgeCount += 1;
		p.flow_left = 1;
	}
	if (IsBorder(x - 1, y - 1) == Material::WATER)
	{
		p.fluidEdgeCount += 1;
		p.flow_down = 1;
	}
	if (IsBorder(x + 1, y + 0) == Material::WATER)
	{
		p.fluidEdgeCount += 1;
		p.flow_right = 1;
	}
	if (IsBorder(x + 0, y + 1) == Material::WATER)
	{
		p.fluidEdgeCount += 1;
		p.flow_up = 1;
	}

	return p;
}



void FIELD::VirtualLagrangian()
{
	for (int x = 1; x < resx;x++)
	{
		for (int y = 1; y < resy;y++)
		{
			float v_left = GetVelocityAtPoint(x + 0, y + 0, DIRECTION::HORIZONTAL);
			float v_right = GetVelocityAtPoint(x + 1, y + 0, DIRECTION::HORIZONTAL);
			float v_down = GetVelocityAtPoint(x + 0, y + 1, DIRECTION::VERTICAL);
			float v_up = GetVelocityAtPoint(x + 1, y + 1, DIRECTION::VERTICAL);

			int index = Index(x, y);

			if (IsBorder(x, y) == Material::WATER)
			{	
				float back_x = x - (v_left * dt);
				float back_y = y - (v_down * dt);

				if (back_x < 0)
				{
					TVelocityX[index] = in_v;
					TVelocityY[index] = 0.0f;
				}
				else
				{
					TVelocityX[index] = v_left;
					TVelocityY[index] = v_down;
				}
			}
			else if (IsBorder(x, y) == Material::OBJECT)
			{
			}
		}
	}
	SwapBuffers(VelocityX, TVelocityX);
	SwapBuffers(VelocityY, TVelocityY);
}

void FIELD::UpdateVelocities()
{
	for (int x = 1; x < resx;x++)
	{
		for (int y = 1; y < resy;y++)
		{	
			if (vector_field[x][y].material == Material::WATER)
			{
				VelocityX[Index(x, y)] -=  dt * (GetPressureAtPoint(x + 1, y) - GetPressureAtPoint(x, y)) / size;
				VelocityY[Index(x, y)] -= dt * (GetPressureAtPoint(x, y + 1) - GetPressureAtPoint(x, y)) / size;
			}
		}
	}
	
}



void FIELD::Gauss_Seidel()
{
	float r_rel = 1.0;			//relative residual

	int i = 0;

	while (r_rel > u_p)
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
	float p_left = GetPressureAtPoint(x - 1,y + 0);
	float p_right = GetPressureAtPoint(x + 1,y + 0);
	float p_down = GetPressureAtPoint(x + 0,y - 1);
	float p_up = GetPressureAtPoint(x + 0,y + 1);

	float this_divergence = vector_field[x][y].divergence;
	float this_p = GetPressureAtPoint(x,y);

	float SumPressures = p_left + p_right + p_down + p_up - 4 * this_p;

	return (SumPressures) / (size * size) - CalcR_p(x,y); //sum of pressures divided by area - sourceterm of particle
}

float FIELD::CalcS_p(int x, int y)
{
	float divergence = vector_field[x][y].divergence;

	return rho * divergence / dt;
}

void FIELD::Divergence_Vitesse_Cell(int x, int y)		//divergence at cell index x and y
{
	CELL& cell = vector_field[x][y];
	float gradient_x = (GetVelocityAtPoint(x + 1, y, DIRECTION::HORIZONTAL) - GetVelocityAtPoint(x, y, DIRECTION::HORIZONTAL)) / size;
	float gradient_y = (GetVelocityAtPoint(x, y + 1, DIRECTION::VERTICAL) - GetVelocityAtPoint(x, y, DIRECTION::VERTICAL)) / size;

	cell.divergence = gradient_x + gradient_y;
}

float& FIELD::GetVelocityAtPoint(int x, int y, DIRECTION direction)
{
	int index = Index(x,y);

	if (direction == DIRECTION::HORIZONTAL)
		return VelocityX[index];
	else
		return VelocityY[index];
}

float& FIELD::GetPressureAtPoint(int x, int y)
{
	return vector_field[x][y].pressure;
}

Material FIELD::IsBorder(int x, int y)
{
	return vector_field[x][y].material;
}
void FIELD::InFlow(float speed_ms)
{
	for (int i = 1; i < resy - 1;i++)
	{
		GetVelocityAtPoint(0, i, DIRECTION::HORIZONTAL) = speed_ms;
	}

	in_v = speed_ms;

	dt = floor(5 * size / speed_ms);						//time interval calculations forward Euler
}

void FIELD::EnableBorder() 
{
	for (int i = 0; i<resx;i++) 
	{
		vector_field[i][0].material = Material::FREE_WALL;
		vector_field[i][resy - 2].material = Material::FREE_WALL;
	}
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

	for (int i = 0; i < resy; i++)
		delete[] vector_field[i];
	delete[] vector_field;
}