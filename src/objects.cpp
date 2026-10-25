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
{
}
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

FIELD::FIELD(int pos1x, int pos1y, int Resx, int Resy, float Size, float rho, float U_p)
	:InPosx(pos1x),InPosy(pos1y), resx(Resx), resy(Resy), size(Size), rho(rho), u_p(U_p)
{ 
	CreateField();
}

void FIELD::CreateField() 
{
	vector_field = new CELL * [resx];
	for (int i = 0; i < resx; i++)
	{
		vector_field[i] = new CELL[resy];
		for (int j = 0; j < resy; j++)
		{
			vector_field[i][j] = CELL(i, j);
		}
	}



	VelocityX = new float[(resx + 1) * resy]();									
	VelocityY = new float[resx * (resy + 1)]();
	
	TVelocityX = new float[(resx + 1) * resy]();
	TVelocityY = new float[resx * (resy + 1)]();
}

void FIELD::Solve()
{
	auto allFinite = [&](int step)
	{
		for (int y = 0; y < resy; y++)
		{
			for (int x = 0; x <= resx; x++)
			{	
				if (!std::isfinite(VelocityX[IndexU(x, y)]))
				{
					std::cerr << "step " << step << ": non-finite u at (" << x << "," << y << ")\n";
					return false;
				}
			}
		}
		for (int y = 0; y <= resy; y++)
		{
			for (int x = 0; x < resx; x++)
			{
				if (!std::isfinite(VelocityY[IndexV(x, y)]))
				{
					std::cerr << "step " << step << ": non-finite v at (" << x << "," << y << ")\n";
					return false;
				}
			}
		}
		for (int x = 0; x < resx; x++)
			for (int y = 0; y < resy; y++)
				if (!std::isfinite(vector_field[x][y].pressure))
				{
					std::cerr << "step " << step << ": non-finite p at (" << x << "," << y << ")\n";
					return false;
				}
		return true;
	};

	for (int i = 0; i < steps; i++)
	{
		VirtualLagrangian();
		Diffuse();
		Gauss_Seidel();
		UpdateVelocities();
		
		if (i % 10 == 0 && !allFinite(i))
			return;
	}

	if (!allFinite(steps))
		return;

	static bool report = true;
	if (report)
	{
		float umax = 0.0f, pmax = 0.0f;
		for (int y = 0; y < resy; y++)
			for (int x = 0; x <= resx; x++)
				umax = std::max(umax, std::fabs(VelocityX[IndexU(x, y)]));
		for (int x = 0; x < resx; x++)
			for (int y = 0; y < resy; y++)
				pmax = std::max(pmax, std::fabs(vector_field[x][y].pressure));
		std::cerr << "max|u|=" << umax << " max|p|=" << pmax << "\n";
		report = false;
	}
}

float FIELD::IntegrateFriction()			//Integrate friction in x-component
{
	float Fx = 0;
	for (int x = 1;x < resx-1; x++)
	{
		for (int y = 1; y < resy-1; y++)
		{
			if (IsBorder(x, y) != Material::OBJECT)
				continue;
			Pressure p = PrepPressure(x, y);
			if (p.flow_left == 1)
			{
				Fx += GetPressureAtPoint(x-1,y) * size;
			}
			if (p.flow_right == 1)
			{
				Fx -= GetPressureAtPoint(x + 1, y) * size;
			}
			if (p.flow_up == 1)
			{
				float u = 0.5f * (VelocityX[IndexU(x, y + 1)] + VelocityX[IndexU(x + 1, y + 1)]);
				Fx += 2.0f * mu * u;
			}
			if (p.flow_down == 1)
			{
				float u = 0.5f * (VelocityX[IndexU(x, y - 1)] + VelocityX[IndexU(x + 1, y - 1)]);
				Fx += 2.0f * mu * u;
			}
		}
	}
	return Fx;
}

void FIELD::SolvePressureCell(int x, int y)
{
	if (vector_field[x][y].material == Material::WATER)
	{
		Pressure p = PrepPressure(x, y);

		float p_left = vector_field[x - 1][y + 0].pressure * p.flow_left;
		float p_right = vector_field[x + 1][y + 0].pressure * p.flow_right;
		float p_down = vector_field[x + 0][y - 1].pressure * p.flow_down;
		float p_up = vector_field[x + 0][y + 1].pressure * p.flow_up;

		float VelSum = vector_field[x][y].divergence*size;
		float PressureSum = p_up + p_down + p_left + p_right;
		if (p.fluidEdgeCount == 0) { GetPressureAtPoint(x, y) = 0; }
		else
		{
			const float p_gs = (PressureSum - (size * rho / dt) * VelSum) / p.fluidEdgeCount;
			float& pc = GetPressureAtPoint(x, y);
			pc += 1.8f * (p_gs - pc);          // omega = 1.8; stay below 2
		}
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
	if (IsBorder(x + 0, y - 1) == Material::WATER)
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

float FIELD::BilinearPos(float f[], int w, int h, float x, float y)
{

	x = std::clamp(x, 0.0f, float(w - 1));			//stay in boundary
	y = std::clamp(y, 0.0f, float(h - 1));

	int i = std::min(int(x), w - 2);
	int j = std::min(int(y), h - 2);

	float left = x - i;
	float down = y - j;

	float a = f[i + (j * w)]	 * (1 - left) + f[i + 1 + j * w] * left;
	float b = f[i + (j + 1) * w] * (1 - left) + f[i + 1 + (j + 1) * w] * left;

	return a * (1 - down) + b * down;
}

float FIELD::SampleX(float x, float y) { return BilinearPos(VelocityX, resx + 1, resy, x, y - 0.5f); }
float FIELD::SampleY(float x, float y) { return BilinearPos(VelocityY, resx, resy + 1, x - 0.5f, y); }

void FIELD::VirtualLagrangian()
{
	std::copy(VelocityX, VelocityX + (resx + 1) * resy, TVelocityX);
	std::copy(VelocityY, VelocityY + resx * (resy + 1), TVelocityY);

	float K = dt / size;
	for (int x = 1; x < resx;x++)
	{
		for (int y = 0; y < resy;y++)
		{

			if (!IsSolidU(x, y))
			{	
				float px = x, py = y + 0.5f;

				float v_left = GetVelocityAtPoint(x + 0, y + 0, DIRECTION::HORIZONTAL);
				float v = SampleY(px, py);

				float back_x = px - K * v_left;
				float back_y =  py - K * v;

				TVelocityX[IndexU(x, y)] = SampleX(back_x, back_y);
			}
			else
			{
				TVelocityX[IndexU(x, y)] = 0.0f;
			}
		}
	}
	for (int x = 0; x < resx; x++)
	{
		for (int y = 1; y < resy;y++)
		{
			if (!IsSolidV(x, y))
			{
				float px = x + 0.5f, py = y;

				float u = SampleX(px, py);
				float v_down = GetVelocityAtPoint(x + 0, y  + 0, DIRECTION::VERTICAL);

				float back_x = px - K * u;
				float back_y = py - K * v_down;

				TVelocityY[IndexV(x, y)] = SampleY(back_x, back_y);
			}
			else
			{
				TVelocityY[IndexV(x, y)] = 0.0f;
			}
		}
	}
	SwapBuffers(VelocityX, TVelocityX);
	SwapBuffers(VelocityY, TVelocityY);
}

void FIELD::Diffuse()
{
	std::copy(VelocityX, VelocityX + (resx + 1) * resy, TVelocityX);
	std::copy(VelocityY, VelocityY + resx * (resy + 1), TVelocityY);

	const float nu = mu / rho;
	const float a = nu * dt / (size * size);		//diffusion term

	//  a <= 0.25 for explicit stability

	for (int y = 1; y < resy - 1; y++)
	{
		for (int x = 1; x < resx; x++)
		{
			if (IsSolidU(x,y))
				continue;                       

			float centr = VelocityX[IndexU(x, y)];
			float left = VelocityX[IndexU(x - 1, y)];
			float right = VelocityX[IndexU(x + 1, y)];
			float down = VelocityX[IndexU(x, y - 1)];
			float up = VelocityX[IndexU(x, y + 1)];

			TVelocityX[IndexU(x, y)] = centr + a * (left + right + down + up - 4.0f * centr);
		}
	}

	for (int y = 1; y < resy; y++)
	{
		for (int x = 1; x < resx - 1; x++)
		{
			if (IsSolidV(x,y))
				continue;

			float centr = VelocityY[IndexV(x, y)];
			float left = VelocityY[IndexV(x - 1, y)];
			float right = VelocityY[IndexV(x + 1, y)];
			float down = VelocityY[IndexV(x, y - 1)];
			float up = VelocityY[IndexV(x, y + 1)];

			TVelocityY[IndexV(x, y)] = centr + a * (left + right + down + up - 4.0f * centr);
		}
	}

	SwapBuffers(VelocityX,TVelocityX);
	SwapBuffers(VelocityY, TVelocityY);
}

void FIELD::UpdateVelocities()
{
	for (int x = 1; x < resx;x++)
	{
		for (int y = 0; y < resy;y++)
		{	
			if (!IsSolidU(x,y))
			{
				VelocityX[IndexU(x, y)] -=  (dt/rho) * (GetPressureAtPoint(x, y) - GetPressureAtPoint(x-1, y)) / size;
			}
		}
	}
	for (int x = 0; x < resx;x++)
	{
		for (int y = 1; y < resy;y++)
		{
			if (!IsSolidV(x, y))
			{
				VelocityY[IndexV(x, y)] -= (dt / rho) * (GetPressureAtPoint(x, y) - GetPressureAtPoint(x, y - 1)) / size;
			}
		}
	}

}



void FIELD::Gauss_Seidel()
{
	float r_rel = 1.0;			//relative residual

	for (int x = 0; x < resx; x++)
	{
		for (int y = 0; y < resy; y++)
		{
			Divergence_Vitesse_Cell(x, y);
		}
	}

	for (int i = 0; i < 50; i++)
	{
		for (int x = 1; x < resx - 1; x++)		// if x or y in boundary: solve pressure neighbour cells don't exist
		{
			for (int y = 1; y < resy - 1; y++)
			{
				SolvePressureCell(x, y);
			}
		}
		
		if (i % 5 == 0)
		{
			double r_2 = 0.0;			// residual norm (sqrt sum of r^2)
			double s_2 = 0.0;			//source term

			for (int x = 1; x < resx-1; x++)
			{
				for (int y = 1; y < resy-1; y++)
				{
					if (IsBorder(x, y) == Material::WATER)
					{
						double r_p = CalcR_p(x, y);
						double s_p = CalcS_p(x, y);

						r_2 += r_p * r_p;
						s_2 += s_p * s_p;
					}
				}
			}
			if (s_2 == 0.0)
			{
				break;
			}
			if (std::sqrt(r_2 / s_2) < u_p)
				break;
		}
	}
}

float FIELD::CalcR_p(int x, int y)
{
	Pressure p = PrepPressure(x, y);

	float sum = vector_field[x - 1][y].pressure * p.flow_left
		+ vector_field[x + 1][y].pressure * p.flow_right
		+ vector_field[x][y - 1].pressure * p.flow_down
		+ vector_field[x][y + 1].pressure * p.flow_up;

	return (sum - p.fluidEdgeCount * GetPressureAtPoint(x,y)) / (size * size) - CalcS_p(x,y); //sum of pressures divided by area - sourceterm of particle
}

float FIELD::CalcS_p(int x, int y)
{
	float divergence = vector_field[x][y].divergence;

	return rho * divergence / dt;
}

void FIELD::Divergence_Vitesse_Cell(int x, int y)		//divergence at cell index x and y
{
	CELL& cell = vector_field[x][y];
	if (cell.material != Material::WATER)
		cell.divergence = 0.0f;
	else
	{
		float gradient_x = (GetVelocityAtPoint(x + 1, y, DIRECTION::HORIZONTAL) - GetVelocityAtPoint(x, y, DIRECTION::HORIZONTAL)) / size;
		float gradient_y = (GetVelocityAtPoint(x, y + 1, DIRECTION::VERTICAL) - GetVelocityAtPoint(x, y, DIRECTION::VERTICAL)) / size;

		cell.divergence = gradient_x + gradient_y;
	}
}

float& FIELD::GetVelocityAtPoint(int x, int y, DIRECTION direction)
{

	return (direction == DIRECTION::HORIZONTAL) ? VelocityX[IndexU(x, y)]
		: VelocityY[IndexV(x, y)];
}

float& FIELD::GetPressureAtPoint(int x, int y)
{
	return vector_field[x][y].pressure;
}

Material FIELD::IsBorder(int x, int y)
{
	return vector_field[x][y].material;
}

bool FIELD::IsSolidU(int x, int y)
{
	bool left = (x > 0) && IsBorder(x - 1, y) != Material::WATER;
	bool right = (x < resx) && IsBorder(x, y) != Material::WATER;
	return right || left;
}

bool FIELD::IsSolidV(int x, int y)
{
	bool left = (y > 0) && IsBorder(x, y - 1) != Material::WATER;
	bool right = (y < resy) && IsBorder(x, y) != Material::WATER;
	return right || left;
}

void FIELD::InFlow(float speed_ms)
{
	if (speed_ms <= 0.0f)
		throw std::invalid_argument("Vitesse Inflow dout être positive");

	for (int x = 0; x <= resx; x++)
	{
		for (int y = 0; y < resy; y++)
		{
			VelocityX[IndexU(x, y)] = IsSolidU(x, y) ? 0.0f : speed_ms;
		}
	}

	in_v = speed_ms;

	dt = CFL * size / speed_ms;						//time interval calculations forward Euler

	const float a = (mu / rho) * dt / (size * size);
	if (a > 0.25f)
		std::cerr << "Avertissement: nombre de diffusion " << a << " > 0.25, Diffuse() est instable!" << std::endl;
}

void FIELD::FindBoard(ObjectMask object, int posx, int posy)
{
	for (int x = 0; x < (int)object.size(); x++)
	{
		for (int y = 0; y < (int)object[x].size(); y++)
		{
			if (!object[x][y])
				continue;

			const int fx = x + posx, fy = y + posy;
			if (fx <= 0 || fx >= resx || fy <= 0 || fy >= resy)
				continue;											// write on outside field

			vector_field[fx][fy].material = Material::OBJECT;
		}
	}
}

void FIELD::EnableBorder() 
{
	for (int i = 0; i<resx;i++) 
	{
		vector_field[i][0].material = Material::FREE_WALL;
		vector_field[i][resy - 1].material = Material::FREE_WALL;
	}
}

void FIELD::ViewVelocity(int x, int y, DIRECTION direction)
{
	std::cout << GetVelocityAtPoint(x, y, direction) << std::endl;
}

int FIELD::IndexV(int x, int y)
{
	return x + y * resx;
}
int FIELD::IndexU(int x, int y)
{
	return x + y * (resx+1);
}

FIELD::~FIELD()
{
	delete[] VelocityX;
	delete[] VelocityY;

	delete[] TVelocityX;
	delete[] TVelocityY;

	for (int i = 0; i < resx; i++)
		delete[] vector_field[i];
	delete[] vector_field;
}
//----------------------------------3D Field (approximately)------------------------------

FIELD3D::FIELD3D(int pos1x, int pos1y, int Resx, int Resy, float Size, std::string Path, float v_in)
	:pos1x(pos1x), pos1y(pos1y), resx(Resx), resy(Resy), size(Size), path(Path), v_in(v_in)
{}

float FIELD3D::SolveAll()
{
	float total_force = 0;

	rapidobj::Result result = rapidobj::ParseFile(path, rapidobj::MaterialLibrary::Ignore());;

	if (result.error)
	{
		std::cerr << "OBJ load failed: " << result.error.code.message()
			<< " (line " << result.error.line_num << ": " << result.error.line << ")\n";
		return 0.0f;
	}

	
	Bounds3D b0 = ComputeBounds(result);
	const int nx0 = 100;
	const double cell0 = b0.width() / nx0;
	const int ny0 = (int)std::ceil(b0.height() / cell0);
	ObjectMask m0 = FindCrossSection(result, b0, 0.5 * (b0.minz + b0.maxz), nx0, ny0, cell0);

	const double tilt = EstimateTilt(m0);          // about -0.5 rad (-30 degrees) for your slice
	std::cerr << "tilt: " << tilt * 180.0 / 3.14159265 << " deg\n";
	RotateAboutZ(result, -tilt);

	// From here on, everything uses the levelled mesh
	Bounds3D b = ComputeBounds(result);
	float to_m = ComputeScale(b, 0.177f);

	if (!ScaleMatches(b, to_m, 0.039f)) { /* axes are not what you think */ }
	
	GridDims obj = ComputeObjectGrid(b, to_m, size);   // e.g. 89 x 20
	int resx = CellsForLength(0.8f, size);
	int resy = CellsForLength(0.3f, size);
	int posx, posy;
	ComputePlacement(resy, obj, 0.2f, size, posx, posy);
	if (!ObjectFits(resx, resy, obj, posx, posy)) 
	{
		std::cerr << "Object (" << obj.nx << "x" << obj.ny << ") does not fit in "
			<< resx << "x" << resy << " at " << posx << "," << posy << "\n";
		return 0.0f;
	}
	
	float dz = 0.005f;
	int slices = ComputeSliceCount(b, to_m, dz);

	ObjectMask mid = FindCrossSection(result, b, SliceZ(b, to_m, dz, slices / 2),
		obj.nx, obj.ny, size / to_m);
	std::cerr << "mid-slice thickness: " << MaxThickness(mid) * size
		<< " m (expected about 0.039)\n";
	//PrintMask(mid);
	
	std::cerr << "v_in=" << v_in << " size=" << size << " to_m=" << to_m
		<< " resx=" << resx << " resy=" << resy
		<< " posx=" << posx << " posy=" << posy << "\n";


	for (int k = 0; k < slices; k++) {

		float z = SliceZ(b, to_m, dz, k);
		ObjectMask m = FindCrossSection(result, b, z, obj.nx, obj.ny, size / to_m);
		if (!AnyTrue(m) || MaxThickness(m) < 3) continue;
		FIELD field(0, 0, resx, resy, size, rho, 1e-3f);
		field.EnableBorder();
		field.FindBoard(m, posx, posy);
		field.InFlow(v_in);
		field.Solve();
		float t_force = field.IntegrateFriction() * dz;  // N/m * m = N
		std::cerr << "slice " << k << "/" << slices << "  F=" << t_force << std::endl;
		//std::cout << "Slice " << k << ": " << t_force << std::endl;
		total_force += t_force;
	}
	return total_force;
}

