// euler_method.cpp : Defines the entry point for the application.
//
#include <objects.h>

#include <filesystem>
#include <iostream>



enum class PLANCHE
{
	BELLY,
	CONCAVE,
	DOUBLE_CONCAVE,
};

int main()
{
	
	std::string path = "../../../resources/Surf1.obj";
	
	FIELD3D field(1, 1, 940, 490, 0.002f, path,1.224f);
	try
	{
		float F = field.SolveAll();       // use your own variable names
		std::cout << F << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << "exception: " << e.what() << std::endl;
	}

	return 0;
}
