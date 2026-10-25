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
	
	std::string path1 = "../../../resources/BellyUtilisée";
	std::string path2 = "../../../resources/Belly3BigFin";

	int n = 0;

	std::string path;

	std::cout << "Quelle planche voulez vous calculer?/nEcrivez 1 pour celle utilisée pour le TM et 2 pour une alternative avec une fin centrale plus longue: ";
	while (n == 0)
	{
		if (!(std::cin >> n))
		{
			std::cerr << "Vous devez fournir un nombre!!!"<<std::endl;
			std::cin.clear();
			std::cin.ignore(10000, '\n');

			n = 0;
		}
	}

	std::cout << std::endl;

	if (n == 1)
	{
		path = path1;
	}
	else
		path = path2;
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
