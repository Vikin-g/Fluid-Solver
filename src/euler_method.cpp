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
	
	std::string path = "../../../resources/surf_board.obj";
	/*std::cerr << "Working dir: " << std::filesystem::current_path() << "\n";
	std::cerr << "Path:        " << path << "\n";
	std::cerr << "Exists:      " << std::filesystem::exists(path) << "\n";
	*/

	FIELD3D field(1, 1, 940, 490, 0.002f, path,1);
	std::cout << field.SolveAll() << std::endl;

	return 0;
}
