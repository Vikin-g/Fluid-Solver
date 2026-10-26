#define NOMINMAX                 // windows.h defines min/max macros that break std::min/std::max
#define WIN32_LEAN_AND_MEAN

#include <objects.h>

#include <filesystem>
#include <iostream>

#include <windows.h>
#include <filesystem>

std::filesystem::path ExeDir()
{
	wchar_t buf[MAX_PATH];
	DWORD n = GetModuleFileNameW(nullptr, buf, MAX_PATH);
	return std::filesystem::path(std::wstring(buf, n)).parent_path();
}

// in main:
std::string path = (ExeDir() / "resources" / "Belly3BigFin.obj").string();



enum class PLANCHE
{
	BELLY,
	CONCAVE,
	DOUBLE_CONCAVE,
};

int main()
{
	std::filesystem::current_path(ExeDir());

	std::string path1 = "../../../resources/BellyUtilisee";
	std::string path2 = "../../../resources/Belly3BigFin";

	int n = 0;

	std::cout << "Quelle planche voulez vous calculer?\nEcrivez 1 pour celle utilisee pour le TM et 2 pour une alternative avec une fin centrale plus longue: ";
	while (n != 1 && n != 2)
	{
		std::cout << "Quelle planche voulez vous calculer?\n"
			"Ecrivez 1 pour celle utilisée pour le TM et 2 pour une alternative avec une fin centrale plus longue: ";
		if (!(std::cin >> n))
		{
			std::cerr << "Vous devez fournir un nombre!!!" << std::endl;
			std::cin.clear();
			std::cin.ignore(10000, '\n');
			n = 0;
		}
	}

	std::cout << std::endl;
	
	const char* files[] = { "BellyUtilisee.obj", "Belly3BigFin.obj" };
	std::string path = (ExeDir() / "resources" / files[n - 1]).string();

	std::cerr << "path:   [" << path << "]\n"
		<< "exists: " << std::filesystem::exists(path) << "\n";
	for (const auto& e : std::filesystem::directory_iterator("../../../resources"))
		std::cerr << "found:  [" << e.path().filename().string() << "]\n";

	FIELD3D field(1, 1, 940, 490, 0.002f, path,1.224f);

	try
	{
		float F = field.SolveAll();       // use your own variable names
		std::cout << F << std::endl;
	}
	catch (const std::exception& e) { std::cerr << "exception: " << e.what() << std::endl; return 1; }
	catch (...) { std::cerr << "unknown exception" << std::endl; return 1; }


	return 0;
}
