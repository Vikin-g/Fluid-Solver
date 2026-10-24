// euler_method.cpp : Defines the entry point for the application.
//
#include <objects.h>

enum class PLANCHE
{
	BELLY,
	CONCAVE,
	DOUBLE_CONCAVE,
};

int main()
{
	FIELD3D field(1, 1, 940, 490, 0.02f, "../resources/surf+board.obj",1);
	std::cout << field.SolveAll() << std::endl;
	return 0;
}
