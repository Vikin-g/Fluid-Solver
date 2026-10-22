// euler_method.cpp : Defines the entry point for the application.
//
#include <objects.h>

int main()
{
	FIELD Field(3, 4, 960, 540, 2, 1, 1e-3);
	Field.EnableBorder();
	Field.InFlow(0.25);

	std::cout << Field.vector_field[0][0].is_border << std::endl;
	Field.ViewVelocity(0, 1, DIRECTION::HORIZONTAL);
	return 0;
}
