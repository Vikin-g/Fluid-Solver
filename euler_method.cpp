// euler_method.cpp : Defines the entry point for the application.
//
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "src/objects.h"

struct ECRAN
{
	GLFWwindow* window;
	int LARGEUR = 960;
	int HAUTEUR = 540;

	ECRAN() 
	{
		if (!glfwInit())
		{
			std::cerr << "Error" << std::endl;
			exit(EXIT_FAILURE);
		}
		window = glfwCreateWindow(LARGEUR, HAUTEUR, "CFD avec methode d'euler", NULL, NULL);

		if (!window)
		{
			glfwTerminate();
			exit(EXIT_FAILURE);
		}

		glfwMakeContextCurrent(window);
		glewExperimental = GL_TRUE;

		if (glewInit() != GLEW_OK) {
			std::cerr << "Failed to initialize GLEW" << std::endl;
			glfwDestroyWindow(window);
			glfwTerminate();
			exit(EXIT_FAILURE);
		}
		glViewport(0, 0, LARGEUR, HAUTEUR);
	}
};

ECRAN ecran;



int main()
{
	FIELD Field(3, 4, 960, 540, 2);

	Field.CreateField();
	Field.EnableBorder();
	Field.InFlow(0.25);

	std::cout << Field.vector_field[0][0].is_border << std::endl;
	std::cout << Field.GetVelocityAtPoint(0, 1, true);
	return 0;
}
