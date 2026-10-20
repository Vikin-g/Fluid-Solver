#include <visualize.h>



ECRAN::ECRAN()
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
