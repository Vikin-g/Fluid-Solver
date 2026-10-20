#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>

struct ECRAN
{
	GLFWwindow* window;
	int LARGEUR = 960;
	int HAUTEUR = 540;

	ECRAN();
};
