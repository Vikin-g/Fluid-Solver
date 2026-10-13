# CMakeList.txt : CMake project for euler_method, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("euler_method")


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#-------Dependencies----------
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)

add_executable(euler_method 
	"euler_method.cpp"
)
target_include_directories(euler_method PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set(DEPS glfw GLEW::GLEW OpenGL::GL)
target_link_libraries(euler_method PRIVATE ${DEPS})

# TODO: Add tests and install targets if needed.
