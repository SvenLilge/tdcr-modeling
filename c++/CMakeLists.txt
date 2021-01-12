cmake_minimum_required(VERSION 3.2)

project(tdcr-modeling)

#Important CMAKE stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)  

#include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR})

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(tdcr-modeling main.cpp ${H_FILES} ${CPP_FILES})
target_link_libraries(tdcr-modeling Eigen3::Eigen GSL::gsl)
