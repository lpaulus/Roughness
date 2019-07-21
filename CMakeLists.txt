cmake_minimum_required (VERSION 2.8.11)
project (Roughness)

find_package(CGAL REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(Polyhedron)

file(GLOB SOURCES "*.cpp")

add_executable(Roughness ${SOURCES})

if(CGAL_FOUND)
  #required to use basic_viewer
  #add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
  #create the executable of the application
  #add_executable(test_ test.cpp)
  #link it with the required CGAL libraries
  target_link_libraries(Roughness CGAL::CGAL)
else()
  message("ERROR: this program requires CGAL and will not be compiled.")
endif()
