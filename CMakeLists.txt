# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project(Approx_MVBB)


# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5 )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for target: main.out
# ############################

add_executable( main main.cpp Grid.cpp PointsBox.cpp)
if(NOT CGAL_Qt5_FOUND)
	message(STATUS "CGAL_Qt5 not found")
	return()
endif()

# Link the executable to CGAL and third-party libraries
target_link_libraries(main PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Basic_viewer)