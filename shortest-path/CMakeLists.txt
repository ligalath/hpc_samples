cmake_minimum_required(VERSION 3.2)

project(mandelbrotSet)

#find mpi
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
#${MPI_LIBRARIES}

#find opencv
set(OpenCV_DIR /usr/local/opencv/share/OpenCV)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
message(${OpenCV_INCLUDE_DIRS})

# compile hybrid 
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
add_executable(floyd floyd.cpp)
target_link_libraries(floyd ${OpenCV_LIBS} ${MPI_LIBRARIES})
