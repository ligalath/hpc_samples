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

## compile mpi static
# add_executable(mandelbrot_mpi_static mpi_static.cpp)
# target_link_libraries(mandelbrot_mpi_static ${MPI_LIBRARIES} ${OpenCV_LIBS})
add_executable(mandelbrot_mpi_dynamic mpi_dynamic.cpp)
target_link_libraries(mandelbrot_mpi_dynamic ${MPI_LIBRARIES} ${OpenCV_LIBS})