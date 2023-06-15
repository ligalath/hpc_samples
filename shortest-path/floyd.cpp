#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include <math.h>
#include "mpi.h"

typedef struct {
    int begin;
    int destination;
    int weight;
}triple_s;

/**
 * @brief 
 * 
 * @param argc 
 * @param argv 
 * in_file:
 *    [type]        [decimal value] [description]
   32 bit integer         3           # vertices
   32 bit integer         6           # edges
   32 bit integer         0           Src id for edge 0
   32 bit integer         1           Dst id for edge 0
   32 bit integer         3           Weight on edge 0
   32 bit integer                     Src id for edge 1
   …                      …           …
   32 bit integer                     Weight on edge 5
 * @return int 
 */
int main(int argc, char* argv[])
{
    MPI_Status status;
    MPI_Init(&argc, &argv);
	int rank, size;
    //get current rank and size
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    std::string in_file_path(argv[0]);
    std::string out_file_path(argv[1]);
    //parallelly read file
    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, in_file_path.c_str(), MPI_MODE_RDONLY,MPI_INFO_NULL, &mpi_file);
    MPI_File_set_view(mpi_file, 0, MPI_INT,MPI_INT, "internal", MPI_INFO_NULL);
    int num_of_vertices = 0;
    MPI_File_read_all(mpi_file, &num_of_vertices, 1, MPI_INT, &status);
    int num_of_edges = 0;
    MPI_File_read_all(mpi_file, &num_of_edges, 1, MPI_INT, &status);
    //read 3 int(src id, dst id, weight) in one read
    const int read_unit = 3;
    //read area partition
    int read_length = num_of_edges*read_unit;
    int offset = num_of_edges/size * rank * read_unit;
    //read partition of edges
    MPI_File_set_view(mpi_file, 2, MPI_INT, MPI_INT, "internal", MPI_INFO_NULL);
    int file_content[read_length];
    MPI_File_read_all(mpi_file, &file_content, read_length, MPI_INT, &status);
    //parse file content
    int distance[num_of_vertices][num_of_vertices];
    memset(distance, INT64_MAX, num_of_vertices*num_of_vertices*sizeof(int));
    for(int i= 0; i < read_length; i+=3)
    {
        int src_id = file_content[i];
        int dst_id = file_content[i+1];
        int weight = file_content[i+2];
        distance[src_id][dst_id] = weight;
    }
    //read roi, block partition
    int num_of_vertices_to_process = num_of_vertices/size;
    int roi_row_start = num_of_vertices_to_process * rank;
    if(size-1 == rank)
        num_of_vertices_to_process = num_of_vertices - num_of_vertices_to_process*(size-1);
    int roi_row_end = roi_src_start + num_of_vertices_to_process;
    //execute floyd algorithym algorithym
    for (int k = 0; k < num_of_vertices; k++)
    {
        //if Kth row in current node
        int from_k_to[num_of_vertices];
        if(k>=roi_row_start && k < roi_row_end)
        {
            memcpy(from_k_to,&(distance[roi_row_start][0]), num_of_vertices);
        }
        MPI_Bcast(from_k_to,)
        for (int i = roi_row_start; i < roi_row_end; i++)
        {
            for (int j = 0; j < num_of_vertices; j++)
            {
                if(distance[i][j] > distance[i][k] + distance[k][j])
                    distance[i][j] = distance[i][k] + distance[k][j];
            }
            
        }
        //inform other procs
    }
    //write result to file
    

    
    MPI_Finalize();
}