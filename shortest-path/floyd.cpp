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
 * @note mpiexec -n 4 ./program in_file_path out_file_path
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
    int offset = num_of_edges/size * rank * read_unit;
    int read_length = (num_of_edges/size + num_of_edges%size*(rank==size-1?1:0))*read_unit;
    //read partition of edges
    MPI_File_set_view(mpi_file, 2+offset, MPI_INT, MPI_INT, "internal", MPI_INFO_NULL);
    int file_content[read_length];
    MPI_File_read_all(mpi_file, &file_content, read_length, MPI_INT, &status);
    MPI_File_close(&mpi_file);
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
    int roi_row_end = roi_row_start + num_of_vertices_to_process;
    //execute floyd algorithym algorithym
    for (int k = 0; k < num_of_vertices; k++)
    {
        //if Kth row in current node
        int from_k_to[num_of_vertices];
        const int owner = k/(num_of_vertices/size);
        if(owner == rank)
        {
            memcpy(from_k_to,&(distance[roi_row_start][0]), num_of_vertices);
        }
        MPI_Bcast(from_k_to,num_of_vertices,MPI_INT, owner,MPI_COMM_WORLD);
        //use openmp
        #pragma openmp parallel for collapse(2), schedule(dynamic)
        for (int i = roi_row_start; i < roi_row_end; i++)
        {
            for (int j = 0; j < num_of_vertices; j++)
            {
                if(distance[i][j] > distance[i][k] + from_k_to[j])
                    distance[i][j] = distance[i][k] + from_k_to[j];
            }
            
        }
    }
    MPI_File mpi_out_file;
    MPI_File_open(MPI_COMM_WORLD, out_file_path.c_str(), MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &mpi_out_file);
    std::vector<int> out_data;
    //write result to file
    for (int i = roi_row_start; i < roi_row_end; i++)
    {
        for (int j = 0; j < num_of_vertices; j++)
        {
            if(distance[i][j] < INT64_MAX)
            {
                out_data.push_back(i);
                out_data.push_back(j);
                out_data.push_back(distance[i][j]);
            }
        }
    }
    int write_len = out_data.size();
    int distance_num[size -1] = {0};
    MPI_Gather(&write_len, 1, MPI_INT, &distance_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_File_seek_shared(mpi_out_file, sizeof(int), MPI_SEEK_SET);
    MPI_File_write_ordered(mpi_out_file, out_data.data(), write_len, MPI_INT, &status);
    if(rank == 0)
    {
        int distance_num_sum = 0;
        for(auto i:distance_num)
            distance_num_sum += i;
        MPI_File_write_at(mpi_out_file, 0, &distance_num, 1, MPI_INT, &status);
    }
    MPI_File_close(&mpi_out_file);

    MPI_Finalize();
}