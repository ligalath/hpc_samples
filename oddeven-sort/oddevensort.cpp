#include "mpi.h"
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <sstream>

void merge(float* block1, int block1_len, float* block2, int block2_len, float* block3)
{
    int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::ostringstream oss;
    oss << "recv block: ";
    for(int i = 0; i < block1_len; ++i)
    {
        oss << block1[i] << ",";
    }
    oss << "cur block: ";
    for(int i = 0; i < block2_len; ++i)
    {
        oss << block2[i] << ",";
    }
    printf("rank: %d, merging: %s\n", rank, oss.str().c_str());
    int block1_pos = 0, block2_pos = 0, block3_pos = 0;
    while(block1_pos < block1_len && block2_pos < block2_len)
    {
        if(block1[block1_pos] > block2[block2_pos])
        {
            block3[block3_pos++] = block2[block2_pos++];
        }
        else
        {
            block3[block3_pos++] = block1[block1_pos++];
        }
    }
    if(block1_pos == block1_len)
    {
        while(block2_pos < block2_len)
            block3[block3_pos++] = block2[block2_pos++];
    }
    else
    {
        while(block1_pos < block1_len)
            block3[block3_pos++] = block1[block1_pos++];
    }

    for(int i = 0; i < block1_len; ++i)
    {
        block1[i] = block3[i];
    }
    for(int i = 0; i < block2_len; ++i)
    {
        block2[i] = block3[i + block1_len];
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    std::vector<float> array = {10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0,2.0,1.0};
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(0 == rank)
    {
	    printf("rank: %d, size: %d, input vec:\n", rank, size);
        for (auto i : array)
        {
            printf("%f, ", i);
        }
        printf("\n");
    }
    int sort_count = 0;
    //initialize block
    int block_len = array.size()/size + 1;
    int offset = rank * block_len;
    int cur_block_len = block_len;
    if(offset + block_len > array.size())
        cur_block_len = array.size() - offset;
    //sort block
    auto begin_ite = array.begin() + offset;
    auto end_ite = offset + block_len > array.size() ? array.end() : begin_ite + block_len;
    std::sort(begin_ite, end_ite);
    //create block
    float block[cur_block_len];
    for (int i = 0; i < cur_block_len; ++i)
    {
        block[i] = array[offset + i];
    }

    while(sort_count < size)
    {
        // if current node is out of range
        if (offset >= array.size())
            return 0;
        MPI_Status status;
        MPI_Request request[4];
        // exchange data
        bool is_odd = rank % 2 == 0 ? true : false, is_edge = end_ite == array.end() ? true : false;
        // odd sort
        if (sort_count % 2 == 0)
        {
            if (is_odd && !is_edge)
            {
                MPI_Send(&block, block_len, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                float recv_block[block_len];
                MPI_Irecv(&recv_block, block_len, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, &request[0]);
                MPI_Wait(&request[0], MPI_STATUS_IGNORE);
                // update data
                for(int i = 0; i < cur_block_len; ++i)
                {
                    block[i] = recv_block[i];
                }
            }
            else if (!is_odd)
            {
                float recv_block[block_len];
                MPI_Irecv(&recv_block, block_len, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &request[1]);
                MPI_Wait(&request[1], MPI_STATUS_IGNORE);
                float res_block[block_len + cur_block_len];
                merge(recv_block, block_len, block, cur_block_len, res_block);
                MPI_Send(&res_block ,block_len, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);
            }
        }
        // even sort
        else
        {
            if (!is_odd && !is_edge)
            {
                MPI_Send(&block, block_len, MPI_FLOAT, rank + 1, 2, MPI_COMM_WORLD);
                float recv_block[block_len];
                MPI_Irecv(&recv_block, block_len, MPI_FLOAT, rank + 1, 3, MPI_COMM_WORLD, &request[2]);
                MPI_Wait(&request[2], MPI_STATUS_IGNORE);
                // update data
                for(int i = 0; i < cur_block_len; ++i)
                {
                    block[i] = recv_block[i];
                }
            }
            else if (is_odd && rank != 0)
            {
                float recv_block[block_len];
                MPI_Irecv(&recv_block, block_len, MPI_FLOAT, rank - 1, 2, MPI_COMM_WORLD, &request[3]);
                MPI_Wait(&request[3], MPI_STATUS_IGNORE);
                float res_block[block_len + cur_block_len];
                merge(recv_block, block_len, block, cur_block_len, res_block);
                MPI_Send(&res_block, block_len, MPI_FLOAT, rank - 1, 3, MPI_COMM_WORLD);
            }
        }
        sort_count++;
    }
    std::ostringstream oss;
    for(int i = 0; i < cur_block_len;  ++i)
    {
        oss << block[i] << ",";
    }
	printf("rank: %d,  finalize with cur block: %s\n", rank, oss.str().c_str());
    
    MPI_Finalize();
}