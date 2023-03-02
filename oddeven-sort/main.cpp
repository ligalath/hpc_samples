#include <stdio.h>
#include <vector>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[])
{
	
	MPI_Init(&argc, &argv);
	std::vector<float> float_vec = {1.0, 2.0, 5.0, 3.0, 6.0, 4.0};
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("rank: %d, size: %d, input vec:\n", rank, size);
	for (auto i : float_vec)
	{
		printf("%f, ", i);
	}
	printf("\n");
	bool last_block = false, finish =false;
	int offset = 0;
	while (!finish)
	{
		//odd sort
		int local_swap = 0, sum_swap = 0, recv_swap = 0;
		int i = rank * 2 + offset;
		float eles[2];
		if (i < float_vec.size() - 1)
		{
			if (float_vec[i] > float_vec[i + 1])
			{
				float tmp = float_vec[i];
				float_vec[i] = float_vec[i + 1];
				float_vec[i + 1] = tmp;
				local_swap++;
			}
			eles[0] = float_vec[i];
			eles[2] = float_vec[i+1];
		}
		else
			last_block = true;
		
		MPI_Status status;
	printf("rank: %d, size: %d, 1\n", rank, size);
		//send element to last process
		if (rank != 0 && i < float_vec.size() - 1)
			MPI_Send(&eles, 2, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
		//receive element from next process
		if (i < float_vec.size() - 2 )
			MPI_Recv(&eles, 2, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &status);
		//update element
		if (i < float_vec.size() - 2 )
		{
			float_vec[i+2] = eles[0];
			float_vec[i+3] = eles[1];
		}

	printf("rank: %d, size: %d, 2\n", rank, size);
		int unit_finish = 1;
		float even_eles[2];
		//even sort
		i = rank * 2 + 1 + offset;
		if (i < float_vec.size() - 1)
		{
			if (float_vec[i] > float_vec[i + 1])
			{
				float tmp = float_vec[i];
				float_vec[i] = float_vec[i + 1];
				float_vec[i + 1] = tmp;
				local_swap++;
				even_eles[0] = float_vec[i];
				even_eles[1] = float_vec[i+1];
			}
		}
	printf("rank: %d, size: %d, 3\n", rank, size);
		//send element to next process
		if (i < float_vec.size() - 2 )
			MPI_Send(&even_eles, 2, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD);
		//receive element from next process
		if (rank != 0 && i < float_vec.size() - 1)
			MPI_Recv(&even_eles, 2, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status);
		if(i > 2 && i < float_vec.size() -1)
		{
			float_vec[i-2] = even_eles[0];
			float_vec[i-1] = even_eles[1];
		}

		if (offset + size * 2 > float_vec.size())
			offset = 0;
		else
			offset += size * 2;

	printf("rank: %d, size: %d, 4\n", rank, size);
		if (rank != 0)
		{
			printf("rank: %d, size: %d, 4.1\n", rank, size);
			MPI_Send(&local_swap, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
			printf("rank: %d, size: %d, 4.2\n", rank, size);
			MPI_Recv(&recv_swap, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
			printf("rank: %d, size: %d, 4.3\n", rank, size);
			if(recv_swap == 0)
			{
				MPI_Send(&local_swap, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
				MPI_Send()
				break;
			}
		}
		else
		{
			printf("rank: %d, size: %d, 4.1\n", rank, size);
			for(int i =1; i < size ; ++i)
			{
				MPI_Recv(&recv_swap, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
				sum_swap += recv_swap;
			}
			printf("rank: %d, size: %d, 4.2\n", rank, size);
			for(int i =1; i < size; ++i)
			{
				int tmp_sum_swap = sum_swap;
			}
			printf("rank: %d, size: %d, sum_swap: %d \n", rank, size, sum_swap);
			if(sum_swap == 0)
			{
				printf("after sorted: ");
				printf("%d,%d", float_vec[i],float_vec[i+1]);
				break;
			}
		}
	printf("rank: %d, size: %d, 5\n", rank, size);
	}
	printf("rank: %d, size: %d, 6\n", rank, size);
	MPI_Finalize();
	for (auto i : float_vec)
	{
		printf("%f, ", i);
	}
}
