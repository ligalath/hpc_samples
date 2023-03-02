#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include "mpi.h"
int main(int argc, char* argv[])
{
    
    int a = 1;
	printf("input vec: %d", a++);
	
    MPI_Status status;
    char message[100];
    MPI_Init(&argc, &argv);
    int numprocs, myid, source;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (myid != 0) {  //非0号进程发送消息
        strcpy(message, "Hello World!");
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,
            MPI_COMM_WORLD);
    }
    else {   // myid == 0，即0号进程接收消息
        for (source = 1; source < numprocs; source++) {
            MPI_Recv(message, 100, MPI_CHAR, source, 99,
                MPI_COMM_WORLD, &status);
            printf("接收到第%d号进程发送的消息：%s\n", source, message);
        }
    }
    MPI_Finalize();
    return 0;
} /* end main */