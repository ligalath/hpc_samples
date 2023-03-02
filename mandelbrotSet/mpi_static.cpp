//./exe $procs $left $right $lower $upper $w $h $out_path

#include "mpi.h"
#include <stdlib.h>
#include <opencv2/opencv.hpp>
//set max iteration num is 300 
#define MAX_ITE_NUM 300

//check whether plural = a + bi is in mandelbrot set
bool IsMandelbrot(double a, double b, int* ite)
{
    double init_a = a;
    double init_b = b;
    for(int i = 0; i< MAX_ITE_NUM; ++i)
    {
        if(a*a + b*b <= 4.0)
        {
            double tmp_a = a*a - b*b + init_a;
            double tmp_b = 2*a*b + init_b; 
            a = tmp_a; 
            b = tmp_b;
        }
        else
        {
            *ite = i;
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    MPI_Status status;
    MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(rank == 0)
    {
        if (argc < 9)
            return 1;
        // read command params
        int procs = atoi(argv[1]);
        double left = atof(argv[2]);
        double right = atof(argv[3]);
        double lower = atof(argv[4]);
        double upper = atof(argv[5]);
        double w = atof(argv[6]);
        double h = atof(argv[7]);
        std::string output_path(argv[8]);
        cv::Mat img(h, w, CV_8UC3);
        //allocate tasks statically
        double w_step = (right - left)/w;
        printf("w_step:%.30f", w_step);
        double h_step = (upper - lower)/h;
        double offset_h[size-1] = {0};
        int cell_h = h/(size-1);
        for(int i = 1; i<size; ++i)
        {
            offset_h[i] = cell_h * (i);
        }
        cv::Mat output_img(h, w, CV_8UC3);
        for(int i = 1; i < size - 1; ++i)
        {
            double params[7] = {offset_h[i-1], w_step, h_step, cell_h, left, lower, w};
            MPI_Send(params, 7, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        double params[7] = {offset_h[size -2], w_step, h_step, h - cell_h * (size - 2), left, lower, w};
        MPI_Send(params, 7, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD);
        printf("mpi proc: %d, line: %d, left: %f, right: %f, lower:%f\n",rank, __LINE__, left, right, lower);
        //recv img data
        MPI_Status status;
        for(int i = 1; i < size -1; ++i)
        {
            unsigned char data_ptr[cell_h*(int)w*3];
            MPI_Recv(data_ptr, cell_h*w*3, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
            int offset = w*3*offset_h[i-1];
            printf("mpi proc: %d, line: %d\n", rank, __LINE__);
            memcpy(output_img.data + offset, data_ptr, cell_h * w * 3);
            printf("mpi proc: %d, line: %d\n", rank, __LINE__);
        }
        unsigned char data_ptr[((int)h-cell_h *(size -2))*(int)w*3];
        printf("mpi proc: %d, line: %d\n", rank, __LINE__);
        MPI_Recv(data_ptr, ((int)h-cell_h*(size-2)) * (int)w * 3, MPI_CHAR, size-1, 1, MPI_COMM_WORLD, &status);
        int offset = w * 3 * offset_h[size-2];
        printf("mpi proc: %d, line: %d\n", rank, __LINE__);
        memcpy(output_img.data + offset, data_ptr, ((int)h-cell_h*(size-2))* w * 3);
        printf("mpi proc: %d, line: %d\n", rank, __LINE__);
        cv::imwrite(output_path, output_img);
    }
    else
    {
        MPI_Status status;
        double params[7] = {0};
        //recv info from process 0
		MPI_Recv(params, 7, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        double offset_h = params[0];
        double w_step = params[1];
        double h_step = params[2];
        int cell_h = (int)params[3];
        int left = (int)params[4];
        int lower = (int)params[5];
        int w = (int)params[6];
        printf("mpi proc: %d, line: %d offset_h: %f,w_step: %f, h_step: %f, cell_h: %d, w: %d\n",
            rank, __LINE__, offset_h, w_step, h_step, cell_h, w);
        unsigned char* data_ptr = new unsigned char[cell_h*w*3];
        printf("mpi proc: %d, line: %d\n", rank, __LINE__);
        //construct plural
        for(int col = 0; col<(int)w; ++col)
        {
            for(int row = 0; row<(int)cell_h; ++row)
            {
                double a = left + (double)(col)*w_step;
                double b = lower + h_step*(offset_h+(double)row);
                int i = -1;
                int pixel_value = 0;
                if(!IsMandelbrot(a, b, &i))
                {
                    pixel_value = i>300?200+(i-300)/(MAX_ITE_NUM-300)*40.0: i/300.0*200;
                }
                data_ptr[row*w*3 + col*3+0] = pixel_value;
                data_ptr[row*w*3 + col*3+1] = pixel_value<<4;
                data_ptr[row*w*3 + col*3+2] = pixel_value<<8;
            }
        }
        printf("mpi proc: %d, line: %d offset_h: %f,w_step: %.30f, h_step: %f, cell_h: %d\n",
            rank, __LINE__, offset_h, w_step, h_step, cell_h);
        MPI_Send(data_ptr, cell_h*w*3, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }
    
    
    MPI_Finalize();
}