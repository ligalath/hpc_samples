#include "mpi.h"
#include <stdlib.h>
#include <opencv2/opencv.hpp>

#include "helper.h"
//set max iteration num is 300 
#define MAX_ITE_NUM 300
#define OCCUPY_FLAG 0
#define TASK_DATA_FLAG 1
#define RESULT_DATA_FLAG 2
#define H_STEP 10

typedef struct 
{
    int w_start;
    int w_end;
    int h_start;
    int h_end;
    double left;
    double lower;
    double scale_real;
    double scale_imaginary;
}task_param_s;


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
    //define MPI custom type
    int block_length[8] = {1,1,1,1,1,1,1,1};
    MPI_Aint displacements[8] = {offsetof(task_param_s, w_start), offsetof(task_param_s, w_end),
                                 offsetof(task_param_s, h_start), offsetof(task_param_s, h_end),
                                 offsetof(task_param_s, left), offsetof(task_param_s, lower),
                                 offsetof(task_param_s, scale_real), offsetof(task_param_s, scale_imaginary)};
    MPI_Datatype types[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_task_param_s;
    int ret = MPI_Type_create_struct(8, block_length, displacements, types, &mpi_task_param_s);
    MPI_Type_commit(&mpi_task_param_s);

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
        unsigned char* data_ptr = img.data;
        /// @brief divide tasks into 2*num of procs
        /// othersize num of procs,
        /// otherwise biggest num of procs we could use
        const int total_task_num_w = w/(size*2)>1?(size*2):(w/size>1?size:w);
        const int total_task_num_h = h/(size*2)>1?(size*2):(h/size>1?size:h);
        // const int total_task_num = total_task_num_w > total_task_num_h?total_task_num_h:total_task_num_w;
        const int total_task_num = ceil((float)h/H_STEP);

        LOGD("total_task_num: %d", total_task_num);
        const int cell_w = w;
        const int cell_h = H_STEP;
        const double scale_real = (right - left)/w;
        const double scale_imaginary = (upper-lower)/h;
        int task_allocated = 0;
        int process_task_num = 0;
        // recv proc signal asynchronously
        MPI_Request requests[size - 1];
        for (int i = 1; i < size; i++)
        {
            char unoccupied = 0;
            MPI_Irecv(&unoccupied, 1, MPI_CHAR, i, OCCUPY_FLAG, MPI_COMM_WORLD, &requests[i - 1]);
        }
        // send tasks
        while(true)
        {
            //check whether any one of Irecv requests is done
            int indexs[size -1] = {0};
            MPI_Status status[size-1] = {0};
            int out_count = 0;
            MPI_Waitsome(size - 1, requests, &out_count, indexs,status);
            char unoccupied = 0;
            MPI_Request res_requests[out_count];
            int ite_array[out_count][int(w)*H_STEP+ 4] = {0};
            for(int i = 0; i < out_count; ++i)
            {
                // send tasks
                int w_start = 0;
                int w_end = w;
                int h_start = 0 + task_allocated * cell_h;
                int h_end = h_start + cell_h;
                h_end = h_end > h ? h : h_end;
                // condition to stop allocate tasks 
                if (task_allocated >= total_task_num)
                    break;
                task_param_s task_param = {
                    0, w, h_start, h_end, left, lower, scale_real, scale_imaginary};
                MPI_Request request;
                //send task
                int index = indexs[i];
                MPI_Isend(&task_param, 1, mpi_task_param_s, index+1, TASK_DATA_FLAG, MPI_COMM_WORLD, &request);
                task_allocated++;
                // recv result data
                // MPI_Irecv(ite_array[i], 4+(w_end - w_start) * H_STEP, MPI_INT, index+1, RESULT_DATA_FLAG, MPI_COMM_WORLD, &res_requests[i]);
                MPI_Status status;
                MPI_Recv(ite_array[i], 4+(w_end - w_start) * H_STEP, MPI_INT, index+1, RESULT_DATA_FLAG, MPI_COMM_WORLD, &status);
                LOGD("recv result data from rank: %d", index+1);
                // recv occupy flag
                MPI_Irecv(&unoccupied, 1, MPI_CHAR, index+1, OCCUPY_FLAG, MPI_COMM_WORLD, &requests[i]);
            }
            // process result
            // MPI_Status res_statuses[out_count];
            // MPI_Waitall(out_count, res_requests, res_statuses);
            int *matrix[H_STEP];
            for(int i = 0; i<out_count; ++i)
            {
                int w_start = ite_array[i][0];
                int w_end = ite_array[i][1];
                int h_start = ite_array[i][2];
                int h_end = ite_array[i][3];
                LOGE("rank: 0,process result w_start:%d, w_end:%d, h_start:%d, h_end:%d", w_start, w_end, h_start, h_end);
                int arr_end = ite_array[i][4+(H_STEP)*(w_end-w_start)];
                LOGD("arr_end: %d", arr_end);
                #pragma omp parallel for schedule(dynamic) collapse(2)
                for (int col = h_start; col < h_end; ++col)
                {
                    for (int row = w_start; row < w_end; ++row)
                    {
                        // convert 1d array from 5th element to 2d array
                        matrix[col - h_start] = &(ite_array[i][4 + (col - h_start) * (w_end - w_start)]);
                        // LOGD("process result col:%d, row: %d", col, row);
                        int ite = matrix[col - h_start][row - w_start];
                        unsigned char pixel_value = ite > 300 ? 200 + (ite - 300) / (MAX_ITE_NUM - 300) * 40.0 : ite / 300.0 * 200;
                        data_ptr[row * (int)w * 3 + col * 3 + 0] = pixel_value;
                        data_ptr[row * (int)w * 3 + col * 3 + 1] = pixel_value << 4;
                        data_ptr[row * (int)w * 3 + col * 3 + 2] = pixel_value << 8;
                    }
                }
                process_task_num++;
                if(process_task_num == total_task_num)
                    break;
            }
            //terminate all processes
            if(process_task_num == total_task_num)
            {
                for(int i = 1; i<size; ++i)
                {
                    task_param_s task_param;
                    task_param.w_start = 0;
                    task_param.w_end = 0;
                    MPI_Request request;
                    MPI_Isend(&task_param, 1, mpi_task_param_s, i, TASK_DATA_FLAG, MPI_COMM_WORLD, &request);
                    MPI_Status terminate_status;
                    MPI_Wait(&request, &terminate_status);
                }
                cv::imwrite(output_path, img);
                LOGD("terminate all processes");
                return 0;
            }
            
        }

    }
    else
    {
        while(true)
        {
            // send unoccupied notice
            char occupied = 0;
            MPI_Request occupy_request;
            MPI_Isend(&occupied, 1, MPI_CHAR, 0, OCCUPY_FLAG, MPI_COMM_WORLD,&occupy_request);
            LOGD("rank: %d send unoccupyed", rank);
            // recv task
            MPI_Request task_request;
            task_param_s task_param;
            MPI_Irecv(&task_param, 1, mpi_task_param_s, 0, TASK_DATA_FLAG, MPI_COMM_WORLD, &task_request);
            LOGD("rank %d, recv task", rank);
            MPI_Status task_status;
            MPI_Wait(&task_request, &task_status);
            //terminate condition
            if(task_param.w_start == task_param.w_end)
            {
                LOGD("rank %d, task terminate", rank);
                break;
            }
            LOGD("rank %d, process task", rank);
            //unserialize task param
            int w_start = task_param.w_start;
            int w_end = task_param.w_end;
            int h_start = task_param.h_start;
            int h_end = task_param.h_end;
            double left = task_param.left;
            double lower =  task_param.lower;
            double scale_imaginary = task_param.scale_imaginary;
            double scale_real = task_param.scale_real;
            int ite_array[(w_end-w_start)*H_STEP + 4] = {0};
            ite_array[0] = w_start;
            ite_array[1] = w_end;
            ite_array[2] = h_start;
            ite_array[3] = h_end;
            LOGD("rank: %d receive task w_start:%d, w_end:%d,h_start:%d, h_end:%d, scale_real: %f, scale_imaginary:%f", rank,w_start, w_end, h_start, h_end, scale_real, scale_imaginary);
            int* matrix[h_end-h_start];
            #pragma omp parallel for schedule(dynamic) collapse(2)
            for (int j = h_start; j < h_end; ++j)
            {
                for (int i = w_start; i < w_end; ++i)
                {
                    // convert 1d array from 5th element to 2d array
                    matrix[j - h_start] = &ite_array[4 + (w_end - w_start) * (j - h_start)];
                    double a = left + i * scale_real;
                    double b = lower + j * scale_imaginary;
                    int ite = 0;
                    IsMandelbrot(a, b, &ite);
                    matrix[j-h_start][i-w_start] = ite;
                }
            }
            // send data back
            // MPI_Request send_back_request;
            // MPI_Isend(ite_array, 4+(w_end - w_start)*H_STEP, MPI_INT, 0, RESULT_DATA_FLAG, MPI_COMM_WORLD, &send_back_request);
            MPI_Send(ite_array, 4+(w_end - w_start)*H_STEP, MPI_INT, 0, RESULT_DATA_FLAG, MPI_COMM_WORLD);
            LOGD("rank:%d send back res, w_start:%d, w_end:%d, h_start:%d, h_end:%d, ite end: %d", rank, ite_array[0], ite_array[1], ite_array[2], ite_array[3], ite_array[4+(H_STEP)*(w_end-w_start)]);
        }
    }
    MPI_Finalize();
    return 0;
}