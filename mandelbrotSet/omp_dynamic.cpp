#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "helper.h"
//set max iteration num is 300 
#define MAX_ITE_NUM (300)
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
    if(argc < 9)
        return 0;
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
    unsigned char *data_ptr = img.data;
    const double scale_real = (right - left) / w;
    const double scale_imaginary = (upper - lower) / h;

    #pragma openmp parallel for collapse(2), schedule(dynamic)
    for(int row = 0; row < w; ++row)
    {
        for(int col = 0; col < h; ++col)
        {
            double a = left + row * scale_real;
            double b = lower + col * scale_imaginary;
            int ite = 0;
            IsMandelbrot(a, b, &ite);
            unsigned char pixel_value = ite > 300 ? 200 + (ite - 300) / (MAX_ITE_NUM - 300) * 40.0 : ite / 300.0 * 200;
            data_ptr[row * (int)w * 3 + col * 3 + 0] = pixel_value;
            data_ptr[row * (int)w * 3 + col * 3 + 1] = pixel_value << 4;
            data_ptr[row * (int)w * 3 + col * 3 + 2] = pixel_value << 8;
        }
    }
    cv::imwrite(output_path, img);
    return 0;
}