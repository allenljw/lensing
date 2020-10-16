/* 
   159735 Parallel Programming

   Startup program for sequential implementation of simulation by ray
   tracing of gravitational lensing.
 */
#include <ctime>

#include <iostream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <cuda.h>

#include "lenses.h"
#include "arrayff.hxx"


// Global variables! Not nice style, but we'll get away with it here.

// Boundaries in physical units on the lens plane
const float WL  = 2.0;
const float XL1 = -WL;
const float XL2 =  WL;
const float YL1 = -WL;
const float YL2 =  WL;


// Implement lens equation, given the lens position (xl, yl) and the
// lens system configuration, shoot a ray back to the source position
// (xs, ys)
__device__ void shoot(float& xs, float& ys, float xl, float yl, 
  float* xlens, float* ylens, float* eps, int nlenses)
{
  float dx, dy, dr;
  xs = xl;
  ys = yl;
  for (int p = 0; p < nlenses; ++p) {
    dx = xl - xlens[p];
    dy = yl - ylens[p];
    dr = dx * dx + dy * dy;
    xs -= eps[p] * dx / dr;
    ys -= eps[p] * dy / dr;
  }

}


// Used to time code. OK for single threaded programs but not for
// multithreaded programs. See other demos for hints at timing CUDA
// code.
double diffclock(clock_t clock1,clock_t clock2)
{
  double diffticks = clock1 - clock2;
  double diffms = (diffticks * 1000) / CLOCKS_PER_SEC;
  return diffms; // Time difference in milliseconds
}

__global__ void mx_shoot(float* xlens, float* ylens, float* eps, float* d_lensim, float XL1, float YL1, int nlenses, float lens_scale, int cols, int rows) 
{
    // Source star parameters. You can adjust these if you like - it is
    // interesting to look at the different lens images that result
    const float rsrc = 0.1;      // radius
    const float ldc = 0.5;      // limb darkening coefficient
    const float xsrc = 0.0;      // x and y centre on the map
    const float ysrc = 0.0;
    float rsrc2 = rsrc * rsrc;

    float xl, yl, xs, ys, sep2, mu;
    float xd, yd;

    int n = blockDim.x * blockIdx.x + threadIdx.x;

    if (n < (cols * rows)) {
      int ix = n % cols;
      int iy = (n - ix) / cols;

      yl = YL1 + iy * lens_scale;
      xl = XL1 + ix * lens_scale;

      shoot(xs, ys, xl, yl, xlens, ylens, eps, nlenses);
      
      xd = xs - xsrc;
      yd = ys - ysrc;
      sep2 = xd * xd + yd * yd;
      if (sep2 < rsrc2) {
          mu = sqrt(1 - sep2 / rsrc2);
          //lensim(row, col) = 1.0 - ldc * (1 - mu);
          d_lensim[n] = 1.0 - ldc * (1 - mu);
      }
    }

}

int main(int argc, char* argv[]) 
{
  // Set up lensing system configuration - call example_1, _2, _3 or
  // _n as you wish. The positions and mass fractions of the lenses
  // are stored in these arrays
  float* xlens;
  float* ylens;
  float* eps;
  const int nlenses = set_example_3(&xlens, &ylens, &eps);
  std::cout << "# Simulating " << nlenses << " lens system" << std::endl;

  // Pixel size in physical units of the lens image. You can try finer
  // lens scale which will result in larger images (and take more
  // time).
  const float lens_scale = 0.005;

  // Size of the lens image
  const int npixx = static_cast<int>(floor((XL2 - XL1) / lens_scale)) + 1;
  const int npixy = static_cast<int>(floor((YL2 - YL1) / lens_scale)) + 1;
  std::cout << "# Building " << npixx << "X" << npixy << " lens image" << std::endl;

  // Put the lens image in this array
  Array<float, 2> lensim(npixy, npixx);

  //declare the variables for device function here
  // copy the host variables to device variables
  float *d_xlens, *d_ylens, *d_eps, *d_lensim;
  size_t size = nlenses * sizeof(float);
  size_t size_img = npixx * npixy * sizeof(float);
  size_t pitch;
  //float* h_lensim = (float*)malloc(size_img);

  cudaMalloc(&d_xlens, size);
  cudaMalloc(&d_ylens, size);
  cudaMalloc(&d_eps, size);
  //cudaMalloc(&d_lensim, size_img);

  cudaMemcpy(d_xlens, xlens, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_ylens, ylens, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_eps, eps, size, cudaMemcpyHostToDevice);

  int total = npixx * npixy;
  int threadsPerBlock = 512;
  int blocksPerGrid = (total + threadsPerBlock - 1) / threadsPerBlock;

  clock_t tstart = clock();

  mx_shoot<<<blocksPerGrid, threadsPerBlock>>>(d_xlens, d_ylens, d_eps, d_lensim, XL1, YL1, nlenses, lens_scale, npixx, npixy);

  clock_t tend = clock();
  double tms = diffclock(tend, tstart);
  std::cout << "# Time elapsed: " << tms << " ms " << std::endl;

  cudaMemcpy(lensim.buffer, d_lensim, size_img, cudaMemcpyDeviceToHost);

  // Draw the lensing image map here. For each pixel, shoot a ray back
  // to the source plane, then test whether or or not it hits the
  // source star
  // replace with device function
  /*const float rsrc2 = rsrc * rsrc;
  float xl, yl, xs, ys, sep2, mu;
  float xd, yd;
  int numuse = 0;
  for (int iy = 0; iy < npixy; ++iy) 
  for (int ix = 0; ix < npixx; ++ix) { 

    yl = YL1 + iy * lens_scale;
    xl = XL1 + ix * lens_scale;

    shoot(xs, ys, xl, yl, xlens, ylens, eps, nlenses);

    xd = xs - xsrc;
    yd = ys - ysrc;
    sep2 = xd * xd + yd * yd;
    if (sep2 < rsrc2) {
      mu = sqrt(1 - sep2 / rsrc2);
      lensim(iy, ix) = 1.0 - ldc * (1 - mu);
    }
  }*/

  // for (int iy = 0; iy < npixy; ++iy) 
  // for (int ix = 0; ix < npixx; ++ix) { 
  //   lensim(iy, ix) = h_lensim[iy * npixy + ix];
  // }

  // Write the lens image to a FITS formatted file. You can view this
  // image file using ds9
  dump_array<float, 2>(lensim, "lens.fit");

  delete[] xlens;
  delete[] ylens;
  delete[] eps;

  cudaFree( d_xlens );
  cudaFree( d_ylens );
  cudaFree( d_eps );
  cudaFree( d_lensim );
  
}

