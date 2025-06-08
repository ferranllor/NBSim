// Baseline version
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

typedef float real;
#define SOFTENING_SQUARED  0.01f
#define BLOCKSIZE 1024

// Data structures real3 and real4
typedef struct { real x, y, z; }    real3;
typedef struct { real x, y, z, w; } real4;

typedef struct { real* __restrict x, * __restrict y, * __restrict z; } real3array;
typedef struct { real* __restrict x, * __restrict y, * __restrict z, * __restrict w; } real4array;

bool debug = true;
double cp_Wtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

__global__ void integrate(real4array out, real4array in,
    real3array vel, real3array force,
    real    dt, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x, j;

    if (i < n)
    {
        real fx = 0, fy = 0, fz = 0;

        for (j = 0; j < n; j++)
        {
            real rx, ry, rz, distSqr, si;

            rx = in.x[j] - in.x[i];  ry = in.y[j] - in.y[i];  rz = in.z[j] - in.z[i];

            distSqr = rx * rx + ry * ry + rz * rz;

            if (distSqr < SOFTENING_SQUARED)
                si = in.w[j] / powf(SOFTENING_SQUARED, 1.5f);
            else
                si = in.w[j] / powf(distSqr, 1.5f);

            fx += rx * si;  fy += ry * si; fz += rz * si;
        }

        force.x[i] = fx;  force.y[i] = fy;  force.z[i] = fz;

        real px = in.x[i], py = in.y[i], pz = in.z[i], invMass = 1.0f / in.w[i];
        real vx = vel.x[i], vy = vel.y[i], vz = vel.z[i];

        // acceleration = force / mass; 
        // new velocity = old velocity + acceleration * deltaTime
        vx += (fx * invMass) * dt;
        vy += (fy * invMass) * dt;
        vz += (fz * invMass) * dt;

        // new position = old position + velocity * deltaTime
        px += vx * dt;
        py += vy * dt;
        pz += vz * dt;

        out.x[i] = px;
        out.y[i] = py;
        out.z[i] = pz;
        out.w[i] = invMass;

        vel.x[i] = vx;
        vel.y[i] = vy;
        vel.z[i] = vz;
    }

}


real dot(real v0[3], real v1[3])
{
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}


real normalize(real vector[3])
{
    float dist = sqrt(dot(vector, vector));
    if (dist > 1e-6)
    {
        vector[0] /= dist;
        vector[1] /= dist;
        vector[2] /= dist;
    }
    return dist;
}


void cross(real out[3], real v0[3], real v1[3])
{
    out[0] = v0[1] * v1[2] - v0[2] * v1[1];
    out[1] = v0[2] * v1[0] - v0[0] * v1[2];
    out[2] = v0[0] * v1[1] - v0[1] * v1[0];
}


void randomizeBodies(real4array pos,
    real3array vel,
    float clusterScale,
    float velocityScale,
    int   n)
{
    srand(42);
    float scale = clusterScale;
    float vscale = scale * velocityScale;
    float inner = 2.5f * scale;
    float outer = 4.0f * scale;

    int i = 0;

    while (i < n)
    {
        real x, y, z;
        // generate real numbers between -1.0 and +1.0
        x = rand() / (float)RAND_MAX * 2 - 1;
        y = rand() / (float)RAND_MAX * 2 - 1;
        z = rand() / (float)RAND_MAX * 2 - 1;

        real point[3] = { x, y, z };
        real len = normalize(point);
        if (len > 1) // discard position and generate new one
            continue;

        pos.x[i] = point[0] * (inner + (outer - inner) * rand() / (real)RAND_MAX);
        pos.y[i] = point[1] * (inner + (outer - inner) * rand() / (real)RAND_MAX);
        pos.z[i] = point[2] * (inner + (outer - inner) * rand() / (real)RAND_MAX);
        pos.w[i] = 1.0f;

        real axis[3] = { 0.0f, 0.0f, 1.0f };

        if (1 - dot(point, axis) < 1e-6)
        {
            axis[0] = point[1];
            axis[1] = point[0];
            normalize(axis);
        }
        real vv[3] = { (real)pos.x[i], (real)pos.y[i], (real)pos.z[i] };
        real vv0[3];

        cross(vv0, vv, axis);
        vel.x[i] = vv0[0] * vscale;
        vel.y[i] = vv0[1] * vscale;
        vel.z[i] = vv0[2] * vscale;

        i++;
    }
}


real3 average(real4array p, int n)
{
    int i;
    real3 av = { 0.0, 0.0, 0.0 };
    for (i = 0; i < n; i++)
    {
        av.x += p.x[i];
        av.y += p.y[i];
        av.z += p.z[i];
    }
    av.x /= n;
    av.y /= n;
    av.z /= n;
    return av;
}

cudaError_t cudaSim(int n, int iterations, real dt, real4array h_pin, real4array h_pout, real3array h_v, real3array h_f)
{
    cudaError_t cudaStatus;
    dim3 blockDim(BLOCKSIZE);
    dim3 gridDim((n + BLOCKSIZE - 1) / BLOCKSIZE);
    double ini, time, avgt = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    if (debug) {
        ini = cp_Wtime();
        cudaDeviceSynchronize();
        time = cp_Wtime() - ini;

        printf("Time to do initial synchronize: %lf sec\n", time);

        ini = cp_Wtime();
    }

    real4array d_pin;
    d_pin.x = 0;
    d_pin.y = 0;
    d_pin.z = 0;
    d_pin.w = 0;

    cudaStatus = cudaMalloc((void**)&d_pin.x, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_pin.y, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_pin.z, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_pin.w, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    real4array d_pout;
    d_pout.x = 0;
    d_pout.y = 0;
    d_pout.z = 0;
    d_pout.w = 0;

    cudaStatus = cudaMalloc((void**)&d_pout.x, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_pout.y, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_pout.z, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_pout.w, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    real3array d_v;
    d_v.x = 0;
    d_v.y = 0;
    d_v.z = 0;

    cudaStatus = cudaMalloc((void**)&d_v.x, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_v.y, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_v.z, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    real3array d_f;
    d_f.x = 0;
    d_f.y = 0;
    d_f.z = 0;

    cudaStatus = cudaMalloc((void**)&d_f.x, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_f.y, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&d_f.z, n * sizeof(real));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    if (debug) {
        cudaDeviceSynchronize();
        time = cp_Wtime() - ini;
        printf("Time to do cudaMallocs: %lf sec\n", time);
        ini = cp_Wtime();
    }

    // Copy input vectors from host memory to GPU buffers.

    cudaStatus = cudaMemcpy(d_pin.x, h_pin.x, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_pin.y, h_pin.y, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_pin.z, h_pin.z, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_pin.w, h_pin.w, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(d_pout.x, h_pout.x, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_pout.y, h_pout.y, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_pout.z, h_pout.z, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_pout.w, h_pout.w, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(d_v.x, h_v.x, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_v.y, h_v.y, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_v.z, h_v.z, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(d_f.x, h_f.x, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_f.y, h_f.y, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(d_f.z, h_f.z, n * sizeof(real), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    if (debug) {
        cudaDeviceSynchronize();
        time = cp_Wtime() - ini;
        printf("Time to do cudaMemcpys: %lf sec\n", time);
    }

    for (int i = 0; i < iterations; i++)
    {
        if (debug)
            ini = cp_Wtime();

        integrate<<<gridDim, blockDim>>>(d_pout, d_pin, d_v, d_f, dt, n);

        real* tmpx, * tmpy, * tmpz, * tmpw;

        tmpx = d_pout.x;
        tmpy = d_pout.y;
        tmpz = d_pout.z;
        tmpw = d_pout.w;

        d_pout.x = d_pin.x;
        d_pout.y = d_pin.y;
        d_pout.z = d_pin.z;
        d_pout.w = d_pin.w;

        d_pin.x = tmpx;
        d_pin.y = tmpy;
        d_pin.z = tmpz;
        d_pin.w = tmpw;

        if (debug)
        {
            time = cp_Wtime() - ini;
            if (i > 0) {
                printf("Time to do step %d: %lf sec\n", i, time);
                avgt += time;
            }
            else
                printf("Time to launch first kernel: %lf sec\n", time);
        }
    }

    if (debug) {
        printf("Average step time: %lf sec\n", avgt/(iterations-1));
        ini = cp_Wtime();
    }

    cudaStatus = cudaMemcpy(h_pin.x, d_pin.x, n * sizeof(real), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(h_pin.y, d_pin.y, n * sizeof(real), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(h_pin.z, d_pin.z, n * sizeof(real), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
    cudaStatus = cudaMemcpy(h_pin.w, d_pin.w, n * sizeof(real), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    if (debug) {
        time = cp_Wtime() - ini;
        printf("Time to memcpy results to host: %lf sec\n", avgt);
    }

Error:
    cudaFree(d_pin.x);  cudaFree(d_pin.y);  cudaFree(d_pin.z);  cudaFree(d_pin.w);
    cudaFree(d_pout.x); cudaFree(d_pout.y); cudaFree(d_pout.z); cudaFree(d_pout.w);
    cudaFree(d_f.x); cudaFree(d_f.y); cudaFree(d_f.z);
    cudaFree(d_v.x); cudaFree(d_v.y); cudaFree(d_v.z);

    return cudaStatus;
}


int main(int argc, char** argv)
{
    int n = 20000;
    int iterations = 10;
    real dt = 0.001667;
    double ini;

    if (argc >= 2) n = atoi(argv[1]);
    if (argc >= 3) iterations = atoi(argv[2]);

    if (debug) 
        ini = cp_Wtime();


    real4array pin;
    pin.x = (real*)malloc(n * sizeof(real));
    pin.y = (real*)malloc(n * sizeof(real));
    pin.z = (real*)malloc(n * sizeof(real));
    pin.w = (real*)malloc(n * sizeof(real));

    real4array pout;
    pout.x = (real*)malloc(n * sizeof(real));
    pout.y = (real*)malloc(n * sizeof(real));
    pout.z = (real*)malloc(n * sizeof(real));
    pout.w = (real*)malloc(n * sizeof(real));

    real3array v;
    v.x = (real*)malloc(n * sizeof(real));
    v.y = (real*)malloc(n * sizeof(real));
    v.z = (real*)malloc(n * sizeof(real));

    real3array f;
    f.x = (real*)malloc(n * sizeof(real));
    f.y = (real*)malloc(n * sizeof(real));
    f.z = (real*)malloc(n * sizeof(real));

    randomizeBodies(pin, v, 1.54f, 8.0f, n);

    printf("n=%d bodies for %d iterations:\n", n, iterations);

    real3 p_av = average(pin, n);
    printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);
    printf("Body-0  position: (%f,%f,%f)\n", pin.x[0], pin.y[0], pin.z[0]);

    cudaError_t cudaStatus = cudaSim(n, iterations, dt, pin, pout, v, f);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSim failed!");
        return 1;
    }

    p_av = average(pin, n);
    printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);
    printf("Body-0  position: (%f,%f,%f)\n", pin.x[0], pin.y[0], pin.z[0]);

    free(pin.x);  free(pout.x);  free(v.x);  free(f.x);
    free(pin.y);  free(pout.y);  free(v.y);  free(f.y);
    free(pin.z);  free(pout.z);  free(v.z);  free(f.z);
    free(pin.w);  free(pout.w);

    if (debug) {
        double time = cp_Wtime() - ini;
        printf("Total elapsed time: %lf sec\n", time);
    }

    return 0;
}