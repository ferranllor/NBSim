// GPU version
// module add nvhpc/21.2
// nvc -fast -tp=nehalem -ta=tesla:managed -acc -Minfo=accel NBSimGPU.c -o NBGPU
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef float real;
#define SOFTENING_SQUARED  0.01f

// Data structures real3 and real4
typedef struct { real x, y, z; }    real3;
typedef struct { real x, y, z, w; } real4;

typedef struct { real* __restrict x, * __restrict y, * __restrict z; } real3array;
typedef struct { real* __restrict x, * __restrict y, * __restrict z, * __restrict w; } real4array;

real3 bodyBodyInteraction(real4 iPos, real4 jPos)
{
    real rx, ry, rz, distSqr, s;

    rx = jPos.x - iPos.x;  ry = jPos.y - iPos.y;  rz = jPos.z - iPos.z;

    distSqr = rx * rx + ry * ry + rz * rz;

    if (distSqr < SOFTENING_SQUARED) s = jPos.w / powf(SOFTENING_SQUARED, 1.5f);
    else                             s = jPos.w / powf(distSqr, 1.5f);

    real3 f;
    f.x = rx * s;  f.y = ry * s; f.z = rz * s;

    return f;
}


void integrate(real4array out, real4array in,
    real3array vel, real3array force,
    real    dt, int n)
{
    int i, j;
    #pragma acc data present(in.x, in.y, in.z, in.w, out.x, out.y, out.z, out.w, force.x, force.y, force.z, vel.x, vel.y, vel.z)
    {
        #pragma acc parallel loop
        for (i = 0; i < n; i++)
        {
            real fx = 0, fy = 0, fz = 0;

            #pragma acc loop reduction(+:fx,fy,fz)
            for (j = 0; j < i; j++)
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

            #pragma acc loop reduction(+:fx,fy,fz)
            for (j = i + 1; j < n; j++)
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
        }

        #pragma acc parallel loop
        for (i = 0; i < n; i++)
        {
            real fx = force.x[i], fy = force.y[i], fz = force.z[i];
            real px = in.x[i], py = in.y[i], pz = in.z[i], invMass = in.w[i];
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


int main(int argc, char** argv)
{
    int i, n = 20000;
    int iterations = 10;
    real dt = 0.01667;

    if (argc >= 2) n = atoi(argv[1]);
    if (argc >= 3) iterations = atoi(argv[2]);


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

    #pragma acc data copy(pin.x[0:n], pin.y[0:n], pin.z[0:n], pin.w[0:n], pout.x[0:n], pout.y[0:n], pout.z[0:n], pout.w[0:n], f.x[0:n] ,f.y[0:n], f.z[0:n], v.x[0:n], v.y[0:n], v.z[0:n])
    for (i = 0; i < iterations; i++)
    {
        integrate(pout, pin, v, f, dt, n);

        real4array tmp = pin;
        pin = pout;
        pout = tmp;
    }

    real3 p_av = average(pin, n);
    printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);
    printf("Body-0  position: (%f,%f,%f)\n", pin.x[0], pin.y[0], pin.z[0]);

    free(pin.x);  free(pout.x);  free(v.x);  free(f.x);
    free(pin.y);  free(pout.y);  free(v.y);  free(f.y);
    free(pin.z);  free(pout.z);  free(v.z);  free(f.z);
    free(pin.w);  free(pout.w);

    return 0;
}
