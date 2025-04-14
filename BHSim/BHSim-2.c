// BASELINE VERSION
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef float real;
#define SOFTENING_SQUARED  0.01

// Data structures real3 and real4
typedef struct { real x, y, z; }    real3;
typedef struct { real x, y, z, w; } real4;

typedef struct octreeNode octreeNode;

struct octreeNode { real x, y, z, w; octreeNode** p; unsigned int id; };

void insert(octreeNode* root, octreeNode node, real size, real3 relPos)
{
    _Bool trobat = 0;
    octreeNode* act = root;
    real x = node.x, y = node.y, z = node.z;

    while (!trobat)
    {
        if (act->p == NULL)
        {
            act->p = (octreeNode**)malloc(8 * sizeof(octreeNode*));
            for (int i = 0; i < 8; i++)
                act->p[i] = NULL;

            octreeNode tmp = (octreeNode){ act->x, act->y, act->z, act->w, NULL, act->id };
            act->id = 0;
            act->x = 0;
            act->y = 0;
            act->z = 0;
            act->w = 0;

            insert(act, tmp, size, relPos);

        }

        size /= 2;
        int quad;

        if (x > relPos.x)
        {
            relPos.x += size;
            if (y > relPos.y)
            {
                relPos.y += size;
                if (z > relPos.z)
                {
                    relPos.z += size;
                    quad = 0;
                }
                else
                {
                    relPos.z -= size;
                    quad = 1;
                }
            }
            else
            {
                relPos.y -= size;
                if (z > relPos.z)
                {
                    relPos.z += size;
                    quad = 2;
                }
                else
                {
                    relPos.z -= size;
                    quad = 3;
                }
            }
        }
        else
        {
            relPos.x -= size;
            if (y > relPos.y)
            {
                relPos.y += size;
                if (z > relPos.z)
                {
                    relPos.z += size;
                    quad = 4;
                }
                else
                {
                    relPos.z -= size;
                    quad = 5;
                }
            }
            else
            {
                relPos.y -= size;
                if (z > relPos.z)
                {
                    relPos.z += size;
                    quad = 6;
                }
                else
                {
                    relPos.z -= size;
                    quad = 7;
                }
            }
        }

        act->x = act->x * act->w;
        act->y = act->y * act->w;
        act->z = act->z * act->w;

        act->x += node.x * node.w;
        act->y += node.y * node.w;
        act->z += node.z * node.w;

        act->x /= (act->w + node.w);
        act->y /= (act->w + node.w);
        act->z /= (act->w + node.w);

        act->w = (act->w + node.w);

        if (act->p[quad] == NULL)
        {
            act->p[quad] = (octreeNode*)malloc(sizeof(octreeNode));
            act->p[quad]->x = node.x;
            act->p[quad]->y = node.y;
            act->p[quad]->z = node.z;
            act->p[quad]->w = node.w;
            act->p[quad]->p = NULL;
            act->p[quad]->id = node.id;

            return;
        }
        else
        {
            act = act->p[quad];
        }
    }
}


real buildOctree(real4* in, octreeNode* root, int n)
{
    real maxX = fabs(in[0].x), maxY = fabs(in[0].y), maxZ = fabs(in[0].z);

    for (int i = 1; i < n; i++)
    {
        real x = fabs(in[i].x), y = fabs(in[i].y), z = fabs(in[i].z);
        if (x > maxX)
            maxX = x;
        if (y > maxY)
            maxY = y;
        if (z > maxZ)
            maxZ = z;
    }

    int X = log2(maxX);
    int Y = log2(maxY);
    int Z = log2(maxZ);

    real size = X;

    if (size < Y)
    {
        if (Y < Z)
            size = Z;
        else
            size = Y;
    }
    else if (size < Z)
    {
        if (Z < Y)
            size = Y;
        else
            size = Z;
    }

    size = powf(2, size);

    root->x = in[0].x;
    root->y = in[0].y;
    root->z = in[0].z;
    root->w = in[0].w;
    root->id = 1;
    root->p = NULL;

    for (int i = 1; i < n; i++)
        insert(root, (octreeNode) { in[i].x, in[i].y, in[i].z, in[i].w, NULL, i + 1 }, size, (real3) { 0, 0, 0 });

    return size;
}

void freeOctree(octreeNode* node)
{
    if (node->p != NULL)
    {
        for (int i = 0; i < 8; i++)
            if (node->p[i] != NULL)
                freeOctree(node->p[i]);
        free(node->p);
    }

    free(node);

    return;
}

real3 bodyBodyInteraction(real4 iPos, real4 jPos)
{
    real rx, ry, rz, distSqr, s;

    rx = jPos.x - iPos.x;  ry = jPos.y - iPos.y;  rz = jPos.z - iPos.z;

    distSqr = rx * rx + ry * ry + rz * rz;

    if (distSqr < SOFTENING_SQUARED) s = jPos.w / powf(SOFTENING_SQUARED, 1.5);
    else                             s = jPos.w / powf(distSqr, 1.5);

    real3 f;
    f.x = rx * s;  f.y = ry * s; f.z = rz * s;

    return f;
}

real3 integrateNode(real4 node, octreeNode* root, unsigned int n, unsigned int id, real theta, real size)
{
    int sp = 1;
    // theta = d/r

    octreeNode** stack = (octreeNode**)malloc(n * sizeof(octreeNode*));
    stack[0] = root;

    real* stackP = (real*)malloc(n * sizeof(real));
    stackP[0] = size;

    real fx = 0, fy = 0, fz = 0;

    while (sp > 0)
    {
        sp = sp - 1;
        octreeNode* act = stack[sp];

        real rx, ry, rz, distSqr;

        rx = act->x - node.x;  ry = act->y - node.y;  rz = act->z - node.z;

        distSqr = rx * rx + ry * ry + rz * rz;

        if (act->id > 0)
        {
            real s;

            if (distSqr < SOFTENING_SQUARED) s = act->w / powf(SOFTENING_SQUARED, 1.5);
            else                             s = act->w / powf(distSqr, 1.5);

            fx += rx * s;  fy += ry * s; fz += rz * s;
        }
        else
        {
            real p = stackP[sp];

            real dist = sqrtf(distSqr);

            if (p / dist < theta)
            {
                real s;

                if (distSqr < SOFTENING_SQUARED) s = act->w / powf(SOFTENING_SQUARED, 1.5);
                else                             s = act->w / powf(distSqr, 1.5);

                fx += rx * s;  fy += ry * s; fz += rz * s;
            }
            else if (act->id != id)
            {
                for (int child = 0; child < 8; child++)
                {
                    if (act->p[child] != NULL)
                    {
                        stack[sp] = act->p[child];
                        stackP[sp] = p / 2;
                        sp++;
                    }
                }
            }
        }

    }

    free(stack);
    free(stackP);

    return (real3) { fx, fy, fz };
}


void integrateOctree(real4* in, real4* out, octreeNode* root, real3* force, real3* vel, real dt, unsigned int n, real theta, real size)
{
    int i;

    for (i = 0; i < n; i++)
    {
        real3 f = integrateNode(in[i], root, n, i + 1, theta, size);

        force[i].x = f.x;  force[i].y = f.y;  force[i].z = f.z;

        real fx = force[i].x, fy = force[i].y, fz = force[i].z;
        real px = in[i].x, py = in[i].y, pz = in[i].z, invMass = in[i].w;
        real vx = vel[i].x, vy = vel[i].y, vz = vel[i].z;

        // acceleration = force / mass; 
        // new velocity = old velocity + acceleration * deltaTime
        vx += (fx * invMass) * dt;
        vy += (fy * invMass) * dt;
        vz += (fz * invMass) * dt;

        // new position = old position + velocity * deltaTime
        px += vx * dt;
        py += vy * dt;
        pz += vz * dt;

        out[i].x = px;
        out[i].y = py;
        out[i].z = pz;
        out[i].w = invMass;

        vel[i].x = vx;
        vel[i].y = vy;
        vel[i].z = vz;
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


void randomizeBodies(real4* pos,
    real3* vel,
    float clusterScale,
    float velocityScale,
    int   n)
{
    srand(42);
    float scale = clusterScale;
    float vscale = scale * velocityScale;
    float inner = 2.5f * scale;
    float outer = 4.0f * scale;

    int p = 0, v = 0;
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

        pos[i].x = point[0] * (inner + (outer - inner) * rand() / (real)RAND_MAX);
        pos[i].y = point[1] * (inner + (outer - inner) * rand() / (real)RAND_MAX);
        pos[i].z = point[2] * (inner + (outer - inner) * rand() / (real)RAND_MAX);
        pos[i].w = 1.0f;

        real axis[3] = { 0.0f, 0.0f, 1.0f };

        if (1 - dot(point, axis) < 1e-6)
        {
            axis[0] = point[1];
            axis[1] = point[0];
            normalize(axis);
        }
        real vv[3] = { (real)pos[i].x, (real)pos[i].y, (real)pos[i].z };
        real vv0[3];

        cross(vv0, vv, axis);
        vel[i].x = vv0[0] * vscale;
        vel[i].y = vv0[1] * vscale;
        vel[i].z = vv0[2] * vscale;

        i++;
    }
}


real3 average(real4* p, int n)
{
    int i;
    real3 av = { 0.0, 0.0, 0.0 };
    for (i = 0; i < n; i++)
    {
        av.x += p[i].x;
        av.y += p[i].y;
        av.z += p[i].z;
    }
    av.x /= n;
    av.y /= n;
    av.z /= n;
    return av;
}


int main(int argc, char** argv)
{
    int i, j, n = 1000000;
    int iterations = 10;
    real dt = 0.001667;
    real theta = 0.5;

    if (argc >= 2) theta = atof(argv[1]);
    if (argc >= 3) n = atoi(argv[2]);
    if (argc >= 4) iterations = atoi(argv[3]);

    real4* pin = (real4*)malloc(n * sizeof(real4));
    real4* pout = (real4*)malloc(n * sizeof(real4));
    real4* tmp;
    real3* v = (real3*)malloc(n * sizeof(real3));
    real3* f = (real3*)malloc(n * sizeof(real3));

    randomizeBodies(pin, v, 1.54f, 8.0f, n);

    printf("n=%d bodies for %d iterations:\n", n, iterations);

    for (i = 0; i < iterations; i++)
    {
        octreeNode* root = (octreeNode*)malloc(sizeof(octreeNode));
        real size = buildOctree(pin, root, n);
        integrateOctree(pin, pout, root, f, v, dt, n, theta, size);
        tmp = pout;
        pout = pin;
        pin = tmp;
        freeOctree(root);
    }

    real3 p_av = average(pin, n);
    printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);
    printf("Body-0  position: (%f,%f,%f)\n", pin[0].x, pin[0].y, pin[0].z);

    free(pin);  free(pout);  free(v);  free(f);

    return 0;
}