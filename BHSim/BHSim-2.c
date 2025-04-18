#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
//#include <GL/glew.h>
//#include <GLFW/glfw3.h>

typedef float real;
#define SOFTENING_SQUARED  0.01f

// Data structures real3 and real4
typedef struct { real* __restrict x, * __restrict y, * __restrict z; }    real3array;
typedef struct { real* __restrict x, * __restrict y, * __restrict z, * __restrict w; } real4array;

typedef struct { real x, y, z; }    real3;
typedef struct { real x, y, z, w; } real4;

typedef struct {
    real x, y, z, w;    // Pos i massa
    int children[8];    // Pos dels fills
    unsigned id;        // ID del node, li sumem 8 i el fem servir tamb per el nº de fills
} octreeNode;

typedef struct {
    octreeNode* __restrict array;   // Array amb els nodes
    int nBodies;                    // Array amb els nodes
    int nNodes;                     // nº de nodes (cossos i centres de masses)
    int arraySize;                  // Mida del array de nodes
    real spaceSize;                 // Mida del espai des de root
    real despl;                     // Desplaçament respecte el origen (0,0,0), si el sistema es desplaça gaire el corregirem restant a tots els cossos aquest, podem obtenir la pos correcta un altre cop sumant aizo a tots als cosses
} octreeArray;

typedef struct {
    octreeArray* octree;            // Arbre octal
    real4array bodies;              // Arbre octal
    real3array force;               // Força dels nodes
    real3array vel;                 // Velocitat dels nodes
    int nBodies;                    // Nº de cossos
    real theta;                     // Valor per determinar si fem la aprox amb un centre de masses o no
    real dt;                        // Pas de temps
} BHSim;

void insert(octreeArray* octree, octreeNode node, int pos, real3 relPos, real size)
{
    octreeNode act = octree->array[pos];
    real x = node.x, y = node.y, z = node.z;

    while (true)
    {
        if (act.id > 8)
        {
            octree->array[pos].id = 0;
            insert(octree, act, pos, relPos, size);
            act = octree->array[pos];
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

        octreeNode* tmp = &octree->array[pos];
        tmp->x = tmp->x * tmp->w;
        tmp->y = tmp->y * tmp->w;
        tmp->z = tmp->z * tmp->w;

        tmp->x += node.x * node.w;
        tmp->y += node.y * node.w;
        tmp->z += node.z * node.w;

        tmp->x /= (tmp->w + node.w);
        tmp->y /= (tmp->w + node.w);
        tmp->z /= (tmp->w + node.w);

        tmp->w = (tmp->w + node.w);

        if (act.children[quad] == 0)
        {
            octree->array[pos].id++;
            octree->array[pos].children[quad] = octree->nNodes;
            octree->array[octree->nNodes].x = node.x;
            octree->array[octree->nNodes].y = node.y;
            octree->array[octree->nNodes].z = node.z;
            octree->array[octree->nNodes].w = node.w;
            octree->array[octree->nNodes].id = node.id;

            for (int i = 0; i < 8; i++)
                octree->array[octree->nNodes].children[i] = 0;

            octree->nNodes++;

            return;
        }
        else
        {
            pos = act.children[quad];
            act = octree->array[pos];
        }
    }
}

real4 propagateMassTask(octreeNode* array, int pos, int prof)
{
    if (array[pos].id > 8)
        return (real4) { array[pos].w * array[pos].x, array[pos].w * array[pos].y, array[pos].w * array[pos].z, array[pos].w };
    else if (prof < 5)
    {
        array[pos].x = 0;
        array[pos].y = 0;
        array[pos].z = 0;
        array[pos].w = 0;

        for (int i = 0; i < 8; i++)
        {
            if (array[pos].children[i] > 0)
            {
                #pragma omp task
                {
                    real4 res = propagateMassTask(array, array[pos].children[i], prof + 1);

                    #pragma omp critical
                    {
                        array[pos].x += res.x;
                        array[pos].y += res.y;
                        array[pos].z += res.z;
                        array[pos].w += res.w;
                    }
                }
            }
        }

        #pragma omp taskwait

        array[pos].x /= array[pos].w;
        array[pos].y /= array[pos].w;
        array[pos].z /= array[pos].w;
    }
    else
    {
        array[pos].x = 0;
        array[pos].y = 0;
        array[pos].z = 0;
        array[pos].w = 0;

        for (int i = 0; i < 8; i++)
        {
            if (array[pos].children[i] > 0)
            {
                real4 res = propagateMassTask(array, array[pos].children[i], prof + 1);

                array[pos].x += res.x;
                array[pos].y += res.y;
                array[pos].z += res.z;
                array[pos].w += res.w;
            }
        }

        array[pos].x /= array[pos].w;
        array[pos].y /= array[pos].w;
        array[pos].z /= array[pos].w;
    }

    return (real4) { array[pos].w* array[pos].x, array[pos].w* array[pos].y, array[pos].w* array[pos].z, array[pos].w };
}

void propagateMass(octreeArray* octree)
{
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            propagateMassTask(octree->array, 0, 0);
        }
    }
}

void buildOctreeArray(BHSim* simulation)
{
    octreeArray* octree = simulation->octree;

    real maxX = fabs(simulation->bodies.x[0]), maxY = fabs(simulation->bodies.y[0]), maxZ = fabs(simulation->bodies.z[0]);

    for (int i = 1; i < octree->nBodies; i++)
    {
        real x = fabs(simulation->bodies.x[i]), y = fabs(simulation->bodies.y[i]), z = fabs(simulation->bodies.z[i]);
        if (x > maxX)
            maxX = x;
        if (y > maxY)
            maxY = y;
        if (z > maxZ)
            maxZ = z;
    }

    int X = ceil(log2(maxX));
    int Y = ceil(log2(maxY));
    int Z = ceil(log2(maxZ));

    int size = X;

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

    size = pow(2, size);
    octree->spaceSize = size;
    octree->nNodes = 1;
    octree->array[0] = (octreeNode){ simulation->bodies.x[0], simulation->bodies.y[0], simulation->bodies.z[0], simulation->bodies.w[0], { 0, 0, 0, 0, 0, 0, 0, 0 }, 9 };

    for (int i = 1; i < octree->nBodies; i++)
        insert(octree, (octreeNode) { simulation->bodies.x[i], simulation->bodies.y[i], simulation->bodies.z[i], simulation->bodies.w[i], { 0, 0, 0, 0, 0, 0, 0, 0 }, i + 9 }, 0, (real3) { 0, 0, 0 }, octree->spaceSize);

    //propagateMass(octree);
}

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

real3 integrateArrayNode(real4 node, octreeArray* octree, octreeNode* stackNodes, real* stackP, real theta, real size, int id)
{
    stackNodes[0] = octree->array[0];
    stackP[0] = size;

    real fx = 0, fy = 0, fz = 0;

    // real blockDist = p / (theta * blocksize)

    int sp = 1;

    while (sp > 0)
    {
        sp = sp - 1;
        octreeNode act = stackNodes[sp];

        if (act.id != id)
        {
            if (act.id > 8)
            {
                real3 ff = bodyBodyInteraction(node, (real4) { act.x, act.y, act.z, act.w });
                fx += ff.x;  fy += ff.y; fz += ff.z;
            }
            else
            {
                real p = stackP[sp];

                real rx = act.x - node.x, ry = act.y - node.y, rz = act.z - node.z;

                real dist = sqrtf(rx * rx + ry * ry + rz * rz);

                if (p / dist < theta)
                {
                    real3 ff = bodyBodyInteraction(node, (real4) { act.x, act.y, act.z, act.w });
                    fx += ff.x;  fy += ff.y; fz += ff.z;
                }
                else
                {
                    for (int child = 0; child < 8; child++)
                    {
                        if (act.children[child] > 0)
                        {
                            stackNodes[sp] = octree->array[act.children[child]];
                            stackP[sp] = p / 2;
                            sp++;
                        }
                    }
                }
            }
        }
    }

    return (real3) { fx, fy, fz };
}

void integrateOctreeArray(BHSim* simulation)
{
    #pragma omp parallel
    {
        octreeNode* stack = (octreeNode*)malloc(simulation->nBodies * sizeof(octreeNode));
        real* stackP = (real*)malloc(simulation->nBodies * sizeof(real));

        int i;
        
        #pragma omp for
        for (i = 0; i < simulation->nBodies; i++)
        {
            real4 body = (real4){ simulation->bodies.x[i], simulation->bodies.y[i], simulation->bodies.z[i], simulation->bodies.w[i] };
            real3 f = integrateArrayNode(body, simulation->octree, stack, stackP, simulation->theta, simulation->octree->spaceSize, i + 9);

            real fx = f.x, fy = f.y, fz = f.z;
            real px = simulation->bodies.x[i], py = simulation->bodies.y[i], pz = simulation->bodies.z[i], invMass = simulation->bodies.w[i];
            real vx = simulation->vel.x[i], vy = simulation->vel.y[i], vz = simulation->vel.z[i];

            // EULER STEP
            // acceleration = force / mass; 
            // new velocity = old velocity + acceleration * deltaTime
            vx += (fx * invMass) * simulation->dt;
            vy += (fy * invMass) * simulation->dt;
            vz += (fz * invMass) * simulation->dt;

            // new position = old position + velocity * deltaTime
            px += vx * simulation->dt;
            py += vy * simulation->dt;
            pz += vz * simulation->dt;

            simulation->bodies.x[i] = px;
            simulation->bodies.y[i] = py;
            simulation->bodies.z[i] = pz;

            simulation->vel.x[i] = vx;
            simulation->vel.y[i] = vy;
            simulation->vel.z[i] = vz;
        }

        free(stack);
        free(stackP);
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
void freeAll(BHSim* simulation)
{
    free(simulation->octree->array);
    free(simulation->octree);

    free(simulation->vel.x);
    free(simulation->vel.y);
    free(simulation->vel.z);

    free(simulation->force.x);
    free(simulation->force.y);
    free(simulation->force.z);

    free(simulation->bodies.x);
    free(simulation->bodies.y);
    free(simulation->bodies.z);
    free(simulation->bodies.w);

    free(simulation);
}

int main(int argc, char** argv)
{
    int iterations = 10;

    BHSim* simulation = (BHSim*)malloc(sizeof(BHSim));

    simulation->nBodies = 1000000;
    simulation->dt = 0.001667;
    simulation->theta = 0.5;

    if (argc >= 2) simulation->theta = atof(argv[1]);
    if (argc >= 3) simulation->nBodies = atoi(argv[2]);
    if (argc >= 4) iterations = atoi(argv[3]);

    simulation->bodies.x = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->bodies.y = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->bodies.z = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->bodies.w = (real*)malloc(simulation->nBodies * sizeof(real));

    simulation->vel.x = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->vel.y = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->vel.z = (real*)malloc(simulation->nBodies * sizeof(real));

    simulation->force.x = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->force.y = (real*)malloc(simulation->nBodies * sizeof(real));
    simulation->force.z = (real*)malloc(simulation->nBodies * sizeof(real));

    simulation->octree = (octreeArray*)malloc(sizeof(octreeArray));

    simulation->octree->nBodies = simulation->nBodies;
    simulation->octree->arraySize = simulation->nBodies * 2;
    simulation->octree->array = (octreeNode*)malloc(simulation->octree->arraySize * sizeof(octreeNode));

    randomizeBodies(simulation->bodies, simulation->vel, 1.54f, 8.0f, simulation->nBodies);

    printf("n=%d bodies for %d iterations:\n", simulation->nBodies, iterations);

    for (int it = 0; it < iterations; it++)
    {
        buildOctreeArray(simulation);
        integrateOctreeArray(simulation);
    }

    real3 p_av = average(simulation->bodies, simulation->nBodies);
    printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);
    printf("Body-0  position: (%f,%f,%f)\n", simulation->bodies.x[0], simulation->bodies.y[0], simulation->bodies.z[0]);

    freeAll(simulation);

    return 0;
}


