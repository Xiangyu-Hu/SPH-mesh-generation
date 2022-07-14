#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "polygonise.h"

namespace Polygonise_MQ{

  int Polygonise (GRIDCELL grid,double isolevel,std::vector<TRIANGLE> &triangles)
  {
    int i,ntriang;
    int cubeindex;
    XYZ vertlist[12];

     /*
        Determine the index into the edge table which
        tells us which vertices are inside of the surface
     */
     cubeindex = 0;
     if (grid.val[0] < isolevel) cubeindex |= 1;
     if (grid.val[1] < isolevel) cubeindex |= 2;
     if (grid.val[2] < isolevel) cubeindex |= 4;
     if (grid.val[3] < isolevel) cubeindex |= 8;
     if (grid.val[4] < isolevel) cubeindex |= 16;
     if (grid.val[5] < isolevel) cubeindex |= 32;
     if (grid.val[6] < isolevel) cubeindex |= 64;
     if (grid.val[7] < isolevel) cubeindex |= 128;

     /* Cube is entirely in/out of the surface */
     if (edgeTable[cubeindex] == 0)
        return(0);

     /* Find the vertices where the surface intersects the cube */
     if (edgeTable[cubeindex] & 1)
        vertlist[0] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
     if (edgeTable[cubeindex] & 2)
        vertlist[1] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
     if (edgeTable[cubeindex] & 4)
        vertlist[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
     if (edgeTable[cubeindex] & 8)
        vertlist[3] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
     if (edgeTable[cubeindex] & 16)
        vertlist[4] = VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
     if (edgeTable[cubeindex] & 32)
        vertlist[5] = VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
     if (edgeTable[cubeindex] & 64)
        vertlist[6] = VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
     if (edgeTable[cubeindex] & 128)
        vertlist[7] = VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
     if (edgeTable[cubeindex] & 256)
        vertlist[8] = VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
     if (edgeTable[cubeindex] & 512)
        vertlist[9] = VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
     if (edgeTable[cubeindex] & 1024)
        vertlist[10] = VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
     if (edgeTable[cubeindex] & 2048)
        vertlist[11] = VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

     /* Create the triangle */
     ntriang = 0;
     for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
      TRIANGLE tri;
      tri.p[0] = vertlist[triTable[cubeindex][i  ]];
      tri.p[1] = vertlist[triTable[cubeindex][i+1]];
      tri.p[2] = vertlist[triTable[cubeindex][i+2]];
      triangles.push_back(tri);
      ntriang++;
     }

     return(ntriang);
  }

  /*
     Linearly interpolate the position where an isosurface cuts
     an edge between two vertices, each with their own scalar value
  */
  XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2)
  {
     double mu;
     XYZ p;

     if (fabs(isolevel-valp1) < 1.0e-15)
        return(p1);
     if (fabs(isolevel-valp2) < 1.0e-15)
        return(p2);
     if (fabs(valp1-valp2) < 1.0e-15)
        return(p1);
     mu = (isolevel - valp1) / (valp2 - valp1);
     p.x = p1.x + mu * (p2.x - p1.x);
     p.y = p1.y + mu * (p2.y - p1.y);
     p.z = p1.z + mu * (p2.z - p1.z);

     return(p);
  }
}