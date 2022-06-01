/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}


class TriangleList
{
   public:
                   TriangleList() { maxTriangles = 1000000; triangleIdx = 0; pts = new float[9*maxTriangles]; };
     virtual      ~TriangleList() { delete [] pts; };

     void          AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxTriangles;
     int           triangleIdx;
};

void
TriangleList::AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
{
    pts[9*triangleIdx+0] = X1;
    pts[9*triangleIdx+1] = Y1;
    pts[9*triangleIdx+2] = Z1;
    pts[9*triangleIdx+3] = X2;
    pts[9*triangleIdx+4] = Y2;
    pts[9*triangleIdx+5] = Z2;
    pts[9*triangleIdx+6] = X3;
    pts[9*triangleIdx+7] = Y3;
    pts[9*triangleIdx+8] = Z3;
    triangleIdx++;
}

vtkPolyData *
TriangleList::MakePolyData(void)
{
    int ntriangles = triangleIdx;
    int numPoints = 3*(ntriangles);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *tris = vtkCellArray::New();
    tris->EstimateSize(numPoints,4);
    for (int i = 0 ; i < ntriangles ; i++)
    {
        double pt[3];
        pt[0] = pts[9*i];
        pt[1] = pts[9*i+1];
        pt[2] = pts[9*i+2];
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[9*i+3];
        pt[1] = pts[9*i+4];
        pt[2] = pts[9*i+5];
        vtk_pts->SetPoint(ptIdx+1, pt);
        pt[0] = pts[9*i+6];
        pt[1] = pts[9*i+7];
        pt[2] = pts[9*i+8];
        vtk_pts->SetPoint(ptIdx+2, pt);
        vtkIdType ids[3] = { ptIdx, ptIdx+1, ptIdx+2 };
        tris->InsertNextCell(3, ids);
        ptIdx += 3;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetPolys(tris);
    tris->Delete();
    vtk_pts->Delete();

    return pd;
}

/* BEGINNING OF THE IMPLEMENTATION OF MY ALGORITHM FOR THIS PROJECT */
class Tetrahedron
{
  public:
    float X[4];
    float Y[4];
    float Z[4];
    float F[4];
    void PrintSelf()
    {
        for (int i = 0 ; i < 4 ; i++)
            printf("\tV%d: (%f, %f, %f) = %f\n", i, X[i], Y[i], Z[i], F[i]);
        printf("\n");
    };
};
float LinearEQ(float FA, float FB, float A, float B, float x)
{   
    float t = (x - A) / (B - A);
    float FX = FA + (t * (FB - FA));
    return FX;
}
void IsosurfaceTet(Tetrahedron &tet, TriangleList &tl, float isoval)
{
    int numSegments[16];
    numSegments[0]  = 0;
    numSegments[1]  = 1;
    numSegments[2]  = 1;
    numSegments[3]  = 2;
    numSegments[4]  = 1;
    numSegments[5]  = 2;
    numSegments[6]  = 2;
    numSegments[7]  = 1;
    numSegments[8]  = 1;
    numSegments[9]  = 2;
    numSegments[10] = 2;
    numSegments[11] = 1;
    numSegments[12] = 2;
    numSegments[13] = 1;
    numSegments[14] = 1;
    numSegments[15] = 0;

    int lup[16][6];
    lup[0][0] = -1;
    lup[0][1] = -1;
    lup[0][2] = -1;
    lup[0][3] = -1;
    lup[0][4] = -1;
    lup[0][5] = -1;
   
    lup[1][0] = 0;
    lup[1][1] = 2;
    lup[1][2] = 5;
    lup[1][3] = -1;
    lup[1][4] = -1;
    lup[1][5] = -1;
   
    lup[2][0] = 0;
    lup[2][1] = 1;
    lup[2][2] = 3;
    lup[2][3] = -1;
    lup[2][4] = -1;
    lup[2][5] = -1;

    lup[3][0] = 1;
    lup[3][1] = 2;
    lup[3][2] = 3;
    lup[3][3] = 5;
    lup[3][4] = -1;
    lup[3][5] = -1;

    lup[4][0] = 1;
    lup[4][1] = 2;
    lup[4][2] = 4;
    lup[4][3] = -1;
    lup[4][4] = -1;
    lup[4][5] = -1;

    lup[5][0] = 0;
    lup[5][1] = 1;
    lup[5][2] = 4;
    lup[5][3] = 5;
    lup[5][4] = -1;
    lup[5][5] = -1;

    lup[6][0] = 0; 
    lup[6][1] = 2;
    lup[6][2] = 3;
    lup[6][3] = 4;
    lup[6][4] = -1;
    lup[6][5] = -1;

    lup[7][0] = 3;
    lup[7][1] = 4;
    lup[7][2] = 5;
    lup[7][3] = -1;
    lup[7][4] = -1;
    lup[7][5] = -1;

    lup[8][0] = 3;
    lup[8][1] = 4;
    lup[8][2] = 5;
    lup[8][3] = -1;
    lup[8][4] = -1;
    lup[8][5] = -1;

    lup[9][0] = 0;
    lup[9][1] = 2;
    lup[9][2] = 3;
    lup[9][3] = 4;
    lup[9][4] = -1;
    lup[9][5] = -1;

    lup[10][0] = 0;
    lup[10][1] = 1;
    lup[10][2] = 4;
    lup[10][3] = 5;
    lup[10][4] = -1;
    lup[10][5] = -1;

    lup[11][0] = 1;
    lup[11][1] = 2;
    lup[11][2] = 4;
    lup[11][3] = -1;
    lup[11][4] = -1;
    lup[11][5] = -1;

    lup[12][0] = 1;
    lup[12][1] = 2;
    lup[12][2] = 3;
    lup[12][3] = 5;
    lup[12][4] = -1;
    lup[12][5] = -1;

    lup[13][0] = 0;
    lup[13][1] = 1;
    lup[13][2] = 3;
    lup[13][3] = -1;
    lup[13][4] = -1;
    lup[13][5] = -1;

    lup[14][0] = 0;
    lup[14][1] = 2;
    lup[14][2] = 5;
    lup[14][3] = -1;
    lup[14][4] = -1;
    lup[14][5] = -1;

    lup[15][0] = -1;
    lup[15][1] = -1;
    lup[15][2] = -1;
    lup[15][3] = -1;
    lup[15][4] = -1;
    lup[15][5] = -1;

    int icase;
    int V3, V2, V1, V0;
    int edge1, edge2, edge3, edge4;
    float x[3];
    float y[3];
    float z[3];
    float p[3];
    for(int i = 0; i < 4; i++){
    	if(i == 0){
	    if(tet.F[0] < isoval){
		V0 = 0;
	    }else{
		V0 = 1;
	    }
	}
    	if(i == 1){
	    if(tet.F[1] < isoval){
		V1 = 0;
	    }else{
		V1 = 2;
	    }
	}
    	if(i == 2){
	    if(tet.F[2] < isoval){
		V2 = 0;
	    }else{
		V2 = 4;
	    }
	}
    	if(i == 3){
	    if(tet.F[3] < isoval){
		V3 = 0;
	    }else{
		V3 = 8;
	    }
	}
     }
     icase = V3 + V2 + V1 + V0;
     int nsegments = numSegments[icase];
     if(nsegments > 0){
         int countx = 0;
    	 int county = 0;
    	 int countz = 0;
	 float x1[4];
         float y1[4];
	 float z1[4];
	
	 int countxx = 0;	
	 int countyy = 0;	
	 int countzz = 0;	
	 if(nsegments == 1){	
	     edge1 = lup[icase][0];
	     edge2 = lup[icase][1];
	     edge3 = lup[icase][2];
	     if( (edge1 == 0) || (edge2 == 0) || (edge3 == 0) ){
	         y[county++] = LinearEQ(tet.Y[0], tet.Y[1], tet.F[0], tet.F[1], isoval);
		 z[countz++] = LinearEQ(tet.Z[0], tet.Z[1], tet.F[0], tet.F[1], isoval);
		 x[countx++] = LinearEQ(tet.X[0], tet.X[1], tet.F[0], tet.F[1], isoval);			
	     } 
	     if( (edge1 == 1) || (edge2 == 1) || (edge3 == 1) ){
		 x[countx++] = LinearEQ(tet.X[1], tet.X[2], tet.F[1], tet.F[2], isoval);
		 y[county++] = LinearEQ(tet.Y[1], tet.Y[2], tet.F[1], tet.F[2], isoval);
		 z[countz++] = LinearEQ(tet.Z[1], tet.Z[2], tet.F[1], tet.F[2], isoval);			
	     } 
	     if( (edge1 == 2) || (edge2 == 2) || (edge3 == 2) ){
		 x[countx++] = LinearEQ(tet.X[0], tet.X[2], tet.F[0], tet.F[2], isoval);
		 z[countz++] = LinearEQ(tet.Z[0], tet.Z[2], tet.F[0], tet.F[2], isoval);
		 y[county++] = LinearEQ(tet.Y[0], tet.Y[2], tet.F[0], tet.F[2], isoval);
	     } 
	     if( (edge1 == 3) || (edge2 == 3) || (edge3 == 3) ){
		 y[county++] = LinearEQ(tet.Y[1], tet.Y[3], tet.F[1], tet.F[3], isoval);
		 x[countx++] = LinearEQ(tet.X[1], tet.X[3], tet.F[1], tet.F[3], isoval);
		 z[countz++] = LinearEQ(tet.Z[1], tet.Z[3], tet.F[1], tet.F[3], isoval);
	     }
	     if( (edge1 == 4) || (edge2 == 4) || (edge3 == 4) ){
		 x[countx++] = LinearEQ(tet.X[2], tet.X[3], tet.F[2], tet.F[3], isoval);
		 y[county++] = LinearEQ(tet.Y[2], tet.Y[3], tet.F[2], tet.F[3], isoval);
		 z[countz++] = LinearEQ(tet.Z[2], tet.Z[3], tet.F[2], tet.F[3], isoval);
	     }
	     if( (edge1 == 5) || (edge2 == 5) || (edge3 == 5) ){
		 z[countz++] = LinearEQ(tet.Z[0], tet.Z[3], tet.F[0], tet.F[3], isoval);
		 x[countx++] = LinearEQ(tet.X[0], tet.X[3], tet.F[0], tet.F[3], isoval);
		 y[county++] = LinearEQ(tet.Y[0], tet.Y[3], tet.F[0], tet.F[3], isoval);
	     }
	     tl.AddTriangle(x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2]);
	 }
	 if(nsegments == 2){
             edge1 = lup[icase][0];
             edge2 = lup[icase][1];
             edge3 = lup[icase][2];
             edge4 = lup[icase][3];
             if( (edge1 == 0) || (edge2 == 0) || (edge3 == 0) || (edge4 == 0) ){
                 y1[countyy++] = LinearEQ(tet.Y[0], tet.Y[1], tet.F[0], tet.F[1], isoval);
                 z1[countzz++] = LinearEQ(tet.Z[0], tet.Z[1], tet.F[0], tet.F[1], isoval);
		 x1[countxx++] = LinearEQ(tet.X[0], tet.X[1], tet.F[0], tet.F[1], isoval);			
             }
             if( (edge1 == 1) || (edge2 == 1) || (edge3 == 1) || (edge4 == 1) ){
                 x1[countxx++] = LinearEQ(tet.X[1], tet.X[2], tet.F[1], tet.F[2], isoval);
                 y1[countyy++] = LinearEQ(tet.Y[1], tet.Y[2], tet.F[1], tet.F[2], isoval);
                 z1[countzz++] = LinearEQ(tet.Z[1], tet.Z[2], tet.F[1], tet.F[2], isoval);
             }
             if( (edge1 == 2) || (edge2 == 2) || (edge3 == 2) || (edge4 == 2) ){
                 x1[countxx++] = LinearEQ(tet.X[0], tet.X[2], tet.F[0], tet.F[2], isoval);
                 z1[countzz++] = LinearEQ(tet.Z[0], tet.Z[2], tet.F[0], tet.F[2], isoval);
                 y1[countyy++] = LinearEQ(tet.Y[0], tet.Y[2], tet.F[0], tet.F[2], isoval);
             }
             if( (edge1 == 3) || (edge2 == 3) || (edge3 == 3) || (edge4 == 3) ){
                 y1[countyy++] = LinearEQ(tet.Y[1], tet.Y[3], tet.F[1], tet.F[3], isoval);
                 x1[countxx++] = LinearEQ(tet.X[1], tet.X[3], tet.F[1], tet.F[3], isoval);
                 z1[countzz++] = LinearEQ(tet.Z[1], tet.Z[3], tet.F[1], tet.F[3], isoval);
             }
             if( (edge1 == 4) || (edge2 == 4) || (edge3 == 4) || (edge4 == 4) ){
                 x1[countxx++] = LinearEQ(tet.X[2], tet.X[3], tet.F[2], tet.F[3], isoval);
                 y1[countyy++] = LinearEQ(tet.Y[2], tet.Y[3], tet.F[2], tet.F[3], isoval);
                 z1[countzz++] = LinearEQ(tet.Z[2], tet.Z[3], tet.F[2], tet.F[3], isoval);
             }
             if( (edge1 == 5) || (edge2 == 5) || (edge3 == 5) || (edge4 == 5) ){
                 z1[countzz++] = LinearEQ(tet.Z[0], tet.Z[3], tet.F[0], tet.F[3], isoval);
                 x1[countxx++] = LinearEQ(tet.X[0], tet.X[3], tet.F[0], tet.F[3], isoval);
                 y1[countyy++] = LinearEQ(tet.Y[0], tet.Y[3], tet.F[0], tet.F[3], isoval);
             }
	     tl.AddTriangle(x1[0], y1[0], z1[0], x1[1], y1[1], z1[1], x1[2], y1[2], z1[2]);
	     if( (icase == 5) || (icase == 10) ){
	         tl.AddTriangle(x1[0], y1[0], z1[0], x1[2], y1[2], z1[2], x1[3], y1[3], z1[3]);
	     }else if( (icase == 3) || (icase == 12) || (icase == 6) || (icase == 9) ){
	    	 tl.AddTriangle(x1[1], y1[1], z1[1], x1[2], y1[2], z1[2], x1[3], y1[3], z1[3]);
	     }
	 }
    }
}

int main()
{
    int  i, j;

/*
    Tetrahedron t;
    t.X[0] = .51;
    t.Y[0] = .52;
    t.Z[0] = .21;
    t.F[0] = 0;
    t.X[1] = .51;
    t.Y[1] = .61;
    t.Z[1] = .22;
    t.F[1] = 1;
    t.X[2] = .52;
    t.Y[2] = .52;
    t.Z[2] = .22;
    t.F[2] = 1;
    t.X[3] = .51;
    t.Y[3] = .52;
    t.Z[3] = .22;
    t.F[3] = 0;
    TriangleList ti;
  //  IsosurfaceTet(t, ti, 0.5);
*/
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Could not find input file." << endl;
        exit(EXIT_FAILURE);
    }

    vtkUnstructuredGrid *ugrid = (vtkUnstructuredGrid *) rdr->GetOutput();
    float *pts = (float *) ugrid->GetPoints()->GetVoidPointer(0);
    float *F = (float *) ugrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    vtkCellArray *cell_array = ugrid->GetCells();

    TriangleList tl;
    cell_array->InitTraversal();
    int ncells = cell_array->GetNumberOfCells();
    cerr << "Number of cells to tetrahedralize is " << ncells << endl;
    int cnt = 0;
    float isoval = 12.2;
    for (int i = 0 ; i < ncells ; i++)
    {
        vtkIdType npts;
        vtkIdType *ids;
        cell_array->GetNextCell(npts, ids);
        if (npts == 4)
        {
            Tetrahedron tet;
            for (int j = 0 ; j < 4 ; j++)
            {
                // This data set is in a huge bounding box.  Normalize as we go.
                tet.X[j] = (pts[3*ids[j]]+3e+7)/6e+7;
                tet.Y[j] = (pts[3*ids[j]+1]+3e+7)/6e+7;
                tet.Z[j] = (pts[3*ids[j]+2]+3e+7)/6e+7;
                tet.F[j] = F[ids[j]];
            }
            IsosurfaceTet(tet, tl, isoval);
        }
        else
        {
            cerr << "Input was non-tetrahedron!!  Ignoring..." << endl;
            cerr << "Type is " << npts << endl;
            cerr << "Cell is " << i << endl;
            cnt++;
            continue;
        }
    }


    vtkPolyData *pd = tl.MakePolyData();

/*
    //This can be useful for debugging
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("proj6_out.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkCleanPolyData *cpd = vtkCleanPolyData::New();
    cpd->SetInputData(pd);
    cpd->SetAbsoluteTolerance(0);
    cpd->PointMergingOn();
    cpd->Update();
    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInputData(cpd->GetOutput());
    //pdn->SetInputData(pd);
    pdn->Update();

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pdn->GetOutput());
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0.5);
    ren1->GetActiveCamera()->SetPosition(0,0,-2);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(0.01, 4);
    //ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
