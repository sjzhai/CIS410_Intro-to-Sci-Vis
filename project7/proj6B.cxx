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
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "TriangleList.h"
#include "tricase.cxx"
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
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
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


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

void triPos(int edge, float **pos, int *ptidx, float *result, float *F, float isovalue){
    if (edge == 0) {
        result[0] = pos[0][0] + (abs(isovalue - F[ptidx[0]]) / abs(F[ptidx[1]] - F[ptidx[0]])) * (pos[1][0] - pos[0][0]);
        result[1] = pos[0][1];
        result[2] = pos[0][2];
    } else if (edge == 1) {
        result[0] = pos[1][0];
        result[1] = pos[1][1] + (abs(isovalue - F[ptidx[1]]) / abs(F[ptidx[3]] - F[ptidx[1]])) * (pos[3][1] - pos[1][1]);
        result[2] = pos[1][2];
    } else if (edge == 2) {
        result[0] = pos[2][0] + (abs(isovalue - F[ptidx[2]]) / abs(F[ptidx[3]] - F[ptidx[2]])) * (pos[3][0] - pos[2][0]);
        result[1] = pos[2][1];
        result[2] = pos[2][2];
    } else if (edge == 3) {
        result[0] = pos[0][0];
        result[1] = pos[0][1] + (abs(isovalue - F[ptidx[0]]) / abs(F[ptidx[2]] - F[ptidx[0]])) * (pos[2][1] - pos[0][1]);
        result[2] = pos[0][2];
    } else if (edge == 4) {
        result[0] = pos[4][0] + (abs(isovalue - F[ptidx[4]]) / abs(F[ptidx[5]] - F[ptidx[4]])) * (pos[5][0] - pos[4][0]);
        result[1] = pos[4][1];
        result[2] = pos[4][2];
    } else if (edge == 5) {
        result[0] = pos[5][0];
        result[1] = pos[5][1] + (abs(isovalue - F[ptidx[5]]) / abs(F[ptidx[7]] - F[ptidx[5]])) * (pos[7][1] - pos[5][1]);
        result[2] = pos[5][2];
    } else if (edge == 6) {
        result[0] = pos[6][0] + (abs(isovalue - F[ptidx[6]]) / abs(F[ptidx[7]] - F[ptidx[6]])) * (pos[7][0] - pos[6][0]);
        result[1] = pos[6][1];
        result[2] = pos[6][2];
    } else if (edge == 7) {
        result[0] = pos[4][0];
        result[1] = pos[4][1] + (abs(isovalue - F[ptidx[4]]) / abs(F[ptidx[6]] - F[ptidx[4]])) * (pos[6][1] - pos[4][1]);
        result[2] = pos[4][2];
    } else if (edge == 8) {
        result[0] = pos[0][0];
        result[1] = pos[0][1];
        result[2] = pos[0][2] + (abs(isovalue - F[ptidx[0]]) / abs(F[ptidx[4]] - F[ptidx[0]])) * (pos[4][2] - pos[0][2]);
    } else if (edge == 9) {
        result[0] = pos[1][0];
        result[1] = pos[1][1];
        result[2] = pos[1][2] + (abs(isovalue - F[ptidx[1]]) / abs(F[ptidx[5]] - F[ptidx[1]])) * (pos[5][2] - pos[1][2]);
    } else if (edge == 10) {
        result[0] = pos[2][0];
        result[1] = pos[2][1];
        result[2] = pos[2][2] + (abs(isovalue - F[ptidx[2]]) / abs(F[ptidx[6]] - F[ptidx[2]])) * (pos[6][2] - pos[2][2]);
    } else if (edge == 11) {
        result[0] = pos[3][0];
        result[1] = pos[3][1];
        result[2] = pos[3][2] + (abs(isovalue - F[ptidx[3]]) / abs(F[ptidx[7]] - F[ptidx[3]])) * (pos[7][2] - pos[3][2]);
    }
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    TriangleList tl;
    int x, y, z, tcase;
    float isovalue = 3.2;
    int idx[3], posidx[3];
    float **pos = new float *[8];
    for (int p = 0; p < 8; p++) { pos[p] = new float[3]; }
    int *ptidx = new int[8];
    int bin0, bin1, bin2, bin3, bin4, bin5, bin6, bin7;

    for (x = 0; x < dims[0]-1; x++) {
        for (y = 0; y < dims[1]-1; y++) {
            for (z = 0; z < dims[2]-1; z++) {
                pos[0][0] = X[x]; pos[0][1] = Y[y]; pos[0][2] = Z[z];
                idx[0] = x; idx[1] = y; idx[2] = z;
                ptidx[0] = GetPointIndex(idx, dims);

                pos[1][0] = X[x+1]; pos[1][1] = Y[y]; pos[1][2] = Z[z];
                idx[0] = x+1; idx[1] = y; idx[2] = z;
                ptidx[1] = GetPointIndex(idx, dims);

                pos[2][0] = X[x]; pos[2][1] = Y[y+1]; pos[2][2] = Z[z];
                idx[0] = x; idx[1] = y+1; idx[2] = z;
                ptidx[2] = GetPointIndex(idx, dims);

                pos[3][0] = X[x+1]; pos[3][1] = Y[y+1]; pos[3][2] = Z[z];
                idx[0] = x+1; idx[1] = y+1; idx[2] = z;
                ptidx[3] = GetPointIndex(idx, dims);

                pos[4][0] = X[x]; pos[4][1] = Y[y]; pos[4][2] = Z[z+1];
                idx[0] = x; idx[1] = y; idx[2] = z+1;
                ptidx[4] = GetPointIndex(idx, dims);

                pos[5][0] = X[x+1]; pos[5][1] = Y[y]; pos[5][2] = Z[z+1];
                idx[0] = x+1; idx[1] = y; idx[2] = z+1;
                ptidx[5] = GetPointIndex(idx, dims);

                pos[6][0] = X[x]; pos[6][1] = Y[y+1]; pos[6][2] = Z[z+1];
                idx[0] = x; idx[1] = y+1; idx[2] = z+1;
                ptidx[6] = GetPointIndex(idx, dims);

                pos[7][0] = X[x+1]; pos[7][1] = Y[y+1]; pos[7][2] = Z[z+1];
                idx[0] = x+1; idx[1] = y+1; idx[2] = z+1;
                ptidx[7] = GetPointIndex(idx, dims);

                if (F[ptidx[0]] < isovalue) {bin0 = 0;}
                else if (F[ptidx[0]] > isovalue) {bin0 = 1;}
                if (F[ptidx[1]] < isovalue) {bin1 = 0;}
                else if (F[ptidx[1]] > isovalue) {bin1 = 1;}
                if (F[ptidx[2]] < isovalue) {bin2 = 0;}
                else if (F[ptidx[2]] > isovalue) {bin2 = 1;}
                if (F[ptidx[3]] < isovalue) {bin3 = 0;}
                else if (F[ptidx[3]] > isovalue) {bin3 = 1;}
                if (F[ptidx[4]] < isovalue) {bin4 = 0;}
                else if (F[ptidx[4]] > isovalue) {bin4 = 1;}
                if (F[ptidx[5]] < isovalue) {bin5 = 0;}
                else if (F[ptidx[5]] > isovalue) {bin5 = 1;}
                if (F[ptidx[6]] < isovalue) {bin6 = 0;}
                else if (F[ptidx[6]] > isovalue) {bin6 = 1;}
                if (F[ptidx[7]] < isovalue) {bin7 = 0;}
                else if (F[ptidx[7]] > isovalue) {bin7 = 1;}

                tcase = bin0 * pow(2, 0) + bin1 * pow(2, 1) + bin2 * pow(2, 2) + bin3 * pow(2, 3)
                + bin4 * pow(2, 4) + bin5 * pow(2, 5) + bin6 * pow(2, 6) + bin7 * pow(2, 7);

                int edge1, edge2, edge3; // find which edge has point
                float pt1[3], pt2[3], pt3[3];
                float result[10];
                for(int e = 0; e < 15; e=e+3) {
                    edge1 = triCase[tcase][e];
                    edge2 = triCase[tcase][e+1];
                    edge3 = triCase[tcase][e+2];
                    if(edge1 != -1 && edge2 != -1 && edge3 != -1){
                        triPos(edge1, pos, ptidx, pt1, F, isovalue);
                        triPos(edge2, pos, ptidx, pt2, F, isovalue);
                        triPos(edge3, pos, ptidx, pt3, F, isovalue);
                        tl.AddTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
                    }
                }
            }
        }
    }
    vtkPolyData *pd = tl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
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

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
