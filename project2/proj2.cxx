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
#include <vtkFloatArray.h>


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
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
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
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
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
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
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
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
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
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.
//              The first number is the size of the array in argument X,
//              the second the size of Y.
//     X: an array (size is specified by dims).
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims,
                        const float *X, const float *Y, const float *F)
{
    int idx1[3], idx2[3], idx3[3], idx4[3];
    int newidx[3];
    int x, y;
    float sf1, sf2, sf3, sf4;
    float posx1, posy1, posx2, posy2, posx3, posy3, posx4, posy4;
    float t, ttop, tpos, fx, fxtop, fxpos;

    for (x = 0; x < dims[0]-1; x++) {
    //for (y = 0; y < dims[1]-1; y++) {
        for (y = 0; y < dims[1]-1; y++) {
        //for (x = 0; x < dims[0]-1; x++) {

            idx1[0] = x;
            idx1[1] = y;

            idx2[0] = x+1;
            idx2[1] = y;

            idx3[0] = x;
            idx3[1] = y+1;

            idx4[0] = x+1;
            idx4[1] = y+1;

            sf1 = F[GetPointIndex(idx1, dims)];
            posx1 = X[idx1[0]];
            posy1 = Y[idx1[1]];

            sf2 = F[GetPointIndex(idx2, dims)];
            posx2 = X[idx2[0]];
            posy2 = Y[idx2[1]];

            sf3 = F[GetPointIndex(idx3, dims)];
            posx3 = X[idx3[0]];
            posy3 = Y[idx3[1]];

            sf4 = F[GetPointIndex(idx4, dims)];
            posx4 = X[idx4[0]];
            posy4 = Y[idx4[1]];

            if (pt[0] >= posx1 && pt[0] <= posx2) {
                if (pt[1] >= posy1 && pt[1] <= posy3) {

                    //get the bottom linear interpolation for scalar field
                    t = (pt[0] - posx1) / (posx2 - posx1);
                    fx = sf1 + t * (sf2 - sf1);
                    //get the top linear interpolation for scalar field
                    ttop = (pt[0] - posx3) / (posx4 - posx3);
                    fxtop = sf3 + ttop * (sf4 - sf3);
                    //get the interpolated field value
                    tpos = (pt[1] - posy1) / (posy3 - posy1);
                    fxpos = fx + tpos * (fxtop - fx);
                    return fxpos;
                }
            }
        }
    }
    return 0;
}

// ****************************************************************************
//  Function: BoundingBoxForCell
//
//  Arguments:
//     X: an array (size is specified by dims).
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.
//              The first number is the size of the array in argument X,
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//     bbox (output): the bounding box of cellId.  Format should be
//                     bbox[0]: the minimum X value in cellId.
//                     bbox[1]: the maximum X value in cellId.
//                     bbox[2]: the minimum Y value in cellId.
//                     bbox[3]: the maximum Y value in cellId.
//
//  Returns:  None (argument bbox is output)
//
// ****************************************************************************

void
BoundingBoxForCell(const float *X, const float *Y, const int *dims,
                   int cellId, float *bbox)
{
    bbox[0] = -100;
    bbox[1] = +100;
    bbox[2] = -100;
    bbox[3] = +100;
    
    int idx[3];
    GetLogicalCellIndex(idx, cellId, dims);
    if(cellId == 0 || idx[0] != 0){
        bbox[0] = X[idx[0]];
        bbox[1] = X[idx[0] + 1];
        bbox[2] = Y[idx[1]];
        bbox[3] = Y[idx[1] + 1];
    }
}

// ****************************************************************************
//  Function: CountNumberOfStraddlingCells
//
//  Arguments:
//     X: an array (size is specified by dims).
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.
//              The first number is the size of the array in argument X,
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//  Returns:  the number of cells that straddle 0, i.e., the number of cells
//            that contains points who have F>0 and also have points with F<0.
//
// ****************************************************************************

int
CountNumberOfStraddlingCells(const float *X, const float *Y, const int *dims,
                             const float *F)
{

    int idx1[3], idx2[3], idx3[3], idx4[3], newidx[3];
    int x, y;
    int pt1idx, pt2idx, pt3idx, pt4idx;
    int pos = 0, neg = 0, count = 0;

    for (y = 0; y < dims[1]-1; y++) {
        for (x = 0; x < dims[0]-1; x++) {

            idx1[0] = x;
            idx1[1] = y;

            idx2[0] = x+1;
            idx2[1] = y;

            idx3[0] = x;
            idx3[1] = y+1;

            idx4[0] = x+1;
            idx4[1] = y+1;

            pt1idx = GetPointIndex(idx1, dims);
            pt2idx = GetPointIndex(idx2, dims);
            pt3idx = GetPointIndex(idx3, dims);
            pt4idx = GetPointIndex(idx4, dims);

            if (F[pt1idx] > 0 || F[pt2idx] > 0 ||
                F[pt3idx] > 0 || F[pt4idx] > 0) {
                pos++;
            }
            if (F[pt1idx] < 0 || F[pt2idx] < 0 ||
                F[pt3idx] < 0 || F[pt4idx] < 0) {
                neg++;
            }
            if (pos != 0 && neg != 0) {
                count++;
            }
            pos = 0, neg = 0;
        }
    }
    return count;
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    int numCells = CountNumberOfStraddlingCells(X, Y, dims, F);
    cerr << "The number of cells straddling zero is " << numCells << endl;

    float bbox[4];
    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        BoundingBoxForCell(X, Y, dims, cellIds[i], bbox);
        cerr << "The bounding box for cell " << cellIds[i] << " is "
             << bbox[0] << "->" << bbox[1] << ", " << bbox[2] << "->" << bbox[3]
             << endl;
    }

    const int npts = 10;
    float pt[npts][3] =
         {
            {1.01119, 0.122062, 0},
            {0.862376, 1.33839, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.138024, 0},
            {0.384128, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };



    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }


    // cerr << "Infinite loop here, else Windows people may have the terminal "
    //      << "disappear before they see the output."
    //      << " Remove these lines if they annoy you." << endl;
    // cerr << "(press Ctrl-C to exit program)" << endl;
    // while (1) ;
}
