
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


class TriangleList
{
   public:
                   TriangleList() { maxTriangles = 1000000; triangleIdx = 0; pts = new float[4*maxTriangles]; };
     virtual      ~TriangleList() { delete [] pts; };

     inline void          AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
     {
         if (triangleIdx >= maxTriangles)
         {
             cerr << "No room for more triangles!" << endl;
             return;
         }
     
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
     };

     inline vtkPolyData  *MakePolyData(void)
     {
         int ntriangles = triangleIdx;
         int numPoints = 3*(ntriangles);
         vtkPoints *vtk_pts = vtkPoints::New();
         vtk_pts->SetNumberOfPoints(numPoints);
         int ptIdx = 0;
         vtkCellArray *tris = vtkCellArray::New();
         tris->EstimateSize(numPoints,3);
         for (int i = 0 ; i < ntriangles ; i++)
         {
             double pt[3];
             pt[0] = pts[9*i+0];
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
     };
     
   protected:
     float        *pts;
     int           maxTriangles;
     int           triangleIdx;
};


