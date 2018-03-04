
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
#include <vtkColorTransferFunction.h>
#include <vtkPlaneSource.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkPointData.h>
#include <vtkArrowSource.h>

#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPointSource.h>
#include <vtkMaskPoints.h>
#include <vtkTransform.h>

#include <vtkRungeKutta4.h>
#include <vtkStreamTracer.h>
#include <vtkLineSource.h>

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


int main(int, char*[])
{
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj8.vtk");
    rdr->Update();

    // One render window
    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->SetSize(800, 800);

    // One interactor
    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);

    // Define viewport ranges
    double Vport1[4] = {0.0, 0.0, 0.5, 0.5};
    double Vport2[4] = {0.0, 0.5, 0.5, 1.0};
    double Vport3[4] = {0.5, 0.0, 1.0, 0.5};
    double Vport4[4] = {0.5, 0.5, 1.0, 1.0};

    // Set four renderers
//window1 (bottomleft window)
//--------------------------------------------------------//
    vtkSmartPointer<vtkContourFilter> iso =
      vtkSmartPointer<vtkContourFilter>::New();
    iso->SetNumberOfContours(2);
    iso->SetValue(0, 2.5);
    iso->SetValue(1, 5.0);
    iso->SetInputConnection(rdr->GetOutputPort());
    iso->Update();

    vtkSmartPointer<vtkPolyDataMapper> win1Mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    win1Mapper->SetInputConnection(iso->GetOutputPort());
    win1Mapper->ScalarVisibilityOn();

    //coloring isosurface
    vtkSmartPointer<vtkColorTransferFunction> lut =
    vtkSmartPointer<vtkColorTransferFunction>::New();
    lut->SetColorSpaceToRGB();
    lut->AddRGBPoint(2.5,1,0,0);//red
    lut->AddRGBPoint(5.0,1,1,1);//white
    lut->SetScaleToLinear();
    win1Mapper->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> Renderer1 =
      vtkSmartPointer<vtkRenderer>::New();

    renWin->AddRenderer(Renderer1);
    Renderer1->SetViewport(Vport1);
    Renderer1->AddActor(win1Actor);
    Renderer1->SetBackground(0.0, 0.0, 0.0);
//--------------------------------------------------------//

//window2 (upperleft window)
//--------------------------------------------------------//

     // Create a plane
     // xz=(1,0,0), XY=(0,0,1),YZ=(0,1,0)
     vtkSmartPointer<vtkPlane> plane1 =
       vtkSmartPointer<vtkPlane>::New();
     plane1->SetNormal(1.0,0.0,0.0);

     vtkSmartPointer<vtkPlane> plane2 =
       vtkSmartPointer<vtkPlane>::New();
     plane2->SetNormal(0.0,0.0,1.0);

     vtkSmartPointer<vtkPlane> plane3 =
       vtkSmartPointer<vtkPlane>::New();
     plane3->SetNormal(0.0,1.0,0.0);

    //create cutter
    vtkSmartPointer<vtkCutter> cutter1 =
      vtkSmartPointer<vtkCutter>::New();
    cutter1->SetCutFunction(plane1);
    cutter1->SetInputConnection(rdr->GetOutputPort());
    cutter1->Update();

    vtkSmartPointer<vtkCutter> cutter2 =
      vtkSmartPointer<vtkCutter>::New();
    cutter2->SetCutFunction(plane2);
    cutter2->SetInputConnection(rdr->GetOutputPort());
    cutter2->Update();

    vtkSmartPointer<vtkCutter> cutter3 =
      vtkSmartPointer<vtkCutter>::New();
    cutter3->SetCutFunction(plane3);
    cutter3->SetInputConnection(rdr->GetOutputPort());
    cutter3->Update();

    //get colors from reader
    double range[2];
    rdr->GetOutput()->GetPointData()->GetScalars()->GetRange(range);

    vtkSmartPointer<vtkPolyDataMapper> win2cutMapper1 =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    win2cutMapper1->SetInputConnection(cutter1->GetOutputPort());
    win2cutMapper1->SetScalarRange(range);

    vtkSmartPointer<vtkPolyDataMapper> win2cutMapper2 =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    win2cutMapper2->SetInputConnection(cutter2->GetOutputPort());
    win2cutMapper2->SetScalarRange(range);

    vtkSmartPointer<vtkPolyDataMapper> win2cutMapper3 =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    win2cutMapper3->SetInputConnection(cutter3->GetOutputPort());
    win2cutMapper3->SetScalarRange(range);

    vtkSmartPointer<vtkActor> planeActor1 =
      vtkSmartPointer<vtkActor>::New();
    planeActor1->SetMapper(win2cutMapper1);

    vtkSmartPointer<vtkActor> planeActor2 =
      vtkSmartPointer<vtkActor>::New();
    planeActor2->SetMapper(win2cutMapper2);

    vtkSmartPointer<vtkActor> planeActor3 =
      vtkSmartPointer<vtkActor>::New();
    planeActor3->SetMapper(win2cutMapper3);

    vtkSmartPointer<vtkRenderer> Renderer2 =
      vtkSmartPointer<vtkRenderer>::New();

    renWin->AddRenderer(Renderer2);
    Renderer2->SetViewport(Vport2);
    Renderer2->AddActor(planeActor1);
    Renderer2->AddActor(planeActor2);
    Renderer2->AddActor(planeActor3);
    Renderer2->SetBackground(0.0, 0.0, 0.0);
//--------------------------------------------------------//


//window3 (bottomright window)
//--------------------------------------------------------//
    rdr->GetOutput()->GetPointData()->SetActiveAttribute("grad", vtkDataSetAttributes::VECTORS);

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    vtkSmartPointer<vtkPointSource> points =
      vtkSmartPointer<vtkPointSource>::New();
    points->SetNumberOfPoints(GetNumberOfPoints(dims));
    points->Update();

    vtkSmartPointer<vtkArrowSource> arrowSource =
      vtkSmartPointer<vtkArrowSource>::New();
    arrowSource->Update();

    vtkSmartPointer<vtkMaskPoints> ptMask =
      vtkSmartPointer<vtkMaskPoints>::New();
    ptMask->SetInputConnection(rdr->GetOutputPort());
    ptMask->SetOnRatio(203);
    ptMask->Update();

    vtkSmartPointer<vtkGlyph3D> glyph =
      vtkSmartPointer<vtkGlyph3D>::New();
    glyph->SetInputConnection(ptMask->GetOutputPort());
    glyph->SetSourceConnection(arrowSource->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> win3Mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    win3Mapper->SetInputConnection(glyph->GetOutputPort());
    win3Mapper->SetScalarRange(range);

    vtkSmartPointer<vtkActor> win3Actor =
      vtkSmartPointer<vtkActor>::New();
    win3Actor->SetMapper(win3Mapper);

    vtkSmartPointer<vtkRenderer> Renderer3 =
      vtkSmartPointer<vtkRenderer>::New();

    renWin->AddRenderer(Renderer3);
    Renderer3->SetViewport(Vport3);
    Renderer3->AddActor(win3Actor);
    Renderer3->SetBackground(0.0, 0.0, 0.0);

//--------------------------------------------------------//
//window4 (upperright window)
//--------------------------------------------------------//
    vtkSmartPointer<vtkLineSource> lines =
      vtkSmartPointer<vtkLineSource>::New();
    lines->SetResolution(18);
    lines->SetPoint1(-9.0, 0.0, 0.0);
    lines->SetPoint2(9.0, 0.0, 0.0);

    vtkRungeKutta4 *rk4 = vtkRungeKutta4::New();

    vtkSmartPointer<vtkStreamTracer> streamer =
      vtkSmartPointer<vtkStreamTracer>::New();
    streamer->SetInputConnection(rdr->GetOutputPort());
    streamer->SetSourceConnection(lines->GetOutputPort());
    streamer->SetIntegrator(rk4);
    streamer->SetMaximumPropagation(GetNumberOfPoints(dims));
    streamer->SetInitialIntegrationStep(.2);
    streamer->SetMinimumIntegrationStep(.01);
    streamer->Update();

    vtkSmartPointer<vtkDataSetMapper> win4Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win4Mapper->SetInputConnection(streamer->GetOutputPort());
    win4Mapper->SetScalarRange(range);

    vtkSmartPointer<vtkActor> win4Actor =
      vtkSmartPointer<vtkActor>::New();
    win4Actor->SetMapper(win4Mapper);

    vtkSmartPointer<vtkRenderer> Renderer4 =
      vtkSmartPointer<vtkRenderer>::New();

    renWin->AddRenderer(Renderer4);
    Renderer4->SetViewport(Vport4);
    Renderer4->AddActor(win4Actor);
    Renderer4->SetBackground(0.0, 0.0, 0.0);
//--------------------------------------------------------//

//set active cameras for four renderers
    Renderer1->GetActiveCamera()->SetFocalPoint(0,0,0);
    Renderer1->GetActiveCamera()->SetPosition(0,0,60);
    Renderer1->GetActiveCamera()->SetViewUp(0,1,0);
    Renderer1->GetActiveCamera()->SetClippingRange(20, 120);
    Renderer1->GetActiveCamera()->SetDistance(30);

    Renderer2->GetActiveCamera()->SetPosition(60,70,40);
    Renderer2->GetActiveCamera()->SetViewUp(1,0,0);
    Renderer2->GetActiveCamera()->Zoom(1.5);

    iren->SetRenderWindow(renWin);
    renWin->Render();
    iren->Start();

}
