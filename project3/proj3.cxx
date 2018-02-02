#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
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

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X,
                              const float *Y, const float *F)
{
    int idx1[3], idx2[3], idx3[3], idx4[3];
    int newidx[3];
    int x, y;
    float sf1, sf2, sf3, sf4;
    float posx1, posy1, posx2, posy2, posx3, posy3, posx4, posy4;
    float t, ttop, tpos, fx, fxtop, fxpos;

    for (x = 0; x < dims[0]-1; x++) {
        for (y = 0; y < dims[1]-1; y++) {

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


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose:
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128)
//        F=1: (255,255,255)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    // if (F == 0) {
    //     RGB[0] = 0;
    //     RGB[1] = 0;
    //     RGB[2] = 128;
    // }
    // else if (F == 1) {
    //     RGB[0] = 255;
    //     RGB[1] = 255;
    //     RGB[2] = 255;
    // }
    if (F > 0 && F < 1){
        RGB[0] = 255 * ((F - 0) / (1 - 0));
        RGB[1] = 255 * ((F - 0) / (1 - 0));
        RGB[2] = (255 - 128) * ((F - 0) / (1 - 0)) + 128;
    }
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose:
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128)
//        F=0.5: (255,255,255)
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    // if (F == 0) {
    //     RGB[0] = 0;
    //     RGB[1] = 0;
    //     RGB[2] = 128;
    // }
    // else if (F == 0.5) {
    //     RGB[0] = 255;
    //     RGB[1] = 255;
    //     RGB[2] = 255;
    // }
    // else if (F == 1) {
    //     RGB[0] = 128;
    //     RGB[1] = 0;
    //     RGB[2] = 0;
    // }
    if (F > 0 && F < 0.5){
        RGB[0] = 255 * (F / 0.5);
        RGB[1] = 255 * (F / 0.5);
        RGB[2] = (255 - 128) * (F / 0.5) + 128;
    }
    else if (F > 0.5 && F < 1){
        RGB[0] = (128 - 255) * ((F - 0.5) / (1 - 0.5)) + 255;
        RGB[1] = (0 - 255) * ((F - 0.5) / (1 - 0.5)) + 255;
        RGB[2] = (0 - 255) * ((F - 0.5) / (1 - 0.5)) + 255;
    }

}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose:
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0,
//     and interpolates hue from 0 to 360 degrees
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float Saturation = 1.0, Value = 1.0, Hue = 0.0, f = 0.0;
    float r, g, b, p, q, t;
    int i;

    Hue = (360 * F) / 60.f;
    i = floor(Hue);
    f = Hue - i;

    p = Value * (1.0 - Saturation);
    q = Value * (1.0 - Saturation * f);
    t = Value * (1.0 - Saturation * (1.0 - f));

    if (i == 0) {
        r = Value; g = t; b = p;
    }
    else if (i == 1) {
        r = q; g = Value; b = p;
    }
    else if (i == 2) {
        r = p; g = Value; b = t;
    }
    else if (i == 3) {
        r = p; g = q; b = Value;
    }
    else if (i == 4) {
        r = t; g = p; b = Value;
    }
    else if (i == 5) {
        r = Value; g = p; b = q;
    }

    RGB[0] = r * 255;
    RGB[1] = g * 255;
    RGB[2] = b * 255;

}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = i / (500.0 - 1.0) * 18.0 - 9.0;
            pt[1] = j / (500.0 - 1.0) * 18.0 - 9.0;

            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            float normalizedF = (f - 1.2) / (5.02 - 1.2); //...; see step 5 re 1.2->5.02

            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
