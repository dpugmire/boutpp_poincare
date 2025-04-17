#include "NetCDFLoader.h"
#include "SplineInterpolation.h"
#include "RK4.h"
#include "RK4_vtkm.h"
#include "Evaluation.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include <vtkm/Particle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/CellLocatorUniformGrid.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/filter/flow/worklet/CellInterpolationHelper.h>
#include <vtkm/filter/resampling/Probe.h>
#include <vtkm/exec/CellInterpolate.h>


std::ofstream trajspline("/Users/dpn/trajspline.v.txt", std::ofstream::out);
std::ofstream puncFid("/Users/dpn/punc.v.txt");
std::ofstream punc_ip_Fid("/Users/dpn/punc_ip.v.txt");
std::ofstream puncFid2("/Users/dpn/punc2.v.txt");
std::ofstream puncSplineFid("/Users/dpn/puncspline.v.txt");
std::ofstream rawPunc("/Users/dpn/rawpunc.v.txt", std::ofstream::out);
std::ofstream trajOut("/Users/dpn/traj.v.txt", std::ofstream::out);
std::ofstream stepOut("/Users/dpn/steps.v.txt", std::ofstream::out);
std::ofstream rk4Out("/Users/dpn/rk4.v.txt", std::ofstream::out);

double
double_mod(double val, double mod_base)
{
    double result = std::fmod(val, mod_base);
    if (result < 0)
        result += mod_base;
    return result;
}

void
min_max_values(const std::vector<double>& arr, double& vmin, double& vmax)
{
    auto [minIt, maxIt] = std::minmax_element(arr.begin(), arr.end());
    vmin = *minIt;
    vmax = *maxIt;
}

class Point
{
    public:
    Point() = default;
#if 0
    template <typename T, typename U, typename V>
    Point(T _x, U _y, V _z, int _xi, int _yi, int _iter)
    : x(static_cast<double>(_x)), y(static_cast<double>(_y)), z(static_cast<double>(_z)), xi(_xi), yi(_yi), iter(_iter) {}
    Point(const Point& pt) : x(pt.x), y(pt.y), z(pt.z), xi(pt.xi), yi(pt.yi), iter(pt.iter) {}

    Point& operator=(const Point& pt)
    {
        if (this != &pt)
        {
            x = pt.x;
            y = pt.y;
            z = pt.z;
            xi = pt.xi;
            yi = pt.yi;
            iter = pt.iter;
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& pt)
    {
        os << "(" << pt.x << ", " << pt.y << ", " << pt.z << ")";
        return os;
    }
    #endif

    double traj1 = 0.0, traj2 = 0.0, traj3 = 0.0, traj4=0.0, traj5=0.0, traj6=0.0, traj7=0.0;
    /*
    int xi = 0;
    int yi = 0;
    int iter = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    */
};

class Options
{
    public:
    Options(const std::string& fileName)
    : loader(fileName)
    {
        this->nx = this->loader.getDim("nx");
        this->ny = this->loader.getDim("ny");
        this->nz = this->loader.getDim("nz");
        this->nx_cfr = this->loader.getDim("nx_cfr");
        this->ny_cfr = this->loader.getDim("ny_cfr");

        this->rxy = this->loader.read2DVariable("rxy");
        this->zxy = this->loader.read2DVariable("zxy");
        this->rxy_cfr = this->loader.read2DVariable("rxy_cfr");
        this->zxy_cfr = this->loader.read2DVariable("zxy_cfr");
        this->psixy = this->loader.read2DVariable("psixy");
        this->zShift = this->loader.read2DVariable("zShift");
        this->zShift_cfr = this->loader.read2DVariable("zShift_cfr");
        this->shiftAngle = this->loader.read1DVariable("shiftAngle");
        this->dxdy = this->loader.read3DVariable("dxdy");
        this->dzdy = this->loader.read3DVariable("dzdy");
        this->dxdy_m1 = this->loader.read2DVariable("dxdy_m1");
        this->dxdy_p1 = this->loader.read2DVariable("dxdy_p1");
        this->dzdy_m1 = this->loader.read2DVariable("dzdy_m1");
        this->dzdy_p1 = this->loader.read2DVariable("dzdy_p1");
        this->nzG = this->nz * this->zperiod;

        std::cout<<"*******************************************************************************************"<<std::endl;
        std::cout<<"*******************************************************************************************"<<std::endl;
        std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
        std::cout<<" Fix me. zarray length issue..."<<std::endl;
        std::cout<<"      Right now using dz = (zmax - zmin) / (nzG-1)"<<std::endl;
        std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
        std::cout<<"*******************************************************************************************"<<std::endl;
        std::cout<<"*******************************************************************************************"<<std::endl;
        /*
        this->dz = (this->zmax - this->zmin) / this->nzG;
        this->init_array(this->ziarray, this->nzG+1);
        this->zarray.resize(this->nzG+1);
        for (int i = 0; i <= this->nzG; i++)
            this->zarray[i] = this->ziarray[i] * this->dz;
        */
       this->init_array(this->ziarray, this->nzG);
       this->dz = (this->zmax - this->zmin) / (this->nzG-1);
       this->zarray.resize(this->nzG);
       for (int i = 0; i < this->nzG; i++)
           this->zarray[i] = this->ziarray[i] * this->dz;

        this->init_array(this->xiarray, this->nx);
        this->init_array(this->yiarray, this->ny);
        this->init_array(this->xiarray_cfr, this->ixsep);
        this->init_array(this->yiarray_cfr, this->nypf1+1, this->nypf2+2);

        this->init_array(this->yiarray_cfr, this->nypf1+1-1, this->nypf2+2-1);

        writeArray1DToFile(this->xiarray_cfr, "xiarray_cfr");
        writeArray1DToFile(this->yiarray_cfr, "yiarray_cfr");
        writeArray2DToFile(this->zShift_cfr, "zshift_cfr");
        writeArray2DToFile(this->rxy_cfr, "rxy_cfr");
        writeArray2DToFile(this->zxy_cfr, "zxy_cfr");
        writeArray3DToFile(this->dxdy, "dxdy_0");
        writeArray3DToFile(this->dzdy, "dzdy_0");

        //zs_cfr.
        #if 0
        size_t numRows = this->ixsep, numCols = this->nypf2 - this->nypf1 + 1;
        this->zs_cfr.resize(numRows, std::vector<double>(numCols, 0.0));
        for (size_t i = 0; i < numRows; ++i)
            for (size_t j = 0; j < numCols; ++j)
                this->zs_cfr[i][j] = this->zShift[i][this->nypf1 + j];

        // Add an additional column to zs_cfr
        std::cout<<"Compute nu!!! "<<std::endl;
        for (size_t i = 0; i < numRows; ++i)
        {
            // Compute the value for the additional column
            double nuVal1 = 0.0; /* Replace with nu(xiarray_cfr[i], nypf1 + 1) */
            double nuVal2 = 0.0; /* Replace with nu(xiarray_cfr[i], nypf2) */
            double newVal = 0.5 * (nuVal1 + nuVal2) * dy0 + zs_cfr[i][numCols - 1];

            // Append the value to the row
            this->zs_cfr[i].push_back(newVal);
        }
        #endif

        // Find the index of the maximum value in the last row of rxy
        auto& lastRow = this->rxy.back();
        auto maxIt = std::max_element(lastRow.begin(), lastRow.end());
        this->jyomp = std::distance(lastRow.begin(), maxIt);

        for (const auto& row : psixy)
            this->xarray.push_back(row[this->jyomp]);

        auto minmaxIt = std::minmax_element(this->xarray.begin(), this->xarray.end());
        this->xMin = *minmaxIt.first;
        this->xMax = *minmaxIt.second;

        //VTKm stuff.
        vtkm::cont::DataSetBuilderUniform builder;
        vtkm::Id3 dims(this->nx, this->ny, 1);
        vtkm::Vec3f origin(0.0, 0.0, 0.0), spacing(1.0, 1.0, 1.0);
        this->Grid2D = builder.Create(dims, origin, spacing);
        dims = { this->nx_cfr, this->ny_cfr, 1 };
        this->Grid2D_cfr = builder.Create(dims, origin, spacing);
        dims = { this->nx, this->nz, 1 };
        this->Grid2D_xz = builder.Create(dims, origin, spacing);

        //Add fields.
        this->AddField(this->psixy, "psixy", this->Grid2D);
        this->AddField(this->zShift, "zShift", this->Grid2D);
        this->AddField(this->rxy, "rxy", this->Grid2D);
        this->AddField(this->zxy, "zxy", this->Grid2D);

        this->AddField(this->zShift_cfr, "zShift_cfr", this->Grid2D_cfr);
        this->AddField(this->rxy_cfr, "rxy_cfr", this->Grid2D_cfr);
        this->AddField(this->zxy_cfr, "zxy_cfr", this->Grid2D_cfr);


        this->AddField(this->dxdy_m1, "dxdy_m1", this->Grid2D_xz);
        this->AddField(this->dxdy_p1, "dxdy_p1", this->Grid2D_xz);
        this->AddField(this->dzdy_m1, "dzdy_m1", this->Grid2D_xz);
        this->AddField(this->dzdy_p1, "dzdy_p1", this->Grid2D_xz);

        vtkm::cont::DataSetBuilderRectilinear builderRect;
        std::vector<vtkm::FloatDefault> yarray;
        yarray.reserve(this->ny);
        for (int i = 0; i < this->ny; i++)
            yarray.push_back(static_cast<vtkm::FloatDefault>(i));

        std::cout<<"Create rectilinear: "<<this->xarray.size()<<" "<<yarray.size()<<" "<<this->zarray.size()<<std::endl;
        this->Grid3D = builderRect.Create(this->xarray, yarray, this->zarray);
        std::cout<<"Create rectilinear: "<<this->xarray.size()<<" "<<yarray.size()<<" "<<this->zarray.size()<<std::endl;
        this->AddField(this->dxdy, "dxdy", this->Grid3D);
        this->AddField(this->dzdy, "dzdy", this->Grid3D);
        std::cout<<std::setprecision(15)<<"Grid3d bounds: "<<this->Grid3D.GetCoordinateSystem().GetBounds()<<std::endl;
        std::cout<<"  dz= "<<this->dz<<" z0: "<<this->zarray[0]<<" "<<this->zarray[1]<<" ... "<<this->zarray[this->zarray.size()-2]<<" "<<this->zarray[this->zarray.size()-1]<<std::endl;


        #if 0
        this->Grid3D.PrintSummary(std::cout);
        std::cout<<"dims: "<<this->dxdy[0][0].size()<<" "<<this->dxdy[0].size()<<" "<<this->dxdy.size()<<std::endl;
        std::cout<<" npts / ncells: "<<this->Grid3D.GetNumberOfPoints()<<" "<<this->Grid3D.GetNumberOfCells()<<std::endl;

        vtkm::io::VTKDataSetWriter writer("/Users/dpn/grid2D.vtk");
        writer.WriteDataSet(this->Grid2D);
        writer = vtkm::io::VTKDataSetWriter("/Users/dpn/grid2D_cfr.vtk");
        writer.WriteDataSet(this->Grid2D_cfr);
        writer = vtkm::io::VTKDataSetWriter("/Users/dpn/grid2D_xz.vtk");
        writer.WriteDataSet(this->Grid2D_xz);
        writer = vtkm::io::VTKDataSetWriter("/Users/dpn/grid3D.vtk");
        writer.WriteDataSet(this->Grid3D);
        #endif

        this->XArray = vtkm::cont::make_ArrayHandle(this->xarray, vtkm::CopyFlag::On);
        this->XiArray = vtkm::cont::make_ArrayHandle(this->xiarray, vtkm::CopyFlag::On);
        this->ZArray = vtkm::cont::make_ArrayHandle(this->zarray, vtkm::CopyFlag::On);
        this->ZiArray = vtkm::cont::make_ArrayHandle(this->ziarray, vtkm::CopyFlag::On);
        this->ShiftAngle = vtkm::cont::make_ArrayHandle(this->shiftAngle, vtkm::CopyFlag::On);
    }

    void AddField(const std::vector<std::vector<double>>& vals, const std::string& fieldName, vtkm::cont::DataSet& ds)
    {
        vtkm::Id n_y = vals[0].size();
        vtkm::Id n_x = vals.size();

        std::vector<double> field;
        field.reserve(n_x * n_y);
        for (vtkm::Id j = 0; j < n_y; ++j)
            for (vtkm::Id i = 0; i < n_x; ++i)
                field.push_back(vals[i][j]);
        ds.AddPointField(fieldName, field);
    }

    void AddField(const std::vector<std::vector<std::vector<double>>>& vals, const std::string& fieldName, vtkm::cont::DataSet& ds)
    {
        vtkm::Id n_z = vals[0][0].size();
        vtkm::Id n_y = vals[0].size();
        vtkm::Id n_x = vals.size();
        std::cout<<"Add Field3D: "<<fieldName<<" "<<n_x<<" "<<n_y<<" "<<n_z<<std::endl;


        std::vector<double> field;
        field.reserve(n_x * n_y * n_z);
        for (vtkm::Id k = 0; k < n_z; ++k)
            for (vtkm::Id j = 0; j < n_y; ++j)
                for (vtkm::Id i = 0; i < n_x; ++i)
                    field.push_back(vals[i][j][k]);

        std::cout<<"  **** field size: "<<field.size()<<std::endl;
        ds.AddPointField(fieldName, field);
    }

    void init_array(std::vector<double>& arr, size_t n)
    {
        this->init_array(arr, 0, n);
    }

    void
    init_array(std::vector<double>& arr, size_t n0, size_t n1)
    {
        size_t n = n1-n0;
        arr.resize(n);
        for (size_t i = 0; i < n; i++)
            arr[i] = static_cast<double>(i+n0);
    }

    int GetRegion(const vtkm::Vec3f& p) const
    {
        int region = -1;
        if (p[0] < static_cast<double>(this->ixsep + 0.5))
        {
            region = 0;  //Closed flux surface
            std::cout<<"Check this +1 stuff.."<<std::endl;
            if (p[1] < this->nypf1 + 1 || p[1] > nypf2-1)
            {
                region = 2;  //PFR
                std::cout<<" We hit the +1 stuff..."<<std::endl;
            }
        }
        else
            region = 1; //SOL

        // Check for divertor starting points
        if (this->direction == 1 && p[1] == this->ny - 1)
          region = 14;
        else if (this->direction == -1 && p[1] == 0)
            region = 13;

        return region;
    }

    int GetRegion(double xind, int yStart) const
    {
        int region = -1;
        if (xind < static_cast<double>(this->ixsep + 0.5))
        {
            region = 0;  //Closed flux surface
            std::cout<<"Check this +1 stuff.."<<std::endl;
            if (yStart < this->nypf1 + 1 || yStart > nypf2-1)
            {
                region = 2;  //PFR
                std::cout<<" We hit the +1 stuff..."<<std::endl;
            }
        }
        else
            region = 1; //SOL

        // Check for divertor starting points
        if (this->direction == 1 && yStart == this->ny-1)
            region = 14;
        else if (this->direction == -1 && yStart == 0)
            region = 13;

        return region;
    }

    //members.
    int nx = 260;
    int ny = 128;
    int nz = 256;
    int nx_cfr = 195;
    int ny_cfr = 97;
    double xMin, xMax;
    int jyomp;
    double zmin = 0.0;
    double zmax = 2.0 * M_PI;
    double dz;
    int zperiod = 1;
    int jyseps1_1 = 16;
    int jyseps2_2 = 112;
    int ixsep1 = 195;
    int ixsep2 = 260;
    int ixsep = 195;
    int nypf1 = 16;
    int nypf2 = 112;

    int direction = 1;

    int nzG;
    std::vector<double> xiarray, xarray, xiarray_cfr;
    std::vector<double> yiarray, yiarray_cfr;
    std::vector<double> ziarray, zarray;

    std::vector<double> dy;
    double dy0 = 0.0;
    std::vector<std::vector<std::vector<double>>> dxdy, dzdy;
    std::vector<std::vector<double>> rxy, zxy, psixy, zShift;
    std::vector<std::vector<double>> rxy_cfr, zxy_cfr, zShift_cfr;
    std::vector<std::vector<double>> dxdy_p1, dzdy_p1, dxdy_m1, dzdy_m1;
    std::vector<double> shiftAngle;

    NetCDFLoader loader;

    //vtkm data.
    vtkm::cont::DataSet Grid2D, Grid2D_cfr, Grid2D_xz;
    vtkm::cont::DataSet Grid3D;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> XArray, ZArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> XiArray, YiArray, ZiArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> ShiftAngle;
};



double INTERP(const std::vector<double>& X, const std::vector<double>& Y, double val)
{
    // Ensure input vectors have the same size
    /*
    if (X.size() != Y.size())
        throw std::invalid_argument("X and Y must have the same size.");
    */

    // Find the appropriate interval for the value
    int idx = -1;
    for (size_t i = 0; i < X.size() - 1; ++i)
    {
        if (X[i] <= val && val <= X[i + 1])
        {
            idx = i;
            break;
        }
    }

    // Handle out-of-bounds value
    if (idx == -1)
        throw std::out_of_range("val is out of the interpolation range.");

    // Perform linear interpolation
    double x1 = X[idx];
    double x2 = X[idx + 1];
    double y1 = Y[idx];
    double y2 = Y[idx + 1];
    double y_i = y1 + (y2 - y1) * (val - x1) / (x2 - x1);

    return y_i;
}

void
dumpTrajSamples(int iline, const std::vector<vtkm::Vec3f>& points)
{
    int n = points.size();
    std::vector<double> tivals, xvals, yvals, zvals;

    for (int i = 0; i < n; i++)
    {
        tivals.push_back(double(i));
        xvals.push_back(points[i][0]);
        yvals.push_back(points[i][1]);
        zvals.push_back(points[i][2]);
    }

    //Replace with Hermite...
    SplineInterpolation splineX(tivals, xvals), splineY(tivals, yvals), splineZ(tivals, zvals);

    int region = -1;
    std::vector<double> tvals;
    double tmin = 0.0, tmax = double(n-1);
    int id = 0;
    double x0 = -1.0;
    for (double t = tmin; t < tmax; t+= 0.0005, id++)
    {
        double x = splineX.evaluate(t);
        double y = splineY.evaluate(t);
        double z = splineZ.evaluate(t);
        trajspline<<iline<<", "<<t<<", "<<x<<", "<<y<<", "<<z<<", "<<region<<std::endl;
        if (id > 0 && x * x0 < 0.0 && y > 0.0)
            puncFid<<iline<<", "<<t<<", "<<x<<", "<<y<<", "<<z<<std::endl;
        x0 = x;
    }
}

bool signChange(const double& x0, const double& x1)
{
    return (x0 * x1 < 0 && x0 != 0 && x1 != 0);
}

std::vector<vtkm::Vec3f>
ConvertToXYZSpace(const Options& opts, int iline, int id, int region, const std::vector<Point>& pts)
{
    std::vector<vtkm::Vec3f> ptsXYZ;

    auto _fid = fopen("/Users/dpn/problem.c.txt", "w");
    auto trajvals_fid = fopen("/Users/dpn/trajvals.c.txt", "w");
    fprintf(trajvals_fid, "iter, xind, yend, zind, zend, x3d, y3d, z3d\n");

    for (const auto& pt : pts)
    {
        //double x = pt.x, y = pt.y, z = pt.z;
        //int xi = static_cast<int>(x), yi = static_cast<int>(y), zi = static_cast<int>(z);
        int yi = int(pt.traj3);

        // Get the column vector rxy[:, yi]
        std::vector<double> rxy_column(opts.rxy.size());
        for (size_t i = 0; i < opts.rxy.size(); ++i)
            rxy_column[i] = opts.rxy[i][yi];

        // Interpolate values
        double rxyvalue = INTERP(opts.xiarray, rxy_column, pt.traj2);

        std::vector<double> zShift_column(opts.zShift.size());
        for (size_t i = 0; i < opts.zShift.size(); ++i)
            zShift_column[i] = opts.zShift[i][yi];
        const auto& __xi = opts.xiarray;
        double zsvalue = INTERP(opts.xiarray, zShift_column, pt.traj2);

        double zvalue = INTERP(opts.ziarray, opts.zarray, pt.traj4);

        // Compute x3d_tmp and y3d_tmp
        double x3d_tmp = rxyvalue * std::cos(zsvalue);
        double y3d_tmp = rxyvalue * std::sin(zsvalue);

        // Compute x3d and y3d
        double x3d = x3d_tmp * std::cos(zvalue) - y3d_tmp * std::sin(zvalue);
        double y3d = x3d_tmp * std::sin(zvalue) + y3d_tmp * std::cos(zvalue);

        // Get the column vector zxy[:, yi]
        std::vector<double> zxy_column(opts.zxy.size());
        for (size_t i = 0; i < opts.zxy.size(); ++i)
            zxy_column[i] = opts.zxy[i][static_cast<size_t>(yi)];
        double z3d = INTERP(opts.xiarray, zxy_column, pt.traj2);

        trajOut<<iline<<", "<<id<<", "<<x3d<<", "<<y3d<<", "<<z3d<<", "<<region<<", "<<yi<<", "<<zsvalue<<", "<<zvalue<<std::endl;
        //fprintf(_fid, "*** istep= %d traj= %d %12.10e %d %12.10e\n", id, (int)pt.traj1, pt.traj2, (int)pt.traj3, pt.traj4);
        fprintf(_fid, "%d  traj4, zvalue= %12.10e %12.10e\n", id, pt.traj4, zvalue);
        //printf("*** istep= %d r= %12.10e zs= %12.10e z= %12.10e t23= %12.10e %d\n", id, rxyvalue, zsvalue, zvalue, pt.traj2, int(pt.traj3));
        ptsXYZ.push_back({x3d, y3d, z3d});

        fprintf(trajvals_fid, "%d, %10.8f, %d, %10.8f, %10.8f, %10.8f, %10.8f, %10.8f\n", id, pt.traj2, (int)pt.traj3, pt.traj4, pt.traj7, x3d, y3d, z3d);

        id++;
    }

    fclose(trajvals_fid);
    dumpTrajSamples(iline, ptsXYZ);

    return ptsXYZ;
}

std::vector<vtkm::Vec3f>
FindPunctures(int iline, const std::vector<vtkm::Vec3f>& ptsXYZ)
{
    std::size_t n = ptsXYZ.size();
    vtkm::FloatDefault t = 0.0;
    std::vector<vtkm::FloatDefault> xVals, yVals, zVals, tVals;
    for (const auto& pt : ptsXYZ)
    {
        xVals.push_back(pt[0]);
        yVals.push_back(pt[1]);
        zVals.push_back(pt[2]);
        tVals.push_back(t);
        t += 1.0;
    }

    SplineInterpolation splineX(tVals, xVals), splineY(tVals, yVals), splineZ(tVals, zVals);
    std::vector<vtkm::Vec3f> punctures;

    for (std::size_t i = 1; i < n; i++)
    {
        const auto& p0 = ptsXYZ[i-1];
        const auto& p1 = ptsXYZ[i];
        if (p0[0] * p1[0] > 0.0) // check for change sign in X.
            continue;

        if (p0[1] < 0.0 || p1[1] < 0.0)
            continue;

          //crosses the x=0 plane between i-1 and i.
          vtkm::FloatDefault t0 = static_cast<vtkm::FloatDefault>(i - 1);
        vtkm::FloatDefault t1 = static_cast<vtkm::FloatDefault>(i);

        auto x0 = splineX.evaluate(t0);
        auto x1 = splineX.evaluate(t1);
        for (int cnt = 0; cnt < 100; cnt++)
        {
            vtkm::FloatDefault tMid = (t0 + t1) / 2.0;
            auto xMid = splineX.evaluate(tMid);
            if (xMid * x0 < 0.0)
            {
                t1 = tMid;
                x1 = xMid;
            }
            else
            {
                t0 = tMid;
                x0 = xMid;
            }
            vtkm::FloatDefault diff = vtkm::Abs(t0-t1);
            if (vtkm::Abs(t0-t1) < vtkm::Epsilon<vtkm::FloatDefault>())
                break;
        }
        vtkm::FloatDefault tVal = (t0+t1)/2.0;
        vtkm::Vec3f pt(splineX.evaluate(tVal), splineY.evaluate(tVal), splineZ.evaluate(tVal));
        VTKM_ASSERT(vtkm::Abs(pt[0]) <= vtkm::Epsilon<vtkm::FloatDefault>());
        puncSplineFid<<iline<<", "<<tVal<<", "<<pt[0]<<", "<<pt[1]<<", "<<pt[2]<<std::endl;
        punctures.push_back(pt);
    }

    return punctures;
}

int main()
{
    trajspline  <<"ID, STEP, X, Y, Z, REGION\n";
    trajOut<<"ID, STEP, X, Y, Z, REGION, YI, ZSVALUE, ZVALUE"<<std::endl;
    rawPunc<<"ID, STEP, X, Y, Z\n";
    puncFid<<"ID, STEP, X, Y, Z\n";
    puncFid2<<"ID, STEP, X, Y, Z\n";
    puncSplineFid<<"ID, STEP, X, Y, Z\n";
    punc_ip_Fid<<"ID, STEP, X, Y, Z\n";
    stepOut<<"ID, STEP, X, Y, Z\n";
    rk4Out<<"ID, STEP, X, Y, Z\n"<<std::scientific<<std::setprecision(6);
    auto TRAJ_FID = fopen("/Users/dpn/pt_traj.c.txt", "w");
    fprintf(TRAJ_FID, "IT, xind, yEnd, zind, REG, zEnd\n");

    std::string fname = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/poincare_clean/stuff.nc";
    Options opts(fname);

    int divertor = 1; //single null
    double xind = 0.0f;

    std::vector<int> LINES;// = {0, 50, 100, 150, 200, 250};
    for (int i = 0; i < 250; i+= 5)
        LINES.push_back(i);
    int nturns = 50;
    nturns = 300;
    //LINES = {149};
    //LINES = {0,50,100,150,200,250};
    LINES = {150};

    std::vector<Point> Points;

    //std::ofstream trajOut("traj.c.txt", std::ofstream::out);
    //trajOut<<"XI, ITER, X, Y, Z"<<std::endl;

/*
    double _xStart = -0.135524, _yStart = 55, _zStart = 0;
    auto step0 = RK4_FLT1(_xStart, _yStart, _zStart, opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, 0, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2);

    _yStart--;
    auto step2 = RK4_FLT1(_xStart, _yStart, _zStart, opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, 0, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2);
*/

    auto start = std::chrono::high_resolution_clock::now();

    for (const auto& iline : LINES)
    {
        xind = static_cast<double>(iline);
        int yyy = opts.jyomp;
        int yStart = opts.jyomp;
        int zzz = 0;
        auto zStart = opts.zarray[zzz];
        std::cout<<"**** xind= "<<xind<<" opts.jyomp= "<<opts.jyomp<<" zStart= "<<zStart<<std::endl;
        double xStart = opts.psixy[static_cast<int>(xind)][opts.jyomp];
        double xStart2 = scalarField2DEval(opts.Grid2D, "psixy", vtkm::Vec3f(xind, (double)yStart, 0));

        //int yind = yStart;
        double zind = INTERP(opts.zarray, opts.ziarray, zStart);
        auto zind2 = scalarField1DEval(opts.ZiArray, opts.ZArray, {zStart});

        //vtkm::Particle p({ xStart, static_cast<vtkm::FloatDefault>(yStart), zStart }, 0);
        vtkm::Vec3f p0(xStart, static_cast<vtkm::FloatDefault>(yStart), zStart);

        std::cout<<std::setprecision(12);
        std::cout<<"   xStart=  "<<xStart<<std::endl;
        std::cout<<"   xStart2= "<<xStart2<<std::endl;
        std::cout<<"   zind= "<<zind<<std::endl;
        std::cout<<"   zind2= "<<zind<<std::endl;

        int region = opts.GetRegion(xind, p0[1]);
        int region_ = opts.GetRegion(p0);
        int iturn = 0, it = 0;

        auto zindFID = fopen("/Users/dpn/zind.c.txt", "w");
        std::cout<<"Region= "<<region<<std::endl;
        while (region < 10 && iturn < nturns)
        {
            // Start field-line tracing.
            // iy is not used -- just a loop...
            for (int iy_ = 0; iy_ < opts.ny-1; iy_++)
            {
                //if (iturn > 5) break;

                //trajOut<<iline<<", "<<iy<<", "<<it<<", "<<iturn<<", "<<xStart<<", "<<p0[1]<<", "<<p0[2]<<std::endl;
                if (it == 0)
                {
                    Point _p;
                    _p.traj1 = it; _p.traj2 = xind; _p.traj3 = p0[1];
                    _p.traj4 = zind; _p.traj5 = region; _p.traj7 = p0[2];
                    fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", (int)_p.traj1, _p.traj2, (int)_p.traj3, _p.traj4, (int)_p.traj5, _p.traj7);
                    Points.push_back(_p);
                    //Points.push_back({xind, p0[1], zind, iline, iy, it});
                }

                if (p0[1]+1 == opts.dxdy[0].size())
                {
                    std::cout<<"Overflow of some kind... Need to track this down."<<std::endl;
                    break;
                }

                //double xEnd, zEnd, yEnd;
                vtkm::Vec3f p1;
                if (region == 0 && p0[1] >= opts.nypf1 && p0[1] < opts.nypf2+1)
                {
                    bool dumpFiles = false;
                    //auto step = RK4_FLT1(xStart, p0[1], p0[2], opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, region, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2, rk4Out, iline, it, dumpFiles);
                    std::cout<<"Begin Step: iturn= "<<iturn<<" iy= "<<iy_<<std::endl;
                    p1 = RK4_FLT1_vtkm(p0, opts.Grid2D, opts.Grid2D_cfr, opts.Grid2D_xz, opts.Grid3D, opts.XArray, opts.ZArray, region, 1, opts.nypf1, opts.nypf2, rk4Out, iline, it, dumpFiles);
                    p1[1] += 1.0;
                }
                std::cout<<"  *** vRK4: end= "<<p1<<std::endl;
                stepOut<<iline<<", "<<it<<", "<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<std::endl;

                // Check where the field line ends
                if (p1[0] > opts.xMax)
                {
                    std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" reaches outer boundary"<<std::endl;
                    region = 12;
                }
                else if (p1[0] < opts.xMin)
                {
                    std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" reaches inner boundary"<<std::endl;
                    region = 11;
                }
                else
                {
                    //xind = INTERP(opts.xarray, opts.xiarray, xEnd);
                    xind = scalarField1DEval(opts.XArray, opts.XiArray, p1[0]);
                    std::cout<<"   vINTERP:  xind "<<xind<<std::endl;
                    if (xind > static_cast<double>(opts.ixsep1) + 0.5)
                    {
                        region = 1;
                        std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" enters the SOL."<<std::endl;
                    }
                }

                //Twist-shift at branch cut.
                if (p0[1] == opts.nypf2-1 && region == 0)
                {
                    std::cout<<"Branch cut: "<<p0[1]<<" "<<opts.nypf2<<std::endl;
                    //double shiftAngle = INTERP(opts.xiarray, opts.shiftAngle, xind);
                    double shiftAngle = scalarField1DEval(opts.XiArray, opts.ShiftAngle, xind);
                    p1[2] = p1[2] + shiftAngle;
                    p1[1] = opts.nypf1;
                }

                Point _p;
                _p.traj1 = it;
                _p.traj2 = xind;
                _p.traj3 = p1[1];
                _p.traj5 = region;
                _p.traj7 = p1[2];

                double zEnd_no_mod = p1[2];
                //Relabel toroidal location.
                if (p1[2] < opts.zmin || p1[2] > opts.zmax)
                    p1[2] = double_mod(p1[2], opts.zmax);
                //zind = INTERP(opts.zarray, opts.ziarray, zEnd);
                zind = scalarField1DEval(opts.ZArray, opts.ZiArray, {p1[2]});
                _p.traj4 = zind;
                Points.push_back(_p);

                //std::cout<<"********** it= "<<it<<" pt1= "<<xEnd<<" "<<yEnd<<" "<<zEnd<<" zind "<<zind<<std::endl;
                fprintf(zindFID, "%d %12.10f --> %12.10f\n", it, p1[2], zind);
                //Points.push_back({xind, yEnd, zind, iline, iy, it});
                if (p1[2] > opts.zmax)
                {
                    std::cout<<"We have a problem now..."<<std::endl;
                    double diff = std::abs(p1[2] - opts.zmax);
                    std::cout<<"diff is "<<diff<<std::endl;
                    std::cout<<"*******"<<std::endl;
                    throw std::runtime_error("Meow");
                }

                it = it+1;
                p0 = p1;
                //throw std::runtime_error("Meow");

                fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", it, _p.traj2, (int)_p.traj3, _p.traj4, (int)_p.traj5, _p.traj7);
            }
            iturn++;
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
	 
        //Convert to XYZ space.
        int id = 0;
        auto PointsXYZ = ConvertToXYZSpace(opts, iline, id, region, Points);
        auto punctures = FindPunctures(iline, PointsXYZ);
        return 0;
    }
}

