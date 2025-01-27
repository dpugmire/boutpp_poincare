#include "NetCDFLoader.h"
#include "SplineInterpolation.h"
#include "RK4.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

std::ofstream trajspline("/Users/dpn/trajspline.c.txt", std::ofstream::out);
std::ofstream puncFid("/Users/dpn/punc.c.txt");
std::ofstream rawPunc("/Users/dpn/rawpunc.c.txt", std::ofstream::out);
std::ofstream trajOut("/Users/dpn/traj.c.txt", std::ofstream::out);
std::ofstream stepOut("/Users/dpn/steps.c.txt", std::ofstream::out);
std::ofstream rk4Out("/Users/dpn/rk4.c.txt", std::ofstream::out);

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
        this->dz = (this->zmax - this->zmin) / this->nzG;

        this->init_array(this->ziarray, this->nzG+1);
        this->zarray.resize(this->nzG+1);
        for (int i = 0; i <= this->nzG; i++)
            this->zarray[i] = this->ziarray[i] * this->dz;
        this->init_array(this->xiarray, this->nx);
        this->init_array(this->yiarray, this->ny);
        this->init_array(this->xiarray_cfr, this->ixsep);
        this->init_array(this->yiarray_cfr, this->nypf1+1, this->nypf2+2);

        writeArray1DToFile(this->xiarray_cfr, "xiarray_cfr");
        writeArray1DToFile(this->yiarray_cfr, "yiarray_cfr");
        writeArray2DToFile(this->zShift_cfr, "zshift_cfr");
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
    }


    void
    init_array(std::vector<double>& arr, size_t n)
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
};

double INTERP(const std::vector<double>& X, const std::vector<double>& Y, double val)
{
    // Ensure input vectors have the same size
    if (X.size() != Y.size())
        throw std::invalid_argument("X and Y must have the same size.");

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

// Function to find zero crossings
std::tuple<std::vector<double>, std::vector<size_t>, std::vector<double>>
find_zero_crossings(const std::vector<double>& itarray, const std::vector<double>& fl_x3d, double dx)
{
    //x = itarray, y = fl_x3d;
    double x0 = 0.0, x1 = static_cast<double>(itarray.size());
    //min_max_values(x, x0, x1);

    SplineInterpolation spline(itarray, fl_x3d);

    std::vector<double> ffl_x3d, fit;
    for (double xi = x0; xi < x1-1; xi += dx)
    {
        ffl_x3d.push_back(spline.evaluate(xi));
        fit.push_back(xi);
    }

    std::cout<<" spline: "<<ffl_x3d[190]<<" xi= "<<fit[190]<<std::endl;
    std::vector<size_t> crossings;
    std::vector<double> crossing_vals;
    for (size_t i = 1; i < ffl_x3d.size(); ++i)
    {
        if (ffl_x3d[i - 1] * ffl_x3d[i] <= 0 && ffl_x3d[i] != 0.0 && ffl_x3d[i-1]) //sign change
        {
            crossing_vals.push_back(ffl_x3d[i]);
            crossings.push_back(i);
        }
    }
    return std::tuple(fit, crossings, crossing_vals);
}

void
dumpTrajSamples(int iline, std::vector<std::vector<double>> points)
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
    SplineInterpolation splineX(tivals, xvals), splineY(tivals, yvals), splineZ(tivals, zvals);

    int region = -1;
    std::vector<double> tvals;
    double tmin = 0.0, tmax = double(n-1);
    int id = 0;
    double x0 = -1.0;
    for (double t = tmin; t < tmax; t+= 0.01, id++)
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


void processPunctures(
    const std::vector<double>& itarray,
    const std::vector<double>& fl_x3d,
    const std::vector<std::vector<double>>& traj,
    const SplineInterpolation& xiSpline,
    const SplineInterpolation& yiSpline,
    const SplineInterpolation& rxySpline,
    const SplineInterpolation& zxySpline,
    const SplineInterpolation& zsSpline,
    const SplineInterpolation& thetaSpline,
    double zmax, double ixsep, double nypf1, double nypf2, int direction)
{
    std::vector<double> px, py, pz, ptheta, ppsi;
    size_t itmax = itarray.size();

    // Fit spline over itarray and fl_x3d
    SplineInterpolation flSpline(itarray, fl_x3d);

    // Compute zero crossings
    auto findZeroCrossings = [](const std::vector<double>& data) {
        std::vector<size_t> crossings;
        for (size_t i = 1; i < data.size(); ++i)
        {
            if (data[i - 1] * data[i] <= 0)
                crossings.push_back(i - 1);
        }
        return crossings;
    };

    auto iit = findZeroCrossings(fl_x3d);
    size_t nc = iit.size();

    if (nc > 0)
    {
        for (size_t i = 0; i < nc; ++i)
        {
            size_t it = static_cast<size_t>(std::floor(itarray[iit[i]]));
            double a = itarray[iit[i]] - it;
            double b = 1.0 - a;

            // Linear interpolation along field-line
            double xind_tmp = b * traj[1][it] + a * traj[1][it + 1];
            double yind_tmp = b * traj[2][it] + a * traj[2][it + 1];
            double zvalue = b * traj[6][it] + a * traj[6][it + 1];

            if (std::fabs(traj[6][it] - traj[6][it + 1]) > 1.0)
                zvalue = b * std::fmod(traj[6][it], zmax) + a * std::fmod(traj[6][it + 1], zmax);

            if (traj[2][it] == nypf2 && direction == 1 && xind_tmp < ixsep + 0.5)
                yind_tmp = b * traj[2][it] + a * (nypf2 + 1);
            else if (traj[2][it] == nypf1 + 1 && direction == -1 && xind_tmp < ixsep + 0.5)
            {
                yind_tmp = b * (nypf2 + 1) + a * traj[2][it + 1];
                double shiftangle = xiSpline.evaluate(xind_tmp);
                zvalue = std::fmod(zvalue - shiftangle, zmax);
            }

            double rxyvalue, zxyvalue, zsvalue;
            if (xind_tmp < ixsep + 0.5)
            {
                rxyvalue = rxySpline.evaluate(xind_tmp, yind_tmp);
                zxyvalue = zxySpline.evaluate(xind_tmp, yind_tmp);
                zsvalue = zsSpline.evaluate(xind_tmp, yind_tmp);
            }
            else
            {
                rxyvalue = rxySpline.evaluate(xind_tmp, yind_tmp);
                zxyvalue = zxySpline.evaluate(xind_tmp, yind_tmp);
                zsvalue = zsSpline.evaluate(xind_tmp, yind_tmp);
            }

            double ipx3d_tmp = rxyvalue * std::cos(zsvalue);
            double ipy3d_tmp = rxyvalue * std::sin(zsvalue);
            double ipx = ipx3d_tmp * std::cos(zvalue) - ipy3d_tmp * std::sin(zvalue);
            double ipy = ipx3d_tmp * std::sin(zvalue) + ipy3d_tmp * std::cos(zvalue);
            double ipz = zxyvalue;

            if (ipy > 0)
            {
                px.push_back(ipx);
                py.push_back(ipy);
                pz.push_back(ipz);
                ptheta.push_back(thetaSpline.evaluate(yind_tmp));
                ppsi.push_back(xiSpline.evaluate(xind_tmp));
                //puncFid << i + 1 << ", " << ipx << ", " << ipy << ", " << ipz << "\n";
            }
        }
    }
}

int main()
{
    trajspline  <<"ID, STEP, X, Y, Z, REGION\n";
    trajOut<<"ID, STEP, X, Y, Z, REGION\n";
    rawPunc<<"ID, STEP, X, Y, Z\n";
    puncFid<<"ID, STEP, X, Y, Z\n";
    stepOut<<"ID, STEP, X, Y, Z\n";
    rk4Out<<"ID, STEP, X, Y, Z\n"<<std::scientific<<std::setprecision(6);

    std::string fname = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/poincare_clean/stuff.nc";
    Options opts(fname);

    int divertor = 1; //single null
    double xind = 0.0f;

    std::vector<int> LINES;// = {0, 50, 100, 150, 200, 250};
    for (int i = 0; i < 250; i+= 5)
        LINES.push_back(i);
    int nturns = 5;
    LINES = {149};

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
        double xStart = opts.psixy[static_cast<int>(xind)][opts.jyomp];
        int yind = yStart;
        double zind = 0;

        int region = opts.GetRegion(xind, yStart);
        int iturn = 0, it = 0;

        std::cout<<"Region= "<<region<<std::endl;
        while (region < 10 && iturn < nturns)
        {
            // Start field-line tracing.
            for (int iy = 0; iy < opts.ny-1; iy++)
            {
                //trajOut<<iline<<", "<<iy<<", "<<it<<", "<<iturn<<", "<<xStart<<", "<<yStart<<", "<<zStart<<std::endl;
                if (it == 0)
                {
                    Point _p;
                    _p.traj1 = it; _p.traj2 = xind; _p.traj3 = yStart;
                    _p.traj4 = zind; _p.traj5 = region; _p.traj7 = zStart;
                    Points.push_back(_p);
                    //Points.push_back({xind, yStart, zind, iline, iy, it});
                }

                if (yStart+1 == opts.dxdy[0].size())
                {
                    std::cout<<"Overflow of some kind... Need to track this down."<<std::endl;
                    break;
                }

                double xEnd, zEnd, yEnd;
                if (region == 0 && yStart >= opts.nypf1 && yStart < opts.nypf2+1)
                {
                    bool dumpFiles = false;
                    if (it == 343)
                    {
                        std::cout<<"Problem coming...."<<std::endl;
                        dumpFiles = true;
                    }
                    auto step = RK4_FLT1(xStart, yStart, zStart, opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, region, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2, rk4Out, iline, it, dumpFiles);
                    xEnd = step.first;
                    zEnd = step.second;
                    yEnd = yStart+1;
                }
                stepOut<<iline<<", "<<it<<", "<<xEnd<<", "<<yEnd<<", "<<zEnd<<std::endl;

                // Check where the field line ends
                if (xEnd > opts.xMax)
                {
                    std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" reaches outer boundary"<<std::endl;
                    region = 12;
                }
                else if (xEnd < opts.xMin)
                {
                    std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" reaches inner boundary"<<std::endl;
                    region = 11;
                }
                else
                {
                    xind = INTERP(opts.xarray, opts.xiarray, xEnd);
                    if (xind > static_cast<double>(opts.ixsep1) + 0.5)
                    {
                        region = 1;
                        std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" enters the SOL."<<std::endl;
                    }
                }

                //Twist-shift at branch cut.
                if (yStart == opts.nypf2-1 && region == 0)
                {
                    double shiftAngle = INTERP(opts.xiarray, opts.shiftAngle, xind);
                    zEnd = zEnd + shiftAngle;
                    yEnd = opts.nypf1;
                }

                //Relabel toroidal location.
                if (zEnd < opts.zmin || zEnd > opts.zmax)
                    zEnd = double_mod(zEnd, opts.zmax);
                zind = INTERP(opts.zarray, opts.ziarray, zEnd);

                //Points.push_back({xind, yEnd, zind, iline, iy, it});
                Point _p;
                _p.traj1 = it; _p.traj2 = xind; _p.traj3 = yEnd; _p.traj4 = zind;
                _p.traj5 = region; _p.traj7 = zEnd;
                Points.push_back(_p);

                it = it+1;
                xStart = xEnd;
                yStart = yEnd;
                zStart = zEnd;
            }
            iturn++;
        }
	 auto end = std::chrono::high_resolution_clock::now();
	 std::chrono::duration<double> elapsed = end - start;
	 std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
	 
        //Convert to XYZ space.
        std::vector<std::vector<double>> PointsXYZ;

        int id = 0;
        for (const auto& pt : Points)
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
            double zsvalue = INTERP(opts.xiarray, zShift_column, pt.traj4);

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

            trajOut<<iline<<", "<<id<<", "<<x3d<<", "<<y3d<<", "<<z3d<<", "<<region<<std::endl;
            PointsXYZ.push_back({x3d, y3d, z3d});
            id++;
        }
        dumpTrajSamples(iline, PointsXYZ);

        //find the intersections.
        std::vector<double> fl_x3d, fl_y3d, fl_z3d, itarray;
        double xi = 0.0;
        for (const auto pt : PointsXYZ)
        {
            fl_x3d.push_back(pt[0]);
            fl_y3d.push_back(pt[1]);
            fl_z3d.push_back(pt[2]);
            itarray.push_back(xi);
            xi = xi+1.0;
        }

        SplineInterpolation ffl_x3d(itarray, fl_x3d), ffl_y3d(itarray, fl_y3d), ffl_z3d(itarray, fl_z3d);
        auto [fit, iit, iit_vals] = find_zero_crossings(itarray, fl_x3d, 0.0001);
        int nc = iit.size();

        for (int i = 0; i < nc; i++)
        {
            int iit_i = iit[i];
            double tval = fit[iit_i];
            double valX = ffl_x3d.evaluate(tval);
            double valY = ffl_y3d.evaluate(tval);
            double valZ = ffl_z3d.evaluate(tval);
            rawPunc<<i<<", "<<valX<<", "<<valY<<", "<<valZ<<std::endl;
        }

        for (int i = 0; i < nc; i++)
        {
            int iit_i = iit[i];
            double fit_i = fit[iit_i];
            int _it = std::floor(fit_i);
            double a = fit[iit[i]] - static_cast<double>(_it);
            double b = 1.0-a;
            auto pt = Points[_it];
            auto pt_1 = Points[_it+1];

            auto xind_tmp = b*pt.traj2 + a*pt_1.traj2;
            auto yind_tmp = b*pt.traj3 + a*pt_1.traj3;
            auto zvalue = b*pt.traj7 + a*pt_1.traj7;

            auto pt_m1 = Points[_it-1];
            auto _nypf2 = opts.nypf2, _nypf1 = opts.nypf1;
            if (std::abs(pt.traj7 - pt_1.traj7) > 1.0)
                zvalue = b*double_mod(pt.traj7, opts.zmax) + a*double_mod(pt_1.traj7, opts.zmax);
            if (pt.traj3 == double(opts.nypf2) && xind_tmp < double(opts.ixsep)+0.5)
                yind_tmp = b*pt.traj3 + a*double(opts.nypf2+1);
            else if (_it > 0 && (Points[_it-1].traj3 == double(opts.nypf2) || (Points[_it-1].traj3 == double(opts.nypf1+1))))
            {
                std::cout<<"******** Need to support the update of zvalue "<<__LINE__<<std::endl;
                //throw std::string("Need to support the update of zvalue");
                //zvalue = b*
            }

            double rxyvalue, zxyvalue, zsvalue;
            if (xind_tmp < double(opts.ixsep)+0.5)
            {
                //rxyvalue = interpolate2D(opts.xiarray_cfr, opts.yiarray_cfr, opts.rxy_cfr,   xind_tmp, yind_tmp);
                //zxyvalue =   interpolate2D(opts.xiarray_cfr, opts.yiarray_cfr, opts.zxy_cfr,  xind_tmp, yind_tmp);
                //zsvalue =   interpolate2D(opts.xiarray_cfr, opts.yiarray_cfr, opts.zShift_cfr, xind_tmp, yind_tmp);
                rxyvalue = bilinear_interp2(opts.xiarray_cfr, opts.yiarray_cfr, opts.rxy_cfr,   xind_tmp, yind_tmp);
                zxyvalue = bilinear_interp2(opts.xiarray_cfr, opts.yiarray_cfr, opts.zxy_cfr,  xind_tmp, yind_tmp);
                zsvalue =  bilinear_interp2(opts.xiarray_cfr, opts.yiarray_cfr, opts.zShift_cfr, xind_tmp, yind_tmp);
            }
            else
            {
                rxyvalue = bilinear_interp2(opts.xiarray, opts.yiarray, opts.rxy, xind_tmp, yind_tmp);
                zxyvalue = bilinear_interp2(opts.xiarray, opts.yiarray, opts.zxy, xind_tmp, yind_tmp);
                zsvalue = bilinear_interp2(opts.xiarray, opts.yiarray, opts.zShift, xind_tmp, yind_tmp);
                //rxyvalue = interpolate2D(opts.xiarray, opts.yiarray, opts.rxy, xind_tmp, yind_tmp);
                //zxyvalue = interpolate2D(opts.xiarray, opts.yiarray, opts.zxy, xind_tmp, yind_tmp);
                //zsvalue = interpolate2D(opts.xiarray, opts.yiarray, opts.zShift, xind_tmp, yind_tmp);
            }

            double ipx3d_tmp = rxyvalue*cos(zsvalue);
            double ipy3d_tmp = rxyvalue*sin(zsvalue);
            double ipx = ipx3d_tmp*cos(zvalue)-ipy3d_tmp*sin(zvalue);
            double ipy = ipx3d_tmp*sin(zvalue)+ipy3d_tmp*cos(zvalue);
            double ipz = zxyvalue;

            if (ipy > 0.0)
            {
                // do the rxy/zxy interpolation
                //if (i > 0) puncFid<<iit[i]-1<<", "<<fl_x3d[iit[i]-1]<<", "<<fl_y3d[iit[i]-1]<<", "<<fl_z3d[iit[i]-1]<<std::endl;
                //puncFid<<i<<", "<<ipx<<", "<<ipy<<", "<<ipz<<std::endl;
                //puncFid<<iline<<", "<<i<<", "<<ipx<<", "<<ipy<<", "<<ipz<<std::endl;
                //puncFid<<iit[i]+1<<", "<<fl_x3d[iit[i]+1]<<", "<<fl_y3d[iit[i]+1]<<", "<<fl_z3d[iit[i]+1]<<std::endl;
            }
        }
        std::cout<<"All done"<<std::endl;
        std::cout<<" Exiting now."<<std::endl;
    }

    return 0;
}
