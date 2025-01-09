#include "NetCDFLoader.h"
#include "SplineInterpolation.h"
#include "RK4.h"
#include <iostream>
#include <fstream>
#include <cmath>

float
float_mod(float val, float mod_base)
{
    float result = std::fmod(val, mod_base);
    if (result < 0)
        result += mod_base;
    return result;
}

class Point
{
    public:
    Point() = default;
    //Point(float _x, float _y, float _z) : x(_x), y(_y), z(_z){}
    template <typename T, typename U, typename V>
    Point(T _x, U _y, V _z, int _id, int _iter) : x(static_cast<float>(_x)), y(static_cast<float>(_y)), z(static_cast<float>(_z)), id(_id), iter(_iter) {}
    Point(const Point& pt) : x(pt.x), y(pt.y), z(pt.z), id(pt.id), iter(pt.iter) {}

    Point& operator=(const Point& pt)
    {
        if (this != &pt)
        {
            x = pt.x;
            y = pt.y;
            z = pt.z;
            id = pt.id;
            iter = pt.iter;
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& pt)
    {
        os << "(" << pt.x << ", " << pt.y << ", " << pt.z << ")";
        return os;
    }

    int id = 0;
    int iter = 0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
};

class Options
{
    public:
    Options(const std::string& fileName)
    : loader(fileName)
    {
        this->rxy = this->loader.read2DVariable("rxy");
        this->zxy = this->loader.read2DVariable("zxy");
        this->psixy = this->loader.read2DVariable("psixy");
        this->zShift = this->loader.read2DVariable("zShift");
        this->shiftAngle = this->loader.read1DVariable("shiftAngle");
        this->dxdy = this->loader.read3DVariable("dxdy");
        this->dzdy = this->loader.read3DVariable("dzdy");
        this->dxdy_m1 = this->loader.read2DVariable("dxdy_m1");
        this->dxdy_p1 = this->loader.read2DVariable("dxdy_p1");
        this->dzdy_m1 = this->loader.read2DVariable("dzdy_m1");
        this->dzdy_p1 = this->loader.read2DVariable("dzdy_p1");
        this->nzG = this->nz * this->zperiod;
        this->dz = (this->zmax - this->zmin) / this->nzG;

        this->ziarray.resize(this->nzG+1);
        this->zarray.resize(this->nzG+1);

        for (int i = 0; i <= this->nzG; i++)
        {
            this->ziarray[i] = static_cast<float>(i);
            this->zarray[i] = this->ziarray[i] * this->dz;
        }
        this->xiarray.resize(this->nx);
        for (int i = 0; i < this->nx; i++)
            xiarray[i] = static_cast<float>(i);

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

    int GetRegion(float xind, int yStart) const
    {
        int region = -1;
        if (xind < static_cast<float>(this->ixsep + 0.5))
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
    float xMin, xMax;
    int jyomp;
    float zmin = 0.0;
    float zmax = 2.0 * M_PI;
    float dz;
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
    std::vector<float> xiarray, xarray;
    std::vector<float> ziarray, zarray;

    std::vector<std::vector<std::vector<float>>> dxdy, dzdy;
    std::vector<std::vector<float>> rxy, zxy, psixy, zShift;
    std::vector<std::vector<float>> dxdy_p1, dzdy_p1, dxdy_m1, dzdy_m1;
    std::vector<float> shiftAngle;

    NetCDFLoader loader;
};

float INTERP(const std::vector<float>& X, const std::vector<float>& Y, float val)
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
    float x1 = X[idx];
    float x2 = X[idx + 1];
    float y1 = Y[idx];
    float y2 = Y[idx + 1];
    float y_i = y1 + (y2 - y1) * (val - x1) / (x2 - x1);

    return y_i;
}

int main()
{
    std::string fname = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/poincare_clean/stuff.nc";
    Options opts(fname);

    int divertor = 1; //single null
    float xind = 0.0f;

    std::vector<int> LINES = {149};
    int nturns = 5;

    std::vector<Point> Points;

    for (const auto& iline : LINES)
    {
        xind = static_cast<float>(iline);
        int yyy = opts.jyomp;
        int yStart = opts.jyomp;
        int zzz = 0;
        auto zStart = opts.zarray[zzz];
        float xStart = opts.psixy[static_cast<int>(xind)][opts.jyomp];
        int yind = yStart;
        float zind = 0;


        int region = opts.GetRegion(xind, yStart);
        int iturn = 0, it = 0;

        std::cout<<"Region= "<<region<<std::endl;
        Points.push_back({xind, yStart, zind, iline, it});
        while (region < 10 && iturn < nturns)
        {

            // Start field-line tracing.
            for (int iy = 0; iy < opts.ny; iy++)
            {
                if (iy == 15)
                    std::cout<<"Funky town coming up."<<std::endl;
                if (yStart+1 == opts.dxdy[0].size())
                {
                    std::cout<<"Overflow of some kind... Need to track this down."<<std::endl;
                    break;
                }

                float xEnd, zEnd, yEnd;
                if (region == 0 && yStart >= opts.nypf1 && yStart < opts.nypf2+1)
                {
                    auto step = RK4_FLT1(xStart, yStart, zStart, opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, region, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2);
                    xEnd = step.first;
                    zEnd = step.second;
                    yEnd = yStart+1;
                }

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
                    if (xind > static_cast<float>(opts.ixsep1) + 0.5)
                    {
                        region = 1;
                        std::cout<<"  Starting xind= "<<xind<<" line= "<<iline<<" enters the SOL."<<std::endl;
                    }
                }

                //Twist-shift at branch cut.
                if (yStart == opts.nypf2-1 && region == 0)
                {
                    float shiftAngle = INTERP(opts.xiarray, opts.shiftAngle, xind);
                    zEnd = zEnd + shiftAngle;
                    yEnd = opts.nypf1;
                }

                //Relabel toroidal location.
                if (zEnd < opts.zmin || zEnd > opts.zmax)
                    zEnd = float_mod(zEnd, opts.zmax);
                zind = INTERP(opts.zarray, opts.ziarray, zEnd);

                Points.push_back({xind, yStart, zind, iline, it});

                it = it+1;
                xStart = xEnd;
                yStart = yEnd;
                zStart = zEnd;
            }
        }

        //Convert to XYZ space.
        std::vector<Point> PointsXYZ;

        std::ofstream trajOut("traj.c.txt", std::ofstream::out);
        trajOut<<"ID, ITER, X, Y, Z"<<std::endl;
        for (const auto& pt : Points)
        {
            float x = pt.x, y = pt.y, z = pt.z;
            int xi = static_cast<int>(x), yi = static_cast<int>(y), zi = static_cast<int>(z);

            // Get the column vector rxy[:, yi]
            std::vector<float> rxy_column(opts.rxy.size());
            for (size_t i = 0; i < opts.rxy.size(); ++i)
                rxy_column[i] = opts.rxy[i][static_cast<size_t>(yi)];

            // Interpolate values
            float rxyvalue = INTERP(opts.xiarray, rxy_column, xi);

            std::vector<float> zShift_column(opts.zShift.size());
            for (size_t i = 0; i < opts.zShift.size(); ++i)
                zShift_column[i] = opts.zShift[i][static_cast<size_t>(yi)];
            float zsvalue = INTERP(opts.xiarray, zShift_column, xi);

            float zvalue = INTERP(opts.ziarray, opts.zarray, zi);

            // Compute x3d_tmp and y3d_tmp
            float x3d_tmp = rxyvalue * std::cos(zsvalue);
            float y3d_tmp = rxyvalue * std::sin(zsvalue);

            // Compute x3d and y3d
            float x3d = x3d_tmp * std::cos(zvalue) - y3d_tmp * std::sin(zvalue);
            float y3d = x3d_tmp * std::sin(zvalue) + y3d_tmp * std::cos(zvalue);

            // Get the column vector zxy[:, yi]
            std::vector<float> zxy_column(opts.zxy.size());
            for (size_t i = 0; i < opts.zxy.size(); ++i)
                zxy_column[i] = opts.zxy[i][static_cast<size_t>(yi)];
            float z3d = INTERP(opts.xiarray, zxy_column, xi);

            trajOut<<pt.id<<", "<<pt.iter<<", "<<x3d<<", "<<y3d<<", "<<z3d<<std::endl;
        }

    }


#if 0
    // Example data
    std::vector<float> xarray = {0, 1, 2, 3, 4};
    std::vector<float> zarray = {0, 1, 2, 3};
    std::vector<std::vector<float>> dxdyp = {
        {0, 1, 4, 9},
        {1, 2, 5, 10},
        {4, 5, 8, 13},
        {9, 10, 13, 18},
        {16, 17, 20, 25}
    };

    float xStart = 2.5, zStart = 1.5;
    try {
        float result = interp2D(xarray, zarray, dxdyp, xStart, zStart);
        std::cout << "Interpolated value: " << result << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }



    std::vector<std::vector<float>> rxy = opts.rxy;
    std::cout << "Read 2D variable 'rxy' with dimensions: " << rxy.size() << " x " << rxy[0].size() << "\n";

    // Print a few values
    std::cout << "First value in rxy: " << rxy[0][0] << "\n";
#endif
    return 0;
}
