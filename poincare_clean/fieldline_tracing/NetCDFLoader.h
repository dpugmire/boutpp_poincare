#pragma once

#include <netcdf.h>
#include <string>
#include <vector>


class NetCDFLoader
{
public:
    explicit NetCDFLoader(const std::string& filePath)
    {
        int status = nc_open(filePath.c_str(), NC_NOWRITE, &this->ncid);
        if (status != NC_NOERR)
            throw std::runtime_error("Error opening NetCDF file: " + std::string(nc_strerror(status)));
    }

    ~NetCDFLoader()
    {
        nc_close(this->ncid);
    }

    std::vector<double> read1DVariable(const std::string& varName);
    std::vector<std::vector<double>> read2DVariable(const std::string& varName, bool transpose=true);
    std::vector<std::vector<std::vector<double>>> read3DVariable(const std::string& varName, bool permute=true);

private:
    int ncid;

    void getVariableDimensions(int varid, size_t* dimSizes, int numDims)
    {
        int ndims;
        nc_inq_varndims(this->ncid, varid, &ndims);
        if (ndims != numDims)
            throw std::runtime_error("Variable does not have " + std::to_string(numDims) + " dimensions.");

        int dimids[numDims];
        nc_inq_vardimid(this->ncid, varid, dimids);

        for (int i = 0; i < numDims; ++i)
        {
            size_t dimSize;
            nc_inq_dimlen(this->ncid, dimids[i], &dimSize);
            dimSizes[i] = dimSize;
        }
    }

    // Transpose a 2D vector
    std::vector<std::vector<double>> transpose2D(const std::vector<std::vector<double>>& data)
    {
        size_t rows = data.size();
        size_t cols = data[0].size();
        std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows, 0.0));

        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < cols; ++j)
                transposed[j][i] = data[i][j];
        }

        return transposed;
    }
};
std::vector<std::vector<double>>
make2Darray(size_t nx, size_t ny);
void
copy2darray(const std::vector<double>& flatArray, std::vector<std::vector<double>>& data);

std::vector<std::vector<std::vector<double>>>
make3Darray(size_t nx, size_t ny, size_t nz);
void
copy3darray(const std::vector<double>& flatArray, std::vector<std::vector<std::vector<double>>>& data);
