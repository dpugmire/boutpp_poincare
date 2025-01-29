#include "NetCDFLoader.h"
#include <iostream>
#include <stdexcept>

std::vector<std::vector<double>>
make2Darray(size_t nx, size_t ny)
{
    std::vector<std::vector<double>> arr;
    arr.resize(nx);
    for (size_t i = 0; i < nx; i++)
        arr[i].resize(ny, 0.0);

    return arr;
}

std::vector<std::vector<std::vector<double>>>
make3Darray(size_t nx, size_t ny, size_t nz)
{
    std::vector<std::vector<std::vector<double>>> arr;
    arr.resize(nx);
    for (size_t i = 0; i < nx; i++)
    {
        arr[i].resize(ny);
        for (size_t j = 0; j < ny; j++)
            arr[i][j].resize(nz, 0.0);
    }

    return arr;
}

void
copy2darray(const std::vector<double>& flatArray, std::vector<std::vector<double>>& data)
{
    size_t nx = data.size();
    size_t ny = data[0].size();

    if (nx * ny != flatArray.size())
        throw std::runtime_error("Error. 2D flat array is the wrong size.");

    size_t idx = 0;
    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            data[i][j] = flatArray[idx++];
}

void
copy3darray(const std::vector<double>& flatArray, std::vector<std::vector<std::vector<double>>>& data)
{
    size_t nx = data.size();
    size_t ny = data[0].size();
    size_t nz = data[0][0].size();

    if (nx * ny * nz != flatArray.size())
        throw std::runtime_error("Error. 3D flat array is the wrong size.");

    size_t idx = 0;
    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
                data[i][j][k] = flatArray[idx++];
}

// Read 1D variable
std::vector<double>
NetCDFLoader::read1DVariable(const std::string& varName)
{
    int varid;
    int status = nc_inq_varid(this->ncid, varName.c_str(), &varid);
    if (status != NC_NOERR)
        throw std::runtime_error("Error getting variable ID: " + std::string(nc_strerror(status)));

    // Get dimensions
    size_t dimSize;
    getVariableDimensions(varid, &dimSize, 1);

    // Allocate storage
    std::vector<double> data(dimSize);

    // Read data
    status = nc_get_var_double(this->ncid, varid, data.data());
    if (status != NC_NOERR)
        throw std::runtime_error("Error reading 1D variable: " + std::string(nc_strerror(status)));

    return data;
}

    // Read 2D variable
    std::vector<std::vector<double>>
    NetCDFLoader::read2DVariable(const std::string& varName, bool transpose)
    {
        int varid;
        nc_type varType;

        // Get the variable ID
        if (nc_inq_varid(ncid, varName.c_str(), &varid) != NC_NOERR)
        {
            throw std::runtime_error("Error getting variable ID for: " + varName);
        }

        // Get the variable type
        if (nc_inq_var(ncid, varid, nullptr, &varType, nullptr, nullptr, nullptr) != NC_NOERR)
        {
            throw std::runtime_error("Error inquiring variable type for: " + varName);
        }

        // Get the number of dimensions
        int numDims;
        if (nc_inq_varndims(ncid, varid, &numDims) != NC_NOERR)
        {
            throw std::runtime_error("Error getting number of dimensions for: " + varName);
        }

        if (numDims != 2)
        {
            throw std::runtime_error("Variable " + varName + " is not 2D.");
        }

        // Get dimension sizes
        int dimIds[2];
        if (nc_inq_vardimid(ncid, varid, dimIds) != NC_NOERR)
        {
            throw std::runtime_error("Error getting dimension IDs for: " + varName);
        }

        size_t dimSizes[2];
        for (int i = 0; i < 2; ++i)
        {
            if (nc_inq_dimlen(ncid, dimIds[i], &dimSizes[i]) != NC_NOERR)
            {
                throw std::runtime_error("Error getting dimension size for: " + varName);
            }
        }

        // Allocate the 2D vector
        auto data = make2Darray(dimSizes[0], dimSizes[1]);

        // Read the variable data
        if (varType == NC_DOUBLE)
        {
            std::vector<double> flatData(dimSizes[0] * dimSizes[1], 0.0);
            if (nc_get_var_double(ncid, varid, flatData.data()) != NC_NOERR)
                throw std::runtime_error("Error reading variable: " + varName);
            copy2darray(flatData, data);
        }
        else
            throw std::runtime_error("Unsupported variable type for: " + varName);

        // Transpose if needed
        if (transpose)
            data = this->transpose2D(data);

        return data;
    }

    // Read 3D variable
    std::vector<std::vector<std::vector<double>>>
    NetCDFLoader::read3DVariable(const std::string& varName, bool permute)
     {
        int varid;
        nc_type varType;

        // Get the variable ID
        if (nc_inq_varid(ncid, varName.c_str(), &varid) != NC_NOERR)
        {
            throw std::runtime_error("Error getting variable ID for: " + varName);
        }

        // Get the variable type
        if (nc_inq_var(ncid, varid, nullptr, &varType, nullptr, nullptr, nullptr) != NC_NOERR)
        {
            throw std::runtime_error("Error inquiring variable type for: " + varName);
        }

        // Get the number of dimensions
        int numDims;
        if (nc_inq_varndims(ncid, varid, &numDims) != NC_NOERR)
        {
            throw std::runtime_error("Error getting number of dimensions for: " + varName);
        }

        if (numDims != 3)
        {
            throw std::runtime_error("Variable " + varName + " is not 3D.");
        }

        // Get dimension sizes
        int dimIds[3];
        if (nc_inq_vardimid(ncid, varid, dimIds) != NC_NOERR)
        {
            throw std::runtime_error("Error getting dimension IDs for: " + varName);
        }

        size_t dimSizes[3];
        for (int i = 0; i < 3; ++i)
        {
            if (nc_inq_dimlen(ncid, dimIds[i], &dimSizes[i]) != NC_NOERR)
            {
                throw std::runtime_error("Error getting dimension size for: " + varName);
            }
        }

        // Allocate the 3D vector
        auto data = make3Darray(dimSizes[0], dimSizes[1], dimSizes[2]);
        double val0 = data[123][79][101];
        double val1 = data[192][47][200];

        // Read the variable data
        if (varType == NC_DOUBLE)
        {
            std::vector<double> flatData(dimSizes[0] * dimSizes[1] * dimSizes[2], 0.0);
            if (nc_get_var_double(ncid, varid, flatData.data()) != NC_NOERR)
                throw std::runtime_error("Error reading variable: " + varName);
            copy3darray(flatData, data);
        }
        else
            throw std::runtime_error("Unsupported variable type for: " + varName);

        // Permute if needed
        if (permute)
        {
            std::vector<std::vector<std::vector<double>>> tmp(dimSizes[2],
                std::vector<std::vector<double>>(dimSizes[1], std::vector<double>(dimSizes[0], 0.0)));

            for (size_t i = 0; i < dimSizes[0]; ++i)
            {
                for (size_t j = 0; j < dimSizes[1]; ++j)
                {
                    for (size_t k = 0; k < dimSizes[2]; ++k)
                    {
                        tmp[k][j][i] = data[i][j][k];
                    }
                }
            }
            data = tmp;
        }
        val0 = data[123][79][101];
        val1 = data[192][47][200];

        std::cout<<"Read: "<<varName<<" ("<<dimSizes[0]<<" "<<dimSizes[1]<<" "<<dimSizes[2]<<")";
        std::cout<<" --> "<<data.size()<<" "<<data[0].size()<<" "<<" "<<data[0][0].size()<<std::endl;
        std::cout<<" val= "<<val0<<std::endl;
        std::cout<<" val= "<<val1<<std::endl;
        return data;
    }
