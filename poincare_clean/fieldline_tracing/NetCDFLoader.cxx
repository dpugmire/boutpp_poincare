#include "NetCDFLoader.h"
#include <iostream>
#include <stdexcept>

// Read 1D variable
std::vector<float>
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
    std::vector<float> data(dimSize);

    // Read data
    status = nc_get_var_float(this->ncid, varid, data.data());
    if (status != NC_NOERR)
        throw std::runtime_error("Error reading 1D variable: " + std::string(nc_strerror(status)));

    return data;
}

    // Read 2D variable
    std::vector<std::vector<float>>
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
        std::vector<std::vector<float>> data(dimSizes[0], std::vector<float>(dimSizes[1], 0.0f));

        // Read the variable data
        if (varType == NC_FLOAT)
        {
            if (nc_get_var_float(ncid, varid, &data[0][0]) != NC_NOERR)
            {
                throw std::runtime_error("Error reading variable: " + varName);
            }
        }
        else if (varType == NC_DOUBLE)
        {
            // Handle conversion from double to float
            std::vector<double> tempData(dimSizes[0] * dimSizes[1]);
            if (nc_get_var_double(ncid, varid, tempData.data()) != NC_NOERR)
            {
                throw std::runtime_error("Error reading variable: " + varName);
            }

            size_t idx = 0;
            for (size_t i = 0; i < dimSizes[0]; ++i)
            {
                for (size_t j = 0; j < dimSizes[1]; ++j)
                {
                    data[i][j] = static_cast<float>(tempData[idx++]);
                }
            }
        }
        else
        {
            throw std::runtime_error("Unsupported variable type for: " + varName);
        }

        // Transpose if needed
        if (transpose)
            data = this->transpose2D(data);

        return data;
    }

    // Read 3D variable
    std::vector<std::vector<std::vector<float>>>
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
        std::vector<std::vector<std::vector<float>>> data(dimSizes[0],
            std::vector<std::vector<float>>(dimSizes[1], std::vector<float>(dimSizes[2], 0.0f)));

        // Read the variable data
        if (varType == NC_FLOAT)
        {
            if (nc_get_var_float(ncid, varid, &data[0][0][0]) != NC_NOERR)
            {
                throw std::runtime_error("Error reading variable: " + varName);
            }
        }
        else if (varType == NC_DOUBLE)
        {
            // Handle conversion from double to float
            std::vector<double> tempData(dimSizes[0] * dimSizes[1] * dimSizes[2]);
            if (nc_get_var_double(ncid, varid, tempData.data()) != NC_NOERR)
            {
                throw std::runtime_error("Error reading variable: " + varName);
            }

            size_t idx = 0;
            for (size_t i = 0; i < dimSizes[0]; ++i)
            {
                for (size_t j = 0; j < dimSizes[1]; ++j)
                {
                    for (size_t k = 0; k < dimSizes[2]; ++k)
                    {
                        data[i][j][k] = static_cast<float>(tempData[idx++]);
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error("Unsupported variable type for: " + varName);
        }


        // Permute if needed
        if (permute)
        {
            std::vector<std::vector<std::vector<float>>> tmp(dimSizes[2],
                std::vector<std::vector<float>>(dimSizes[1], std::vector<float>(dimSizes[0], 0.0f)));

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
        std::cout<<"Read: "<<varName<<" ("<<dimSizes[0]<<" "<<dimSizes[1]<<" "<<dimSizes[2]<<")";
        std::cout<<" --> "<<data.size()<<" "<<data[0].size()<<" "<<" "<<data[0][0].size()<<std::endl;
        std::cout<<" val= "<<data[49][99][49]<<std::endl;
        std::cout<<" val= "<<data[149][47][249]<<std::endl;
        return data;
    }