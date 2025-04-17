#pragma once

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>

#include <string>
#include <vector>

vtkm::FloatDefault scalarField1DEval(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& x,
                                     const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& y,
                                     const std::vector<vtkm::FloatDefault>& vals);

vtkm::FloatDefault scalarField1DEval(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& x,
                                     const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& y,
                                     const vtkm::FloatDefault& val);
/*
vtkm::FloatDefault scalarField2DEval(const vtkm::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const std::vector<vtkm::Vec3f>& pts);
*/

vtkm::FloatDefault scalarField2DEval(const vtkm::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const vtkm::Vec3f& pt);

vtkm::FloatDefault scalarField3DEval(const vtkm::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const vtkm::Vec3f& pt);

vtkm::FloatDefault scalarField3DEvalFwdAvg(const vtkm::cont::DataSet& dataset,
                                            const std::string& fieldname,
                                            const vtkm::Vec3f& pt);

/*
vtkm::FloatDefault scalarField3DEval(const vtkm::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const std::vector<vtkm::Vec3f>& pts);
*/
