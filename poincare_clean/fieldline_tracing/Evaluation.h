#pragma once

#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/DataSet.h>

#include <string>
#include <vector>

viskores::FloatDefault scalarField1DEval(const viskores::cont::ArrayHandle<viskores::FloatDefault>& x,
                                     const viskores::cont::ArrayHandle<viskores::FloatDefault>& y,
                                     const std::vector<viskores::FloatDefault>& vals);

viskores::FloatDefault scalarField1DEval(const viskores::cont::ArrayHandle<viskores::FloatDefault>& x,
                                     const viskores::cont::ArrayHandle<viskores::FloatDefault>& y,
                                     const viskores::FloatDefault& val);
/*
viskores::FloatDefault scalarField2DEval(const viskores::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const std::vector<viskores::Vec3f>& pts);
*/

viskores::FloatDefault scalarField2DEval(const viskores::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const viskores::Vec3f& pt);

viskores::FloatDefault scalarField3DEval(const viskores::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const viskores::Vec3f& pt);

viskores::FloatDefault scalarField3DEvalFwdAvg(const viskores::cont::DataSet& dataset,
                                            const std::string& fieldname,
                                            const viskores::Vec3f& pt);

/*
viskores::FloatDefault scalarField3DEval(const viskores::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const std::vector<viskores::Vec3f>& pts);
*/
