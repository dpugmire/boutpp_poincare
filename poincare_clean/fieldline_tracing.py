import numpy as np
from scipy.interpolate import interp1d, interp2d
from netCDF4 import Dataset
import utils
from scipy.interpolate import RectBivariateSpline


# RK4 Field-Line Tracing Function
# the interpoloation for steps 2,3,4 are a little off the values in matlab.
def RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_pm1, dzdy_pm1, direction, nypf1, nypf2):
    hh = 0.5
    h6 = 1 / 6.0

    # Need half-step and full-step info
    if direction == 1:
        dxdyp = dxdy[:, yStart, :]
        dzdyp = dzdy[:, yStart, :]
        utils.write_array_to_file(dxdyp, 'dxdyp_1')
        utils.write_array_to_file(dzdyp, 'dzdyp_1')

        if region == 0 and yStart == nypf2:
            dxdyn = dxdy_pm1
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyn = dzdy_pm1
            dzdyh = 0.5 * (dzdyp + dzdyn)
        else:
            dxdyn = dxdy[:, yStart + 1, :]
            dxdyh = 0.5 * (dxdy[:, yStart, :] + dxdy[:, yStart + 1, :])
            dzdyn = dzdy[:, yStart + 1, :]
            dzdyh = 0.5 * (dzdy[:, yStart, :] + dzdy[:, yStart + 1, :])
    elif direction == -1:
        dxdyp = dxdy[:, yStart, :]
        dzdyp = dzdy[:, yStart, :]
        if region == 0 and yStart == nypf1 + 1:
            dxdyn = dxdy_pm1
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyn = dzdy_pm1
            dzdyh = 0.5 * (dzdyp + dzdyn)
        else:
            dxdyn = dxdy[:, yStart - 1, :]
            dxdyh = 0.5 * (dxdy[:, yStart, :] + dxdy[:, yStart - 1, :])
            dzdyn = dzdy[:, yStart - 1, :]
            dzdyh = 0.5 * (dzdy[:, yStart, :] + dzdy[:, yStart - 1, :])
    else:
        raise ValueError("Check direction parameter setting!")

    utils.write_array_to_file(dxdyn, 'dxdyn')
    utils.write_array_to_file(dxdyh, 'dxdyh')
    utils.write_array_to_file(dzdyn, 'dzdyn')
    utils.write_array_to_file(dzdyh, 'dzdyh')

    # Pathing last z point
    dxdyp = np.pad(dxdyp, ((0, 0), (0, 1)), 'edge')
    dxdyn = np.pad(dxdyn, ((0, 0), (0, 1)), 'edge')
    dxdyh = np.pad(dxdyh, ((0, 0), (0, 1)), 'edge')
    dzdyp = np.pad(dzdyp, ((0, 0), (0, 1)), 'edge')
    dzdyn = np.pad(dzdyn, ((0, 0), (0, 1)), 'edge')
    dzdyh = np.pad(dzdyh, ((0, 0), (0, 1)), 'edge')

    utils.write_array_to_file(dxdyp, 'dxdyp_')
    utils.write_array_to_file(dxdyn, 'dxdyn_')
    utils.write_array_to_file(dxdyh, 'dxdyh_')
    utils.write_array_to_file(dzdyp, 'dzdyp_')
    utils.write_array_to_file(dzdyn, 'dzdyn_')
    utils.write_array_to_file(dzdyh, 'dzdyh_')


    print(xarray.shape, zarray.shape, dxdyp.T.shape, dxdyp.shape)
    spline_dxdyp = RectBivariateSpline(xarray, zarray, dxdyp)
    spline_dzdyp = RectBivariateSpline(xarray, zarray, dzdyp)
    spline_dxdyh = RectBivariateSpline(xarray, zarray, dxdyh)
    spline_dzdyh = RectBivariateSpline(xarray, zarray, dzdyh)
    spline_dxdyn = RectBivariateSpline(xarray, zarray, dxdyn)
    spline_dzdyn = RectBivariateSpline(xarray, zarray, dzdyn)

    # First step
    print('xzStart= %12.10e %12.10e\n' %(xStart, zStart))
    #dxdy1 = interp2d(xarray, zarray, dxdyp.T, kind='cubic')(xStart, zStart)
    #dzdy1 = interp2d(xarray, zarray, dzdyp.T, kind='cubic')(xStart, zStart)
    dxdy1 = spline_dxdyp(xStart, zStart)
    dzdy1 = spline_dzdyp(xStart, zStart)
    x1 = xStart + direction * hh * dxdy1
    z1 = zStart + direction * hh * dzdy1
    print('dxdy1/dzdy1= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy1, dzdy1, x1,z1))

    # Second step
    #dxdy2 = interp2d(xarray, zarray, dxdyh.T, kind='cubic')(x1, np.mod(z1, 2 * np.pi))
    #dzdy2 = interp2d(xarray, zarray, dzdyh.T, kind='cubic')(x1, np.mod(z1, 2 * np.pi))
    dxdy2 = spline_dxdyh(x1, np.mod(z1, 2 * np.pi))
    dzdy2 = spline_dzdyh(x1, np.mod(z1, 2 * np.pi))
    x2 = xStart + direction * hh * dxdy2
    z2 = zStart + direction * hh * dzdy2
    print('dxdy2/dzdy2= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy2, dzdy2, x2,z2))

    # Third step
    #dxdy3 = interp2d(xarray, zarray, dxdyh.T, kind='cubic')(x2, np.mod(z2, 2 * np.pi))
    #dzdy3 = interp2d(xarray, zarray, dzdyh.T, kind='cubic')(x2, np.mod(z2, 2 * np.pi))
    dxdy3 = spline_dxdyh(x2, np.mod(z2, 2 * np.pi))
    dzdy3 = spline_dzdyh(x2, np.mod(z2, 2 * np.pi))
    x3 = xStart + direction * dxdy3
    z3 = zStart + direction * dzdy3
    print('dxdy3/dzdy3= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy3, dzdy3, x3,z3))

    # Fourth step
    #dxdy4 = interp2d(xarray, zarray, dxdyn.T, kind='cubic')(x3, np.mod(z3, 2 * np.pi))
    #dzdy4 = interp2d(xarray, zarray, dzdyn.T, kind='cubic')(x3, np.mod(z3, 2 * np.pi))
    dxdy4 = spline_dxdyn(x3, np.mod(z3, 2 * np.pi))
    dzdy4 = spline_dzdyn(x3, np.mod(z3, 2 * np.pi))
    # Accumulate increments with proper weights
    xEnd = xStart + direction * h6 * (dxdy1 + 2 * dxdy2 + 2 * dxdy3 + dxdy4)
    zEnd = zStart + direction * h6 * (dzdy1 + 2 * dzdy2 + 2 * dzdy3 + dzdy4)
    print('dxdy4/dzdy4= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy4, dzdy4, xEnd,zEnd))

    return xEnd, zEnd

def printEval(arr, i,j,k) :
    i = i-1
    j = j-1
    k = k-1
    print('arr[%d,%d,%d]= %12.10e\n'%(i,j,k, arr[i,j,k]))

# Reading NetCDF data
def read_netcdf_data(filename, varname):
    with Dataset(filename, 'r') as nc:
        data = nc.variables[varname][:]
    return data

# Parameters
xind = 175-1


nx = 260
ny = 128
nz = 256
zperiod = 1

jyseps1_1 = 16
jyseps2_2 = 112
ixsep1 = 195
ixsep2 = 260
ixsep = 195
nypf1 = 16
nypf2 = 112

direction = 1
nzG = nz * zperiod  # Full torus
zmin = 0.0
zmax = 2 * np.pi
dz = (zmax - zmin) / nzG
print('***********************\n   drp: might be some +1 issues here.\n\n')
ziarray = np.arange(0, nzG + 1)  # MATLAB's 1:nzG+1
zarray = ziarray * dz
xiarray = np.arange(1, nx + 1)

# Example: Load data from NetCDF
filename = 'stuff.nc'
rxy = np.transpose(read_netcdf_data(filename, 'rxy'))
zxy = np.transpose(read_netcdf_data(filename, 'zxy'))
psixy = np.transpose(read_netcdf_data(filename, 'psixy'))
print('psixy: ', psixy.shape)

dxdy = read_netcdf_data(filename, 'dxdy')
utils.write_array_to_file(dxdy, 'dxdy_0')
dxdy = np.transpose(dxdy)
utils.write_array_to_file(dxdy, 'dxdy_0T')
print('dxdy: ', dxdy.shape)
printEval(dxdy, 124, 80, 102)
printEval(dxdy, 102, 80, 124)
printEval(dxdy, 80, 124, 102)
printEval(dxdy, 14, 100, 28)

dzdy = np.transpose(read_netcdf_data(filename, 'dzdy'))
utils.write_array_to_file(dzdy, 'dzdy_0')
dxdy_p1 = np.transpose(read_netcdf_data(filename, 'dxdy_p1'))
dzdy_p1 = np.transpose(read_netcdf_data(filename, 'dzdy_p1'))
dxdy_m1 = np.transpose(read_netcdf_data(filename, 'dxdy_m1'))
dzdy_m1 = np.transpose(read_netcdf_data(filename, 'dzdy_m1'))

if ixsep2 < nx:
    divertor = 2  # Double null
    print("\tDouble null configuration")
elif ixsep1 < nx:
    divertor = 1  # Single null
    print("\tSingle null configuration")
    ixsep = ixsep1
    nypf1 = jyseps1_1
    nypf2 = jyseps2_2
else:
    divertor = 0  # Circular
    ixsep = nx
    nypf1 = 0
    nypf2 = ny
    print("\tCircular configuration")

# Determine maximum in the last row of rxy
jyomp = np.argmax(rxy[-1, :])
xarray = psixy[:, jyomp]  # Extract the jyomp-th column from psixy
xMin = xarray.min()
xMax = xarray.max()
iturn = 1
it = 0

print('psixy= ', psixy.shape, ' xi/jy=', xind, jyomp)
xStart = psixy[xind, jyomp]
yyy = jyomp
yStart = jyomp
zzz = 0  # Adjust for Python indexing
zStart = zarray[zzz]
print('zarray= ', zarray)
print(f"FieldLine: startXYZ({xind} {jyomp}) = {xStart} {yStart} {zStart}")

if xind < ixsep + 0.5:
    region = 0  # Closed flux surface
    if yStart < nypf1 + 1 or yStart > nypf2:
        region = 2  # PFR
else:
    region = 1  # SOL

yind = yStart
zind = interp1d(zarray, ziarray, kind='linear')(zStart)
print('\tline 0 started at indices (%12.10f,%12.10f,%12.10f), region= %d\n' %(xind,yind,zind, region))

for iy in range(1, ny):
    print(f"iy = {iy}")

    if region == 0 and nypf1 < yStart < nypf2 + 1:  # In CFR
        if direction == 1:
            xEnd, zEnd = RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_p1, dzdy_p1, 1, nypf1, nypf2)
            yEnd = yStart + 1
        elif direction == -1:
            xEnd, zEnd = RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_m1, dzdy_m1, -1, nypf1, nypf2)
            yEnd = yStart - 1

        print('  step: xyz: %12.10e %12.10e %12.10e --> xz: %12.10e %12.10e' % (xStart, yStart, zStart, xEnd, zEnd))

        # Check the end of the fieldline
        if xEnd > xMax:  # Outer boundary
            print(f"\tstarting xind={xind}, line {iline} reaches outer boundary")
            region = 12
        elif xEnd < xMin:  # Inner boundary
            print(f"\tstarting xind={xind}, line {iline} reaches inner boundary")
            region = 11
        else:  # Remains in simulation domain
            xind = interp1d(xarray, xiarray, kind='linear')(xEnd)
            if xind > ixsep1 + 0.5:
                region = 1
                print(f"\tending xind={xind}, line {iline} enters the SOL")

        # Twist-shift at the branch cut
        if direction == 1 and yStart == nypf2 and region == 0:
            shiftangle = interp1d(xiarray, sa, kind='linear')(xind)
            zEnd += shiftangle
            yEnd = nypf1 + 1

        if direction == -1 and yStart == nypf1 + 1 and region == 0:
            shiftangle = interp1d(xiarray, sa, kind='linear')(xind)
            zEnd -= shiftangle
            yEnd = nypf2

        # Relabel toroidal location
        if zEnd < zmin or zEnd > zmax:
            zEnd = zEnd % zmax

        zind = interp1d(zarray, ziarray, kind='linear')(zEnd)
        print(f"    record: {iturn} ({xEnd} {yEnd} {zEnd}) r={region}")

        it += 1

        # Update start point
        xStart = xEnd
        yStart = yEnd
        zStart = zEnd