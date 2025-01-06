import numpy as np
from scipy.interpolate import interp1d, interp2d
from netCDF4 import Dataset
import utils
from scipy.interpolate import RectBivariateSpline

def INTERP(X,Y,val) :
    idx = next((i for i in range(len(X) - 1) if X[i] <= val <= X[i+1]), None)
    
    # Handle out-of-bounds x_i
    if idx is None:
        raise ValueError("x_i is out of the interpolation range.")
    
    # Linear interpolation formula
    x1, x2 = X[idx], X[idx + 1]
    y1, y2 = Y[idx], Y[idx + 1]
    y_i = y1 + (y2 - y1) * (val - x1) / (x2 - x1)
    return y_i


### matlab structure
'''
for iline in range(lines)
  set xyzStart
  it = 0
  region = set_region(xind)
  yind=ystart, zind = interp
  set region from diverter and ystart.

  while region < 10 and it < nturns
    for iy in range(ny) :
        # assume region == 0 and yStart in nypf1,2
        xEnd,zEnd = RK4
        yEnd = yStart+1

        region = set_region(xEnd)

        #update zEnd,yEnd at branch cut. twist-shift
        #relabel zEnd if needed.

        zind = interp1(zarray, ziarray, zEnd)
        it++

        #record xind, yEnd, zind, region

        start = end



'''


# RK4 Field-Line Tracing Function
# the interpoloation for steps 2,3,4 are a little off the values in matlab.
def RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_pm1, dzdy_pm1, direction, nypf1, nypf2):
    hh = 0.5
    h6 = 1 / 6.0

    dumpFiles = False
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
            #dumpFiles = True
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

    if dumpFiles :
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

    if dumpFiles :
        utils.write_array_to_file(dxdyp, 'dxdyp_')
        utils.write_array_to_file(dxdyn, 'dxdyn_')
        utils.write_array_to_file(dxdyh, 'dxdyh_')
        utils.write_array_to_file(dzdyp, 'dzdyp_')
        utils.write_array_to_file(dzdyn, 'dzdyn_')
        utils.write_array_to_file(dzdyh, 'dzdyh_')


    #print(xarray.shape, zarray.shape, dxdyp.T.shape, dxdyp.shape)
    spline_dxdyp = RectBivariateSpline(xarray, zarray, dxdyp)
    spline_dzdyp = RectBivariateSpline(xarray, zarray, dzdyp)
    spline_dxdyh = RectBivariateSpline(xarray, zarray, dxdyh)
    spline_dzdyh = RectBivariateSpline(xarray, zarray, dzdyh)
    spline_dxdyn = RectBivariateSpline(xarray, zarray, dxdyn)
    spline_dzdyn = RectBivariateSpline(xarray, zarray, dzdyn)

    # First step
    #print('xzStart= %12.10e %12.10e\n' %(xStart, zStart))
    #dxdy1 = interp2d(xarray, zarray, dxdyp.T, kind='cubic')(xStart, zStart)
    #dzdy1 = interp2d(xarray, zarray, dzdyp.T, kind='cubic')(xStart, zStart)
    dxdy1 = spline_dxdyp(xStart, zStart)
    dzdy1 = spline_dzdyp(xStart, zStart)
    x1 = xStart + direction * hh * dxdy1
    z1 = zStart + direction * hh * dzdy1
    #print('dxdy1/dzdy1= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy1, dzdy1, x1,z1))

    # Second step
    #dxdy2 = interp2d(xarray, zarray, dxdyh.T, kind='cubic')(x1, np.mod(z1, 2 * np.pi))
    #dzdy2 = interp2d(xarray, zarray, dzdyh.T, kind='cubic')(x1, np.mod(z1, 2 * np.pi))
    dxdy2 = spline_dxdyh(x1, np.mod(z1, 2 * np.pi))
    dzdy2 = spline_dzdyh(x1, np.mod(z1, 2 * np.pi))
    x2 = xStart + direction * hh * dxdy2
    z2 = zStart + direction * hh * dzdy2
    #print('dxdy2/dzdy2= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy2, dzdy2, x2,z2))

    # Third step
    #dxdy3 = interp2d(xarray, zarray, dxdyh.T, kind='cubic')(x2, np.mod(z2, 2 * np.pi))
    #dzdy3 = interp2d(xarray, zarray, dzdyh.T, kind='cubic')(x2, np.mod(z2, 2 * np.pi))
    dxdy3 = spline_dxdyh(x2, np.mod(z2, 2 * np.pi))
    dzdy3 = spline_dzdyh(x2, np.mod(z2, 2 * np.pi))
    x3 = xStart + direction * dxdy3
    z3 = zStart + direction * dzdy3
    #print('dxdy3/dzdy3= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy3, dzdy3, x3,z3))

    # Fourth step
    #dxdy4 = interp2d(xarray, zarray, dxdyn.T, kind='cubic')(x3, np.mod(z3, 2 * np.pi))
    #dzdy4 = interp2d(xarray, zarray, dzdyn.T, kind='cubic')(x3, np.mod(z3, 2 * np.pi))
    dxdy4 = spline_dxdyn(x3, np.mod(z3, 2 * np.pi))
    dzdy4 = spline_dzdyn(x3, np.mod(z3, 2 * np.pi))
    # Accumulate increments with proper weights
    xEnd = xStart + direction * h6 * (dxdy1 + 2 * dxdy2 + 2 * dxdy3 + dxdy4)
    zEnd = zStart + direction * h6 * (dzdy1 + 2 * dzdy2 + 2 * dzdy3 + dzdy4)
    #print('dxdy4/dzdy4= %12.10e %12.10e xz= %12.10e %12.10e' % (dxdy4, dzdy4, xEnd,zEnd))

    return (xEnd[0][0], zEnd[0][0])

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
xind = 195-1
nlines = 1
nturns = 2


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
stepsF = open('steps.p.txt', 'w', buffering=1)
#stepsF.write('ID,YI,iTURN,xStart,yStart,zStart,xind,yind,zind,region\n')
trajF = open('traj.p.txt', 'w', buffering=1)
trajF.write('ID, ITER, X, Y, Z\n')

direction = 1
nzG = nz * zperiod  # Full torus
zmin = 0.0
zmax = 2 * np.pi
dz = (zmax - zmin) / nzG
print('***********************\n   drp: might be some +1 issues here.\n\n')
ziarray = np.arange(0, nzG + 1)  # MATLAB's 1:nzG+1
zarray = ziarray * dz
xiarray = np.arange(0, nx)

# Example: Load data from NetCDF
filename = 'stuff.nc'
rxy = np.transpose(read_netcdf_data(filename, 'rxy'))
zxy = np.transpose(read_netcdf_data(filename, 'zxy'))
psixy = np.transpose(read_netcdf_data(filename, 'psixy'))
zShift = np.transpose(read_netcdf_data(filename, 'zShift'))
shiftAngle = read_netcdf_data(filename, 'shiftAngle')
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
#print('zarray= ', zarray)
print(f"FieldLine: startXYZ({xind} {jyomp}) = {xStart} {yStart} {zStart}")

region = -1
if xind < ixsep + 0.5:
    region = 0  # Closed flux surface
    if yStart < nypf1 + 1 or yStart > nypf2:
        region = 2  # PFR
else:
    region = 1  # SOL

yind = yStart
zind = interp1d(zarray, ziarray, kind='linear')(zStart)
print('\tline 0 started at indices (%12.10f,%12.10f,%12.10f), region= %d\n' %(xind,yind,zind, region))

VALUES = np.round(np.linspace(0, 255, nlines))
VALUES = list(map(int, VALUES))
COUNTER = 0
NUM_ROUND = 3
nturns = 5
LINES = [194]

# Loop over lines
for iline in LINES:
    xind = iline
    xStart = psixy[xind, jyomp]
    yyy = jyomp
    yStart = jyomp
    zzz = 0  # Adjust for Python indexing
    zStart = zarray[zzz]

    if xind < ixsep + 0.5:
        region = 0  # Closed flux surface
        if yStart < nypf1 + 1 or yStart > nypf2:
            region = 2  # PFR
        else:
            region = 1  # SOL
    yind = yStart
    #zind = interp1d(zarray, ziarray, fill_value="extrapolate")(zStart)
    zind = INTERP(zarray, ziarray, zStart)

    # Check for divertor starting points
    if direction == 1 and yStart == ny-1:
        region = 14
    elif direction == -1 and yStart == 0:
        region = 13

    stepsF.write('COUNTER= %d region= %d xyzind= %d %d %d\n' % (COUNTER, region, xind, yind, zind))

    it = 0
    iturn = 1
    
    # Determine the region
    if xind < ixsep1 + 0.5:
        region = 0  # Closed flux surface
        if yStart < nypf1 or yStart > nypf2:
            region = 2  # Private flux region
    else:
        region = 1  # Scrape-off layer


    print(f"Line {iline} started at indices ({xind}, {yind}, {zind})")
    
    # Check for divertor starting points
    if direction == 1 and yStart == ny:
        print(f"Line {iline} starts on the divertor.")
        region = 14
        traj[5, it] = 0.0
    elif direction == -1 and yStart == 1:
        print(f"Line {iline} starts on the divertor.")
        region = 13
        traj[5, it] = 0.0

    # Record field-line location info in trajectory array
    #traj[:, it] = [1, xind, yStart, zind, region, 0.0, zStart]
    XYZVals = [(xind, yStart, zind)]

    while region < 10 and iturn < nturns:
        stepsF.write('COUNTER= %d region= %d iturn= %d\n' % (COUNTER, region, iturn))
        if iturn % 50 == 1:
            print(f"\tLine {iline}, turn {iturn}/{nturns}...")

        # Start field-line tracing
        for iy in range(ny):
            if yStart +1 == dxdy.shape[1] :
                print('overflow of some kind..... Need to track this down.')
                break

            stepsF.write('COUNTER= %d region= %d iy= %d\n   xyzStart= %.8f %.8f %.8f\n   xyzInd= %.8f %.8f %.8f\n' % (COUNTER, region, iy, xStart, yStart, zStart, xind, yind, zind))

            if region == 0 and yStart >= nypf1 and yStart < nypf2+1 :
                # take FWD RK4 step.
                xEnd, zEnd = RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_p1, dzdy_p1, 1, nypf1, nypf2)
                yEnd = yStart + 1

                #traj[6, it + 1] = zEnd

                # Check where the field line ends
                if xEnd > xMax:
                    print(f"\tStarting xind={xind}, line {iline} reaches outer boundary.")
                    region = 12
                elif xEnd < xMin:
                    print(f"\tStarting xind={xind}, line {iline} reaches inner boundary.")
                    region = 11
                else:
                    #xind = interp1d(xarray, xiarray, fill_value="extrapolate")(xEnd)
                    xind = INTERP(xarray, xiarray, xEnd)
                    if xind > float(ixsep1) + 0.5 :
                        region = 1
                        print(f"\tEnding xind={xind}, line {iline} enters the SOL.")

                # Twist-shift at the branch cut
                if yStart == nypf2-1 and region == 0 :
                    #shiftangle = interp1d(xiarray, sa, fill_value="extrapolate")(xind)
                    shiftangle = INTERP(xiarray, shiftAngle, xind)
                    zEnd = zEnd + shiftangle
                    yEnd = nypf1

                # Relabel toroidal location
                if zEnd < zmin or zEnd > zmax:
                    zEnd = zEnd % zmax
                #zind = interp1d(zarray, ziarray, fill_value="extrapolate")(zEnd)
                zind = INTERP(zarray, ziarray, zEnd)
                XYZVals.append((xind, yEnd, zind))

                it = it+1
                COUNTER = COUNTER+1
                #traj[:, it] = [iturn, xind, yEnd, zind, region, 0.0, zEnd]
                
                # Update starting points
                xStart, yStart, zStart = xEnd, yEnd, zEnd

            elif region == 1 or region == 2 :
                print('region= ', region)
                stepsF.write('region 1 or 2\n')
                print('******* in SOL/PFR')
                #if direction == 1
                print('yStart= ', yStart, dxdy.shape)
                xEnd, zEnd = RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_p1, dzdy_p1, 1, nypf1, nypf2)
                yEnd = yStart + 1
                print('xxx ', it, xStart, yStart, zStart, xEnd, zEnd)
                print('  step_%d: xyz: %12.10e %12.10e %12.10e --> xz: %12.10e %12.10e' % (it, xStart, yStart, zStart, xEnd, zEnd))

                #if direction == 1 :
                xEnd, zEnd = RK4_FLT1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray, region, dxdy_p1, dzdy_p1, 1, nypf1, nypf2)
                #traj(7,it+1)=zEnd; % for better interpolation of puncture point

                #correct yEnd for two PFR cases
                if (direction == 1 and yStart == nypf1 and region == 2) :
                    yEnd = nypf2+1
                elif (direction == -1 and yStart == nypf2+1 and region == 2) :
                    yEnd = nypf1

                if xEnd > xMax :
                    print('Starting xind=%f, line %i reaches outer bndry' %(xind, iline))
                    stepsF.write('region 12: %.8f > %.8f\n' % (xEnd, xMax))
                    region = 12
                elif xEnd < xMin :
                    print('Starting xind=%f, line %i reaches inner bndry\n' % (xind,iline))
                    stepsF.write('region 11: %d < %d\n' % (xEnd, xMin))
                    region = 11
                else:
                    print('********** in end of field-line remains in simulation\n');
                    xind = INTERP(xarray, xiarray, xEnd)
                    stepsF.write('xind INTERP(%.8f)= %.8f\n' % (xEnd, xind))

                if xind < float(ixsep1)+0.5 and yEnd > nypf1 and yEnd < nypf2+1 :
                    if region != 0 :
                        print('\tending xind=%f, line %i enters the CFR\n' %(xind,iline))
                    stepsF.write('region 0: %.8f < %d AND %d > %d AND %d < %d\n' % (xind, ixsep1, yEnd, nypf1, yEnd, nypf2+1))
                    region = 0
                elif xind < float(ixsep1)+0.5 and (yEnd > nypf2-1 or yEnd < nypf1) :
                    if region != 2 :
                        print('\tending xind=%f, line %i enters the PFR\n' % (xind,iline))
                    region = 2
                    stepsF.write('region 2: %.8f < %d AND %d > %d OR %d < %d\n' %(xind, ixsep1,  yEnd, nypf2-1, yEnd, nypf1))

                zind = INTERP(zarray, ziarray, zEnd)
                stepsF.write('zind INTERP(%.8f)= %.8f\n' % (zEnd, zind))

                if direction == 1 and yEnd == ny :
                    print('\tstarting xind=%f, line %i reaches divertor\n' %(xind,iline))
                    stepsF.write('region 14: %d == %d\n' % (yEnd, ny))
                    region = 14
                elif direction == -1 and yEnd == 1 :
                    print('\tstarting xind=%f, line %i reaches divertor\n' %(xind,iline))
                    stepsF.write('region 13: %d == 1\n' % (yEnd))
                    region = 13

                zind = INTERP(zarray, ziarray, zEnd)
                XYZVals.append((xind, yEnd, zind))
                it = it+1
                COUNTER = COUNTER+1
                # Update starting points
                xStart, yStart, zStart = xEnd, yEnd, zEnd

            else:
                print("In other region. Need to fix this. region= ", region)
                shit()

        iturn = iturn + 1

    cnt = 0
    for xyz in XYZVals :
        xi = xyz[0]
        yi = xyz[1]
        zi = xyz[2]
        tmp = rxy[:,yi]

        rxyvalue = INTERP(xiarray, rxy[:,xyz[1]], xyz[0])
        zsvalue = INTERP(xiarray, zShift[:, xyz[1]], xyz[0])
        zvalue = INTERP(ziarray, zarray, xyz[2])
        x3d_tmp = rxyvalue * np.cos(zsvalue)
        y3d_tmp = rxyvalue * np.sin(zsvalue)
        x3d = x3d_tmp * np.cos(zvalue) - y3d_tmp * np.sin(zvalue)
        y3d = x3d_tmp * np.sin(zvalue) + y3d_tmp * np.cos(zvalue)
        z3d = INTERP(xiarray, zxy[:,xyz[1]], xyz[0])
        trajF.write('%d, %d, %f, %f, %f\n' % (iline, cnt, x3d, y3d, z3d))
        cnt = cnt+1
    # Record the maximum valid steps
    itmax = it
    #traj = traj[:, :itmax]
    
    # Save trajectory if needed
    #if save_traj:
     #   filename = f"./mat_traj/x{iline}_y{yStart}_z{zStart}_traj.npy"
      #  np.save(filename, traj)

