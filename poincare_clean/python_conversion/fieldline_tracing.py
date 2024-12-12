import netCDF4 as nc
import numpy as np
import scipy.io
from get_apar_sc import *
from get_apar_sn import *
from scipy.interpolate import interp1d
import utils
import os


nx,ny = None,None
rxy,zxy,psixy,rmag = None,None,None,None
zShift = None
ixsep1,ixsep2,nypf11,nypf12,nypf21,nypf22 = None,None,None,None,None,None
nypf1,nypf2 = None,None
bxy,btxy,bpxy,hthe,sinty, bxcvx,bxcvy,bxcvz = None,None,None,None,None,None,None,None
jpar0,dy,dy0,sa = None,None,None,None
diverter = 0
nz = 256
zperiod = 1
deltaix, ixoffset = 1, 1  # Line tracing starts at (deltaix*ilines+ixoffset, iy, iz)

# Fill the torus if needed
nzG = nz * zperiod  # Full torus extension
zmin, zmax = 0.0, 2 * np.pi
dz = (zmax - zmin) / nzG
ziarray = np.arange(1, nzG + 2)  # 1-based index
zarray = (ziarray - 1) * dz





def readData(fname) :
    global nx, ny, rxy, zxy, psixy, rmag, zShift
    global ixsep, ixsep1,ixsep2,nypf11, nypf12, nypf21, nypf22
    global nypf1, nypf2
    global bxy,btxy,bpxy,hthe,sinty, bxcvx,bxcvy,bxcvz
    global divertor
    global jpar0,dy,dy0,sa

    ds = nc.Dataset(fname, 'r')
    nx = ds.variables['nx'][0]
    ny = ds.variables['ny'][0]
    rxy = np.array(ds.variables['Rxy'])
    zxy = np.array(ds.variables['Zxy'])
    psixy = np.array(ds.variables['psixy'])
    rmag = np.array(ds.variables['rmag'])
    zShift = np.array(ds.variables['zShift'])
    ixsep1 = ds.variables['ixseps1'][0]
    ixsep2 = ds.variables['ixseps2'][0]
    if ixsep2 < nx:
        divertor = 2  # Double null
        print("\tDouble null configuration")
        nypf11 = ds.variables['jyseps1_1'][0] + 1
        nypf21 = ds.variables['jyseps2_1'][0] + 1
        nypf12 = ds.variables['jyseps1_2'][0] + 1
        nypf22 = ds.variables['jyseps2_2'][0] + 1
    elif ixsep1 < nx:
        divertor = 1  # Single null
        print("\tSingle null configuration")
        ixsep = ixsep1
        nypf1 = ds.variables['jyseps1_1'][0] + 1
        nypf2 = ds.variables['jyseps2_2'][0] + 1
    else:
        divertor = 0  # Circular configuration
        ixsep = nx
        nypf1, nypf2 = 0, ny
        print("\tCircular configuration")

    bxy = np.array(ds.variables['Bxy'])
    btxy = np.array(ds.variables['Btxy'])
    bpxy = np.array(ds.variables['Bpxy'])
    hthe = np.array(ds.variables['hthe'])
    sinty = np.array(ds.variables['sinty'])

    bxcvx = np.array(ds.variables['bxcvx'])
    bxcvy = np.array(ds.variables['bxcvy'])
    bxcvz = np.array(ds.variables['bxcvz'])

    jpar0 = np.array(ds.variables['Jpar0'])
    dy = np.array(ds.variables['dy'])
    dy0 = dy[0, 0]

    sa = ds.variables['ShiftAngle']



gridfile =  '../data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc'

readData(gridfile)

xiarray = np.arange(1, nx + 1)  # 1-based index
yiarray = np.arange(1, ny + 1)  # 1-based index
jyomp = np.argmax(rxy[-1, :])  # Select the last row with `-1`
xarray = psixy[:, jyomp]  # This selects all rows (:) and column `jyomp`

# Find the minimum and maximum values of `xarray`
xMin = np.min(xarray)
xMax = np.max(xarray)


# Additional calculations
dz = 2 * np.pi / zperiod / nz
nu = btxy * hthe / bpxy / rxy

# Calculate SOL, X-point locations, and poloidal angle theta
theta = np.zeros(ny)  # Initialize theta array

if divertor == 1:  # Single null configuration
    print(ny, nypf1)
    corex = rxy[0, nypf1:ny-nypf1]
    corey = zxy[0, nypf1:ny-nypf1]

    solx = rxy[0, :nypf1].tolist()
    soly = zxy[0, :nypf1].tolist()

    solx.extend(rxy[0, ny-nypf1:ny].tolist())
    soly.extend(zxy[0, ny-nypf1:ny].tolist())

    solx.extend(rxy[:, -1].tolist())
    soly.extend(zxy[:, -1].tolist())

    solx.extend(rxy[-1, ny-1::-1].tolist())
    soly.extend(zxy[-1, ny-1::-1].tolist())

    solx.extend(rxy[::-1, 0].tolist())
    soly.extend(zxy[::-1, 0].tolist())

    ##tmp = list(range(1, nypf1+1)) + list(range(ny-nypf1, nypf1, -1)) + list(range(ny-nypf1, ny))
    tmp = list(range(1, nypf1+1)) + list(range(ny-nypf1, nypf1, -1)) + list(range(ny-nypf1+1, ny+1))
    tmp = [x-1 for x in tmp]
    print('len(tmp)=',len(tmp), 'ixsep=', ixsep, 'rxy=',rxy.shape)
    print('len(rxy)= ', rxy.shape)
    print('tmp: ', tmp[0], tmp[1])
    #sepx = 0.5 * (rxy[ixsep, tmp] + rxy[ixsep+1, tmp])
    sepx = 0.5 * (rxy[ixsep-1, tmp] + rxy[ixsep, tmp])
    #sepy = 0.5 * (zxy[ixsep, tmp] + zxy[ixsep+1, tmp])
    sepy = 0.5 * (zxy[ixsep-1, tmp] + zxy[ixsep, tmp])

    center_x = 0.5 * (np.max(corex) + np.min(corex))
    center_y = 0.5 * (np.max(corey) + np.min(corey))
    print('center: ', center_x, center_y)

    # Reference point: X-point (theta = 0)
    print('ny= ', ny,' nypf1= ', nypf1)
    #print('sepy: ', sepy)
    #print('sepx: %12.10f %12.10f %12.10f' % (sepx[0], sepx[1], sepx[3]))
    #orig xpoint_x = 0.25 * (sepx[nypf1] + sepx[nypf1+1] + sepx[ny-nypf1] + sepx[ny-nypf1+1])
    #print('--xp=  %f %f %f %f' % (sepx[nypf1-1]), sepx[nypf1], sepx[ny-nypf1], sepx[ny-nypf1])
    A,B,C,D = sepx[nypf1-1],  sepx[nypf1],  sepx[ny-nypf1-1],  sepx[ny-nypf1]
    print('xp= ',      sepx[nypf1-1],  sepx[nypf1],  sepx[ny-nypf1-1],  sepx[ny-nypf1])
    print('xp= ', A,B,C,D)
    xpoint_x = 0.25 * (sepx[nypf1-1] + sepx[nypf1] + sepx[ny-nypf1-1] + sepx[ny-nypf1])

    #orig xpoint_y = 0.25 * (sepy[nypf1] + sepy[nypf1+1] + sepy[ny-nypf1] + sepy[ny-nypf1+1])
    xpoint_y = 0.25 * (sepy[nypf1-1] + sepy[nypf1] + sepy[ny-nypf1-1] + sepy[ny-nypf1])
    print('xpoint= %12.10f %12.10f' % (xpoint_x, xpoint_y))
    print('xpoint_x= ', (0.25*(A+B+C+D)))

    u = np.array([center_x - xpoint_x, center_y - xpoint_y, 0])

    # Calculate theta for each `iy`
    for iy in range(ny):
        v = np.array([center_x - rxy[0, iy], center_y - zxy[0, iy], 0])
        theta[iy] = np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))

    theta /= np.pi

    # Adjust theta values
    itheta = np.argmax(theta)
    theta[itheta:ny] = 2 - theta[itheta:ny]

    itheta = np.argmax(theta)
    if itheta != ny - 1:
        theta[itheta:ny] = 4 - theta[itheta:ny]

    # Shift theta to have reference at (1, nypf1+1)
    theta -= theta[nypf1]

elif divertor == 0:  # Shift-circular configuration
    center_x = 0.5 * (np.max(rxy[0, :]) + np.min(rxy[0, :]))
    center_y = 0.5 * (np.max(zxy[0, :]) + np.min(zxy[0, :]))

    # Reference point: (1,1)
    u = np.array([center_x - rxy[0, 0], center_y - zxy[0, 0], 0])

    for iy in range(ny):
        v = np.array([center_x - rxy[0, iy], center_y - zxy[0, iy], 0])
        theta[iy] = np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))

    theta /= np.pi

    # Adjust theta values
    itheta = np.argmax(theta)
    theta[itheta:ny] = 2 - theta[itheta:ny]

    itheta = np.argmax(theta)
    if itheta != ny - 1:
        theta[itheta:ny] = 4 - theta[itheta:ny]

else:
    print("\t\tConfiguration to be implemented!")



# Construct/patch closed flux surface region for a better Poincare plot
xiarray_cfr = np.arange(1, ixsep + 1, dtype=float)
yiarray_cfr = np.arange(nypf1, nypf2 + 1, dtype=float)  # Note: +1 to include nypf2
theta_cfr = theta[yiarray_cfr.astype(int) - 1]  # Adjust for Python 0-based indexing
theta_cfr[-1] = 2.0  # theta is pi-based

rxy_cfr = rxy[:ixsep, nypf1:nypf2]  # Slicing for rxy
rxy_cfr = np.hstack((rxy_cfr, rxy_cfr[:, :1]))  # Append the first column at the end

zxy_cfr = zxy[:ixsep, nypf1:nypf2]  # Slicing for zxy
zxy_cfr = np.hstack((zxy_cfr, zxy_cfr[:, :1]))  # Append the first column at the end

zs_cfr = zShift[:ixsep, nypf1:nypf2]
# Additional zShift info on the last point
zs_cfr = np.hstack((
    zs_cfr,
    0.5 * (nu[xiarray_cfr.astype(int) - 1, nypf1] + nu[xiarray_cfr.astype(int) - 1, nypf2]) * dy0 + zs_cfr[:, -1:]
))

# Calculate geometric coefficients -- refer to the note
A1 = rxy * bpxy * btxy / hthe
A2 = bxy**2
A3 = sinty * A1

JJ = 4 * np.pi * 1.e-7 * bpxy / hthe / (bxy**2) * jpar0
print('A3: ', A1[0,1], sinty[0,1], '= ', A3[0,1])
utils.write_array_to_file(rxy, 'rxy.p.txt')
utils.write_array_to_file(bpxy, 'bpxy.p.txt')
utils.write_array_to_file(btxy, 'btxy.p.txt')
utils.write_array_to_file(hthe, 'hthe.p.txt')
utils.write_array_to_file(sinty, 'sinty.p.txt')
utils.write_array_to_file(jpar0, 'jpar0.p.txt')
utils.write_array_to_file(A1, 'A1.p.txt')
utils.write_array_to_file(A2, 'A2.p.txt')
utils.write_array_to_file(A3, 'A3.p.txt')
utils.write_array_to_file(JJ, 'JJ.p.txt')

print('Loading perturbed field information ...')

# Initialize arrays for BOUT++ output
apar = np.zeros((nx, ny, nzG))
dapardx = np.zeros((nx, ny, nzG))
dapardy = np.zeros((nx, ny, nzG))
dapardz = np.zeros((nx, ny, nzG))

# Read RMP info on BOUT++ grid from a MATLAB file
rmp = scipy.io.loadmat('../data/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat')
utils.write_array_to_file(rmp['apar'], 'apar.p.txt')
print('****** Divertor= ', divertor)

# Handle configurations for divertor
if divertor == 0:
    apar0, dapardx0, dapardy0, dapardz0 = get_apar_sc(
        rmp['apar'], bxy, psixy, zShift, sa, sinty, dy0, dz, zperiod, 0, 1, 1
    )
elif divertor == 1:
    apar0, dapardx0, dapardy0, dapardz0 = get_apar_sn(
        rmp['apar'], bxy, psixy, zShift, sa, sinty, dy0, dz, ixsep, nypf1, nypf2, zperiod, 0, 0, 1
    )
else:
    print('\tConfiguration to be implemented!')

utils.write_array_to_file(apar0, 'apar0.p.txt')
utils.write_array_to_file(dapardx0, 'dapardx0.p.txt')
utils.write_array_to_file(dapardy0, 'dapardy0.p.txt')
utils.write_array_to_file(dapardz0, 'dapardz0.p.txt')

# Populate arrays for the full torus
for zp in range(1, zperiod + 1):
    start_idx = (zp - 1) * nzG // zperiod
    end_idx = zp * nzG // zperiod
    apar[:, :, start_idx:end_idx] = apar0
    dapardx[:, :, start_idx:end_idx] = dapardx0
    dapardy[:, :, start_idx:end_idx] = dapardy0
    dapardz[:, :, start_idx:end_idx] = dapardz0


utils.write_array_to_file(apar, 'apar.p.txt')
utils.write_array_to_file(dapardx, 'dapardx.p.txt')
utils.write_array_to_file(dapardy, 'dapardy.p.txt')
utils.write_array_to_file(dapardz, 'dapardz.p.txt')

# STEP 2: Field-line tracing

# Calculate perturbed field
b0dgy = bpxy / hthe

bdgx = np.zeros_like(apar)
bdgy = np.zeros_like(apar)
bdgz = np.zeros_like(apar)
dxdy = np.zeros_like(apar)
dzdy = np.zeros_like(apar)

for k in range(nzG):
    bdgx[:, :, k] = (1.0 / bxy) * (-A1 * dapardy[:, :, k] - A2 * dapardz[:, :, k]) + apar[:, :, k] * bxcvx
    bdgy[:, :, k] = (1.0 / bxy) * (A1 * dapardx[:, :, k] - A3 * dapardz[:, :, k]) + apar[:, :, k] * (bxcvy + JJ)
    bdgz[:, :, k] = (1.0 / bxy) * (A2 * dapardx[:, :, k] + A3 * dapardz[:, :, k]) + apar[:, :, k] * bxcvz

    dxdy[:, :, k] = bdgx[:, :, k] / (b0dgy + bdgy[:, :, k])
    dzdy[:, :, k] = bdgz[:, :, k] / (b0dgy + bdgy[:, :, k])

utils.write_array_to_file(dxdy, 'dxdy.p.txt')
print('nzG= ', nzG)
shit()


#write_array_to_file(bdgx, 'bdgx.p.txt')
print('bdgx: %12.10e' % bdgx[100,10,10])
shit()

# Additional grid point information near the branch-cut for the closed flux surface region
# Assuming zeroth-order continuity

dxdyt = dxdy[:, nypf1, :]
dxdyt = np.concatenate([dxdyt, dxdyt[:, :1]], axis=1)
dzdyt = dzdy[:, nypf1, :]
dzdyt = np.concatenate([dzdyt, dzdyt[:, :1]], axis=1)

# p1: First point to ny+1, twist-shift from the first point
dxdy_p1 = np.zeros((nx, nzG))
dzdy_p1 = np.zeros((nx, nzG))

for ix in range(ixsep):
    zarray_shift = np.mod(zarray[:nzG] + sa[ix], zmax)
    dxdy_p1[ix, :] = CubicSpline(zarray, dxdyt[ix, :])(zarray_shift)
    dzdy_p1[ix, :] = CubicSpline(zarray, dzdyt[ix, :])(zarray_shift)

# m1: Last point to minus 1, reverse twist-shift from the last point
dxdyt = dxdy[:, nypf2, :]
dxdyt = np.concatenate([dxdyt, dxdyt[:, :1]], axis=1)
dzdyt = dzdy[:, nypf2, :]
dzdyt = np.concatenate([dzdyt, dzdyt[:, :1]], axis=1)

dxdy_m1 = np.zeros((nx, nzG))
dzdy_m1 = np.zeros((nx, nzG))

for ix in range(ixsep):
    zarray_rshift = np.mod(zarray[:nzG] - sa[ix], zmax)
    dxdy_m1[ix, :] = CubicSpline(zarray, dxdyt[ix, :])(zarray_rshift)
    dzdy_m1[ix, :] = CubicSpline(zarray, dzdyt[ix, :])(zarray_rshift)

print('dxdy: ', dxdy[0,0,0])

print("Starting field-line tracing...")
print()

# Colormap for visualization
#cm = plt.get_cmap("jet")(np.linspace(0, 1, nlines))



def RK4_FLT1(*args):
    # Placeholder for the RK4_FLT1 function implementation
    pass

def save_trajectory(filepath, traj):
    """Save trajectory data as a .mat file."""
    sio.savemat(filepath, {'traj': traj})

def main():
    VALUES = np.round(np.linspace(1, 256, nlines))

    for iline in range(1, nlines + 1):
        # Pick starting points
        if testPoint == 1:
            xind = 175
        else:
            xind = VALUES[iline - 1]

        xStart = psixy[int(xind), jyomp]  # Note: jyomp doesn't matter here
        yyy = jyomp
        yStart = jyomp
        zzz = 1
        zStart = zarray[zzz - 1]

        # Initialize trajectory and puncture points arrays
        traj = np.zeros((7, nsteps))
        fl_x3d = np.zeros(nsteps)
        fl_y3d = np.zeros(nsteps)
        fl_z3d = np.zeros(nsteps)
        px = np.zeros(np_)
        py = np.zeros(np_)
        pz = np.zeros(np_)
        ptheta = np.zeros(np_)
        ppsi = np.zeros(np_)

        it = 0
        iturn = 1

        if xind < ixsep + 0.5:
            region = 0
            if yStart < nypf1 + 1 or yStart > nypf2:
                region = 2
        else:
            region = 1

        yind = yStart
        zind = interp1d(zarray, ziarray)(zStart)

        print(f"\tline {iline} started at indices ({xind},{yind},{zind}),")
        print(f"FieldLine: startXYZ({xind} {jyomp}) = {xStart} {yStart} {zStart}")

        # Divertor check
        if divertor == 1:
            if yStart == ny and direction == 1:
                print(f"\tline {iline} starts on the divertor.")
                region = 14
                traj[5, it] = 0.0
            elif yStart == 1 and direction == -1:
                print(f"\tline {iline} starts on the divertor.")
                region = 13
                traj[5, it] = 0.0

        traj[0, it] = 1
        traj[1, it] = xind
        traj[2, it] = yStart
        traj[3, it] = zind
        traj[4, it] = region
        traj[6, it] = zStart

        while region < 10 and iturn < nturns:
            if iturn % 5 == 1:
                print(f"\t\t line{iline}, turn {iturn}/{nturns} ...")

            # Field-line tracing
            for iy in range(ny - 1):
                if region == 0 and nypf1 < yStart < nypf2 + 1:
                    if direction == 1:
                        xEnd, zEnd = RK4_FLT1(
                            xStart, yStart, zStart, dxdy, dzdy, xarray, zarray,
                            region, dxdy_p1, dzdy_p1, 1, nypf1, nypf2
                        )
                        yEnd = yStart + 1
                    elif direction == -1:
                        xEnd, zEnd = RK4_FLT1(
                            xStart, yStart, zStart, dxdy, dzdy, xarray, zarray,
                            region, dxdy_m1, dzdy_m1, -1, nypf1, nypf2
                        )
                        yEnd = yStart - 1

                    traj[6, it + 1] = zEnd

                    if xEnd > xMax:
                        print(f"\tstarting xind={xind}, line {iline} reaches outer bndry")
                        region = 12
                    elif xEnd < xMin:
                        print(f"\tstarting xind={xind}, line {iline} reaches inner bndry")
                        region = 11
                    else:
                        xind = interp1d(xarray, xiarray)(xEnd)
                        if xind > ixsep1 + 0.5:
                            region = 1
                            print(f"\tending xind={xind}, line {iline} enters the SOL")

                    if direction == 1 and yStart == nypf2 and region == 0:
                        shiftangle = interp1d(xiarray, sa)(xind)
                        zEnd += shiftangle
                        yEnd = nypf1 + 1
                    elif direction == -1 and yStart == nypf1 + 1 and region == 0:
                        shiftangle = interp1d(xiarray, sa)(xind)
                        zEnd -= shiftangle
                        yEnd = nypf2

                    if zEnd < zmin or zEnd > zmax:
                        zEnd %= zmax

                    zind = interp1d(zarray, ziarray)(zEnd)

                    it += 1
                    traj[0, it] = iturn
                    traj[1, it] = xind
                    traj[2, it] = yEnd
                    traj[3, it] = zind
                    traj[4, it] = region
                    traj[5, it] = hthe[round(xind), yEnd]

                    xStart = xEnd
                    yStart = yEnd
                    zStart = zEnd

            iturn += 1

        itmax = it
        traj = traj[:, :itmax]

        if save_traj:
            filename = f"./mat_traj/x{iline}y{yyy}z{zzz}_v3lc-01-250{'p' if direction == 1 else 'm'}.mat"
            save_trajectory(filename, traj)

if __name__ == "__main__":
    # Define or load required variables here
    main()
