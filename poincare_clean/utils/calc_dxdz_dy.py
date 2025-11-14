#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import CubicSpline
from scipy.io import loadmat

# ----------------------------
# Utilities
# ----------------------------

def read_xy(fid: Dataset, name: str, nx: int, ny: int, dtype=float) -> np.ndarray:
    arr = np.asarray(fid.variables[name][:], dtype=dtype)
    if arr.ndim != 2:
        raise ValueError(f"{name}: expected 2D array, got shape {arr.shape}")
    if arr.shape == (ny, nx):
        arr = arr.T
    elif arr.shape != (nx, ny):
        raise ValueError(f"{name}: unexpected shape {arr.shape}, expected (nx,ny)=({nx},{ny}) or (ny,nx)=({ny},{nx})")
    return arr

def _ensure_dims(ds: Dataset, dims: dict):
    for name, size in dims.items():
        if name in ds.dimensions:
            if len(ds.dimensions[name]) != size:
                raise ValueError(f"Dimension {name} length mismatch: file={len(ds.dimensions[name])}, expected={size}")
        else:
            ds.createDimension(name, size)

def write_scalar(ds: Dataset, name: str, value: int) :
    var = ds.createVariable(name, "i4")
    var.assignValue(value)

def write_1d(ds: Dataset, name: str, arr_1d: np.ndarray, dim_name: str, dtype="f8"):
    if arr_1d.ndim != 1:
        raise ValueError(f"{name}: expected 1D array, got {arr_1d.shape}")
    if len(ds.dimensions[dim_name]) != arr_1d.size:
        raise ValueError(f"{name}: length {arr_1d.size} != dim {dim_name}={len(ds.dimensions[dim_name])}")
    if name in ds.variables:
        del ds.variables[name]
    var = ds.createVariable(name, dtype, (dim_name,), zlib=True)
    var[:] = arr_1d

def write_2d(ds: Dataset, name: str, arr_nxny: np.ndarray, dtype="f8"):
    if arr_nxny.ndim != 2:
        raise ValueError(f"{name}: expected 2D array (nx,ny), got {arr_nxny.shape}")
    ny = ds.dimensions['ny'].size
    nx = ds.dimensions['nx'].size
    if arr_nxny.shape != (nx, ny):
        raise ValueError(f"{name}: in-memory shape {arr_nxny.shape} != (nx,ny)=({nx},{ny})")
    if name in ds.variables:
        del ds.variables[name]
    var = ds.createVariable(name, dtype, ('ny', 'nx'), zlib=True)
    var[:, :] = arr_nxny.T  # Fortran order (ny,nx)

def write_2d_cfr(ds: Dataset, name: str, arr_nxcfr_nycfr: np.ndarray, dtype="f8"):
    if arr_nxcfr_nycfr.ndim != 2:
        raise ValueError(f"{name}: expected 2D array (nx_cfr,ny_cfr), got {arr_nxcfr_nycfr.shape}")
    ny_cfr = ds.dimensions['ny_cfr'].size
    nx_cfr = ds.dimensions['nx_cfr'].size
    if arr_nxcfr_nycfr.shape != (nx_cfr, ny_cfr):
        raise ValueError(f"{name}: in-memory shape {arr_nxcfr_nycfr.shape} != (nx_cfr,ny_cfr)=({nx_cfr},{ny_cfr})")
    if name in ds.variables:
        del ds.variables[name]
    var = ds.createVariable(name, dtype, ('ny_cfr', 'nx_cfr'), zlib=True)
    var[:, :] = arr_nxcfr_nycfr.T  # (ny_cfr,nx_cfr)

def write_3d(ds: Dataset, name: str, arr_nxny_nz: np.ndarray, dtype="f8"):
    if arr_nxny_nz.ndim != 3:
        raise ValueError(f"{name}: expected 3D array (nx,ny,nz), got {arr_nxny_nz.shape}")
    ny = ds.dimensions['ny'].size
    nx = ds.dimensions['nx'].size
    nz = ds.dimensions['nz'].size
    if arr_nxny_nz.shape != (nx, ny, nz):
        raise ValueError(f"{name}: in-memory shape {arr_nxny_nz.shape} != (nx,ny,nz)=({nx},{ny},{nz})")
    if name in ds.variables:
        del ds.variables[name]
    var = ds.createVariable(name, dtype, ('nz', 'ny', 'nx'), zlib=True)
    var[:, :, :] = np.transpose(arr_nxny_nz, (2, 1, 0))  # (nz,ny,nx)

def write_2d_xz(ds: Dataset, name: str, arr_nxnz: np.ndarray, dtype="f8"):
    if arr_nxnz.ndim != 2:
        raise ValueError(f"{name}: expected 2D array (nx,nz), got {arr_nxnz.shape}")
    nx = ds.dimensions['nx'].size
    nz = ds.dimensions['nz'].size
    if arr_nxnz.shape != (nx, nz):
        raise ValueError(f"{name}: in-memory shape {arr_nxnz.shape} != (nx,nz)=({nx},{nz})")
    if name in ds.variables:
        del ds.variables[name]
    var = ds.createVariable(name, dtype, ('nz', 'nx'), zlib=True)
    var[:, :] = arr_nxnz.T  # (nz,nx)

# -----------------------------------------
# Apar & derivatives (SC and SN)
# -----------------------------------------

def _periodic_spline_dz(values_ijk, dz):
    nx, ny, nz = values_ijk.shape
    z = np.arange(nz, dtype=float) * dz
    zp = np.concatenate([z, z[:1] + z[-1] + dz])
    out = np.empty_like(values_ijk, dtype=float)
    for i in range(nx):
        for j in range(ny):
            f = values_ijk[i, j, :]
            fp = np.concatenate([f, f[:1]])
            cs = CubicSpline(zp, fp, bc_type="periodic")
            out[i, j, :] = cs(z, 1)
    return out

def get_apar_sc(psi_or_apar, bxy, psixy, zs, sa, sinty, dy0, dz,
                zperiod, interp_opt, deriv_opt, true_apar):
    nx, ny, nz = psi_or_apar.shape
    if not true_apar:
        raise NotImplementedError("psi->apar path omitted; true_apar required.")
    apar = psi_or_apar.astype(float)

    kz = np.concatenate([np.arange(0, nz//2 + 1), np.arange(-nz//2 + 1, 0)]).astype(float)[None, None, :]
    apars = np.real(np.fft.ifft(np.fft.fft(apar, axis=2) * np.exp(-1j * zs[:, :, None] * kz), axis=2))

    if deriv_opt == 0:
        dapardpsi = np.zeros_like(apar)
        dapardpsi[1:-1,:,:] = (apars[2:,:,:] - apars[:-2,:,:]) / 2.0
        dapardpsi[0,:,:]    = (-3*apars[0,:,:] + 4*apars[1,:,:] - apars[2,:,:]) / 2.0
        dapardpsi[-1,:,:]   = ( 3*apars[-1,:,:]- 4*apars[-2,:,:]+ apars[-3,:,:]) / 2.0
    else:
        dapardpsi = np.zeros_like(apar)
        for k in range(nz):
            for j in range(ny):
                x = psixy[:, j].astype(float)
                f = np.concatenate([apars[:, j, k], [apars[-1, j, k]]])
                xpad = np.concatenate([x, [2*x[-1]-x[-2]]])
                cs = CubicSpline(xpad, f, bc_type="natural")
                dapardpsi[:, j, k] = cs(x, 1)

    dapardz = _periodic_spline_dz(apar, dz) if deriv_opt == 1 else np.gradient(apar, dz, axis=2, edge_order=2)
    dapardx = dapardpsi + sinty[:, :, None] * dapardz
    dapardy = np.zeros_like(apar)
    dapardy[:,1:-1,:] = (apar[:,2:,:] - apar[:,:-2,:])/(2.0*dy0)

    z = np.arange(nz, dtype=float) * dz
    zmax = dz * nz
    zp = np.concatenate([z, z[:1] + z[-1] + dz])
    for i in range(nx):
        z_shift_rev = (z - sa[i]) % zmax
        fend = np.concatenate([apar[i,-1,:], apar[i,-1,:1]])
        cs = CubicSpline(zp, fend, bc_type="periodic")
        apar_m1 = cs(z_shift_rev)
        dapardy[i,0,:] = (apar[i,1,:] - apar_m1)/(2.0*dy0)

        z_shift_fwd = (z + sa[i]) % zmax
        fbeg = np.concatenate([apar[i,0,:], apar[i,0,:1]])
        cs = CubicSpline(zp, fbeg, bc_type="periodic")
        apar_p1 = cs(z_shift_fwd)
        dapardy[i,-1,:] = (apar_p1 - apar[i,-2,:])/(2.0*dy0)

    return apar, dapardx, dapardy, dapardz

def get_apar_sn(psi_or_apar, bxy, psixy, zs, sa, sinty, dy0, dz,
                ixsep, nypf1, nypf2, zperiod, interp_opt, deriv_opt, true_apar):
    nx, ny, nz = psi_or_apar.shape
    if not true_apar:
        raise NotImplementedError("psi->apar path omitted; true_apar required.")
    apar = psi_or_apar.astype(float)

    kz = np.concatenate([np.arange(0, nz//2 + 1), np.arange(-nz//2 + 1, 0)]).astype(float)[None, None, :]
    apars = np.real(np.fft.ifft(np.fft.fft(apar, axis=2) * np.exp(-1j * zs[:, :, None] * kz), axis=2))

    if deriv_opt == 0:
        dapardpsi = np.zeros_like(apar)
        dapardpsi[1:-1,:,:] = (apars[2:,:,:] - apars[:-2,:,:]) / 2.0
        dapardpsi[0,:,:]    = (-3*apars[0,:,:] + 4*apars[1,:,:] - apars[2,:,:]) / 2.0
        dapardpsi[-1,:,:]   = ( 3*apars[-1,:,:]- 4*apars[-2,:,:]+ apars[-3,:,:]) / 2.0
    else:
        dapardpsi = np.zeros_like(apar)
        for k in range(nz):
            for j in range(ny):
                x = psixy[:, j].astype(float)
                f = np.concatenate([apars[:, j, k], [apars[-1, j, k]]])
                xpad = np.concatenate([x, [2*x[-1]-x[-2]]])
                cs = CubicSpline(xpad, f, bc_type="natural")
                dapardpsi[:, j, k] = cs(x, 1)

    dapardz = _periodic_spline_dz(apar, dz) if deriv_opt == 1 else np.gradient(apar, dz, axis=2, edge_order=2)
    dapardx = -dapardpsi + sinty[:, :, None] * dapardz  # SN convention
    dapardy = np.zeros_like(apar)
    dapardy[:,1:-1,:] = (apar[:,2:,:] - apar[:,:-2,:])/(2.0*dy0)

    z = np.arange(nz, dtype=float) * dz
    zmax = dz * nz
    zp = np.concatenate([z, z[:1] + z[-1] + dz])
    for i in range(nx):
        z_shift_rev = (z - sa[i]) % zmax
        fend = np.concatenate([apar[i,-1,:], apar[i,-1,:1]])
        cs = CubicSpline(zp, fend, bc_type="periodic")
        apar_m1 = cs(z_shift_rev)
        dapardy[i,0,:] = (apar[i,1,:] - apar_m1)/(2.0*dy0)

        z_shift_fwd = (z + sa[i]) % zmax
        fbeg = np.concatenate([apar[i,0,:], apar[i,0,:1]])
        cs = CubicSpline(zp, fbeg, bc_type="periodic")
        apar_p1 = cs(z_shift_fwd)
        dapardy[i,-1,:] = (apar_p1 - apar[i,-2,:])/(2.0*dy0)

    return apar, dapardx, dapardy, dapardz

# ----------------------------
# Main
# ----------------------------
#    gridfile="../../data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc",
#    matfile ="../../data/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat",
#    stuff_file="stuff_python.nc",

def run_minimal(gridFile, aparFile, zperiod, timeStep, outputFile,
    save_fields=True,
    #nx=260, ny=128,
    nz_hint=256,
    interp_opt=0, deriv_opt=0, true_apar=True,
):
    gridfile = Path(gridFile)
    aparfile = Path(aparFile)

    print("Loading grid information ...")
    fid = Dataset(gridfile, "r")
    nx = len(fid.dimensions["x"])
    ny = len(fid.dimensions["y"])

    ixseps1 = np.asarray(fid.variables["ixseps1"][:], dtype=int)
    ixseps2 = np.asarray(fid.variables["ixseps2"][:], dtype=int)
    jyseps1_1 = np.asarray(fid.variables["jyseps1_1"][:], dtype=int)
    jyseps1_2 = np.asarray(fid.variables["jyseps1_2"][:], dtype=int)
    jyseps2_1 = np.asarray(fid.variables["jyseps2_1"][:], dtype=int)
    jyseps2_2 = np.asarray(fid.variables["jyseps2_2"][:], dtype=int)
    ny_inner = np.asarray(fid.variables["ny_inner"][:], dtype=int)

    zShift = read_xy(fid, "zShift", nx, ny, float)
    rxy    = read_xy(fid, "Rxy",    nx, ny, float)
    zxy    = read_xy(fid, "Zxy",    nx, ny, float)
    psixy  = read_xy(fid, "psixy",  nx, ny, float)

    bxy    = read_xy(fid, "Bxy",    nx, ny, float)
    btxy   = read_xy(fid, "Btxy",   nx, ny, float)
    bpxy   = read_xy(fid, "Bpxy",   nx, ny, float)
    hthe   = read_xy(fid, "hthe",   nx, ny, float)
    sinty  = read_xy(fid, "sinty",  nx, ny, float)

    bxcvx  = read_xy(fid, "bxcvx",  nx, ny, float)
    bxcvy  = read_xy(fid, "bxcvy",  nx, ny, float)
    bxcvz  = read_xy(fid, "bxcvz",  nx, ny, float)

    jpar0  = read_xy(fid, "Jpar0",  nx, ny, float)

    ixsep1 = int(fid.variables["ixseps1"][()])
    ixsep2 = int(fid.variables["ixseps2"][()])

    if ixsep2 < nx:
        fid.close()
        raise NotImplementedError("Double-null not implemented.")
    elif ixsep1 < nx:
        divertor = 1
        ixsep = ixsep1
        # treat jyseps as 0-based for slicing
        jy11 = int(fid.variables["jyseps1_1"][()])
        jy21 = int(fid.variables["jyseps2_1"][()])
        jy12 = int(fid.variables["jyseps1_2"][()])
        jy22 = int(fid.variables["jyseps2_2"][()])
    else:
        divertor = 0
        ixsep = nx
        jy11 = 0; jy21 = ny-1
        jy12 = 0; jy22 = -1  # unused

    dy0 = float(fid.variables["dy"][0, 0])
    sa  = np.asarray(fid.variables["ShiftAngle"][:], dtype=float).ravel()
    fid.close()


    if ('.mat' in aparFile) :
        print("Loading A_parallel from .mat ...")
        mat = loadmat(aparfile)
        if "apar" not in mat:
            raise KeyError(f"{aparFile} does not contain 'apar'")
        apar = np.asarray(mat["apar"]).astype(float)
    elif ('.npy' in aparFile) :
        print("Loading A_parallel from .npy ...")
        apar = np.load(aparfile)
        apar = apar[...,timeStep]
        print('removing last z')
        apar = apar[:,:,:-1]
        if zperiod > 1 :
            apar = np.tile(apar, (1,1,zperiod))
        print('apar_t: ', apar.shape)

    if apar.ndim != 3:
        raise ValueError(f"apar has unexpected shape {apar.shape}")
    if apar.shape[:2] != (nx, ny):
        raise ValueError(f"apar[:2] shape {apar.shape[:2]} != (nx,ny)=({nx},{ny})")

    nz = apar.shape[2]
    nzG = nz * zperiod
    dz = 2.0*np.pi / nzG
    zarray = np.arange(nzG, dtype=float) * dz

    # CFR patch (NO wrap column; use unique PF-leg indices only)
    if divertor == 1:
        nx_cfr = ixsep
        leg1 = np.arange(jy11, jy21+1, dtype=int)
        leg2 = np.arange(jy12, jy22+1, dtype=int)
        y_cat = np.concatenate([leg1, leg2])
        # de-duplicate while preserving order
        _, first_idx = np.unique(y_cat, return_index=True)
        y_idx = y_cat[np.sort(first_idx)]
        ny_cfr = y_idx.size  # no wrap
        rxy_cfr      = rxy[:ixsep, :][:, y_idx]
        zxy_cfr      = zxy[:ixsep, :][:, y_idx]
        zShift_cfr   = zShift[:ixsep, :][:, y_idx]
    else:
        nx_cfr = ixsep
        ny_cfr = ny  # no wrap for circular either
        rxy_cfr      = rxy[:ixsep, :]
        zxy_cfr      = zxy[:ixsep, :]
        zShift_cfr   = zShift[:ixsep, :]

    # geometry coefficients (nx,ny)
    mu0 = 4.0e-7*np.pi
    A1 = rxy * bpxy * btxy / hthe
    A2 = bxy**2
    A3 = sinty * A1
    JJ = mu0 * bpxy / hthe / (bxy**2) * jpar0

    # Apar derivatives (nx,ny,nz)
    if divertor == 1:
        apar_, dAx, dAy, dAz = get_apar_sn(
            psi_or_apar=apar, bxy=bxy, psixy=psixy, zs=zShift, sa=sa, sinty=sinty, dy0=dy0, dz=dz,
            ixsep=ixsep-1, nypf1=jy11, nypf2=jy21, zperiod=zperiod,
            interp_opt=interp_opt, deriv_opt=deriv_opt, true_apar=True
        )
    else:
        apar_, dAx, dAy, dAz = get_apar_sc(
            psi_or_apar=apar, bxy=bxy, psixy=psixy, zs=zShift, sa=sa, sinty=sinty, dy0=dy0, dz=dz,
            zperiod=zperiod, interp_opt=interp_opt, deriv_opt=deriv_opt, true_apar=True
        )

    # Field-line slopes (nx,ny,nz)
    b0dgy = bpxy / hthe
    invb = (1.0 / bxy)[:, :, None]
    A1e = A1[:, :, None]; A2e = A2[:, :, None]; A3e = A3[:, :, None]; JJe = JJ[:, :, None]
    bdgx = invb * (-A1e * dAy - A2e * dAz) + apar_ * bxcvx[:, :, None]
    bdgy = invb * ( A1e * dAx - A3e * dAz) + apar_ * (bxcvy[:, :, None] + JJe)
    bdgz = invb * ( A2e * dAx + A3e * dAz) + apar_ * bxcvz[:, :, None]
    denom = (b0dgy[:, :, None] + bdgy)
    dxdy = bdgx / denom
    dzdy = bdgz / denom

    # SPECIAL neighbor fields (nx, nz)
    dxdy_p1 = np.zeros((nx, nz), dtype=float)
    dzdy_p1 = np.zeros((nx, nz), dtype=float)
    dxdy_m1 = np.zeros((nx, nz), dtype=float)
    dzdy_m1 = np.zeros((nx, nz), dtype=float)

    if divertor == 1:
        dxdyt = dxdy[:, jy11, :]
        dzdyt = dzdy[:, jy11, :]
        z = zarray
        zmaxG = dz * z.size
        zp = np.concatenate([z, z[:1] + z[-1] + dz])
        for ix in range(ixsep):
            fdx = np.concatenate([dxdyt[ix, :], dxdyt[ix, :1]])
            fdz = np.concatenate([dzdyt[ix, :], dzdyt[ix, :1]])
            csx = CubicSpline(zp, fdx, bc_type="periodic")
            csz = CubicSpline(zp, fdz, bc_type="periodic")
            z_shift_fwd = (z + sa[ix]) % zmaxG
            z_shift_rev = (z - sa[ix]) % zmaxG
            dxdy_p1[ix, :] = csx(z_shift_fwd)
            dzdy_p1[ix, :] = csz(z_shift_fwd)
            dxdy_m1[ix, :] = csx(z_shift_rev)
            dzdy_m1[ix, :] = csz(z_shift_rev)

    # ---------------------------------------
    # Save (Fortran-order dims) in requested order
    # ---------------------------------------
    if save_fields:
        print("Saving out fields (Fortran order dims) …")
        with Dataset(Path(outputFile), "w") as ds:
            _ensure_dims(ds, {'nx': nx, 'ny': ny, 'nz': nz})
            _ensure_dims(ds, {'nx_cfr': rxy_cfr.shape[0], 'ny_cfr': rxy_cfr.shape[1]})

            write_scalar(ds, "ixsep1", ixseps1)
            write_scalar(ds, "ixsep2", ixseps2)
            write_scalar(ds, "jyseps1_1", jyseps1_1)
            write_scalar(ds, "jyseps1_2", jyseps1_2)
            write_scalar(ds, "jyseps2_1", jyseps2_1)
            write_scalar(ds, "jyseps2_2", jyseps2_2)
            write_scalar(ds, "ny_inner", ny_inner)
            write_scalar(ds, "zperiod", zperiod)

            # Order: psixy, dxdy, dzdy, dxdy_p1, dzdy_p1, dxdy_m1, dzdy_m1,
            #        shiftAngle, zShift, rxy, zxy, rxy_cfr, zxy_cfr, zShift_cfr
            write_2d(ds, "psixy", psixy)
            write_3d(ds, "dxdy", dxdy)
            write_3d(ds, "dzdy", dzdy)

            write_2d_xz(ds, "dxdy_p1", dxdy_p1)
            write_2d_xz(ds, "dzdy_p1", dzdy_p1)
            write_2d_xz(ds, "dxdy_m1", dxdy_m1)
            write_2d_xz(ds, "dzdy_m1", dzdy_m1)

            write_1d(ds, "shiftAngle", np.asarray(sa), "nx")
            write_2d(ds, "zShift", zShift)
            write_2d(ds, "rxy", rxy)
            write_2d(ds, "zxy", zxy)

            write_2d_cfr(ds, "rxy_cfr", rxy_cfr)
            write_2d_cfr(ds, "zxy_cfr", zxy_cfr)
            write_2d_cfr(ds, "zShift_cfr", zShift_cfr)

        return dict(dxdy=dxdy, dzdy=dzdy)

    return dict(dxdy=dxdy, dzdy=dzdy)



if __name__ == "__main__":
    if len(sys.argv) != 6 :
        print('Error. Usage: <grid_file> <apar_file> <zperiod> <time_step> <output_file>')
        sys.exit()

    gridFile = sys.argv[1]
    aparFile = sys.argv[2]
    zperiod = int(sys.argv[3])
    timeStep = int(sys.argv[4])
    outputFile = sys.argv[5]

    run_minimal(gridFile, aparFile, zperiod, timeStep, outputFile)


# case from Ben Z.
## python ./calc_dxdz_dy.py ../../data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc ../../data/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat OUT.nc
#This doesn't work when zperiod is not 1.

# single null case.
#python calc_dxdz_dy.py /Users/dpn/proj/bout++/nersc_data/xpoint_singlenull/d3d_194347.03100_0.85psi1.1_132x64y_bout_exp_withSOL.nc /Users/dpn/proj/bout++/nersc_data/xpoint_singlenull/Apar_194347.03100_0.85psi1.1_132x64y_data-low1.npy  5 100 stuff.100.nc