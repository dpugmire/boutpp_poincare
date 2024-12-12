import numpy as np
from scipy.interpolate import CubicSpline, splrep, splev
from scipy.fft import fft, ifft

def get_apar_sc(psi, bxy, psixy, zs, sa, sinty, dy0, dz, zperiod, interp_opt, deriv_opt, true_apar):
    nx, ny, nz = psi.shape

    apar = np.zeros((nx, ny, nz))
    dapardx = np.zeros((nx, ny, nz))
    dapardy = np.zeros((nx, ny, nz))
    dapardz = np.zeros((nx, ny, nz))
    
    ny_cfr = ny + 1
    iy_cfr = np.arange(ny)
    zarray = np.arange(nz) * dz
    zarrayp = np.arange(nz + 1) * dz
    zmax = dz * nz

    if not true_apar:
        # Step 1: Interpolation of psi back to CELL_CENTER
        psis = np.copy(psi)
        psisp1 = np.zeros((nx, 1, nz))

        # Shift first point across the branch-cut to create ny+1 cell
        for i in range(nx):
            zarray_shift = np.mod(zarray + sa[i], zmax)
            psis_tmp = psi[i, 0, :]
            psis_tmp = np.append(psis_tmp, psis_tmp[0])  # Append first element
            psisp1[i, 0, :] = CubicSpline(zarrayp, psis_tmp)(zarray_shift)
        
        if interp_opt == 0:
            # Linear interpolation
            print('Linear interpolation of psi back to cell center.')
            for i in range(nx):
                psi[i, iy_cfr[:-1], :] = 0.5 * (psis[i, iy_cfr[:-1], :] + psis[i, iy_cfr[1:], :])
                psi[i, ny - 1, :] = 0.5 * (psisp1[i, 0, :] + psis[i, ny - 1, :])
        elif interp_opt == 1:
            # Cubic spline interpolation
            print('Cubic spline interpolation of psi back to cell center.')
            for i in range(nx):
                for k in range(nz):
                    psis_cfr = psis[i, iy_cfr, k]
                    psis_cfr = np.append(psis_cfr, psisp1[i, 0, k])
                    psi[i, :ny, k] = CubicSpline(np.arange(1, ny_cfr + 1), psis_cfr)(np.arange(1, ny_cfr) + 0.5)
        else:
            print('Unknown interpolation method for psi!')

        # Step 2: apar = psi * bxy
        for k in range(nz):
            apar[:, :, k] = psi[:, :, k] * bxy
    else:
        # If the input is a true apar
        apar = np.copy(psi)

    # Step 3: Calculate derivatives of apar
    apars = np.zeros((nx, ny, nz))
    dapardpsi = np.zeros((nx, ny, nz))
    kz = np.concatenate([np.arange(nz // 2 + 1), np.arange(-nz // 2 + 1, 0)]) * zperiod
    ci = 1j

    for i in range(nx):
        for j in range(ny):
            apars[i, j, :] = np.real(ifft(fft(apar[i, j, :]) * np.exp(-ci * zs[i, j] * kz)))

    if deriv_opt == 0:
        # Numerical differentiation (less accurate)
        print('Central finite difference for Apar.')
        dpsi = psixy[101, 39] - psixy[100, 39]
        dapardpsi[1:-1, :, :] = 0.5 * (apars[2:, :, :] - apars[:-2, :, :]) / dpsi

        dapardz[:, :, 0] = 0.5 * (apar[:, :, 1] - apar[:, :, -1]) / dz
        dapardz[:, :, 1:-1] = 0.5 * (apar[:, :, 2:] - apar[:, :, :-2]) / dz
        dapardz[:, :, -1] = 0.5 * (apar[:, :, 0] - apar[:, :, -2]) / dz

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

        # d/dy
        for i in range(nx):
            dapardy[i, iy_cfr[1:-1], :] = 0.5 * (apar[i, iy_cfr[2:], :] - apar[i, iy_cfr[:-2], :]) / dy0
            zarray_shift = np.mod(zarray - sa[i], zmax)
            apar_tmp = apar[i, ny - 1, :]
            apar_tmp = np.append(apar_tmp, apar_tmp[0])
            aparm1 = CubicSpline(zarrayp, apar_tmp)(zarray_shift)
            dapardy[i, 0, :] = 0.5 * (apar[i, 1, :] - aparm1) / dy0
            zarray_shift = np.mod(zarray + sa[i], zmax)
            apar_tmp = apar[i, 0, :]
            apar_tmp = np.append(apar_tmp, apar_tmp[0])
            aparp1 = CubicSpline(zarrayp, apar_tmp)(zarray_shift)
            dapardy[i, ny - 1, :] = 0.5 * (aparp1 - apar[i, ny - 2, :]) / dy0
    elif deriv_opt == 1:
        # Analytical differentiation based on cubic spline
        print('Cubic spline fit then differentiation of Apar.')
        for k in range(nz):
            for j in range(ny):
                apar_tmp = apars[:, j, k]
                apar_tmp = np.append(apar_tmp, apar_tmp[-1])  # Patch last point
                xarray = psixy[:, j]
                xarray = np.append(xarray, 2 * xarray[-1] - xarray[-2])
                pp = CubicSpline(xarray, apar_tmp)
                dapardpsi[:, j, k] = pp(xarray, 1)  # First derivative

        for k in range(nz):
            for i in range(nx):
                apar_tmp = apar[i, :, k]
                apar_tmp = np.append(apar_tmp, apar_tmp[0])
                pp = CubicSpline(zarrayp, apar_tmp)
                dapardy[i, :, k] = pp(np.arange(ny + 1) + 0.5, 1) / dy0

        for j in range(ny):
            for i in range(nx):
                apar_tmp = apar[i, j, :]
                apar_tmp = np.append(apar_tmp, apar_tmp[0])
                pp = CubicSpline(zarrayp, apar_tmp)
                dapardz[i, j, :] = pp(zarrayp, 1)

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

    return apar, dapardx, dapardy, dapardz
