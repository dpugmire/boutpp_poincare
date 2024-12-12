import numpy as np
from scipy.interpolate import CubicSpline
from scipy.fftpack import fft, ifft
import utils

# Define the function to calculate apar and its derivatives
def get_apar_sn(psi, bxy, psixy, zs, sa, sinty, dy0, dz,
                ixsep, nypf1, nypf2, zperiod, interp_opt, deriv_opt, true_apar):
    nx, ny, nz = psi.shape

    apar = np.zeros((nx, ny, nz))
    dapardx = np.zeros((nx, ny, nz))
    dapardy = np.zeros((nx, ny, nz))
    dapardz = np.zeros((nx, ny, nz))

    #ny_cfr = nypf2 - nypf1 + 1
    ny_cfr = nypf2 - nypf1
    print('ny_cfr= ', ny_cfr)
    print('nypf1,2= ', nypf1, nypf2)
    #iy_cfr = np.arange(nypf1, nypf2 + 1)
    iy_cfr = np.arange(nypf1, nypf2)
    #iy_pfr = np.concatenate([np.arange(nypf1), np.arange(nypf2 + 1, ny)])
    iy_pfr = np.concatenate([np.arange(nypf1), np.arange(nypf2, ny)])
    utils.write_array_to_file(iy_cfr, 'iy_cfr')
    utils.write_array_to_file(iy_pfr, 'iy_pfr')

    zarray = np.arange(nz) * dz
    zarrayp = np.arange(nz + 1) * dz
    zmax = dz * nz


    ## current code doesn't do this!!!
    if not true_apar:
        # Step 1: Interpolation of psi back to CELL_CENTER
        psis = psi.copy()
        psisp1 = np.zeros((nx, 1, nz))

        # Shift first point across the branch cut
        for i in range(ixsep):
            zarray_shift = (zarray + sa[i]) % zmax
            psis_tmp = psi[i, nypf1:, :].copy()
            psis_tmp = np.append(psis_tmp, psis_tmp[:, 0:1], axis=1)
            psisp1[i, 0, :] = CubicSpline(zarrayp, psis_tmp[:, 0])(zarray_shift)

        if interp_opt == 0:
            # Linear interpolation
            print("Linear interpolation of psi back to cell center.")
            for i in range(nx):
                if i <= ixsep:
                    # Closed flux surface region
                    psi[i, iy_cfr[:-1], :] = 0.5 * (psis[i, iy_cfr[:-1], :] + psis[i, iy_cfr[1:], :])
                    psi[i, nypf2, :] = 0.5 * (psisp1[i, 0, :] + psis[i, nypf2, :])
                    # Private flux region
                    psi[i, iy_pfr[:-1], :] = 0.5 * (psis[i, iy_pfr[:-1], :] + psis[i, iy_pfr[1:], :])
                    psi[i, iy_pfr[-1], :] = 1.5 * psis[i, iy_pfr[-1], :] - 0.5 * psis[i, iy_pfr[-2], :]
                else:
                    psi[i, :-1, :] = 0.5 * (psis[i, :-1, :] + psis[i, 1:, :])
                    psi[i, -1, :] = 1.5 * psis[i, -1, :] - 0.5 * psis[i, -2, :]

        elif interp_opt == 1:
            # Cubic spline interpolation
            print("Cubic spline interpolation of psi back to cell center.")
            for i in range(nx):
                if i <= ixsep:
                    for k in range(nz):
                        psis_cfr = psis[i, iy_cfr, k]
                        psis_cfr = np.append(psis_cfr, psisp1[i, 0, k])
                        psi[i, nypf1:nypf2, k] = CubicSpline(np.arange(1, ny_cfr + 1), psis_cfr)(np.arange(1, ny_cfr))

                        psis_pfr = psis[i, iy_pfr, k]
                        psis_pfr = np.append(psis_pfr, 2 * psis_pfr[-1] - psis_pfr[-2])
                        psi[i, iy_pfr, k] = CubicSpline(np.arange(1, 2 * nypf1 + 2), psis_pfr)(np.arange(1, 2 * nypf1 + 1))
                else:
                    for k in range(nz):
                        psis_tmp = psis[i, :, k]
                        psis_tmp = np.append(psis_tmp, psis_tmp[-1])
                        psi[i, :, k] = CubicSpline(np.arange(1, ny + 2), psis_tmp)(np.arange(1, ny + 1))

        else:
            raise ValueError("Unknown interpolation method for psi!")

        # Step 2: Calculate apar
        for k in range(nz):
            apar[:, :, k] = psi[:, :, k] * bxy

    else:
        # True apar
        apar = psi.copy()

    # Step 3: Calculate apar derivatives
    apars = np.zeros((nx, ny, nz))
    dapardpsi = np.zeros((nx, ny, nz))
    #kz=[0:nz/2 -nz/2+1:1:-1]*zperiod;
    # kz=kz';   ' is complex conjugate
    kz = np.concatenate([np.arange(0, nz//2 + 1), np.arange(-nz//2 + 1, 0)]) * zperiod
    kz = kz.reshape(-1, 1)  # Equivalent to transposing to a column vector
    utils.write_array_to_file(kz, 'kz')
    utils.write_array_to_file(zs, 'zs')

    ci = 1j

    for i in range(nx):
        for j in range(ny):
            ## apars(i,j,:)=real(ifft(fft(squeeze(apar(i,j,:))).*exp(-ci*zs(i,j)*kz)));
            # Apply FFT, multiply, and inverse FFT
            apar_ij = apar[i, j, :]  # Extract the (i, j, :) slice
            transformed = np.fft.fft(apar_ij) * np.exp(-ci * zs[i, j] * kz.flatten())
            apars[i, j, :] = np.real(np.fft.ifft(transformed))


    utils.write_array_to_file(apars, 'apars')
    shit()


    if deriv_opt == 0:
        # Numerical differentiation (finite difference)
        print("Central finite difference for Apar.")
        dpsi = np.abs(psixy[100, 38] - psixy[99, 38])
        dapardpsi[1:-1, :, :] = 0.5 * (apars[2:, :, :] - apars[:-2, :, :]) / dpsi

        dapardz[:, :, 0] = 0.5 * (apar[:, :, 1] - apar[:, :, -1]) / dz
        dapardz[:, :, 1:-1] = 0.5 * (apar[:, :, 2:] - apar[:, :, :-2]) / dz
        dapardz[:, :, -1] = 0.5 * (apar[:, :, 0] - apar[:, :, -2]) / dz

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

        for i in range(nx):
            if i < ixsep:
                dapardy[i, iy_cfr[1:-1], :] = 0.5 * (apar[i, iy_cfr[2:], :] - apar[i, iy_cfr[:-2], :]) / dy0

                zarray_shift = (zarray - sa[i]) % zmax
                apar_tmp = apar[i, nypf2, :].copy()
                apar_tmp = np.append(apar_tmp, apar_tmp[0])
                aparm1 = CubicSpline(zarrayp, apar_tmp)(zarray_shift)
                dapardy[i, nypf1 + 1, :] = 0.5 * (apar[i, nypf1 + 2, :] - aparm1) / dy0

                zarray_shift = (zarray + sa[i]) % zmax
                apar_tmp = apar[i, nypf1, :].copy()
                apar_tmp = np.append(apar_tmp, apar_tmp[0])
                aparp1 = CubicSpline(zarrayp, apar_tmp)(zarray_shift)
                dapardy[i, nypf2, :] = 0.5 * (aparp1 - apar[i, nypf2 - 1, :]) / dy0

    elif deriv_opt == 1:
        # Analytical differentiation
        print("Cubic spline fit then differentiation of Apar.")
        for k in range(nz):
            for j in range(ny):
                apar_tmp = apars[:, j, k]
                apar_tmp = np.append(apar_tmp, apar_tmp[-1])
                xarray = psixy[:, j]
                xarray = np.append(xarray, 2 * xarray[-1] - xarray[-2])
                pp = CubicSpline(xarray, apar_tmp)
                dapardpsi[:, j, k] = pp(xarray, 1)

    return apar, dapardx, dapardy, dapardz
