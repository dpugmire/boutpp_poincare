

%pkg load netcdf; %% for octave

function dump_fieldline_data(fname, nx, ny, nz, zperiod, ixsep1, ixsep2, jyseps1_1, jyseps1_2, jyseps2_1, jyseps2_2, rxy, zxy, rxy_cfr, zxy_cfr, shiftAngle, zShift, zShift_cfr, psixy, dxdy, dzdy, dxdy_p1, dzdy_p1, dxdy_m1, dzdy_m1)

  fprintf('nxyz= %d %d %d\n', nx,ny,nz);
  drp_print_dims('psixy', psixy);
  drp_print_dims('dxdy', dxdy);
  drp_print_dims('dzdy', dzdy);
  drp_print_dims('dxdy_p1', dxdy_p1);
  drp_print_dims('dzdy_p1', dzdy_p1);
  drp_print_dims('dxdy_m1', dxdy_m1);
  drp_print_dims('dzdy_m1', dzdy_m1);
  drp_print_dims('shiftAngle', shiftAngle);
  drp_print_dims('zShift', zShift);
  nx_cfr = size(rxy_cfr, 1);
  ny_cfr = size(rxy_cfr, 2);

  ncid = netcdf.create(fname, 'CLOBBER');

  % Define dimensions
  nx_dim = netcdf.defDim(ncid, 'nx', nx);
  ny_dim = netcdf.defDim(ncid, 'ny', ny);
  nz_dim = netcdf.defDim(ncid, 'nz', nz);
  zperiod_dim = netcdf.defDim(ncid, 'zperiod', zperiod);
  ixsep1_dim = netcdf.defDim(ncid, 'ixsep1', ixsep1);
  ixsep2_dim = netcdf.defDim(ncid, 'ixsep2', ixsep2);
  jyseps1_1_dim = netcdf.defDim(ncid, 'jyseps1_1', jyseps1_1);
  jyseps1_2_dim = netcdf.defDim(ncid, 'jyseps1_2', jyseps1_2);
  jyseps2_1_dim = netcdf.defDim(ncid, 'jyseps2_1', jyseps2_1);
  jyseps2_2_dim = netcdf.defDim(ncid, 'jyseps2_2', jyseps2_2);
  nx_cfr_dim = netcdf.defDim(ncid, 'nx_cfr', nx_cfr);
  ny_cfr_dim = netcdf.defDim(ncid, 'ny_cfr', ny_cfr);

  % Define the variable
  v_psixy = netcdf.defVar(ncid, 'psixy', 'double', [nx_dim, ny_dim]);
  v_dxdy = netcdf.defVar(ncid, 'dxdy', 'double', [nx_dim ny_dim nz_dim]);
  v_dzdy = netcdf.defVar(ncid, 'dzdy', 'double', [nx_dim ny_dim nz_dim]);
  v_dxdy_p1 = netcdf.defVar(ncid, 'dxdy_p1', 'double', [nx_dim, nz_dim]);
  v_dzdy_p1 = netcdf.defVar(ncid, 'dzdy_p1', 'double', [nx_dim, nz_dim]);
  v_dxdy_m1 = netcdf.defVar(ncid, 'dxdy_m1', 'double', [nx_dim, nz_dim]);
  v_dzdy_m1 = netcdf.defVar(ncid, 'dzdy_m1', 'double', [nx_dim, nz_dim]);
  v_shiftangle = netcdf.defVar(ncid, 'shiftAngle', 'double', [nx_dim]);
  v_zShift = netcdf.defVar(ncid, 'zShift', 'double', [nx_dim, ny_dim]);
  v_rxy = netcdf.defVar(ncid, 'rxy', 'double', [nx_dim, ny_dim]);
  v_zxy = netcdf.defVar(ncid, 'zxy', 'double', [nx_dim, ny_dim]);
  v_rxy_cfr = netcdf.defVar(ncid, 'rxy_cfr', 'double', [nx_cfr_dim, ny_cfr_dim]);
  v_zxy_cfr = netcdf.defVar(ncid, 'zxy_cfr', 'double', [nx_cfr_dim, ny_cfr_dim]);
  v_zShift_cfr = netcdf.defVar(ncid, 'zShift_cfr', 'double', [nx_cfr_dim, ny_cfr_dim]);
  % End define mode

  netcdf.endDef(ncid);
  netcdf.putVar(ncid, v_psixy, psixy);
  netcdf.putVar(ncid, v_dxdy, dxdy);
  netcdf.putVar(ncid, v_dzdy, dzdy);
  netcdf.putVar(ncid, v_dxdy_p1, dxdy_p1);
  netcdf.putVar(ncid, v_dzdy_p1, dzdy_p1);
  netcdf.putVar(ncid, v_dxdy_m1, dxdy_m1);
  netcdf.putVar(ncid, v_dzdy_m1, dzdy_m1);
  netcdf.putVar(ncid, v_shiftangle, shiftAngle);
  netcdf.putVar(ncid, v_zShift, zShift);
  netcdf.putVar(ncid, v_rxy, rxy);
  netcdf.putVar(ncid, v_zxy, zxy);
  netcdf.putVar(ncid, v_rxy_cfr, rxy_cfr);
  netcdf.putVar(ncid, v_zxy_cfr, zxy_cfr);
  netcdf.putVar(ncid, v_zShift_cfr, zShift_cfr);

  netcdf.close(ncid);

end

function drp_print_dims(vname, var)
  szvar = size(var);
  str = mat2str(szvar);
  fprintf('%s : %s\n', vname, str);
end

