%% matlab script to trace field-lines for poloidal Poincare plot and other analysis
%  last updated by B. Zhu 12/2023


%addpath('~/Documents/MATLAB/utilities')
%addpath('~/Documents/MATLAB/colormaps')

clear all;
close all;

set(groot,'DefaultLineLinewidth',2);
set(groot,'DefaultAxesFontSize',20);
set(groot,'DefaultTextFontSize',20);
set(groot,'DefaultLineMarkerSize',12);

% Check if running in Octave
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    disp('Running in Octave');
    pkg load netcdf; %% for octave
else
    disp('Running in MATLAB');
end

doHermite = false

trajFID = fopen('/Users/dpn/traj.m.txt', 'w');
stepsFID = fopen('/Users/dpn/steps.m.txt', 'w');
rk4FID = fopen('/Users/dpn/rk4.m.txt', 'w');
if (doHermite)
  trajSplineFID = fopen('/Users/dpn/trajhermite.m.txt', 'w');
  puncSplineFid = fopen('/Users/dpn/punchermite.m.txt', 'w');
else
  trajSplineFID = fopen('/Users/dpn/trajspline.m.txt', 'w');
  puncSplineFid = fopen('/Users/dpn/puncspline.m.txt', 'w');
end

rawPuncFid = fopen('/Users/dpn/rawpunc.m.txt', 'w');
puncFid = fopen('/Users/dpn/punc.m.txt', 'w');
punc_ip_Fid = fopen('/Users/dpn/punc_ip.m.txt', 'w');
TRAJ_FID = fopen('/Users/dpn/pt_traj.m.txt', 'w');

fprintf(trajFID, 'ID, STEP, X, Y, Z, REGION, YI, ZSVALUE, ZVALUE\n');
fprintf(trajSplineFID, 'ID, STEP, X, Y, Z, REGION\n');
fprintf(rawPuncFid, 'ID, STEP, X, Y, Z\n');
fprintf(puncFid, 'ID, STEP, X, Y, Z\n');
fprintf(puncSplineFid, 'ID, STEP, X, Y, Z, REGION\n');
fprintf(punc_ip_Fid, 'ID, STEP, X, Y, Z\n');
fprintf(stepsFID, 'ID, STEP, X, Y, Z\n');
fprintf(rk4FID, 'ID, STEP, X, Y, Z\n');
fprintf(TRAJ_FID, "IT, xind, yEnd, zind, REG, zEnd\n");

%%% STEP 0: user setup

% BOUT++ grid file
%gridfile =  'kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc';
gridfile =  '../data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc';

% Mesh resolution info
nx = 260; ny = 128; nz = 256; zperiod = 1;
% Field-line tracing direction: 1 (y index increasing); -1 (y index decreasing)
direction = 1;
% Individual field-lines to be traced radially
nlines = 256; % number of field-lines in radial direction
deltaix = 1; ixoffset = 1; % by default, line tracing starts at
                           % index space (deltaix*ilines+ixoffset, iy, iz)
% (Roughly) total poloidal turns
nturns = 250;
nturns = 5+1;
nturns = 15+1;
nsteps = nturns*ny;
np = 1250; % maximum puncture points, rougly nturns*q
% Output option
save_traj = 1;  % save trajectory of each field line
save_pp = 1;    % save puncture point of each field line
saveFields = 0; % save out computed fields for the python version.
saveHighResTraj = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% no user setup in this section %%%%%%%%

% fill the torus if needed, i.e., zperiod>1
nzG=nz*zperiod; % now nz extends full torus
zmin=0.0; zmax=2*pi; dz= (zmax-zmin)/nzG;
ziarray=(1:nzG+1); zarray=(ziarray-1)*dz;

xiarray=(1:nx);
% BZ: currently use index-based selection, could be modified as flux-based
% xset=(0:nlines-1)*(xmax-xmin)/(nlines-1)+xmin;

yiarray = (1:ny);

%%% STEP 1: load grid info and simulation output

    fprintf('Loading grid information ...\n');

    % Import variable from grid file
    fid = netcdf.open(gridfile, 'nc_nowrite');
    vid = netcdf.inqVarID(fid, 'zShift');
    zShift = netcdf.getVar(fid, vid);
    %zShift = double(zShift);
    if ( length(size(zShift)) ~= 2 || sum(size(zShift)) == 2 )
        fprintf('\tPlease, check the grid file. zShift is not 2-dimension variable\n');
        return
    end
    zShift = permute(zShift, [2 1]);

    vid = netcdf.inqVarID(fid, 'Rxy'); rxy = netcdf.getVar(fid, vid); rxy = double(rxy); rxy = permute(rxy, [2 1]);
    vid = netcdf.inqVarID(fid, 'Zxy'); zxy = netcdf.getVar(fid, vid); zxy = double(zxy); zxy = permute(zxy, [2 1]);
    vid = netcdf.inqVarID(fid, 'psixy'); psixy = netcdf.getVar(fid, vid); psixy = double(psixy); psixy = permute(psixy, [2 1]);
    vid = netcdf.inqVarID(fid, 'rmag'); rmag = netcdf.getVar(fid, vid);

    vid = netcdf.inqVarID(fid, 'ixseps1');  ixsep1 = netcdf.getVar(fid, vid);
    vid = netcdf.inqVarID(fid, 'ixseps2');  ixsep2 = netcdf.getVar(fid, vid);
    if (ixsep2 < nx)
        divertor = 2; % double null
        fprintf('\tDouble null configration\n');
        vid = netcdf.inqVarID(fid, 'jyseps1_1'); nypf11 = netcdf.getVar(fid, vid) + 1;
        vid = netcdf.inqVarID(fid, 'jyseps2_1'); nypf21 = netcdf.getVar(fid, vid) + 1;
        vid = netcdf.inqVarID(fid, 'jyseps1_2'); nypf12 = netcdf.getVar(fid, vid) + 1;
        vid = netcdf.inqVarID(fid, 'jyseps2_2'); nypf22 = netcdf.getVar(fid, vid) + 1;
    elseif (ixsep1 < nx)
        divertor = 1;
        fprintf('\tSingle null configration\n');
        ixsep = ixsep1; % index of the LCFS
        vid = netcdf.inqVarID(fid, 'jyseps1_1'); nypf1 = netcdf.getVar(fid, vid) + 1;
        vid = netcdf.inqVarID(fid, 'jyseps2_2'); nypf2 = netcdf.getVar(fid, vid) + 1;
    else
        divertor = 0;
        ixsep = nx; nypf1 = 0; nypf2 = ny;
        fprintf('\tCircular configration\n');
    end

    [~,jyomp]=max(rxy(end,:));
    xarray = psixy(:,jyomp); % in psi
    % note in BOUT++ convention, psixy is always increasing UNLESS ...
    xMin=min(xarray); xMax=max(xarray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% [check 1] plot grid in 3D configuration space
%     x3d = rxy.*cos(zShift);
%     y3d = rxy.*sin(zShift);
%     z3d = zxy;
%     figure(1)
%     if (divertor == 0)
%         mesh(x3d,y3d,z3d);
%     elseif (divertor == 1)
%         mesh(x3d(1:ixsep1,nypf1:nypf2),y3d(1:ixsep1,nypf1:nypf2),z3d(1:ixsep1,nypf1:nypf2))
%         hold on
%         mesh(x3d(ixsep1:end,:),y3d(ixsep1:end,:),z3d(ixsep1:end,:))
%         % missing privite flux region
%     else
%         fprintf('to be implemented!');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for local shear information
    vid = netcdf.inqVarID(fid, 'Bxy');   bxy  = netcdf.getVar(fid, vid); bxy  = permute(bxy,  [2 1]);
    vid = netcdf.inqVarID(fid, 'Btxy');  btxy = netcdf.getVar(fid, vid); btxy = permute(btxy, [2 1]);
    vid = netcdf.inqVarID(fid, 'Bpxy');  bpxy = netcdf.getVar(fid, vid); bpxy = permute(bpxy, [2 1]);
    vid = netcdf.inqVarID(fid, 'hthe');  hthe = netcdf.getVar(fid, vid); hthe = permute(hthe, [2 1]);
    vid = netcdf.inqVarID(fid, 'sinty'); sinty= netcdf.getVar(fid, vid); sinty= permute(sinty,[2 1]);

    vid = netcdf.inqVarID(fid, 'bxcvx'); bxcvx= netcdf.getVar(fid, vid); bxcvx= permute(bxcvx,[2 1]);
    vid = netcdf.inqVarID(fid, 'bxcvy'); bxcvy= netcdf.getVar(fid, vid); bxcvy= permute(bxcvy,[2 1]);
    vid = netcdf.inqVarID(fid, 'bxcvz'); bxcvz= netcdf.getVar(fid, vid); bxcvz= permute(bxcvz,[2 1]);

    vid = netcdf.inqVarID(fid, 'Jpar0'); jpar0= netcdf.getVar(fid, vid); jpar0= permute(jpar0,[2 1]);
    vid = netcdf.inqVarID(fid, 'dy');  dy = netcdf.getVar(fid, vid); dy0 = dy(1,1);

    vid = netcdf.inqVarID(fid, 'ShiftAngle'); sa =netcdf.getVar(fid, vid);

    clear vid;
    netcdf.close(fid);

    dz = 2*pi/zperiod/nz;
    nu = btxy.*hthe./bpxy./rxy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% [check 2] double check zshift and nu
%     xpos = 150;
%     nu_xpos = nu(xpos,[nypf1+1:nypf2 nypf1+1]);
%     zs(1)=0;
%     for iy = 2:nypf2-nypf1+1
%         zs(iy)= 0.5*(nu_xpos(iy-1)+nu_xpos(iy))*dy0+zs(iy-1);
%     end
%     figure(2)
%     %yyaxis left
%     subplot(2,1,1)
%     plot(zs,'s')
%     hold on
%     plot(zShift(xpos,nypf1+1:nypf2)-zShift(xpos,nypf1+1),'+')
%     xlabel('y index (w/o offset) at x index 150')
%     legend('zShift(calculated)','zShift(grid)')
%
%     % calcuate zshift for the last point before the branch cut
%     nu_temp = nu(1:ixsep1,[nypf1+1:nypf2 nypf1+1]);
%     zs_temp = zeros(ixsep1,nypf2-nypf1);
%     for iy = 2:nypf2-nypf1+1
%         zs_temp(:,iy)= 0.5*(nu_temp(:,iy-1)+nu_temp(:,iy))*dy0+zs_temp(:,iy-1);
%     end
%     zs_last=zs_temp(:,end);
%     %yyaxis right
%     subplot(2,1,2)
%     plot(zs_last,'s')
%     hold on
%     plot(sa,'+')
%     xlabel('x index')
%     legend('shiftangle(calculated)','shiftangle(grid)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calculate sol, x-point locations and poloidal angle theta
    if (divertor == 1) % single null
        corex = rxy(1,nypf1+1:ny-nypf1); corey = zxy(1,nypf1+1:ny-nypf1);
        solx = rxy(1,1:nypf1); soly = zxy(1,1:nypf1);
        solx(nypf1+1:2*nypf1) = rxy(1,ny-nypf1+1:ny); soly(nypf1+1:2*nypf1) = zxy(1,ny-nypf1+1:ny);
        solx(2*nypf1:2*nypf1+nx-1) = rxy(:,end); soly(2*nypf1:2*nypf1+nx-1) = zxy(:,end);
        solx(2*nypf1+nx-1:2*nypf1+nx+ny-2) = rxy(end,ny:-1:1); soly(2*nypf1+nx-1:2*nypf1+nx+ny-2) = zxy(end,ny:-1:1);
        solx(2*nypf1+nx+ny-2:2*nypf1+2*nx+ny-3) = rxy(nx:-1:1,1); soly(2*nypf1+nx+ny-2:2*nypf1+2*nx+ny-3) = zxy(nx:-1:1,1);
        tmp = [1:nypf1 ny-nypf1:-1:nypf1+1 ny-nypf1+1:ny];
        sepx = 0.5*(rxy(ixsep,tmp)+rxy(ixsep+1,tmp)); sepy = 0.5*(zxy(ixsep,tmp)+zxy(ixsep+1,tmp));

        center_x=0.5*(max(rxy(1,nypf1+1:ny-nypf1))+min(rxy(1,nypf1+1:ny-nypf1)));
        center_y=0.5*(max(zxy(1,nypf1+1:ny-nypf1))+min(zxy(1,nypf1+1:ny-nypf1)));
        % use x-point as the reference point, i.e., theta=0
        xpoint_x=0.25*(sepx(nypf1)+sepx(nypf1+1)+sepx(ny-nypf1)+sepx(ny-nypf1+1));
        xpoint_y=0.25*(sepy(nypf1)+sepy(nypf1+1)+sepy(ny-nypf1)+sepy(ny-nypf1+1));
        u=[center_x-xpoint_x center_y-xpoint_y 0];
%         % or, use omp as the reference point, i.e., theta=0
%         u=[center_x-3. 0 0];

        for iy=1:ny
            v=[center_x-rxy(1,iy) center_y-zxy(1,iy) 0];
            theta(iy)=atan2(norm(cross(u,v)),dot(u,v));
        end
        theta=theta/pi;
        [c,itheta]=max(theta);
        theta(itheta:ny)=2-theta(itheta:ny);
        [c,itheta]=max(theta);
        if (itheta ~=ny)
            theta(itheta:ny)=4-theta(itheta:ny);
        end
        % for more accuracy, shift to (1,nypf+1) as the reference point
        theta=theta-theta(nypf1+1);

    elseif (divertor == 0) % shift-circular configuration
	    center_x=0.5*(max(rxy(1,:))+min(rxy(1,:)));
        center_y=0.5*(max(zxy(1,:))+min(zxy(1,:)));
	    % use (1,1) as the reference point
	    u=[center_x-rxy(1,1) center_y-zxy(1,1) 0];

        for iy=1:ny
            v=[center_x-rxy(1,iy) center_y-zxy(1,iy) 0];
            theta(iy)=atan2(norm(cross(u,v)),dot(u,v));
        end

        theta=theta/pi;
        [c,itheta]=max(theta);
        theta(itheta:ny)=2-theta(itheta:ny);
        [c,itheta]=max(theta);
        if (itheta ~=ny)
            theta(itheta:ny)=4-theta(itheta:ny);
        end

    else
        fprintf('\t\tConfiguration to be implemented!');
    end

    % construct/patch closed flux surface region for a better Poincare plot
    xiarray_cfr=double(1:ixsep);
    yiarray_cfr=double(nypf1+1:nypf2+1); % note one additional grid point is patched
    theta_cfr=theta(yiarray_cfr); theta_cfr(end)=2.; % theta is pi based

    rxy_cfr=rxy(1:ixsep,nypf1+1:nypf2);   rxy_cfr(:,end+1)=rxy_cfr(:,1);
    zxy_cfr=zxy(1:ixsep,nypf1+1:nypf2);   zxy_cfr(:,end+1)=zxy_cfr(:,1);

    zs_cfr =zShift(1:ixsep,nypf1+1:nypf2);
    % need an additional zshift info on the last point
    zs_cfr(:,end+1)=0.5*(nu(xiarray_cfr,nypf1+1)+nu(xiarray_cfr,nypf2))*dy0+zs_cfr(:,end);
    % may also need to do the same for sinty and bxcvz

    % calculate geometric coefficients -- refer to note
    A1 = rxy.*bpxy.*btxy./hthe;
    A2 = bxy.^2;
    A3 = sinty.*A1;
    JJ = 4.*pi*1.e-7*bpxy./hthe./(bxy.^2).*jpar0;

    fprintf('Loading perturbed field information ...\n');
    % this script use BOUT++ output psi, note apar=psi*B0 -- refer to idl script
    apar    = zeros(nx,ny,nzG);
    dapardx = zeros(nx,ny,nzG);
    dapardy = zeros(nx,ny,nzG);
    dapardz = zeros(nx,ny,nzG);

    write_array_to_file(xiarray_cfr, 'xiarray_cfr');
    write_array_to_file(yiarray_cfr, 'yiarray_cfr');
    write_array_to_file(zs_cfr, 'zs_cfr');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % analytical apar model for test purpose -- only for closed flux region
%     sFactor1=1e-5;  sFactor2=0e-4;
%     alpx=185.637*100;   alpy=2.53303;
%     nz1=5; nz2=60;
%     xMid = xarray(ixsep); yMid = 39;
%
% %    apar=exp(-alpx*(x-xMid).^2)*EXP(-alpy*(y-yMid).^2)*( sFactor1*sin(nz1*z) + sFactor2*sin(nz2*z) );
% %    dAdx=EXP(-alpx*(x-xMid).^2)*EXP(-alpy*(y-yMid).^2)*(-2*(x-xMid)*alpx)*( sFactor1*SIN(nz1*z) + sFactor2*SIN(nz2*z) );
% %    dAdy=EXP(-alpx*(x-xMid).^2)*EXP(-alpy*(y-yMid).^2)*(-2*(y-yMid)*alpy)*( sFactor1*SIN(nz1*z) + sFactor2*SIN(nz2*z) );
% %    dAdz=EXP(-alpx*(x-xMid).^2)*EXP(-alpy*(y-yMid).^2)*( sFactor1*nz1*cos(nz1*z) + sFactor2*nz2*COS(nz2*z) );
%
%     for ix=1:nx
%         for jy=1:ny
%             for kz=1:nz
%                 apar(ix,jy,kz)  =exp(-alpx*(xarray(ix)-xMid).^2)*exp(-alpy*(jy-yMid).^2) ...
%                     *( sFactor1*sin(nz1*zarray(kz)) + sFactor2*sin(nz2*zarray(kz)) );
%                 dapardx(ix,jy,kz)=exp(-alpx*(xarray(ix)-xMid).^2)*exp(-alpy*(jy-yMid).^2) ...
%                     *(-2*(xarray(ix)-xMid)*alpx) ...
%                     *( sFactor1*sin(nz1*zarray(kz)) + sFactor2*sin(nz2*zarray(kz)) );
%                 dapardy(ix,jy,kz)=exp(-alpx*(xarray(ix)-xMid).^2)*exp(-alpy*(jy-yMid).^2) ...
%                     *(-2*(jy-yMid)*alpy) ...
%                     *( sFactor1*sin(nz1*zarray(kz)) + sFactor2*sin(nz2*zarray(kz)) );
%                 dapardz(ix,jy,kz)=exp(-alpx*(xarray(ix)-xMid).^2)*exp(-alpy*(jy-yMid).^2) ...
%                     *( sFactor1*nz1*cos(nz1*zarray(kz)) + sFactor2*nz2*cos(nz2*zarray(kz)) );
%             end
%         end
%     end
%
%     fprintf('\tmodeled apar = %f ...\n',sFactor1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% end "no user setup" section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %parpool('local',8);

    % Read RMP info on BOUT++ grid either from mat file or from grid file
    rmp=load('../data/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat');

    if (divertor == 0)
        [apar0,dapardx0,dapardy0,dapardz0] = ...
            get_apar_sc(rmp.apar,bxy,psixy,zShift,sa,sinty,dy0,dz,...
            zperiod,0,1,1);
    elseif (divertor == 1)
        [apar0,dapardx0,dapardy0,dapardz0] = ...
            get_apar_sn(rmp.apar,bxy,psixy,zShift,sa,sinty,dy0,dz,...
            ixsep,nypf1,nypf2,zperiod,0,0,1);
    else
        fprintf('\tConfiguration to be implemented!');
    end

    for zp=1:zperiod
        apar(:,:,(zp-1)*nzG/zperiod+1:zp*nzG/zperiod)=apar0;
        dapardx(:,:,(zp-1)*nzG/zperiod+1:zp*nzG/zperiod)=dapardx0;
        dapardy(:,:,(zp-1)*nzG/zperiod+1:zp*nzG/zperiod)=dapardy0;
        dapardz(:,:,(zp-1)*nzG/zperiod+1:zp*nzG/zperiod)=dapardz0;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % [check 3] Plot apar and dapardx, dapardy, dapardz
%     figure(3)
%     set(gcf,'Position',[100 100 1200 400])
%     subplot(1,4,1)
%     BOUT_poloidal_snapshot(gridfile,apar,nz,zperiod,0,1)
%     subplot(1,4,2)
%     BOUT_poloidal_snapshot(gridfile,dapardx,nz,zperiod,0,1)
%     subplot(1,4,3)
%     BOUT_poloidal_snapshot(gridfile,dapardy,nz,zperiod,0,1)
%     subplot(1,4,4)
%     BOUT_poloidal_snapshot(gridfile,dapardz,nz,zperiod,0,1)
%     for i=1:4
%         subplot(1,4,i)
%         shading flat; axis equal; axis tight;
%         xlim([1.2,2.3]); ylim([-0.8,1.3]); colorbar;
%     end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 2: field-line tracing

    % calculate perturbed field
    b0dgy = bpxy./hthe;
    for k=1:nzG
        bdgx(:,:,k) = (1./bxy).*(-A1.*dapardy(:,:,k)-A2.*dapardz(:,:,k))+apar(:,:,k).*bxcvx;
        bdgy(:,:,k) = (1./bxy).*( A1.*dapardx(:,:,k)-A3.*dapardz(:,:,k))+apar(:,:,k).*(bxcvy+JJ);
        bdgz(:,:,k) = (1./bxy).*( A2.*dapardx(:,:,k)+A3.*dapardz(:,:,k))+apar(:,:,k).*bxcvz;

        dxdy(:,:,k)=bdgx(:,:,k)./(b0dgy+bdgy(:,:,k));
        dzdy(:,:,k)=bdgz(:,:,k)./(b0dgy+bdgy(:,:,k));
    end

    % may also need additional grid point information near the branch-cut
    % for the closed flux surface region; we assume zeroth order continuity
    % here, may not hold when a large perturbed field is present.
    dxdyt=squeeze(dxdy(:,nypf1+1,:)); dxdyt(:,end+1)=dxdyt(:,1);
    dzdyt=squeeze(dzdy(:,nypf1+1,:)); dzdyt(:,end+1)=dzdyt(:,1);
    % p1, first point to ny plus 1, twist-shift from first point
    dxdy_p1=zeros(nx,nzG); dzdy_p1=zeros(nx,nzG);
    for ix=1:ixsep
        zarray_shift=mod(zarray(1:nzG)+sa(ix),zmax);
        dxdy_p1(ix,:)=spline(zarray,dxdyt(ix,:),zarray_shift);
        dzdy_p1(ix,:)=spline(zarray,dzdyt(ix,:),zarray_shift);
    end
    dxdyt=squeeze(dxdy(:,nypf2,:)); dxdyt(:,end+1)=dxdyt(:,1);
    dzdyt=squeeze(dzdy(:,nypf2,:)); dzdyt(:,end+1)=dzdyt(:,1);
    %m1, last point to minus 1, reverse twist-shift from last point
    dxdy_m1=zeros(nx,nzG); dzdy_m1=zeros(nx,nzG);
    for ix=1:ixsep
        zarray_rshift=mod(zarray(1:nzG)-sa(ix),zmax);
        dxdy_m1(ix,:)=spline(zarray,dxdyt(ix,:),zarray_rshift);
        dzdy_m1(ix,:)=spline(zarray,dzdyt(ix,:),zarray_rshift);
    end

    %% save out tracing data to netcdf.
    if saveFields == 1
        fprintf('Saving out fields....\n');
        %dxdy = rand(nx,ny,nz);
        %n=nx*ny*nz;
        %dxdy = reshape(0:(n-1), nx,ny,nz);
        val0 = dxdy(124,80,102);
        val1 = dxdy(193, 48, 201);
        write_array_to_file(dxdy, 'dxdy_0');
        write_array_to_file(dzdy, 'dzdy_0');
        write_array_to_file(rxy, 'rxy');
        write_array_to_file(zxy, 'zxy');
        write_array_to_file(rxy_cfr, 'rxy_cfr');
        write_array_to_file(zxy_cfr, 'zxy_cfr');

        %tmp = reshape(0:(n-1), nx,ny,nz);
        %printEval(dxdy, 124, 80, 102);
        %printEval(dxdy, 14, 100, 28);

        %fprintf('ixseps: %d  jyseps: %d %d %d %d  nyfp1,2: %d %d\n', ixsep1, ixsep2, jyseps1_1, jyseps1_2, jyseps2_1, jyseps2_2, nypf1, nypf2);
        dump_fieldline_data('stuff.nc', nx, ny, nz, rxy, zxy, rxy_cfr, zxy_cfr, sa, zShift, zs_cfr, psixy, dxdy, dzdy, dxdy_p1, dzdy_p1, dxdy_m1, dzdy_m1);
        return;
    end

    fprintf('Starting field-line tracing ...\n');
    fprintf('\n');

    cm = jet(nlines);

    COUNTER = 0;
    LINES = 1:nlines;
    LINES = [100, 150, 175, 195, 200, 210];
    LINES = [150]; %% this stays in region 0 the whole time.
    %LINES = [195]; %% this from region 0 to 1 to 2.
    %LINES = [130, 140, 150, 160, 170, 180];
    LINES = [150];
    LINES = [150];
    LINES = 1:50:nlines;
    LINES = [150];
    LINES = [151];
    %LINES = [1,51, 101, 151, 201, 251];
    LINES = 1:50:nlines;

    YVALS = 1:ny-1;
    %YVALS = [60];

    %parfor iline = 1:nlines
    zindFID = fopen('/Users/dpn/zind.m.txt', 'w');

    tic; %% Start timer
    rk4_cnt = 0;
    for iline = LINES

        % pick starting points
        xind = iline;
        xStart = psixy(xind,jyomp); % note here jyomp doesn't matter
        yyy = jyomp;
        yStart = jyomp;
        zzz = 1;
	      zStart = zarray(zzz);

        % declare trajectory and puncture points arrays for this field-line
        traj=zeros(7,nsteps);fl_x3d=zeros(nsteps,1);fl_y3d=zeros(nsteps,1);fl_z3d=zeros(nsteps,1);
        px=zeros(np,1);py=zeros(np,1);pz=zeros(np,1);ptheta=zeros(np,1);ppsi=zeros(np,1);

        it = 1; iturn = 1;

        if (xind < double(ixsep)+0.5)
            region = 0; % =0, closed flux surface;
                        % =1 sol; =2 pfr;
                        % = 11/12 inner/outer budry; 13/14 inner/outer divertor
            if (yStart<nypf1+1 || yStart>nypf2)
                region = 2;
            end
        else
            region = 1;
        end

        yind = yStart;
        zind = interp1(zarray, ziarray, zStart);

        fprintf('\n\n\n\t**********         line %i started at idx: (%f,%f,%f), (%f %f %f)\n',iline,xind,yind,zind, xStart, yStart, zStart);

        % rule out the starting points on the divertor targets and go towards
        % the divertor plates
        if (divertor == 1)
            if (yStart == double(ny) && direction == 1)
                fprintf('\tline %i starts on the divertor.\n',iline);
                region=14;
                traj(6,it)=0.; lc=0.;
            elseif (yStart == double(1) && direction == -1)
                fprintf('\tline %i starts on the divertor.\n',iline);
                region=13;
                traj(6,it)=0.; lc=0.;
            end
        end

        % record field-line location info in trajectory array
        traj(1,it) = 1;
        traj(2,it) = xind;
        traj(3,it) = yStart;
        traj(4,it) = zind;
        traj(5,it) = region;
        % zStart info is stored for better interpolation of puncture point
        traj(7,it) = zStart;
%        traj(:,it)=[1;xind;yStart;zind;region];
        fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", traj(1,it)-1, traj(2,it)-1, traj(3,it)-1, traj(4,it)-1, traj(5,it), traj(7,it));

        %fprintf(stepsFID, 'COUNTER= %d region= %d xyzind= %d %d %d\n', COUNTER, region, xind-1, yind-1, zind-1);

        fprintf('region= %d\n', region);
        while (region < 10 && iturn < nturns)

            %fprintf(stepsFID, 'COUNTER= %d region= %d iturn= %d\n', COUNTER, region, iturn);

            if (mod(iturn,50) == 1)
                fprintf('\t\t line%i, turn %i/%i ...\n',iline,iturn,nturns);
            end

            % start field-line tracing
            for iy = YVALS
                %fprintf(trajFID, '%d, %d, %d, %d, %.8f, %d, %.8f\n', iline-1, iy-1, it-1, iturn-1, xStart, yStart-1, zStart);

                if (it == 20)
                  printf('****** start debugging.\n');
                  printf('\n\n');
                endif

                %fprintf('%d: xi,yi= %d %d region: %d xyzStart= %.8f %.8f %.8f  xyzInd= %.8f %.8f %.8f\n', COUNTER, iline, iy, region, xStart,yStart,zStart, xind,yind,zind);
                %fprintf(stepsFID, 'COUNTER= %d region= %d iy= %d\n   xyzStart= %.8f %.8f %.8f\n   xyzInd= %.8f %.8f %.8f\n', COUNTER, region, iy-1, xStart, yStart-1, zStart, xind-1, yind-1, zind-1);
                if (region == 0 && yStart > nypf1 && yStart < nypf2+1) % in CFR

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                % quick test for equilibrium field
  %                %xEnd = xStart;
  %                %zEnd = zStart;
  %                % if perturbed field, i.e, finite dxdy, dzdy
  %                %[xEnd,zEnd]=RK4_FLT0(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                  %disp(size(dxdy));
                  %fprintf('dxdy: %e\n', dxdy(50, 100, 50));
                  %fprintf('dxdy: %e\n', dxdy(150, 48, 250));

                  if (direction == 1)
                      % This is at step 343 where things go wonky in the c++ code.
                      _x = -0.13563311719504134;
                      _y = 69;
                      _z = 6.2769540966961062;

                      %[xTMP,zTMP] = RK4_FLT1(_x,_y,_z,dxdy,dzdy,xarray,zarray,region,dxdy_p1,dzdy_p1,1,nypf1,nypf2, rk4FID, iline, it, true);
                      %shit();

                      if rk4_cnt == 14
                        fprintf('here we are at 14.\n');
                      endif

                      [xEnd,zEnd]=RK4_FLT1(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray,region,dxdy_p1,dzdy_p1,1,nypf1,nypf2, rk4FID, iline, it, false);
                      yEnd = yStart+1;
                      rk4_cnt = rk4_cnt+1;
                  elseif (direction == -1)
                      [xEnd,zEnd]=RK4_FLT1(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray,region,dxdy_m1,dzdy_m1,-1,nypf1,nypf2, rk4FID, iline, it, false);
                      yEnd = yStart-1;
                  end
                  fprintf(stepsFID, '%d, %d, %f, %d, %f\n', iline-1, it-1, xEnd, yEnd-1, zEnd);

                  % zEnd info is needed for better interpolation of puncture point
                  if it+1 == 8
                    printf('need to look at this difference.....\n');
                    printf('stop now\n');

                  endif
                  traj(7,it+1)=zEnd;

                  % check the where is the end of the fieldline
                  if (xEnd > xMax) % field-line hits the outer boundary
                      fprintf('\tstarting xind=%f, line %i reaches outer bndry\n',xind,iline);
                      region = 12;
                  elseif (xEnd < xMin) % field-line hits the inner boundary
                      fprintf('\tstarting xind=%f, line %i reaches inner bndry\n',xind,iline);
                      region = 11;
                  else % end of field-line remains in the simulation domain
                      xind = interp1(xarray, xiarray, xEnd);
                      if (xind > double(ixsep1)+0.5)
                          region = 1;
                          fprintf('\tending xind=%f, line %i enters the SOL\n',xind,iline);
                      end
                  end

                  % twist-shift at the branch cut
                  if (direction == 1 && yStart == nypf2 && region ==0)
                      shiftangle = interp1(xiarray,sa,xind);
                      zEnd = zEnd+shiftangle;
                      yEnd = nypf1+1;
                  end

                  if (direction == -1 && yStart == nypf1+1 && region ==0)
                      shiftangle = interp1(xiarray,sa,xind);
                      zEnd = zEnd-shiftangle;
                      yEnd = nypf2;
                  end

                  % re-label toroidal location (zEtnd) if necessary
                  %fprintf('********** it= %d pt0= %f %f %f\n', it-1, xEnd, yEnd-1, zEnd);
                  if (zEnd < zmin || zEnd > zmax)
                      zEnd = mod(zEnd,zmax);
                  end
                  zind = interp1(zarray, ziarray, zEnd);
                  fprintf(zindFID, '%d %12.10f --> %12.10f\n', it-1, zEnd, zind-1);
                  %fprintf('********** it= %d pt1= %f %d %f zind %f\n', it-1, xEnd, yEnd-1, zEnd, zind);

                  it = it+1;
                  traj(1,it) = iturn;
                  traj(2,it) = xind;
                  traj(3,it) = yEnd;
                  traj(4,it) = zind;
                  traj(5,it) = region;
                  % approximate field-line segment length
                  traj(6,it) = hthe(round(xind),yEnd);
  %                traj(:,it)=[iturn;xind;yEnd;zind;region];
                  fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", it-1, traj(2,it)-1, traj(3,it)-1, traj(4,it)-1, traj(5,it), traj(7,it));

                  % now the end-point becomes new start-point
                  xStart = xEnd;
                  yStart = yEnd;
                  zStart = zEnd;

                elseif (region == 1 || region == 2) % for points at SOL/PFR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    % if equilibrium field
%                    %xEnd = xStart;
%                    %zEnd = zStart;
%                    % if perturbed field, i.e, finite dxdy, dzdy
%                    %[xEnd,zEnd]=RK4_FLT0(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if (direction == 1)
                        [xEnd,zEnd]=RK4_FLT1(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray,region,dxdy_p1,dzdy_p1,1,nypf1,nypf2);
                        yEnd = yStart+1;
                    elseif (direction == -1)
                        [xEnd,zEnd]=RK4_FLT1(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray,region,dxdy_m1,dzdy_m1,-1,nypf1,nypf2);
                        yEnd = yStart-1;
                    end

                    traj(7,it+1)=zEnd; % for better interpolation of puncture point

                    % correct yEnd for two PFR cases
                    if (direction == 1 && yStart == nypf1 && region == 2)
                        yEnd = nypf2+1;
                    elseif (direction == -1 && yStart == nypf2+1 && region == 2)
                        yEnd = nypf1;
                    end

                % check the where is the end of the fieldline
                if (xEnd > xMax)
                    fprintf('\tstarting xind=%f, line %i reaches outer bndry\n',xind,iline);
                    region = 12;
                elseif (xEnd < xMin)
                    fprintf('\tstarting xind=%f, line %i reaches inner bndry\n',xind,iline);
                    region = 11;
                else % end of field-line remains in the simulation domain
                    xind = interp1(xarray, xiarray, xEnd);
                    if (xind < double(ixsep1)+0.5 && yEnd > nypf1 && yEnd < nypf2+1)
                        if (region~=0) % if not already in the CFR
                            fprintf('\tending xind=%f, line %i enters the CFR\n',xind,iline);
                        end
                        region = 0;
                    elseif (xind < double(ixsep1)+0.5 && (yEnd > nypf2-1 || yEnd < nypf1))
                        if (region~=2) % if not already in the PFR
                            fprintf('\tending xind=%f, line %i enters the PFR\n',xind,iline);
                        end
                        region = 2;
                    end
                end

                if (direction == 1 && yEnd == ny) % in the sol/pfr, last point is outer divertor no matter what
                    fprintf('\tstarting xind=%f, line %i reaches divertor\n',xind,iline);
                    region = 14;
                elseif (direction == -1 && yEnd == 1) % in the sol/pfr, last point is inner divertor
                    fprintf('\tstarting xind=%f, line %i reaches divertor\n',xind,iline);
                    region = 13;
                end

                % re-label toroidal location (zEnd) if necessary
                if (zEnd < zmin || zEnd > zmax)
                    zEnd = mod(zEnd,zmax);
                end
                zind = interp1(zarray, ziarray, zEnd);

                it = it+1;
                traj(1,it) = iturn;
                traj(2,it) = xind;
                traj(3,it) = yEnd;
                traj(4,it) = zind;
                traj(5,it) = region;
                traj(6,it) = hthe(round(xind),yStart);
%                traj(:,it)=[iturn;xind;yEnd;zind;region];
                fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", traj(1,it)-1, traj(2,it)-1, traj(3,it)-1, traj(4,it)-1, traj(5,it), traj(7,it));

                % now the end-point becomes new start-point
                xStart = xEnd;
                yStart = yEnd;
                zStart = zEnd;

                end % end single (CFR/SOL/PFR) stepping

                COUNTER = COUNTER + 1;

            end % end apparoximate one turn (ny steps)

            iturn = iturn+1;

        end % end nturns

    % record the maximum valid steps
    itmax=it;
    elapsedTime = toc;
    fprintf('*************** time= %.5f s\n', elapsedTime);

    traj=traj(:,1:itmax);
    if (save_traj)
        if (direction == 1)
            parsave(strcat('./mat_traj/x',num2str(iline),'y',num2str(yyy), ...
                'z',num2str(zzz),'_v3lc-01-250','p.mat'),traj);
        elseif (direction == -1)
            parsave(strcat('./mat_traj/x',num2str(iline),'y',num2str(yyy), ...
                'z',num2str(zzz),'_v3lc-01-250','m.mat'),traj);
        end
    end


    %% DRP: Save out high resolution trajectory.
    if (saveHighResTraj)
      _xi = 1:itmax;
      _samples = 1:0.001:itmax;
      nsamples = length(_samples);
      xVals = zeros(1, itmax);
      yVals = zeros(1, itmax);
      zVals = zeros(1, itmax);
      for _it=1:itmax
        xi = traj(2,_it);
        yi = traj(3,_it);
        zi = traj(4,_it);
        %TMP = rxy(:,traj(3,istep))
        rxyvalue = interp1(xiarray,rxy(:,traj(3,_it)),traj(2,_it));
        zsvalue  = interp1(xiarray,zShift(:,traj(3,_it)),traj(2,_it));
        zvalue   = interp1(ziarray,zarray,traj(4,_it));
        x3d_tmp = rxyvalue*cos(zsvalue);
        y3d_tmp = rxyvalue*sin(zsvalue);

        xVals(_it) = x3d_tmp*cos(zvalue)-y3d_tmp*sin(zvalue);
        yVals(_it) = x3d_tmp*sin(zvalue)+y3d_tmp*cos(zvalue);
        zVals(_it) = interp1(xiarray,zxy(:,traj(3,_it)),traj(2,_it));
      endfor

      splineX = spline(_xi, xVals, _samples);
      splineY = spline(_xi, yVals, _samples);
      splineZ = spline(_xi, zVals, _samples);
      splineXEval = spline(_xi, xVals);
      for _it=1:nsamples
        fprintf(trajSplineFID, '%d, %f, %f, %f, %f, %d\n', iline-1, _samples(_it)-1.0, splineX(_it), splineY(_it), splineZ(_it), region);
        % root finding for more accurate punctures.
        if _it > 1 && _it < nsamples
            x0 = ppval(splineXEval, _samples(_it-1));
            x1 = ppval(splineXEval, _samples(_it));
            if x0*x1 < 0 && x0 != 0 && x1 != 0
                %% Zero crossing exists between x0 and x1.
                t0 = _samples(_it-1);
                t1 = _samples(_it);
                cnt = 0;
                tZero = -1.0;
                for _cnt=1:100
                    val0 = ppval(splineXEval, t0);
                    val1 = ppval(splineXEval, t1);
                    dt = t1-t0;
                    tMid = t0 + dt * 0.5;
                    valMid = ppval(splineXEval, tMid);

                    % zero lies between t0 and tMid
                    if val0*valMid < 0 && val0 != 0 && valMid != 0
                        t0 = t0;
                        t1 = tMid;
                        val0 = val0;
                        val1 = valMid;
                    else  % zero lies between tMid and t1
                        t0 = tMid;
                        t1 = t1;
                        val0 = valMid;
                        val1 = val1;
                    endif
                    dx = abs(val0-val1);
                    %fprintf('  *** %d: zero finding: (%12.10f %12.10f) dX= %15.13f\n', _cnt, val0, val1, dx);

                    if dx < 1e-8
                        tZero = t0 + (t1-t0)*0.5;
                        break;
                    endif
                endfor

                px0 = ppval(splineXEval, tZero);
                if (doHermite)
		  splineYEval = pchip(_xi, yVals);
                  splineZEval = pchip(_xi, zVals);
		else 
		  splineYEval = spline(_xi, yVals);
                  splineZEval = spline(_xi, zVals);
                end
                
                py0 = ppval(splineYEval, tZero);
                pz0 = ppval(splineZEval, tZero);
                fprintf(' FOUND Zero: cnt= %d  %d, %f, %f, %f, %f, %d\n', _cnt, iline-1, tZero-1.0, px0, py0, pz0, region);
                fprintf(puncSplineFid, '%d, %f, %f, %f, %f, %d\n', iline-1, tZero-1.0, px0, py0, pz0, region);
            endif
        endif
      endfor

    endif


% % STEP 3: calculate puncture points and generate Poincare plot and others
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % [check 4] plot in y-z index plane for equilibrium field
% %     figure(4)
% %     %plot(traj(3,:),traj(4,:))
% %     hold on
% %     xlim([zmin,zmax]); ylim([0,65])
% %     for it=1:itmax
% %         zvalue = interp1(ziarray,zarray,traj(4,it));
% %         scatter(zvalue,traj(3,it))
% %         pause(0.1)
% %     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % coarse plot of field-line
    % clear fl_x3d fl_y3d fl_z3d ffl_x3d iit px py pz ptheta ppsi

    %DRP
    _fid = fopen("/Users/dpn/problem.m.txt", "w");
    trajvals_fid = fopen("/Users/dpn/trajvals.m.txt", "w");
    fprintf(trajvals_fid, "iter, xind, yend, zind, zend, x3d, y3d, z3d\n");

    for istep=1:itmax
        if istep == 14
          fprintf('begin debugging.\n');
        endif
        _t1 = traj(1,istep); _t2 = traj(2,istep); _t3 = traj(3,istep); _t4 = traj(4,istep); _t7 = traj(7,istep);
        xi = traj(2,istep);
        yi = traj(3,istep);
        zi = traj(4,istep);
        %TMP = rxy(:,traj(3,istep))
        __zshift = zShift(:, traj(3,istep));
        rxyvalue = interp1(xiarray,rxy(:,traj(3,istep)),traj(2,istep));
        zsvalue  = interp1(xiarray,zShift(:,traj(3,istep)),traj(2,istep));
        zvalue   = interp1(ziarray,zarray,traj(4,istep));
        x3d_tmp = rxyvalue*cos(zsvalue);
        y3d_tmp = rxyvalue*sin(zsvalue);
        fl_x3d(istep) = x3d_tmp*cos(zvalue)-y3d_tmp*sin(zvalue);
        fl_y3d(istep) = x3d_tmp*sin(zvalue)+y3d_tmp*cos(zvalue);
        fl_z3d(istep) = interp1(xiarray,zxy(:,traj(3,istep)),traj(2,istep));
        fprintf(trajFID, '%d, %d, %.8f, %.8f, %.8f, %d, %d, %f, %f\n', iline-1, istep-1, fl_x3d(istep), fl_y3d(istep), fl_z3d(istep), region, yi-1, zsvalue, zvalue);
        _t3 = traj(3,istep);
        _t2 = traj(2,istep);
        %fprintf('*** istep= %d r= %12.10e zs= %12.10e z= %12.10e t23= %12.10e %d\n', istep-1, rxyvalue, zsvalue, zvalue, _t2, _t3-1);
        %fprintf(_fid, '*** istep= %d traj= %d %12.10e %d %12.10e\n', istep-1, int32(traj(1,istep)-1), traj(2,istep), int32(traj(3,istep)-1), traj(4,istep));
        %fprintf(_fid, '%d  traj4, zvalue= %12.10e %12.10e\n', istep-1, traj(4,istep), zvalue);

        if istep < 51
          fprintf(trajvals_fid, "%d, %10.8f, %d, %10.8f, %10.8f, %10.8f, %10.8f, %10.8f\n", istep-1, traj(2,istep)-1, floor(traj(3,istep))-1, traj(4,istep)-1, traj(7,istep), fl_x3d(istep), fl_y3d(istep),fl_z3d(istep));
        end
    end
    fclose(_fid);
    fclose(trajvals_fid);

    fprintf('All done dumping the file....\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % [check 5] plot distance between two adjacent points
%     figure(5)
%     r=sqrt((fl_x3d(2:it)-fl_x3d(1:it-1)).^2+ ...
%            (fl_y3d(2:it)-fl_y3d(1:it-1)).^2+ ...
%            (fl_z3d(2:it)-fl_z3d(1:it-1)).^2);
%     plot(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get punctures info needed for Poincare plot

    if (itmax > 1)

      % first find points intercept with x=0 plane (as in Cartesian coordinate)
      itarray = (1:itmax);
      fit = (1:0.0001:itmax);
      ffl_x3d = spline(itarray,fl_x3d(1:itmax),fit);
      zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
      iit = zci(ffl_x3d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % [check 6] quick check to see the accuracy of zero-crossing points
%     figure(6)
%     hold on
%     plot(fl_x3d(itarray),'-d')
%     plot(fit,ffl_x3d)
%     plot(fit(iit),ffl_x3d(iit),'o')
%     pausera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [nc,~]=size(iit); nc=nc-1; ip=0; id=0;

      if (nc > 0)
        ffl_y3d = spline(itarray,fl_y3d(1:itmax),fit);
        ffl_z3d = spline(itarray,fl_z3d(1:itmax),fit);

        for i=1:nc
          iit_i = iit(i);
          valX = ffl_x3d(iit_i);
          valY = ffl_y3d(iit_i);
          valZ = ffl_z3d(iit_i);
          if i == (439)
            fprintf('here is the problem point. follow it...\n');
          endif
          fprintf(rawPuncFid, '%d, %d, %f, %f, %f\n', iline-1, i-1, valX, valY, valZ);
        endfor
      endif


      % only for the field-lines pass x=0 plane, one calculates the
      % corresponding puncture point information
      fprintf('dumping out the punctures.\n');
      fprintf('dumping punctures.\n');
      if (nc > 0)
        for i=1:nc
            %i = 439;
            if i == (439)
              fprintf('here is the problem point. follow it...\n');
            endif

            iit_i= iit(i);
            fit_i = fit(iit_i);
            it=floor(fit(iit(i)));
            a=fit(iit(i))-it; b=1-a;
            %fprintf('iit_i, fit_i, it = %f %f %d\n', iit_i, fit_i, it);

            % linear interpolation along field-line
            %traj(2) = xind, traj(3) = yend: pointsXYZ.x, .y
            xtraj = traj(2, it); xtraj1 = traj(2, it+1);
            ytraj = traj(3, it); ytraj1 = traj(3, it+1);
            ztraj = traj(7, it); ztraj1 = traj(7, it+1);
            traj4 = traj(7, it); traj41 = traj(4, it+1);
            y_m1 = traj(3,it-1);
            fprintf('pt  : %f %f %f %f\n', traj(2,it), traj(3,it), traj(4,it), traj(7,it));
            fprintf('pt_1: %f %f %f %f\n', traj(2,it+1), traj(3,it+1), traj(4,it+1), traj(7,it+1));
            fprintf('a= %f\n', a);

            drp_debug = false;
            if (drp_debug)
              xtraj = 149.54106292358114 + 1.0;
              ytraj = 73 + 1.0;
              ztraj = 0.00029379732422597326;
              traj4 = 0.011970379882940391;

              xtraj1 = 149.61731548903271 + 1.0;
              ytraj1 = 74 + 1.0;
              ztraj1 = 6.2795725343616526;
              traj41 = 255.852802393026;
              y_m1 = 72 + 1.0;
              a = 0.6585999999675245;
              b = 0.3414000000324755;
            endif


            xind_tmp=b*xtraj+a*xtraj1;
            % default, unless at the branch cut
            yind_tmp=b*ytraj+a*ytraj1;
            % with raw zEnd information, it doesn't matter whether the
            % field-line across the branch cut or not, unless ...
            %traj(7) = zEnd
            zvalue = b*ztraj+a*ztraj1;
            fprintf('xind,yind,zvalue= %f %f %f\n', xind_tmp, yind_tmp, zvalue);

            if (abs(ztraj-ztraj1)>1.)
                zvalue=b*mod(ztraj,zmax)+a*mod(ztraj1,zmax);
            end

            % further ajustments for edge cases
            if (ytraj==double(nypf2) && direction==1 && xind_tmp<double(ixsep)+0.5)
                % here y index is extended to one more point for accuracy
                yind_tmp=b*ytraj+a*double(nypf2+1);
            elseif (ytraj==double(nypf1+1) && direction==-1 && xind_tmp<double(ixsep)+0.5)
                % here y index is extended to one more point for accuracy
                yind_tmp=b*double(nypf2+1)+a*ytraj1;
                % one needs to twist-shift in this case
                shiftangle = interp1(xiarray,sa,xind_tmp);
                zvalue = mod(zvalue-shiftangle,zmax);
            elseif (it>1 && (y_m1==double(nypf2) || y_m1==double(nypf1+1)))
                zvalue = b*interp1(ziarray,zarray,traj4) + ...
                     a*interp1(ziarray,zarray,traj41);
            end

            if (xind_tmp < double(ixsep)+0.5)
                rxyvalue = interp2(xiarray_cfr,yiarray_cfr,rxy_cfr',xind_tmp,yind_tmp,'spline');
                zxyvalue = interp2(xiarray_cfr,yiarray_cfr,zxy_cfr',xind_tmp,yind_tmp,'spline');
                zsvalue  = interp2(xiarray_cfr,yiarray_cfr,zs_cfr', xind_tmp,yind_tmp,'spline');
            else
                rxyvalue = interp2(xiarray,yiarray,rxy',xind_tmp,yind_tmp,'spline');
                zxyvalue = interp2(xiarray,yiarray,zxy',xind_tmp,yind_tmp,'spline');
                zsvalue  = interp2(xiarray,yiarray,zShift',xind_tmp,yind_tmp,'spline');
            end

            ipx3d_tmp = rxyvalue*cos(zsvalue);
            ipy3d_tmp = rxyvalue*sin(zsvalue);
            ipx = ipx3d_tmp*cos(zvalue)-ipy3d_tmp*sin(zvalue);
            ipy = ipx3d_tmp*sin(zvalue)+ipy3d_tmp*cos(zvalue);
            ipz = zxyvalue;

            % one side of tokamak vs. the otherside.
            if (ipy>0)
                %if (abs(ipx)>0.05)
                %    id=id+1;
                %else
                ip=ip+1;
                px(ip)=ipx;
                py(ip)=ipy;
                pz(ip)=ipz;

                ptheta(ip)=interp1(yiarray_cfr,theta_cfr,yind_tmp);
                ppsi(ip)=interp1(xiarray,xarray,xind_tmp);
                %% this is the right one.
                %fprintf(puncFid, '%d, %.8f, %.8f, %.8f\n', i, px(ip), py(ip), pz(ip));
                %%
                fprintf(puncFid, '%d, %d, %.8f, %.8f, %.8f\n', iline-1, i-1, ipx, ipy, ipz);
                %end
            else
                fprintf(punc_ip_Fid, '%d, %d, %.8f, %.8f, %.8f\n', iline-1, i-1, ipx, ipy, ipz);
                fprintf('%d: ipy is < 0:  %f \n', i, ipy);
            end

      end

    fprintf('\t\tline %i has %i(+%i) interception points.\n',iline,ip,id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % quick plot of the Poincare plot
%     figure(31)
%     set(gcf,'Position',[200 500 1600 500])
%     subplot(1,3,1)
%     hold on
%     plot(px(1:ip),'color',cm(iline,:))
%     xlabel('zero-crossing points')
%     ylabel('x')
%     subplot(1,3,2)
%     %plot(rxy(ixsep,5:60),zxy(ixsep,5:60),'k')
%     hold on
%     plot(py(1:ip),pz(1:ip),'.','color',cm(iline,:))
%     axis equal;
%     xlabel('R');ylabel('Z')
%     subplot(1,3,3)
%     hold on
%     plot(ptheta(1:ip),ppsi(1:ip),'.','color',cm(iline,:))
%     xlabel('$\theta/(2\pi)$','interpreter','latex');ylabel('$\psi$','interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (ip > 0)
        pxp=px(1:ip);pyp=py(1:ip);pzp=pz(1:ip);ptp=ptheta(1:ip);ppp=ppsi(1:ip);
    else
        pxp=0.; pyp=0.; pzp=0.; ptp=0.; ppp=0.;
    end

    % if field-line not cross x=0 plane at all ...
    else
        pxp=0.; pyp=0.; pzp=0.; ptp=0.; ppp=0.;
    end % end nc>0

    traj0 = zeros(3,itmax);
    traj0(1,1:itmax) = fl_x3d(1:itmax);
    traj0(2,1:itmax) = fl_y3d(1:itmax);
    traj0(3,1:itmax) = fl_z3d(1:itmax);

    % calculate connection length
    lc = sum(traj(6,1:itmax),2);
    fprintf('\t\tline %i: connection length is %f and ends at region=%i.\n',iline,lc,region);

    % save puncture point info for Poincare plot
    if (save_pp)
        if (direction == 1)
            parsave8(strcat('./mat_pp/x',num2str(iline),'y',num2str(yyy), ...
                'z',num2str(zzz),'_v3lc-00-250','p.mat'), ...
                pxp,pyp,pzp,ptp,ppp,traj0,lc,region);
        elseif (direction == -1)
            parsave8(strcat('./mat_pp/x',num2str(iline),'y',num2str(yyy), ...
                'z',num2str(zzz),'_v3lc-00-250','m.mat'), ...
                pxp,pyp,pzp,ptp,ppp,traj0,lc,region);
        end
    end

    end % end itmax>1

    %clear traj fl_x3d fl_y3d fl_z3d fit ffl_x3d iit px py pz ptheta ppsi pxp pyp pzp ptp ppp traj0 lc region

    end % end iline loop

