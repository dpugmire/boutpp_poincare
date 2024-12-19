% field-line tracing with RK4 integration
%
% RK4 method: for dy/dx=f
%       k1=hf(xn,yn)
%       k2=hf(xn+h/2,yn+k1/2)
%       k3=hf(xn+h/2,yn+k2/2)
%       k4=hf(xn+h,yn+k3)
%       yn+1=yn+(k1+2k2+2k3+k4)/6
%
function [xEnd,zEnd] = RK4_FLT1(xStart,yStart,zStart,dxdy,dzdy,xarray,zarray,region,dxdy_pm1,dzdy_pm1,dir,nypf1,nypf2)
    hh=1/2.; h6=1/6.;

    dumpFiles = 0;
    % need half step and full step info
    if (dir == 1)
        dxdyp=squeeze(dxdy(:,yStart,:)); dzdyp=squeeze(dzdy(:,yStart,:));
        %write_array_to_file(dxdyp, 'dxdyp_1');
        %write_array_to_file(dzdyp, 'dzdyn_1');


        if (region ==0 && yStart==nypf2)
            dumpFiles = 0;
            dxdyn=dxdy_pm1;
            dxdyh=0.5*(dxdyp+dxdyn);
            dzdyn=dzdy_pm1;
            dzdyh=0.5*(dzdyp+dzdyn);
        else
            dumpFiles = 0;
            dxdyn=squeeze(dxdy(:,yStart+1,:));
            dxdyh=0.5*squeeze(dxdy(:,yStart,:)+dxdy(:,yStart+1,:));
            dzdyn=squeeze(dzdy(:,yStart+1,:));
            dzdyh=0.5*squeeze(dzdy(:,yStart,:)+dzdy(:,yStart+1,:));
        end
    elseif (dir == -1)
        dxdyp=squeeze(dxdy(:,yStart,:)); dzdyp=squeeze(dzdy(:,yStart,:));
        if (region ==0 && yStart==nypf1+1)
            dxdyn=dxdy_pm1;
            dxdyh=0.5*(dxdyp+dxdyn);
            dzdyn=dzdy_pm1;
            dzdyh=0.5*(dzdyp+dzdyn);
        else
            dxdyn=squeeze(dxdy(:,yStart-1,:));
            dxdyh=0.5*squeeze(dxdy(:,yStart,:)+dxdy(:,yStart-1,:));
            dzdyn=squeeze(dzdy(:,yStart-1,:));
            dzdyh=0.5*squeeze(dzdy(:,yStart,:)+dzdy(:,yStart-1,:));
        end
    else
        fprintf('\tCheck direction parameter setting!! \n');
    end
    if dumpFiles == 1
        write_array_to_file(dxdyn, 'dxdyn');
        write_array_to_file(dxdyh, 'dxdyh');
        write_array_to_file(dzdyn, 'dzdyn');
        write_array_to_file(dzdyh, 'dzdyh');
    end

    % pathing last z point
    dxdyp(:,end+1)=dxdyp(:,end);
    dxdyn(:,end+1)=dxdyn(:,end);
    dxdyh(:,end+1)=dxdyh(:,end);
    dzdyp(:,end+1)=dzdyp(:,end);
    dzdyn(:,end+1)=dzdyn(:,end);
    dzdyh(:,end+1)=dzdyh(:,end);

    if dumpFiles == 1
        write_array_to_file(dxdyp, 'dxdyp_');
        write_array_to_file(dxdyn, 'dxdyn_');
        write_array_to_file(dxdyh, 'dxdyh_');
        write_array_to_file(dzdyp, 'dzdyp_');
        write_array_to_file(dzdyn, 'dzdyn_');
        write_array_to_file(dzdyh, 'dzdyh_');
    end

    % first step
    %fprintf('xzStart= %12.10f %12.10f\n', xStart, zStart);
    % dxdyp': the ' performs the complex conjugate transpose.
    dxdy1=interp2(xarray,zarray,dxdyp',xStart,zStart,'spline');
    dzdy1=interp2(xarray,zarray,dzdyp',xStart,zStart,'spline');
    x1=xStart+dir*hh*dxdy1;
    z1=zStart+dir*hh*dzdy1;
    %fprintf('dxdy1/dzdy1= %12.10e %12.10e xz= %12.10e %12.10e\n', dxdy1, dzdy1, x1,z1);

    % second step
    dxdy2=interp2(xarray,zarray,dxdyh',x1,mod(z1,2*pi),'spline');
    dzdy2=interp2(xarray,zarray,dzdyh',x1,mod(z1,2*pi),'spline');
    x2=xStart+dir*hh*dxdy2;
    z2=zStart+dir*hh*dzdy2;
    %fprintf('dxdy2/dzdy2= %12.10e %12.10e xz= %12.10e %12.10e\n', dxdy2, dzdy2, x2,z2);

    % third step
    dxdy3=interp2(xarray,zarray,dxdyh',x2,mod(z2,2*pi),'spline');
    dzdy3=interp2(xarray,zarray,dzdyh',x2,mod(z2,2*pi),'spline');
    x3=xStart+dir*dxdy3;
    z3=zStart+dir*dzdy3;
    %fprintf('dxdy3/dzdy3= %12.10e %12.10e xz= %12.10e %12.10e\n', dxdy3, dzdy3, x3,z3);
    
    % forth step
    dxdy4=interp2(xarray,zarray,dxdyn',x3,mod(z3,2*pi),'spline');
    dzdy4=interp2(xarray,zarray,dzdyn',x3,mod(z3,2*pi),'spline');

    % accumulate increments with proper weights 
    xEnd = xStart+dir*h6*(dxdy1+2*dxdy2+2*dxdy3+dxdy4);
    zEnd = zStart+dir*h6*(dzdy1+2*dzdy2+2*dzdy3+dzdy4);
    %fprintf('dxdy4/dzdy4= %12.10e %12.10e xz= %12.10e %12.10e\n', dxdy4, dzdy4,xEnd,zEnd);

end