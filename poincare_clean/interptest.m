% Original data
x = 0:5; % x-coordinates
y = 0:5; % y-coordinates
[X, Y] = meshgrid(x, y);
Z = sin(X) .* cos(Y); % Some sample 2D data

% New points for interpolation
xq = 0:0.1:5; % Query x-coordinates
yq = 0:0.1:5; % Query y-coordinates
[Xq, Yq] = meshgrid(xq, yq);

% Spline interpolation
Zq = interp2(X, Y, Z, Xq, Yq, 'spline');

% Plot original and interpolated data
figure;
subplot(1, 2, 1);
surf(X, Y, Z);
title('Original Data');
subplot(1, 2, 2);
surf(Xq, Yq, Zq);
title('Interpolated Data');
