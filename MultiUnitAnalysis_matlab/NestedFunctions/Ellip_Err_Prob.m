
function [meanX, meanY, stdX, stdY, a0, Xe, Ye] =  ... 
    Ellip_Err_Prob(scatter_x_axis, scatter_y_axis, err_percent)

%% Elliptical error probable

% Number of elements in the ellipse
Ne = 5000;

% Standard deviations 
stdX = std(scatter_x_axis);
stdY = std(scatter_y_axis);

% 2x2 Covariance matrix
McovXY = cov(scatter_x_axis, scatter_y_axis);
VarXY = McovXY(2,1);

% Angle of the data
a0 = 0.5*atan2(2*VarXY,(stdY^2-stdX^2));

% Coordinates in ellipse basis
theta = a0-pi/2;
M = [+cos(theta) -sin(theta);...
     +sin(theta) +cos(theta)];
V = M*[scatter_x_axis; scatter_y_axis];
Xr = V(1,:);
Yr = V(2,:);
% Statistic estimators: means and standard deviations
meanX = mean(scatter_x_axis);
meanY = mean(scatter_y_axis);
stdX = std(Xr);
stdY = std(Yr);
% Ellipse coordinates
a = linspace(0,2*pi,Ne);
theta = pi/2-a0;
M = [+cos(theta) -sin(theta);...
     +sin(theta) +cos(theta)];
V =  M*[stdX*cos(a); stdY*sin(a)];
Xe0 = V(1,:);
Ye0 = V(2,:);
% Ellipses coordinates weighted by EEPs
Xe = zeros(1, Ne);
Ye = zeros(1, Ne);

f = sqrt(-2*log(1-err_percent));
Xe(1,:) = f*Xe0+meanX;
Ye(1,:) = f*Ye0+meanY;















