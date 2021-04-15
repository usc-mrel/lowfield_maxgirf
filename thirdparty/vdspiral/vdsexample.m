%
% Examples using vds.m and vdsmex.m
%

% This is a 16-interleave spiral that achieves
% 1mm resolution, and density decreases linearly
% from supporting 24cm FOV at |k|=0 to 12cm FOV
% at |k|=maximum.
%

% 1e4 G = 1T = 1e3 mT => 10G = mT
% [mT/m/msec] * [10G/mT] * [m/1e2 cm] * [1e3 msec/sec] => 1e2 [G/cm/s]

smax   = 15000;     % maximum slew rate [G/cm/sec] (e.g., 150 [mT/m/ms])
gmax   = 4;         % maximum gradient [G/cm]
T      = 4e-6;      % sampling period [sec]
N      = 16;        % number of interleaves
Fcoeff = [24 -12];  % FOV decreases linearly from 24 to 12cm.
res    = 1;         % resolution [mm]
rmax   = 5 / res;   % [1/cm], corresponds to 1mm resolution.

disp('Calculating Gradient');
[k,g,s,time,r,theta] = vds(smax,gmax,T,N,Fcoeff,rmax);

disp('Plotting Gradient');
g = [real(g(:)) imag(g(:))];
plotgradinfo(g,T);


%
% Here the example is repeated with vdsmex, which
% should be a lot faster!

disp('Press Enter to repeat with vdsmex.');
pause;
figure;

% Alternative stopping condition, if duration is to be limited.
ngmax = 1000000; % max # of points along gradient

disp('Calculating Gradient');
[k,g,s,time] = vdsmex(N,Fcoeff,res,gmax,smax,T,ngmax);

disp('Plotting Gradient');
plotgradinfo(g,T);

