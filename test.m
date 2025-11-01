radvel = 3;
f0 = 4e9;
lambda = physconst('LightSpeed')/f0;
doppler_shift = speed2dop(radvel,lambda);


speed_mps = 3; % Convert km/h to m/s
dopplerShift_Hz = (speed_mps / 299792458) * f0;