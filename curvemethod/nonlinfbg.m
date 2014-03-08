function F = nonlinfbg(x,a,b,c)
%Properties of FBG Sensors
Lam_a = 1535*10^-9; %Meters
Lam_b = 1531*10^-9; %Meters
Lam_c = 1539*10^-9; %Meters
P_e = 0.22;

%Specfications determined by calibration
r_a = 356*10^-6;    %Meters
r_b = 312*10^-6;    %Meters
r_c = 266*10^-6;    %Meters
L_ab = 126.2*pi/180; %Radians
L_bc = 113.6*pi/180; %Radians

%Nonlinear set of equations to be solved.
F = [(a/(Lam_a*(1-P_e))-x(1)*r_a*sin(x(2))-x(3));
    (b/(Lam_b*(1-P_e))-x(1)*r_b*sin(x(2)+L_ab)-x(3));
    (c/(Lam_c*(1-P_e))-x(1)*r_c*sin(x(2)+L_ab+L_bc)-x(3))];
end