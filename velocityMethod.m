function [T,globalPos] = velocityMethod(shiftCell)
%Properties of FBG Sensors
Lam_a = 1535*10^-9; %Meters
Lam_b = 1531*10^-9; %Meters
Lam_c = 1539*10^-9; %Meters
Lam = [1535, 1531, 1539]*10^-9;
P_e = 0.22;

%Specfications determined by calibration
r_a = 356*10^-6;    %Meters
r_b = 312*10^-6;    %Meters
r_c = 266*10^-6;    %Meters
L_ab = 126.2*pi/180; %Radians
L_bc = 113.6*pi/180; %Radians


a=[r_a, 0, 0];
b=[r_b*cos(L_ab), r_b*sin(L_ab), 0];
c=[r_c*cos(L_ab+L_bc), r_c*sin(L_ab+L_bc), 0];

pre=eye(4);
for i=1:length(shiftCell)
    for j=1:3
        e(j) = shiftCell{i}(j)/(Lam(j)*(1-P_e));
    end
    A = [-a(1) a(2) 1;
        -b(1) b(2) 1;
        -c(1) c(2) 1;];
    m = e/A;
    w_y=m(1); w_x=m(2); eta=m(3);
    
    xi_hat = [0 0 w_y 0;
        0 0 -w_x 0;
        -w_y w_x 0 1;
        0 0 0 0;];
    T{i}=pre*expm(xi_hat);
    pre=T{i};
    globalPos{i}=T{i}(:,4);
end



end