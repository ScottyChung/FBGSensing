function [transMatrix,globalPos] = velocitymethod(shiftCell)
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

preTrans = eye(4);
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
    
    transMatrix{i}=preTrans*expm(xi_hat);
    globalPos{i}=transMatrix{i}(:,4);
    
    plotcurve(xi_hat,preTrans)
    
    preTrans=transMatrix{i};
    
    StrainPlot = 1;
    figure(StrainPlot)
    
    hold on
    set(StrainPlot,'Position',[996   408   672   504])
    plot(i,w_y,'*b')
    plot(i,w_x,'*r')
    
   
    
end
end

function plotcurve(xi_hat,preTrans)
    figure(2)
    hold on
    k=linspace(0,1,10);
    for n=1:length(k)
        mediateMatrix{n} = preTrans*expm(xi_hat*k(n));
        pos{n} = mediateMatrix{n}(1:3,4);
        plot3(pos{n}(1),pos{n}(2),pos{n}(3),'.')
    end
    legend('Velocity Method')
end
