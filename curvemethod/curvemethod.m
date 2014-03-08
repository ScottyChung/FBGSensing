function [T,globalPos]=curvemethod(shiftCell,vargin)
%Calls the nonLinSolve function for each wavelength shift. Requires the
%input of the random shift cell array and the optional plotting of the
%curvature. Returns the curvature, orientaion, and bias of each sensor.
%%
%Calls the nonLinSolve function for the number of rows in shiftCell.
for i=1:length(shiftCell)
    [curvature(i), orientation(i), bias(i)] = nonLinSolve(shiftCell{i}(1),shiftCell{i}(2),shiftCell{i}(3));
end

if strcmp(vargin,'PlotOn')
    curvaturePlot(curvature)
    %uiPlotCurve(shiftMatrix)
end
[T,globalPos] = transformation(curvature, orientation,bias);
end

function [curvature,orientation,bias] = nonLinSolve(sLam_a,sLam_b,sLam_c)
%Uses matlab's fsolve function to solve system of nonlinear equations.
%Given the shift of each sensor at a cluster, returns the curvature,
%orientation, bias.
%%
%Required function for fsolve method.
x0 = [1;0;0]; %Inital guess for k,p,and e_0. Require for fsolve function.

%Creates the anonymous function necessary for the matlab fsolve
%function
f = @(x)nonlinfbg(x,sLam_a,sLam_b,sLam_c); 

%Calls the fsolve function with optional display of information.
% options = optimoptions('fsolve','Display','iter');
options = optimset('Display','off');
[x,fval] = fsolve(f,x0,options);

%Stores the solved function as the returned values.
curvature = x(1);
orientation = x(2);
bias = x(3);
end

function curvaturePlot(curvature)
%Plots the curvature of each sensor location versus their sensor number.
%%
CurvePlot = 1;
figure(CurvePlot)
clf
hold on
set(CurvePlot,'Position',[996   408   672   504])
plot(0:length(curvature),[0, curvature],'*-')
ylabel('Curviture of Fiber')
xlabel('FBG Cluster Number')

end

function [T,globalPos] = transformation(curvature,orientation,bias)
%Given the curvature, orientation, bias of the sensors, plots the shape of
%the sensor along with the rotated coordinate for each sensor.
%%
%Determines the number of sensors from the length of input. 
numSensors = length(curvature);

%Creates the empty transformation cell array
T = cell(numSensors,1);

%Loops through the number of sensors and will create their relative
%transformation matrix. Then will apply the global transformation.
for i=1:numSensors
    
    phi = orientation(i);
    k = curvature(i);
    th = k*1;   %theta = k*ds, where ds is the length between sensors

    %Builds the transformation matrix
    T{i} = [cos(phi) -sin(phi) 0 0;
        sin(phi) cos(phi) 0 0;
        0 0 1 0;
        0 0 0 1] * ...
        [cos(th) 0 sin(th) 1/k-1/k*cos(th);
        0 1 0 0;
        -sin(th) 0 cos(th) 1/k*sin(th);
        0 0 0 1];
    %Each transformation matrix is for the global frame. 
    if i>1
        T{i} = T{i-1}*T{i};
    end
    %Determines the global coordinates for the origin.
    globalPos{i}=T{i}(:,4);
    
end



end