%Shape Simulation from randomly created wavelength shifts. This function
%uses the fsolve method in matlab to solve the nonlinear equations for 
%determining the curvature and orientation of the shape. Next is an 
%interpolation of constant curvature in order to recreate the shape.
%%
function shapeSim(shiftCell)
%close all
%Used to create random shift in wavelengths based on an average shift.
addpath('curvemethod','velocitymethod')
avgShift = 1.5*10^-10;
randMatrix = rand(3,4);

if nargin==0
%Generates the random shifts using the random matrix and average shift.
shiftCell = wavelengthShift(randMatrix,avgShift);
end

[T,globalPos] = velocitymethod(shiftCell);

%[T,globalPos]=curvemethod(shiftCell,'PlotOn');


simulationPlot(T,globalPos,'gCoordinates')

figure(3)
end

function [shiftCell] = wavelengthShift(randMatrix,avgShift)
%Given the random shift matrix and the average shift, returns a cell matrix
%with each element being a 3x1 array with the wavelength shifts for a,b,c.
%%
shiftCell=cell(length(randMatrix),1); %Creates empy cell array

%Iterates through the shift matrix multiplying the values by the average
%shift value and storing the results in the cell matrix. 
for i=1:length(shiftCell)
    sLam_a = avgShift*randMatrix(1,i);
    sLam_b = avgShift*randMatrix(2,i);
    sLam_c = avgShift*randMatrix(3,i);
    shiftCell{i} = [sLam_a, sLam_b, sLam_c];
end
end


function simulationPlot(T,globalPos,origin,curvature,orientation)
%Determines how the figure will look and function. Pushbutton allows for
%each new step to be determined.
CoordinateFig = 2;
figure(CoordinateFig)
clf
daspect([1, 1, 1])
hold on
grid on
view(3)
axis([-0.5 0.5 -0.5 0.5 0 4])
set(CoordinateFig,'toolbar','figure','Position',[214   405   672   504])
uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'Position', [20 20 50 20],...
    'Callback', 'uiresume(gcbf)');


%Plots the origin coordinates for a general frame of reference.
if strcmp('gCoordinates',origin)
    origin = [1 0 0 0;
        0 1 0 0;
        0 0 1 0];
    for i = 1:3
        x(i) = origin(1,4);
        y(i) = origin(1,4);
        z(i) = origin(1,4);
    end
    u = origin(1,1:3);
    v = origin(2,1:3);
    w = origin(3,1:3);
    quiver3(x,y,z,u,v,w,0.3)
end

%Plots the coordinates of each reference frame at the global position.
trans=eye(3);
[preP1, preP2, preP3] = deal(0);
for j = 1:length(T)
    %uiwait(gcf)
    m = T{j};
    p = globalPos{j};
    for i=1:3
        x(i) = p(1);
        y(i) = p(2);
        z(i) = p(3);
    end
    u = m(1,1:3);
    v = m(2,1:3);
    w = m(3,1:3);
    
    quiver3(x,y,z,u,v,w,0.3)
    
    %Plots a line between each coordinate frame's origin.
    %     t=50;
    %     q1 = linspace(preP1,p(1),t);
    %     q2 = linspace(preP2,p(2),t);
    %     q3 = linspace(preP3,p(3),t);
    %     preP1 = p(1);
    %     preP2 = p(2);
    %     preP3 = p(3);
    %     plot3(q1,q2,q3)
    
    %%
%     %Plots curve between each coordinate
%     k=curvature(j);
%     phi = orientation(j);
%     
%     %Creates the curve for the length of segment
%     t=linspace(0,k);
%     curveX=1/k-1/k*cos(t);
%     curveY=0;
%     curveZ=1/k*sin(t);
%     
%     %Rotation around the z axis
%     zRot = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
%     
%     %Rotation of the curve into the global orientation
%     if j>1
%         trans = T{j-1}(1:3,1:3);
%     end
%     
%     %Rotating each point into the global orientation
%     for i=1:length(curveX)
%         myMatrix{i}= trans*zRot*[curveX(i); 0; curveZ(i);];
%     end
%     
%     %Separating points into separate coordinates
%     for i=1:length(myMatrix)
%         rotCurveX(i)=myMatrix{i}(1);
%         rotCurveY(i)=myMatrix{i}(2);
%         rotCurveZ(i)=myMatrix{i}(3);
%     end
%     
%     %Adding the previous frame's origin in global coordinates.
%     rotCurveX=rotCurveX+preP1;
%     rotCurveY=rotCurveY+preP2;
%     rotCurveZ=rotCurveZ+preP3;
%     
%     %Plotting of the curve
%     plot3(rotCurveX,rotCurveY,rotCurveZ)
%     
%     %Sets the global coordinate of the current frame's origin.
%     preP1 = p(1);
%     preP2 = p(2);
%     preP3 = p(3);
    
    %%
end
end

%% Slider Control
function aShift(hObj,~,shiftMatrix)
% Called to set zlim of surface in figure axes
% when user moves the slider control
val = get(hObj,'Value');
shiftMatrix(1,1) = val;
curvemethod(shiftMatrix);
end

function bShift(hObj,~,shiftMatrix)
% Called to set zlim of surface in figure axes
% when user moves the slider control
val = get(hObj,'Value');
shiftMatrix(2,1) = val;
curvemethod(shiftMatrix);
end

function cShift(hObj,~,shiftMatrix)
% Called to set zlim of surface in figure axes
% when user moves the slider control
val = get(hObj,'Value');
shiftMatrix(3,1) = val;
curvemethod(shiftMatrix);
end

function uiPlotCurve(shiftMatrix)
uicontrol('Style', 'slider',...
    'Min',0,'Max',1,'Value',0.5,...
    'Position', [500 40 120 20],...
    'Callback', {@aShift,shiftMatrix});
uicontrol('Style', 'slider',...
    'Min',0,'Max',1,'Value',0.5,...
    'Position', [500 20 120 20],...
    'Callback', {@bShift,shiftMatrix});
uicontrol('Style', 'slider',...
    'Min',0,'Max',1,'Value',0.5,...
    'Position', [500 0 120 20],...
    'Callback', {@cShift,shiftMatrix});
uicontrol('Style', 'popup',...
    'String', '1|2|3|4',...
    'Position', [20 400 100 50],...
    'Callback', @setmap);

end