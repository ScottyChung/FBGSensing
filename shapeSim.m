%Shape Simulation from randomly created wavelength shifts. This function
%uses the fsolve method in matlab to solve the nonlinear equations for 
%determining the curvature and orientation of the shape. Next is an 
%interpolation of constant curvature in order to recreate the shape.
%%
function shapeSim(shiftCell)
%close all
%Used to create random shift in wavelengths based on an average shift.
addpath('curvemethod','velocitymethod','plotting')
avgShift = 1.5*10^-9;
randMatrix = rand(3,4);

if nargin==0
%Generates the random shifts using the random matrix and average shift.
shiftCell = wavelengthShift(randMatrix,avgShift);
end

[transMatrix,globalPos] = velocitymethod(shiftCell);
simulationplot(transMatrix,globalPos,'gCoordinates')

[transMatrix,globalPos,curvature, orientation]=curvemethod(shiftCell,'PlotOn');
simulationplot(transMatrix,globalPos,'gCoordinates',curvature, orientation)



%figure(3)
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