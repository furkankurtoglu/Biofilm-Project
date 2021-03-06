close all
clear
clc
%%cd C:\Users\Furkan\Documents\GitHub\Biofilm-Project\PhysiCell_Folder
cd output\

%%
load('output00000013_microenvironment0.mat');
XPos = multiscale_microenvironment(1,:);
YPos = multiscale_microenvironment(2,:);
Oxygen = multiscale_microenvironment(5,:);
ECM = multiscale_microenvironment(6,:);
Glucose = multiscale_microenvironment(7,:);

figure(1)
plot3(XPos,YPos,Oxygen,'.')
figure(2)
plot3(XPos,YPos,ECM,'.g')
figure(3)
plot3(XPos,YPos,Glucose,'.r')

cd ..