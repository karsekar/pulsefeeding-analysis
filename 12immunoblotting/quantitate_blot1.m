%Code to quantitate bands on immunoblot
% by Karthik Sekar (karsekar@gmail.com)
% last updated 21.03.2018
% written for Matlab 2015b
close all;
clear all;

I = imread('blot1.bmp');
imshow(I);

hold on;

width = 90;
height = 45;

%define the coordinates
nc_coor = [550 707];
coordinates = {[650 707], [741 707], [832 707], [923 707], [1014 707]};
totals  = zeros(1,5);

patch([nc_coor(1) nc_coor(1)+width nc_coor(1)+width nc_coor(1)],[nc_coor(2) nc_coor(2) nc_coor(2)+height nc_coor(2)+height],'b','EdgeColor','b','FaceColor','none','LineStyle',':','LineWidth',2);
subimage = I(nc_coor(2):(nc_coor(2)+height),nc_coor(1):(nc_coor(1)+width));
nctotal = sum(sum(255 - subimage));


for i=1:length(coordinates)
    patch([coordinates{i}(1) coordinates{i}(1)+width coordinates{i}(1)+width coordinates{i}(1)],[coordinates{i}(2) coordinates{i}(2) coordinates{i}(2)+height coordinates{i}(2)+height],'b','EdgeColor','r','FaceColor','none','LineStyle',':','LineWidth',2);
    subimage = I(coordinates{i}(2):(coordinates{i}(2)+height),coordinates{i}(1):(coordinates{i}(1)+width));
    totals(i) = sum(sum(255 - subimage));
end

totals = totals - nctotal;

figure;
bar(totals);

%patch([659 731 731 659],[707 707 744 744],'r','FaceAlpha',0.5);
