%Analysis of real-time metabolomics under antibiotics at f=0.18
% by Karthik Sekar (karsekar@gmail.com)
% last updated 25.04.2017
% written for Matlab 2015b


clear all;
close all;

addpath('../common');


%designate the data sets
datasets = {'Chloramphenicol','Rifamycin','AZT'};

%%load the data

%load the ion table
load('intersection_data.mat');

%load the chlorampenicol data
load('outputdata36and37.mat');
data36and37 = ztransform(outputdatamat',1:60);

%load the rifamycin data
load('outputdata39and40.mat');
data39and40 = ztransform(outputdatamat',1:60);

%load the AZT data
load('dataE221615_005.mat');
dataE221615 = ztransform(myvar.data,1:60);


%filter parameters
b = [0.333 0.334 0.333];
a = 1;

%%designate the ions of interest

%glutamate, phe, guanine, and thymine
ofinterest = [1 44 47 70]

%(iso)leucine, glutamine, valine, and hypoxanthien
%ofinterest = [53 62 88 27]

%%Visualize the data

figure('Position', [100, 100, 600, 800]);
        
for i=1:length(ofinterest)
    j=ofinterest(i);
    subplot(2,2,i);
    plot(data36and37(:,intersectionoutput.x36and37(j)),'.','Color',	[1 144/255 1]);
    hold on;
    plot(data39and40(:,intersectionoutput.x39and40(j)),'.','Color',	[135/255 222/255 100/255]);
    plot(dataE221615(:,intersectionoutput.E221615(j)),'.','Color',	[255/255 102/255 0/255]);

    ax = gca;


    plot(data36and37(:,intersectionoutput.x36and37(j)),'.','Color',	[1 144/255 1]);
    plot(data39and40(:,intersectionoutput.x39and40(j)),'.','Color',	[135/255 222/255 100/255]);
    plot(dataE221615(:,intersectionoutput.E221615(j)),'.','Color',	[255/255 102/255 0/255]);

    filtered = filter(b,a,data36and37(:,intersectionoutput.x36and37(j)));
    plot(3:length(filtered),filtered(3:end),'-','Color',[1 144/255 1],'LineWidth',2);
    filtered = filter(b,a,data39and40(:,intersectionoutput.x39and40(j)));
    plot(3:length(filtered),filtered(3:end),'-','Color',[135/255 222/255 100/255],'LineWidth',2);
    filtered = filter(b,a,dataE221615(:,intersectionoutput.E221615(j)));
    plot(3:length(filtered),filtered(3:end),'-','Color',[255/255 102/255 0/255],'LineWidth',2);
    set(gca, 'XTickLabel', [], 'YTickLabel', []);

    title(intersectionoutput.FullName(j))
%
end


