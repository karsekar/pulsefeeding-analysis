%Analysis of real-time metabolomics under antibiotics at f=0.18
% by Karthik Sekar (karsekar@gmail.com)
% last updated 05.02.2018
% written for MATLAB 2015B


clear all;
close all;

addpath('../common');

%set the colors
colors = {[0/255 81/255 134/255], [213/255 155/255 15/255], [181/255 69/255 21/255]};
lightcolors = {[0/255 114/255 189/255],[237/255 177/255 32/255],[217/255 83/255 25/255]};



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

%load the negative control data
load('../04analyzeiondata_varyfreq/dataE221271_005.mat');
dataE221271 = ztransform(myvar.data,1:60);


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
    plot(data36and37(:,intersectionoutput.x36and37(j)),'.','Color',	lightcolors{1});
    hold on;
    plot(data39and40(:,intersectionoutput.x39and40(j)),'.','Color',	lightcolors{2});
    plot(dataE221615(:,intersectionoutput.E221615(j)),'.','Color',	lightcolors{3});

    filtered = filter(b,a,data36and37(:,intersectionoutput.x36and37(j)));
    plot(3:length(filtered),filtered(3:end),'-','Color',colors{1},'LineWidth',2);
    filtered = filter(b,a,data39and40(:,intersectionoutput.x39and40(j)));
    plot(3:length(filtered),filtered(3:end),'-','Color',colors{2},'LineWidth',2);
    filtered = filter(b,a,dataE221615(:,intersectionoutput.E221615(j)));
    plot(3:length(filtered),filtered(3:end),'-','Color',colors{3},'LineWidth',2);
    set(gca, 'XTickLabel', [], 'YTickLabel', []);

    title(intersectionoutput.FullName(j))
%
end

%make the negative control plots
ofinterest = [1 44 47 70];

%provide the same indices in the negative control
nc_indices = [44 60 48 23];

%set the times. the NC data set was sampled every 10 s, but the antibiotics
%were samples every 15 s. Need separate time vectors
time1 = [1:length(data36and37(:,1))]*15;
time2 = [1:length(dataE221271(:,1))]*10;
timeRIF = [1:length(data39and40(:,1))]*15;
timeAZT = [1:length(dataE221615(:,1))]*15;

%make the plots
figure('Position', [200, 200, 1000, 800]);
for i=1:length(ofinterest)
    j=ofinterest(i);

    %first to CHP
    subplot(3,length(ofinterest),i);
    plot(time1,data36and37(:,intersectionoutput.x36and37(j)),'.','Color',	lightcolors{1});
    hold on;

    plot(time2,dataE221271(:,nc_indices(i)),'.','Color',	[0.5 0.5 0.5]);

    filtered = filter(b,a,dataE221271(:,nc_indices(i)));
    plot(time2(3:end),filtered(3:end),'-','Color',[0.3 0.3 0.3],'LineWidth',2);
    filtered = filter(b,a,data36and37(:,intersectionoutput.x36and37(j)));
    plot(time1(3:end),filtered(3:end),'-','Color',colors{1},'LineWidth',2);

    %then do RIF
    subplot(3,length(ofinterest),i+length(ofinterest));
    plot(timeRIF,data39and40(:,intersectionoutput.x39and40(j)),'.','Color',	lightcolors{2});
    hold on;

    plot(time2,dataE221271(:,nc_indices(i)),'.','Color',	[0.5 0.5 0.5]);

    filtered = filter(b,a,dataE221271(:,nc_indices(i)));
    plot(time2(3:end),filtered(3:end),'-','Color',[0.3 0.3 0.3],'LineWidth',2);
    filtered = filter(b,a,data39and40(:,intersectionoutput.x39and40(j)));
    plot(timeRIF(3:end),filtered(3:end),'-','Color',colors{2},'LineWidth',2);

    %then do AZT
    subplot(3,length(ofinterest),i+length(ofinterest)*2);
    plot(timeRIF,dataE221615(:,intersectionoutput.E221615(j)),'.','Color',	lightcolors{3});
    hold on;

    plot(time2,dataE221271(:,nc_indices(i)),'.','Color',	[0.5 0.5 0.5]);


    filtered = filter(b,a,dataE221271(:,nc_indices(i)));
    plot(time2(3:end),filtered(3:end),'-','Color',[0.3 0.3 0.3],'LineWidth',2);
    filtered = filter(b,a,dataE221615(:,intersectionoutput.E221615(j)));
    plot(timeRIF(3:end),filtered(3:end),'-','Color',colors{3},'LineWidth',2);

end


