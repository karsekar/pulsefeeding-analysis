%Analysis of long term experiment
% by Karthik Sekar (karsekar@gmail.com)
% last updated 08.08.2017
% written for Matlab 2015b

close all;
clear all;


%for the long term experiment
colors = {[213/255 213/255 213/255]}
filetoopen = './data/longterm.csv';
WTindex = 1;
feedrate = 0.4; %units mmol/g/h
timemax = 720;

addpath('../common');

data = readtable(filetoopen);
conditions = data.Properties.VariableNames(2:end);
h = [];
time = data.Time; 

%take only first points defined in the time frame

data = table2array(data);
%get the std deviation and mean for outlier removal

stddevs = std(data);
means = mean(data);

timepoints = find(time < timemax);
data = data(timepoints,:);
time = time(timepoints);
 
i=2;

%remove points that are 3 Z above the mean over all of the data
outlierthreshold = means(i) + 3*stddevs(i)
outliers  = find(data(:,i) > outlierthreshold);
numOfOutliers(i-1) = length(outliers);
time_temp = time;
time_temp(outliers) = [];
data_temp = data(:,i);
data_temp(outliers) = [];

beta0_200 = [50, 0.07, 1];
beta0_100 = [25, 0.07, 1];
beta0_0 = [0, 0.07, 1];

mdl_0 = fitnlm(time_temp,data_temp,@threshold,beta0_0);
mdl_100 = fitnlm(time_temp,data_temp,@threshold,beta0_100);
mdl_200 = fitnlm(time_temp,data_temp,@threshold,beta0_200);

[dontneed, index] = max([mdl_0.Rsquared.Ordinary, mdl_100.Rsquared.Ordinary, mdl_200.Rsquared.Ordinary]);

if index == 0
    mdl = mdl_0;
elseif index == 1
    mdl = mdl_100;
else
    mdl = mdl_200;
end


plot(time_temp,data_temp,'.','Color',colors{i-1},'MarkerSize',10);
hold on;


WTlagtime = 288.15*exp(-4.75*(feedrate - 0.19943));
limits = ylim;

fittime = 1:timemax;
plot(fittime,threshold([WTlagtime; 0.0975*feedrate/60*mdl.Coefficients.Estimate(3) ;mdl.Coefficients.Estimate(3)],fittime),'k');





ylabel('OD');
xlabel('Time (min)');

legend(h,conditions);


numOfOutliers




