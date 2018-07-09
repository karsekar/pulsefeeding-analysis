%Analysis to compare perturbations on lag time
% by Karthik Sekar (karsekar@gmail.com)
% last updated 21.03.2018
% written for Matlab 2015b

close all;
clear all;

%choose the experiment and comment as necessary

%for the Protease inhibitor experiment
colors = {[213/255 155/255 15/255],[0/255 114/255 189/255]}
filetoopen = './data/01piaddition.csv';
WTindex = 2;
timemax = 300;

%for the PdhR titration experiment
colors = {[213/255 155/255 15/255],[217/255 83/255 25/255],[0/255 114/255 189/255]}
filetoopen = './data/02pdhRtitration.csv';
WTindex = 3;
timemax = 300;

%for the FtsZ titration experiment
colors = flip({[237/255 177/255 32/255],[217/255 83/255 25/255],[0/255 114/255 189/255]});
filetoopen = './data/03ftsZtitration.csv';
WTindex = 1;
timemax = 300;

%for the ClpX titration experiment
colors = flip({[237/255 177/255 32/255],[217/255 83/255 25/255],[0/255 114/255 189/255]});
filetoopen = './data/04clpXtitration.csv';
WTindex = 1;
timemax = 300;


%for the glycerol experiment
% colors = {[217/255 83/255 25/255],[0/255 114/255 189/255],[237/255 177/255 32/255]};
% filetoopen = './data/05glycerol.csv';
% WTindex = 1;
% timemax = 200;

%for the acetate experiment
% colors = {[237/255 177/255 32/255],[0/255 114/255 189/255]};
% filetoopen = './data/06acetate.csv';
% WTindex = 1;
% timemax = 200;

%for the nitrogen experiment
% colors = {[217/255 83/255 25/255],[237/255 177/255 32/255],[0/255 114/255 189/255]};
% filetoopen = './data/07nitrogen.csv';
% WTindex = 1;
% timemax = 400;

%for the ftsA experiment
% colors = {[0/255 114/255 189/255],[237/255 177/255 32/255]};
% filetoopen = './data/08ftsA.csv';
% WTindex = 1;
% timemax = 360;
% 
% %for the ftsL experiment
% colors = {[0/255 114/255 189/255],[217/255 83/255 25/255],[237/255 177/255 32/255]};
% filetoopen = './data/09ftsL.csv';
% WTindex = 1;
% timemax = 360;
% 
% %for the ftsB experiment
% colors = {[0/255 114/255 189/255],[217/255 83/255 25/255],[237/255 177/255 32/255]};
% filetoopen = './data/10ftsB.csv';
% WTindex = 1;
% timemax = 360;
% 
% %for the ftsN experiment
% colors = {[0/255 114/255 189/255],[217/255 83/255 25/255],[237/255 177/255 32/255]};
% filetoopen = './data/11ftsN.csv';
% WTindex = 1;
% timemax = 360;

% %for the ftsZ experiment to shift the critical rate
% colors = {[0/255 114/255 189/255],[217/255 83/255 25/255],[237/255 177/255 32/255]};
% filetoopen = './data/12ftsZshift.csv';
% WTindex = 1;
% timemax = 360;

%for the ftsZ natural starvation experiments
colors = {[0/255 114/255 189/255],[237/255 177/255 32/255]};
filetoopen = './data/13ftsZnatural1.csv';
WTindex = 1;
timemax = 250;

%for the ftsZ natural starvation experiments
colors = {[0/255 114/255 189/255],[237/255 177/255 32/255]};
filetoopen = './data/14ftsZnatural2.csv';
WTindex = 1;
timemax = 250;


beta0_200 = [50, 0.07, 1];
beta0_100 = [25, 0.07, 1];
beta0_0 = [0, 0.07, 1];

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
 
%make all of the individual plots 
for i=2:length(data(1,:));
   
    %remove points that are 3 Z above the mean over all of the data
    outlierthreshold = means(i) + 3*stddevs(i)
    outliers  = find(data(:,i) > outlierthreshold);
    numOfOutliers(i-1) = length(outliers);
    time_temp = time;
    time_temp(outliers) = [];
    data_temp = data(:,i);
    data_temp(outliers) = [];
    
    
    %try different initial conditions
    mdl_0 = fitnlm(time_temp,data_temp,@threshold,beta0_0);
    mdl_100 = fitnlm(time_temp,data_temp,@threshold,beta0_100);
    mdl_200 = fitnlm(time_temp,data_temp,@threshold,beta0_200);

    %take the best one
    [dontneed, index] = max([mdl_0.Rsquared.Ordinary, mdl_100.Rsquared.Ordinary, mdl_200.Rsquared.Ordinary]);

    if index == 0
        mdl = mdl_0;
    elseif index == 1
        mdl = mdl_100;
    else
        mdl = mdl_200;
    end
    
    fitting = threshold(mdl.Coefficients.Estimate,time);
    
    plot(time_temp,data_temp-fitting(1),'.','Color',colors{i-1});
    hold on;
        
    h(i-1) = plot(time,fitting-fitting(1),'Color',colors{i-1});
       
    feedrate = 0.9*5/2500*1/180*1000*1/(mdl.Coefficients.Estimate(3)*0.4)*60/(mean(diff(time)));
    if WTindex == i-1
        wtf = feedrate;
    end
    WTlagtime = 288.15*exp(-4.75*(feedrate - 0.19943));

    
end

ylabel('OD - ODi');
xlabel('Time (min)');

legend(h,conditions);

limits = ylim;

numOfOutliers

%lag time using the empirical fit of the lag time versus feedrate data.
%NOTE this is only accurate for the glucose-feeding conditions
WTlagtime = 288.15*exp(-4.75*(wtf - 0.19943));
plot([WTlagtime WTlagtime],[limits(1)+0.01 limits(2)-0.01],'--k');


