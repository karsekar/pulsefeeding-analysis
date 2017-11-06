%Analysis to compare perturbations on lag time
% by Karthik Sekar (karsekar@gmail.com)
% last updated 06.09.2017
% written for Matlab 2015b

close all;
clear all;

%parameters
timemax = 250;

beta0_200 = [100, 0.07, 1];
beta0_100 = [50, 0.07, 1];
beta0_0 = [0, 0.07, 1];

colors = {[213/255 155/255 15/255],[0/255 114/255 189/255],[217/255 83/255 25/255]}


files = {'20160914crp.csv', '20161207crp.csv', '20160929pdhR.csv'}
addpath('../common');

feedrates = zeros(1,length(files));
handles = zeros(1,length(files));

figure;

hold on;

for k=1:length(files)
    
    
    
    filetoopen = files{k};
    
    data = readtable(filetoopen);
    
    time = data.Time; 
    data = table2array(data);
    
    %take only first points defined in the time frame
    timepoints = find(time < timemax);
    data = data(timepoints,2);
    time = time(timepoints);

    %do not consider outlier 3Z away
    outlierthreshold = mean(data) + 3*std(data)
    outliers  = find(data > outlierthreshold);
    numOfOutliers(k) = length(outliers);
    time_temp = time;
    time_temp(outliers) = [];
    data_temp = data;
    data_temp(outliers) = [];

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
    
    fitted = threshold(mdl.Coefficients.Estimate,time_temp);
    hold on;
    

    feedrates(k) = 0.9*5/2500*1/180*1000*1/(mdl.Coefficients.Estimate(3)*0.4)*60/(mean(diff(time)));
    handles(k) = plot(time_temp,(data_temp-fitted(1))/fitted(1),'.','Color',colors{k});

    plot(time_temp,(fitted-fitted(1))/fitted(1),'Color',colors{k});
end

ylabel('(OD - ODi)/ODi');
xlabel('Time (min)');

legend(handles, strread(num2str(feedrates),'%s'));
 
 
 