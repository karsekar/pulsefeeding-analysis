%Showing the labeling data
% by Karthik Sekar (karsekar@gmail.com)
% last updated 25.09.2017
% written for Matlab 2015b

clear all;
close all;

data = readtable('labelingdata.csv','Delimiter',';');

feedrates = unique(data.feedrate);

%read the values

threonine_avgs = zeros(1,length(feedrates));
threonine_se = zeros(1,length(feedrates));
deoxy_avgs = zeros(1,length(feedrates));
deoxy_se = zeros(1,length(feedrates));

pvalue = 0;

%first load the data and calculate the avg/standard error
for i=1:length(feedrates)
    indices = find(data.feedrate == feedrates(i));
    threonine_avgs(i) = mean(data.threonine(indices));
    threonine_se(i) = std(data.threonine(indices))/sqrt(length(indices));
    deoxy_avgs(i) = mean(data.deoxyribose(indices));
    deoxy_se(i) = std(data.deoxyribose(indices))/sqrt(length(indices));
    if i>1 %calculate all of the adjacent p values using a Student's t test
        [bah pthr] = ttest2(data.threonine(indices), lastthr, 'Tail', 'right');
        [bah pdeoxy] = ttest2(data.deoxyribose(indices), lastdeoxy, 'Tail', 'right');
        pvalue = max([pvalue pthr pdeoxy]); %keep the highest
    end
    lastthr = data.threonine(indices);
    lastdeoxy = data.deoxyribose(indices);
end

%plot the data

figure('Position', [100, 100, 250, 800]);

subplot(2,1,1);

bar(feedrates,[threonine_avgs' [1-threonine_avgs']],0.4, 'stacked');
title('Threonine (Protein)');
ylabel('[+4] fraction');
ylim([0.0 0.045]);
xlim([-0.03 0.21]);
set(gca,'XTick',[])
hold on;
errorbar(feedrates,threonine_avgs,threonine_se,'.');


subplot(2,1,2);

bar(feedrates,[deoxy_avgs' [1-deoxy_avgs']],0.4, 'stacked');
title('Deoxyribose (DNA)');
ylabel('[+5] fraction');
ylim([0.0 0.045]);
xlim([-0.03 0.21]);
set(gca,'XTick',[])
hold on;
errorbar(feedrates,deoxy_avgs,deoxy_se,'.')

%max P value is...
pvalue