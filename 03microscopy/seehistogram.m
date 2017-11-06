%Histograms of cell lengths
% by Karthik Sekar (karsekar@gmail.com)
% last updated 08.08.2017
% written for Matlab 2015b

clear all;
close all;

ratio = 15.625; %ratio of pixels to microns

bins = [0:5:110] / ratio;

time0table = readtable('time0.xlsx');
time1table = readtable('time1.xlsx');
time2table = readtable('time2.xlsx');

time0counts = histcounts(time0table.length / ratio, bins);
time1counts = histcounts(time1table.length / ratio, bins);
time2counts = histcounts(time2table.length / ratio, bins);

figure; 
hold on;

plot(bins(1:end-1),time0counts / sum(time0counts),'.-','MarkerSize',10);
plot(bins(1:end-1),time1counts / sum(time1counts),'.-','MarkerSize',10);
plot(bins(1:end-1),time2counts / sum(time2counts),'.-','MarkerSize',10);

legend('0 minutes','150 minutes','250 minutes');
xlabel('Cell length (microns)');
ylabel('Norm. count');