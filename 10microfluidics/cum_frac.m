%Generates the cumulative fraction curves from the microfluidic experiments
%written by Roberto Rusconi and Karthik Sekar (karsekar@gmail.com)
%last updated 25.09.2017

clear all;
close all;

files = {'every3minutes.csv','every4minutes.csv','every5minutes.csv'}

figure;
hold on;

for i = 1:length(files)
    inputTable = readtable(strcat('cum_frac_data/',files{i}),'Delimiter',';');
    plot(inputTable.time, inputTable.cumfrac);
end

legend('every 3 minutes','every 4 minutes','every 5 minutes');
xlabel('Time from start of feeding (min)');
ylabel('Cumulative fraction of divided cells');