%Parsing the counts from the flow cytometry data
% by Karthik Sekar (karsekar@gmail.com)
% last updated 08.08.2017
% written for Matlab 2015b
%
% NOTE: data must be downloaded from https://doi.org/10.5281/zenodo.1035825

close all;
clear all;

addpath('../common');


%define experiments

experiments = cell(0);

%first experiment

experiments{1}=struct();
experiments{1}.csvtable = '20170713.csv';
experiments{1}.datafolder = 'data20170713';
experiments{1}.color = [0.6 0.6 0.6];
experiments{1}.period = 4.5;

%second experiment

experiments{2}=struct();
experiments{2}.csvtable = '20170718.csv';
experiments{2}.datafolder = 'data20170718';
experiments{2}.color = [0.6 0.6 0.6];
experiments{2}.period = 2.75;

%define the gate
v = [8.5 8; 14 14; 14 4; 8.5 4];


figure;
hold on;

%iterate through each experiment
for k=1:length(experiments)
    subplot(1,length(experiments),k);
    period = experiments{k}.period;

    %read the summary file for each experiment
    inTable = readtable(experiments{k}.csvtable,'Delimiter',';');
    datafolder = experiments{k}.datafolder;


    odData = table2array(inTable(1:end,1:2));

    prefixes = inTable.prefix;

    %initialize the metrics
    counts = cell(1,length(prefixes));
    totalevents = zeros(1,length(prefixes));
    indices = cell(1,length(prefixes));
    percents = cell(1,length(prefixes));
 

    filelist = cellstr(ls(['./' datafolder '/']));

    for i=1:length(prefixes)
        current = prefixes{i};
        idx = strfind(filelist,current);
        indices{i} = find(~cellfun(@isempty,idx))';
        for j=indices{i}
            x = [exp(v(1,1)):1e4:exp(14)];
            filename = filelist{j};
            filetoopen = ['./' datafolder '/' filename];
            data = readtable(filetoopen);
            inred = inpolygon(log(data.FL1_A), log(data.FSC_A), v(:,1), v(:,2));
            tobepcad = horzcat(data.FL1_A(inred),log(data.FSC_A(find(inred))));
            result = tobepcad*[1; 0];
            ncounts = histcounts(result, x);
            counts{i} = [counts{i} sum(ncounts)];
            percents{i} = [percents{i} sum(ncounts)/height(data)];
            totalevents(i) = height(data);
        end

    end

    
    meancounts = cellfun(@mean, counts);
    firstthree = mean(totalevents(1:3));
    initialOD = mean(inTable.od(1:4));
    feedrate = 2.5*1/180*1000*1/((initialOD)*0.4)*60/(period)* 23/32000;
    
    percent_array = cellfun(@mean, percents);

    errorbar(inTable.time,percent_array*firstthree/initialOD,cellfun(@stderr, percents)*firstthree/initialOD,'.','MarkerSize',2,'Color',experiments{k}.color);

    ylim([1.6e4 2.2e4]);
    ylabel('Counts per 10 nL/ODi');
    xlim([-50 350]);
    WTlagtime = 288.15*exp(-4.75*(feedrate - 0.19943));
    hold on;
    plot([WTlagtime WTlagtime],[1.6e4 1.9e4],'--','Color',[0 0 0]);

    title(num2str(feedrate));
end

