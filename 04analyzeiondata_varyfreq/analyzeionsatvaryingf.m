%Analysis of real-time metabolomics under varied pulse frequencies
% by Karthik Sekar (karsekar@gmail.com)
% last updated 25.04.2017
% written for Matlab 2015b

clear all;
close all;

addpath('../common');


%designate the data sets and enumerate metadata of each
datasets = {'E221423','E221425','E221271'};
bedingungen = {'0.06 mmol/g DCW hr','0.12 mmol/g DCW hr', '0.18 mmol/g DCW hr'};
baselines = {[1:20],[1:20],[1:20]}; %baseline to normalize by
starts = {20,22,23}; %where the pulsing is anticipated to start. (in time index)
deltas = {15, 15, 10}; %time delta between measurements in real-time measurement. (in seconds)
maxtime = {40, 25, 30}; %duration of measurements (in minutes)
spanvector ={[1:7],[9:12],[14:18]}; %for subplotting 
linecolors = {[0/255 77/255 77/255],[0/255 77/255 77/255],[0/255 77/255 77/255],[242/255 208/255 208/255]};
pointcolors = {[185/255 255/255 255/255],[185/255 255/255 255/255],[185/255 255/255 255/255],[242/255 208/255 208/255]};

%filter parameters for the smoothing line
b = [0.333 0.334 0.333];
a = 1;

%%Choose which ions to look at, based off the unionoutputWNC.csv
ofinterest = [67 102 256 134];

%supplemental ccm
%ofinterest = [1 22 168]

%supplemental amino acids/bases
%ofinterest = [142 266 329 284]


%%load the table to correspond metabolites to particular ions
unionoutput = readtable('unionoutputWNC.csv','Delimiter',';');

datacell.maindata =cell(1,length(datasets));
datacell.baseline = cell(1,length(datasets));

%load the data
for i=1:length(datasets)
    datafile = strcat('data',datasets{i},'_005.mat');
    load(datafile);
    datacell.maindata{i} = myvar.data; %(2:59,:);
    datacell.baseline{i} = baselines{i};
    datacell.start{i} = starts{i}
    datacell.shift{i} = max(cell2mat(starts))-starts{i}
    datacell.maindata{i} = ztransform(datacell.maindata{i},datacell.baseline{i}); %z transform the data
end

%load the negative control data

datafile = 'ncdata005.mat';
load(datafile);
ncbaseline = [1:20];
negcondata = data;

%normalize the negative control data
negcondata = ztransform(negcondata,ncbaseline);


%%visualize the data

figure('Position', [100, 100, 400, 200*length(ofinterest)]);
%go through the ions of interest and plot 
for i=1:length(ofinterest)
    j = ofinterest(i);
    for k=1:length(datasets)
        subplot(length(ofinterest),sum(cell2mat(maxtime))/5, spanvector{k}+(i-1)*(sum(cell2mat(maxtime))/5));
        if not(isnan(unionoutput.(datasets{k})(j)))
            ymax = max(datacell.maindata{k}(:,unionoutput.(datasets{k})(j)));
            ymin = min(datacell.maindata{k}(:,unionoutput.(datasets{k})(j)));
            set1 = filter(b,a,datacell.maindata{k}(:,unionoutput.(datasets{k})(j)));

            if not(isnan(unionoutput.E221565(j)))
                ncset = filter(b,a,negcondata(:,unionoutput.E221565(j)));
                hold on;
                ymax = max([max(ncset) ymax]);
                ymin = min([min(ncset) ymin]);
                ylimits(1) = ymin*1.1;
                ylimits(2) = ymax*1.1;
                patch(([length(ncset) 0 3:length(ncset)])*15/60,[ylimits(1) ylimits(1) ncset(3:end)'],pointcolors{4},'LineStyle','none');
                plot(([3:length(ncset)])*15/60,ncset(3:end),'-','Color',[255/255 90/255 90/255],'LineWidth',1);
            end
            
            ylimits(1) = ymin*1.1;
            ylimits(2) = ymax*1.1;
             
            set(gca,'XTickLabel', [], 'YTickLabel', []);
            hold on;
            

            plot(([1:length(datacell.maindata{k}(:,unionoutput.(datasets{k})(j)))]+datacell.shift{k})*deltas{k}/60,datacell.maindata{k}(1:end,unionoutput.(datasets{k})(j)),'.','Color',pointcolors{k});
            plot(([3:length(set1)]+datacell.shift{k})*deltas{k}/60,set1(3:end),'-','Color',linecolors{k});

            set(gca, 'XTickLabel', [], 'YTickLabel', []);
            xlim([0 maxtime{k}+3]);
            ylim(ylimits);
        end
    end
end

