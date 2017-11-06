% Looking at DNA distribution of cells over time and parsing the counts
% from each DNA dist
% by Karthik Sekar (karsekar@gmail.com)
% last updated 08.08.2017
% written for Matlab 2015b
%
% NOTE: Data must be downloaded from https://doi.org/10.5281/zenodo.1035825

close all;
clear all;

addpath('../common');

%filter parameters
b = [0.333 0.334 0.333];
a = 1;

%set the pulsing period
period = 2.75;

%read the summary table
inTable = readtable('20170718.csv','Delimiter',';');

%convert the range to numbers
for i=1:height(inTable)
    inTable.range{i} = str2num(inTable.range{i});
end

odData = table2array(inTable(1:end,1:2));


%make the color vector. red to blue
scalingt = (inTable.time)/330*17+1;
colorVector = [((scalingt)*10.545-9.1458)/255 ((scalingt)*-1.451+95.156)/255 ((scalingt)*-7.9286 + 163.83)/255]*1.4;

%initalize variables
datafolder = 'data20170718';

prefixes = inTable.prefix;
inTable.lowDNADist = cell(height(inTable),1);
inTable.highDNADist = cell(height(inTable),1);

counts = cell(1,length(prefixes));
indices = cell(1,length(prefixes));
percents = cell(1,length(prefixes));
integrateddna = zeros(1,length(prefixes));

ofinterest = [1 13 16];

ratios = [];

%define the gate
v = [8.5 8; 14 14; 14 4; 8.5 4];

%show the blank
subplot(1,2,1);

filetoopen = ['./' datafolder '/blank.csv'];
data = readtable(filetoopen);

scatter(log(data.FL1_A),log(data.FSC_A),'filled','MarkerEdgeColor','none','SizeData',2,'MarkerFaceColor','k');
alpha(0.1);
xlabel('log Fluorescence (DNA)');
ylabel('log Forward Scattering');

axis([4 15 4 14]);

title('blank');

%show the gate with cells

subplot(1,2,2);

filetoopen = ['./' datafolder '/time1200_1.csv'];
data = readtable(filetoopen);


scatter(log(data.FL1_A),log(data.FSC_A),'filled','MarkerEdgeColor','none','SizeData',2,'MarkerFaceColor','k');
inred = inpolygon(log(data.FL1_A), log(data.FSC_A), v(:,1), v(:,2));
alpha(0.1);
patch(v(:,1), v(:,2), 'red', 'EdgeColor', 'none','FaceAlpha',.2);
hold on;
title('Sample with live cells');
xlabel('log Fluorescence (DNA)');
ylabel('log Forward Scattering');

axis([4 15 4 14]);

%initialize

maxes = [];
expected = [];
sizes = [];
distcounts = [];
runningcounts = [];

h = figure;

currentfig = figure;

current2 = figure;

ofintfig = figure;

filelist = cellstr(ls(['./' datafolder '/']));

%go through each time point
for i=1:length(prefixes)
    current = prefixes{i};
    idx = strfind(filelist,current);
    indices{i} = find(~cellfun(@isempty,idx))';
    dnatotal = [];
    x = [0:1e4:exp(v(2,1))];
    
    %go through each file for each time point
    for j=indices{i}
        %open the CSV file and apply the gate to get DNA fluorescence
        %information
        filename = filelist{j};
        filetoopen = ['./' datafolder '/' filename];
        data = readtable(filetoopen);
        inred = inpolygon(log(data.FL1_A), log(data.FSC_A), v(:,1), v(:,2));
        tobepcad = horzcat(data.FL1_A(find(inred)),log(data.FSC_A(find(inred))));
        result = tobepcad*[1; 0];
        dnatotal = [dnatotal result'];
        counts{i} = [counts{i} length(result)];
        percents{i} = [percents{i} length(result)/height(data)];
        totalevents(i) = height(data);

    end
    
    %now split the distributions for each time point
    clf(currentfig);
    figure(currentfig);

    x = 0:1e4:4e5;
    ncounts = histcounts(dnatotal, x) * 1/length(indices{i}); %normalize to the number of technical replicates
    runningcounts(i) = length(dnatotal) * 1/length(indices{i});
    
    rangeOfInt = inTable.range{i}; %get the rough range for the medium and high DNA distribution

    %first find the medium DNA range
    range =  find(x >= rangeOfInt(1) & x <= rangeOfInt(2)); %bounds of where it is
    plot(x(1:end-1),1/sum(ncounts)*ncounts,'.b','MarkerSize',10); 
    hold on;
    queryvals = result(find(result >= x(range(1)) & result <= x(range(end))));
    bParams = [0.1 1e5 4e4];
    mdl = fitnlm(x(range),1/sum(ncounts)*ncounts(range),@normdist,bParams); %fit a normal distribution to it
    fitting = normdist(mdl.Coefficients.Estimate,x);
    inTable.lowDNADist{i} = mdl.Coefficients.Estimate; %save the fitting parameters
    inTable.lowPercent(i) = sum(fitting); %save the fraction of cells that are in medium DNA distribution
    plot(x,fitting ,'-c');

    plot(x(range),1/sum(ncounts)*ncounts(range),'.r','MarkerSize',10);
    title(current);
    
    %find higher DNA distribution
    clf(current2);
     figure(current2);
    
    %subtract the medium DNA distribution
    subtracted = ncounts-fitting(1:end-1)*sum(ncounts);
    
    %get the new range for the high DNA and repeat fitting process
    range2 =  find(x >= rangeOfInt(3) & x <= rangeOfInt(4));
    plot(x(1:end-1),subtracted * 1/sum(ncounts)  ,'m.','MarkerSize',10);
    hold on;
    plot(x(range2),subtracted(range2) * 1/sum(ncounts) ,'.b','MarkerSize',10);

    bParams = [0.1 2.2e5 3e4];
    mdl2 = fitnlm(x(range2),subtracted(range2)  * 1/sum(ncounts) ,@normdist,bParams);
    fitting2 = normdist(mdl2.Coefficients.Estimate,x);
    inTable.highDNADist{i} = mdl2.Coefficients.Estimate;
    inTable.highPercent(i) = sum(fitting2);

   
    plot(x,fitting2 ,'-g')

    title(current);
    hold off;
    
    figure(currentfig);
    plot(x(range2),fitting2(range2)+fitting(range2) ,'-k');
          
    x = 8.5:0.1:14;
    xf = 4:0.1:14;
    ncounts = histcounts(log(dnatotal), x)*1/length(indices{i});
    
    %generate the 3D plot of DNA distribution over time
      figure(h);
      filtered = filter(b,a,ncounts); %use a smoothed histogram
      plot3(x(3:end-1),odData(i,1)*ones(1,length(ncounts)-2),filtered(3:end),'Color',colorVector(i,:));
      hold on;

      integrateddna(i) = sum(dnatotal)*1/length(indices{i});

    %plot DNA distribution of interest
      if length(find(ofinterest == i))
        figure(ofintfig);
        filtered = filter(b,a,ncounts);
        plot(x(3:end-1),filtered(3:end),'Color',colorVector(i,:));

        hold on;
      end
end

xlabel('log DNA fluorescence');
ylabel('Time (min)');
zlabel('Count');

%plot the dist counts

meancounts = cellfun(@mean, counts);
firstthree = mean(totalevents(1:3));
initialOD = mean(inTable.od(1:4));

percent_array = cellfun(@mean, percents);

figure;
plot(inTable.time,inTable.highPercent' .* percent_array*firstthree/initialOD,'.');
title('DNA Dist count');

hold on;;
plot(inTable.time,(inTable.lowPercent)' .* percent_array*firstthree/initialOD,'.');
legend('High DNA','Medium DNA');

