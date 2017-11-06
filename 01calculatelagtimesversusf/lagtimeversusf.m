% Analyzing lag time versus feedrate
% by Karthik Sekar (karsekar@gmail.com)
% last updated: 25.4.2017
%
% Designed to work with MATLAB 2015b

clear all;
close all;

%load the raw data
addpath('../common');
schottlist = readtable('schottlist.csv','Delimiter',';');

platereaderlist = readtable('platereaderlist.csv','Delimiter',';');

%number of files
numOfEntries = height(schottlist);
numOfPlateEntries = height(platereaderlist);

feedrates = [];
delays = [];
slopes = [];
linearslopes = [];
legende = {};
balkenx = [];
balkeny = [];
sugarthres = [];
rsquarevals = [];
lineargr = [];
odi = [];

%from plate reader experiments
feedrates_plates = []; 
delays_plates = []; 
balkenx_plates = [];
balkeny_plates = [];


%initial parameters for the fitting
beta0_200 = [200, 0.07, 1];
beta0_100 = [100, 0.07, 1];
beta0_0 = [0, 0.07, 1];

figure('Position', [100, 100, 1200, 500]);


%parse the data
for i=1:numOfEntries
    %load the data file
    nameToLoad = schottlist.schott{i};
    data = dlmread(nameToLoad,';',3,0);
    
    %get meta information
    fid = fopen(nameToLoad);
    tline = strsplit(fgetl(fid),';');
    fclose(fid);
    data = data(:,1:2); %remove unnecessary info
    compile(i).magnitude = str2num(tline{2}); %parse stats about the specific run
    compile(i).period = str2num(tline{4});
    compile(i).starvation_time = tline{6};
    
    %calculate the feedrate for Schott bottle experiments
    
    % add glucose pulses to plot
    clf;
    minval = round(min(data(:,2)),2)-0.02;
    maxval = round(max(data(:,2)),2)+0.02;
    for j=0:compile(i).period:360   
    %    patch([j j j+1 j+1],[minval maxval maxval minval],[0.85 0.85 0.85],'LineStyle','none');
    end
    hold on;
    %for checking the individual plots
    plot(data(:,1),data(:,2), '.','MarkerSize',15, 'Color',[255/255 158/255 158/255] );
    
    
    %try various initial conditions to fit 
    mdl_0 = fitnlm(data(:,1),data(:,2),@threshold,beta0_0);
    mdl_100 = fitnlm(data(:,1),data(:,2),@threshold,beta0_100);
    mdl_200 = fitnlm(data(:,1),data(:,2),@threshold,beta0_200);

    [dontneed, index] = max([mdl_0.Rsquared.Ordinary, mdl_100.Rsquared.Ordinary, mdl_200.Rsquared.Ordinary]);

    if index == 0
        mdl = mdl_0;
    elseif index == 1
        mdl = mdl_100;
    else
        mdl = mdl_200;
    end

    rsquarevals(end+1) = mdl.Rsquared.Ordinary; 


    delay = mdl.Coefficients.Estimate(1);
    if delay < 0
        delay = 0;
    end
    
    feedrate = compile(i).magnitude*1/180*1000*1/(mdl.Coefficients.Estimate(3)*0.4)*60/compile(i).period * 23/32000;
    compile(i).feedrate = feedrate;
    feedrates(end+1) = feedrate;
    delays = [delays delay];
    sugarthres(end+1) = delay*feedrate;
    slopes(end+1) = mdl.Coefficients.Estimate(2);
    
    x = 0:360;
    plot(x,threshold(mdl.Coefficients.Estimate,x),'r','LineWidth',1 );
    title(strcat(num2str(feedrate),' mmol glc/g DCW hr'));
    xlabel('Time from the start of feeding (minutes)');
    ylabel('Measured OD600');
    hold off;
    totalgrowth = (mdl.Coefficients.Estimate(2)*60)*0.4*0.032;
    
    lineargr(end+1) = mdl.Coefficients.Estimate(2)/mdl.Coefficients.Estimate(3)*60;
    odi(end+1) = mdl.Coefficients.Estimate(3);
    if feedrate < 1
        balkenx(end+1) = mdl.Coefficients.Estimate(2)/mdl.Coefficients.Estimate(3)*60; %mdl.Coefficients.Estimate(2)*60;
        balkeny(end+1) = feedrate; %mdl2.Coefficients.Estimate(2);
    end
    
    title(num2str(feedrate));
end


for i=1:numOfPlateEntries
    filetoopen = platereaderlist.file{i};
    currentTable = readtable(filetoopen);
    
    data = currentTable.OD;
    time = currentTable.Time;

    mdl_0 = fitnlm(time,data,@threshold,beta0_0);
    mdl_100 = fitnlm(time,data,@threshold,beta0_100);
    mdl_200 = fitnlm(time,data,@threshold,beta0_200);
        
    
    [dontneed, index] = max([mdl_0.Rsquared.Ordinary, mdl_100.Rsquared.Ordinary, mdl_200.Rsquared.Ordinary]);

    if index == 0
        mdl = mdl_0;
    elseif index == 1
        mdl = mdl_100;
    else
        mdl = mdl_200;
    end
    
    rsquarevals(end+1) = mdl.Rsquared.Ordinary; 

    
    delay = mdl.Coefficients.Estimate(1);
    if delay < 0
        delay = 0;
    end
    delays_plates(end+1) = delay;
    
    feedrate = 0.9*5/2500*1/180*1000*1/(mdl.Coefficients.Estimate(3)*0.4)*60/(mean(diff(time)));
    feedrates_plates(end+1) = feedrate;

    % add glucose pulses to plot, if you want
    clf;
    minval = round(min(data),2)-0.02;
    maxval = round(max(data),2)+0.02;
    for j=0:mean(diff(time)):360   
    %    patch([j j j+1 j+1],[minval maxval maxval minval],[0.85 0.85 0.85],'LineStyle','none');
    end
    
    hold on;
    plot(time,data,'.', 'MarkerSize', 10, 'Color',[158/255 158/255 255/255] );
    
    
    plot(time,threshold(mdl.Coefficients.Estimate,time),'b-', 'LineWidth', 1 );
    hold off;
     
    totalgrowth = (mdl.Coefficients.Estimate(2)*60)*0.4*0.032;
    
    lineargr(end+1) = mdl.Coefficients.Estimate(2)/mdl.Coefficients.Estimate(3)*60;
    odi(end+1) = mdl.Coefficients.Estimate(3);

    
    if feedrate < 1
        balkenx_plates(end+1) = mdl.Coefficients.Estimate(2)/mdl.Coefficients.Estimate(3)*60;  %mdl.Coefficients.Estimate(2)*60;
        balkeny_plates(end+1) = feedrate; %mdl2.Coefficients.Estimate(2);
    end
    
    title(num2str(feedrate));
end

figure('Position', [100, 100, 900, 650]);

subplot(1,2,1);

mdl = fitnlm([feedrates feedrates_plates]',[delays delays_plates]',@thres_decay,[0.4 2 250])

xrange = [mdl.Coefficients.Estimate(1):0.01:1.5];
plot(xrange,thres_decay(mdl.Coefficients.Estimate,xrange),'--','LineWidth',2,'Color',[0/255 163/255 0/255]);
hold on;
patch([0 0 xrange],[0 mdl.Coefficients.Estimate(3) thres_decay(mdl.Coefficients.Estimate,xrange)'],[1 1 158/255],'LineStyle','none','FaceAlpha',0.4);
patch([1.5 xrange],[mdl.Coefficients.Estimate(3) thres_decay(mdl.Coefficients.Estimate,xrange)'],[158/255 1 158/255],'LineStyle','none','FaceAlpha',0.4);

h1 = plot(feedrates,delays,'.','MarkerSize',20,'Color',[255/255 65/255 65/255]);

h2 = plot(feedrates_plates,delays_plates,'.','MarkerSize',20,'Color',[65/255 65/255 255/255]);
legend([h1 h2],{'\fontsize{10} Spin flask', '\fontsize{10} Plate reader'},'Location','northeast');
%legend('boxoff');
ylim([0 mdl.Coefficients.Estimate(3)]);

xlabel('Time-integrated feedrate, \it f \rm (mmol/g/h)','FontSize',16);
ylabel('Lag time (min)','FontSize',16);
text(0.2,20,['\color{darkGreen} R^{2} = ' num2str(mdl.Rsquared.Ordinary,2)]);



subplot(1,2,2);

plot(balkenx,balkeny,'.','MarkerSize',20,'Color', [255/255 158/255 255/255]);
%legend('boxoff');

hold on;
plot(balkenx_plates,balkeny_plates,'.','MarkerSize',20,'Color', [255/255 218/255 159/255]);
xlabel('Growth rate after lag, \it \mu \rm (1/h)','FontSize',16);
ylabel('Time-integrated feedrate, \it f \rm (mmol/gDCW/h)','FontSize',16);
legend('\fontsize{10} Spin flask','\fontsize{10} Plate reader','Location','northwest');

xrange = 0:0.01:0.08;
axis([xrange(1) xrange(end) 0 1]);
fitflask = fitlm(balkenx,balkeny);
fitplate = fitlm(balkenx_plates,balkeny_plates);

y_schott = polyval(fliplr(fitflask.Coefficients.Estimate'), xrange);
plot(xrange,y_schott,'-','LineWidth',2,'Color',[255/255 65/255 255/255]);
y_plate= polyval(fliplr(fitplate.Coefficients.Estimate'), xrange);
plot(xrange,y_plate,'-','LineWidth',2,'Color', [255/255 183/255 65/255]);

text(0.03,0.12,'$$ f = \frac{1}{Y_{\frac{X}{S}}} \mu + m_{s} $$','Interpreter','latex','FontSize',12);
text(0.005,0.4,['\color{magenta} R^{2} = ' num2str(fitflask.Rsquared.Ordinary,2)]);
text(0.005,0.6,['$$  m_{s} = ' num2str(fitflask.Coefficients.Estimate(1),1) ' \pm' num2str(fitflask.Coefficients.SE(1),1) '$$'],'Interpreter','latex');
text(0.005,0.65,['$$  Y_{\frac{X}{S}} = ' num2str(1/fitflask.Coefficients.Estimate(2),1) ' \pm' num2str(fitflask.Coefficients.SE(2)/(fitflask.Coefficients.Estimate(2)^2),1) '$$'],'Interpreter','latex');

text(0.045,0.45,['\color{orange} R^{2} = ' num2str(fitplate.Rsquared.Ordinary,2)]);
text(0.045,0.6,['$$  m_{s} = ' num2str(fitplate.Coefficients.Estimate(1),1) ' \pm' num2str(fitflask.Coefficients.SE(1),3) '$$'],'Interpreter','latex');
text(0.045,0.65,['$$  Y_{\frac{X}{S}} = ' num2str(1/fitplate.Coefficients.Estimate(2),1) ' \pm' num2str(fitplate.Coefficients.SE(2)/(fitplate.Coefficients.Estimate(2)^2),1) '$$'],'Interpreter','latex');

%generate S6

figure;
plot(feedrates,delays/60.*feedrates,'.','MarkerSize',20,'Color',[255/255 65/255 65/255]);
hold on;
plot(feedrates_plates,feedrates_plates .* delays_plates/60,'.','MarkerSize',20,'Color',[65/255 65/255 255/255]);
xrange = [mdl.Coefficients.Estimate(1):0.01:1.5];
plot(xrange,xrange'.*thres_decay(mdl.Coefficients.Estimate,xrange)/60,'--','LineWidth',2,'Color',[0/255 163/255 0/255]);


%generate the summary table
whichsystem = cell(length([feedrates feedrates_plates]),1);
whichsystem(1:length(feedrates)) = {'Spin flask system'};
whichsystem(length(feedrates)+1:end) = {'Plate reader system'};
summaryT = table(whichsystem,odi',[feedrates feedrates_plates]', [delays delays_plates]', rsquarevals', lineargr', 'VariableNames', {'PulsingSystem','ODi','Feedrates','Lagtime','R2','GrowthRate'});
writetable(summaryT,'summaryinfo.csv');