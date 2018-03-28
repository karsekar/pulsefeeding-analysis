%look at size and protein synthesis of non-dividing and dividing cells from
%microfluidic measurements
% by Karthik Sekar (karsekar@gmail.com) with preprocessing done by Roberto Rusconi
% last updated 26.03.2018

clear all;
close all;

%load all of the data
div_protein = readtable('protein_and_size_data/div_protein.csv');
nondiv_protein = readtable('protein_and_size_data/nondiv_protein.csv');

div_length = readtable('protein_and_size_data/div_length.csv');
nondiv_length = readtable('protein_and_size_data/nondiv_length.csv');

nondividetime = 90; %set a time for where the cells don't divide before


%indicate the time range
range = div_protein.time_min_';
ofinterest =  find(range < nondividetime);


%first make the protein synthesis figure

%make a patch to indicate the error of the protein synthesis in
%non-dividing cells
patchrange = [range flip(range)];
patchys = horzcat([std(nondiv_protein{:,2:end},0,2)/sqrt(5) + mean(nondiv_protein{:,2:end},2)]',flip([mean(nondiv_protein{:,2:end},2) - std(nondiv_protein{:,2:end},0,2)/sqrt(5)]'));

%calculate summary data 

pr_slopes = zeros(1,length(nondiv_protein{1,2:end}));

for i=1:length(nondiv_protein{1,2:end})
    P =  polyfit(range(ofinterest)', nondiv_protein{ofinterest,i+1},1);
    pr_slopes(i) =  P(1);
end

nondiv_mean = mean(pr_slopes);
nondiv_ste =  std(pr_slopes)/sqrt(length(pr_slopes));

pr_slopes = zeros(1,length(nondiv_protein{1,2:end}));

for i=1:length(div_protein{1,2:end})
    P =  polyfit(range(ofinterest)', div_protein{ofinterest,i+1},1);
    pr_slopes(i) =  P(1);
end

div_mean = mean(pr_slopes);
div_ste =  std(pr_slopes)/sqrt(length(pr_slopes));

patch(patchrange,patchys,[250/255 219/255 206/255],'EdgeColor','none');
hold on;

%make a patch to indicate the error of the protein synthesis in
%dividing cells
patchys_dup = horzcat([std(div_protein{:,2:end},0,2)/sqrt(5) + mean(div_protein{:,2:end},2)]',flip([mean(div_protein{:,2:end},2) - std(div_protein{:,2:end},0,2)/sqrt(5)]'))
patch(patchrange,patchys_dup,[222/255 235/255 250/255],'EdgeColor','none');

%plot the averages
h(2) = plot(range,mean(nondiv_protein{:,2:end},2),'LineWidth',1,'Color',[217/255 83/255 25/255]);
h(1) = plot(range,mean(div_protein{:,2:end},2),'LineWidth',1,'Color',[0/255 114/255 189/255]);

xlabel('Time from start of feeding (min)');
ylabel('Normalized GFP activity');
legend(h,'Dividing cells','Non-dividing cells');

figure;

bar([div_mean  nondiv_mean]);
hold on;
errorbar([div_mean  nondiv_mean],[div_ste nondiv_ste],'.');

%make the cell size figures and calculate the extension rate

%first the dividing cells
figure;
hold on;
%drange = length(Dmean_obj(:,1));

slopes = zeros(1,length(div_length{1,2:end}));

drange = div_length.time_min_;
ofinterest =  find(drange < nondividetime);

for i=1:length(div_length{1,2:end})
    plot(drange,div_length{:,i+1},'Color',[0.7 0.7 0.7])
    P = polyfit(drange(ofinterest), div_length{ofinterest,i+1}, 1); %perform a linear regression on the time before division
    slopes(i) =  P(1);
end

avgslope = mean(slopes);
stderr =  std(slopes)/sqrt(length(div_length{1,2:end}));

plot(drange,mean(div_length{:,2:end},2),'LineWidth',2,'Color',[0/255 114/255 189/255])

ylim([1.5 4.5]);
xlabel('Time from start of feeding (min)');
ylabel('Cell length (\mum)');
title('Dividing cells');

text(5,4.2,{'Average extension rate before 90 min: ',[num2str(avgslope,2) ' \mum/min \pm ' num2str(stderr,2) ' \mum/min']});


%then the non-dividing cells
figure;
hold on;
for i=1:length(nondiv_length{1,2:end})
    plot(drange,nondiv_length{:,i+1},'Color',[0.7 0.7 0.7])
    P = polyfit(drange(ofinterest), nondiv_length{ofinterest,i+1}, 1); %perform a linear regression on the time before division
    slopes(i) =  P(1);
end

avgslope = mean(slopes);
stderr =  std(slopes)/sqrt(length(div_length{1,2:end}));

plot(drange,mean(nondiv_length{:,2:end},2),'LineWidth',2,'Color',[217/255 83/255 25/255])
ylim([1.5 4.5]);
xlabel('Time from start of feeding (min)');
ylabel('Cell length (\mum)');

title('Non-dividing cells');

text(5,4.2,{'Average extension rate before 90 min: ',[num2str(avgslope,2) ' \mum/min \pm ' num2str(stderr,2) ' \mum/min']});



