%look at size and protein synthesis of non-dividing and dividing cells from
%microfluidic measurements
% by Karthik Sekar (karsekar@gmail.com) with preprocessing done by Roberto Rusconi
% last updated 24.10.2017

clear all;
close all;

%load all of the data
div_protein = readtable('protein_and_size_data/div_protein.csv');
nondiv_protein = readtable('protein_and_size_data/nondiv_protein.csv');

div_length = readtable('protein_and_size_data/div_length.csv');
nondiv_length = readtable('protein_and_size_data/nondiv_length.csv');

%indicate the time range
range = div_protein.time_min_';

%first make the protein synthesis figure

%make a patch to indicate the error of the protein synthesis in
%non-dividing cells
patchrange = [range flip(range)];
patchys = horzcat([std(nondiv_protein{:,2:end},0,2)/sqrt(5) + mean(nondiv_protein{:,2:end},2)]',flip([mean(nondiv_protein{:,2:end},2) - std(nondiv_protein{:,2:end},0,2)/sqrt(5)]'));
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

%make the cell size figures;

%first the dividing cells
figure;
hold on;
%drange = length(Dmean_obj(:,1));
drange = div_length.time_min_;
for i=1:length(div_length{1,2:end})
    plot(drange,div_length{:,i+1},'Color',[0.7 0.7 0.7])
end

plot(drange,mean(div_length{:,2:end},2),'LineWidth',2,'Color',[0/255 114/255 189/255])

ylim([1.5 4.5]);
xlabel('Time from start of feeding (min)');
ylabel('Cell length (um)');
title('Dividing cells');


%then the non-dividing cells
figure;
hold on;
for i=1:length(nondiv_length{1,2:end})
    plot(drange,nondiv_length{:,i+1},'Color',[0.7 0.7 0.7])
end

plot(drange,mean(nondiv_length{:,2:end},2),'LineWidth',2,'Color',[217/255 83/255 25/255])
ylim([1.5 4.5]);
xlabel('Time from start of feeding (min)');
ylabel('Cell length (um)');

title('Non-dividing cells');


