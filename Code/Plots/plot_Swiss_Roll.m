function plot_Swiss_Roll(matchingMethod)

if nargin<1
    matchingMethod=1;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Data/')));

map1=zeros(5,3);
gr = [0,1,0];
ma = [1,0,1];
cy = [0,1,1];
map1(1,:) = gr;
map1(2,:) = ma;
map1(3,:) = cy;
map1(4,:) = [1,0,0];
map1(5,:) = [0 0 0];
%cmap(6,:) = [0.5 0.5 0];
set(groot,'defaultAxesColorOrder',map1);

for i=1:3
    switch i
        case 1
            fileName=strcat(rootDir,'Data/Results/SwissRoll',num2str(matchingMethod),'Results.mat');
            xla='Number of Training Data';
            figureName='SwissRoll';
            tit='Swiss Roll';
            lgdLoc='NorthWest';
        case 2
            fileName=strcat(rootDir,'Data/Results/SwissRoll',num2str(matchingMethod),'ResultsNoise.mat');
            xla='Noise Level';
            figureName='SwissRollNoise';
            tit='Swiss Roll with Noise';
            lgdLoc='NorthEast';
        case 3
            fileName=strcat(rootDir,'Data/Results/SwissRoll',num2str(matchingMethod),'ResultsOutlier.mat');
            xla='Percentage of Outliers';
            figureName='SwissRollOutlier';
            tit='Swiss Roll with Outliers';
            lgdLoc='NorthEast';
    end
    
    load(fileName);
    x=1:length(numRange);
    
    ax=figure;
    plot(x,accMMSJ,'.-', x,accMDS,'.:', x,accISO,'.--',x,accLLE,'.--',x,accLTSA,'.--', 'LineWidth',4);
    legend('MMSJ', 'MDS','Isomap', 'LLE', 'LTSA','Location',lgdLoc);
    title(tit,'FontSize',32);
    xlabel(xla)
    ylabel('Matching Ratio')
    ylim([0 1])
    xlim([1 length(numRange)])
    set(gca,'FontSize',25);
    set(gca,'XTick',[1,ceil(length(numRange)/2),length(numRange)],'XTickLabel',[numRange(1),numRange(ceil(length(numRange)/2)),numRange(end)],'YTick',[0.5,1]); 
    saveas(ax,fullfile(strcat(rootDir,'Figures/'), strcat(figureName,'Acc',num2str(matchingMethod))),'png')
    
    ax=figure;
    plot(x,powerMMSJ(:,2),'.-', x,powerMDS(:,2),'.:', x,powerISO(:,2),'.--',x,powerLLE(:,2),'.--',x,powerLTSA(:,2),'.--', 'LineWidth',4);
    legend('MMSJ', 'MDS','Isomap', 'LLE', 'LTSA','Location',lgdLoc);
    title(tit,'FontSize',32);
    xlabel(xla)
    ylabel('Testing Power')
    ylim([0 1])
    xlim([1 length(numRange)])
    set(gca,'FontSize',25);
    set(gca,'XTick',[1,ceil(length(numRange)/2),length(numRange)],'XTickLabel',[numRange(1),numRange(ceil(length(numRange)/2)),numRange(end)],'YTick',[0.5,1]); 
    saveas(ax,fullfile(strcat(rootDir,'Figures/'), strcat(figureName,'Power',num2str(matchingMethod))),'png')
end