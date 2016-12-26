function plot_Swiss_Roll(matchingMethod)

if nargin<1
    matchingMethod=2;
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
        case 2
            fileName=strcat(rootDir,'Data/Results/SwissRoll',num2str(matchingMethod),'NoiseResults.mat');
            xla='Noise Level';
            figureName='SwissRollNoise';
        case 3
            fileName=strcat(rootDir,'Data/Results/SwissRoll',num2str(matchingMethod),'OutlierResults.mat');
            xla='Percentage of Outliers';
            figureName='SwissRollOutlier';
    end
    
    load(fileName);
    x=1:length(numRange);
    
    ax=figure;
    plot(x,accMMSJ,'.-', x,accMDS,'.:', x,accIso,'.--',x,accLLE,'.--',x,accLTSA,'.--', 'LineWidth',2);
    legend('MMSJ', 'MDS','Isomap', 'LLE', 'LTSA','Location','SouthEast');
    xlabel(xla)
    ylabel('Matching Ratio')
    ylim([0 1])
    saveas(ax,figureName,'pdf')
end