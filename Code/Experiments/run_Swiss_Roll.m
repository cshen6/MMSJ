function run_Swiss_Roll()
% generate data for figure 1
%%% File path searching
if nargin<1
    type=8;
end
if nargin<2
    n=50;
end
if nargin<3
    dim=1;
end
if nargin<4
    noise=0;
end
if nargin<5
    rep=1000;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
addpath(genpath(strcat(rootDir,'Data/')));
load('SwissRoll.mat');

save(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'dim',num2str(dim),'noise',num2str(noise),'.mat'),'tA','test','testN','type','n','dim','noise','rep','powerMLocal','neighbor','pMLocal','pMGC','optimalInd','k','l','C','D','x','y','A','B','mantelH','mcorrH','A_MGC','B_MGC','C_MGC','RC','RD');