function run_Wiki_Data()
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
load('Wikipedia.mat');

%Wiki Data
clear;
%load ResultMatchingWikiTETF.mat
load Wiki_Data.mat
tran=500;dim=10;tesn=100;K=20;numData=2;reps=100;iter=-1;scale=1;m=50;
dis=[TE TF];
%clear;
clear
load('BrainHippoShape.mat')
dis=[LMLS, LMRS];
tran=70;dim=3;tesn=10;K=12;numData=2;reps=100;iter=-1;scale=1;m=10;
y=squareform(pdist(Label));
%y=(y>0)+1;
y=y+1;
for i=1:n
    y(i,i)=0;
end
dis=[LMLS,y];
%
clear
load('BrainCP.mat')
dis=[distC, distP];
tran=30;dim=5;tesn=4;K=5;numData=2;reps=100;iter=-1;scale=1;m=5;
%
%
% political graphs
clear
load('PoliticalNetwork.mat')
dis=[squareform(pdist(Adj)), squareform(pdist(Label))];
tran=22;dim=2;tesn=4;K=5;numData=2;reps=100;iter=-1;scale=1;m=dim;mm=2;
%
ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
disEuc=[TEEuc TFEuc];
dis=[squareform(pdist(TEEuc')), squareform(pdist(Label))];
%dis=[squareform(pdist(TEEuc')) squareform(pdist(TFEuc'))];
options = struct('nonlinear',0,'match',mm,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[sol,power, dCorr]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test using Original Distance')
options = struct('nonlinear',1,'match',mm,'neighborSize',K,'jointSelection',1,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solIsoJ,powerIsoJ, dCorrIsoJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test using Joint Isomap')
options = struct('nonlinear',1,'match',mm,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solIsoS,powerIsoS, dCorrIsoS]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test using Separate Isomap')
options = struct('nonlinear',2,'match',mm,'neighborSize',K,'jointSelection',1,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLLEJ,powerLLEJ, dCorrLLEJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Joint LLE')
options = struct('nonlinear',2,'match',mm,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLLES,powerLLES, dCorrLLES]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Separate LLE')
%LTSA and Laplacian
%options = struct('nonlinear',2,'match',2,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
%[solLLE,powerLLE]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
options = struct('nonlinear',3,'match',mm,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLTSA,powerLTSA]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
options = struct('nonlinear',4,'match',mm,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLAP,powerLAP]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
%options = struct('nonlinear',5,'match',2,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
%[solHLLE,powerHLLE]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
[accIsoJ,acc,accIsoS,accLLES,accLTSA]=ResultsAcc(solIsoJ,sol,solIsoS,solLLES,solLTSA,numData,tran,tesn,reps);

%three data
clear;
load Wiki_Data.mat
dis=[TE TF GE];
tran=500;dim=10;tesn=100;K=20;numData=3;reps=100;iter=-1;scale=1;m=50;
ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
GEEuc=SMDS(dis(1:ss, 2*ss+1:3*ss),m,0);
disEuc=[TEEuc TFEuc GEEuc];
%title('Wikipedia TE and TF and GE Matching Test Using Original Distance')
%four data
clear;
load Wiki_Data.mat
dis=[TE TF GE GF];
tran=500;dim=10;tesn=100;K=20;numData=4;reps=100;iter=-1;scale=1;m=50;
ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
GEEuc=SMDS(dis(1:ss, 2*ss+1:3*ss),m,0);
GFEuc=SMDS(dis(1:ss, 3*ss+1:4*ss),m,0);
disEuc=[TEEuc TFEuc GEEuc GFEuc];

save('ResultMatchingWikiTETF20', 'sol','power','solIsoJ','powerIsoJ','solIsoS','powerIsoS','solLLEJ','powerLLEJ','solLLES','powerLLES','solLTSA','powerLTSA','solLAP','powerLAP','dCorr','dCorrIsoJ','dCorrIsoS','dCorrLLEJ','dCorrLLES','tran','m','dim','tesn','K','iter','reps','numData')

save(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'dim',num2str(dim),'noise',num2str(noise),'.mat'),'tA','test','testN','type','n','dim','noise','rep','powerMLocal','neighbor','pMLocal','pMGC','optimalInd','k','l','C','D','x','y','A','B','mantelH','mcorrH','A_MGC','B_MGC','C_MGC','RC','RD');