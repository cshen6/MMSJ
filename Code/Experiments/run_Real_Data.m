function run_Real_Data(matchingMethod,reps)
if nargin<1
    matchingMethod=2;
end
if nargin<2
    reps=100;
end

% generate data for figure 1
%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
addpath(genpath(strcat(rootDir,'Data/')));

% 
% load('Wikipedia.mat');
% 
% tran=500;dim=10;tesn=100;K=20;numData=2;reps=100;iter=-1;scale=1;m=50;matchingMethod=2;
% dis=[TE TF];
% ss=size(dis,2)/numData;
% TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
% TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
% disEuc=[TEEuc TFEuc];
% options = struct('nonlinear',0,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
% [~,accMDS,powerMDS]=MatchingTest(dis,dim,tran,tesn,reps,options);
% options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
% [~,accMMSJ,powerMMSJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
% options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
% [~,accISO,powerISO]=MatchingTest(dis,dim,tran,tesn,reps,options);
% options = struct('nonlinear',2,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
% [~,accLLE,powerLLE]=MatchingTest(dis,dim,tran,tesn,reps,options);
% options = struct('nonlinear',3,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
% [~,accLTSA,powerLTSA]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
% load('BrainCP.mat')
% dis=[distC, distP];
% tran=36;dim=3;tesn=2;K=6;numData=2;reps=100;iter=-1;scale=1;m=5;matchingMethod=1;n=42;
% per=zeros(reps,n);


load('Brain.mat')
% dis=[distM2g(ind,ind), squareform(pdist(cci(ind)))];
dim=2;K=7;numData=2;reps=100;iter=-1;scale=1;m=dim;matchingMethod=1;n=109;numRange=1:20;

accMDS=zeros(length(numRange),1);powerMDS=zeros(length(numRange),1);
accMMSJ=zeros(length(numRange),1);powerMMSJ=zeros(length(numRange),1);
accISO=zeros(length(numRange),1);powerISO=zeros(length(numRange),1);
accLLE=zeros(length(numRange),1);powerLLE=zeros(length(numRange),1);
accLTSA=zeros(length(numRange),1);powerLTSA=zeros(length(numRange),1);

dis=[distM2g, distMigrain];
% per=zeros(reps,n);
ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
disEuc=[TEEuc TFEuc];
% for i=1:reps;
% per(i,:)=randperm(n);
% end

for i=1:length(numRange)
tran=109-(21-i)*3;tesn=21-i;
options = struct('nonlinear',0,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accMDS(i),tmp]=MatchingTest(dis,dim,tran,tesn,reps,options);
powerMDS(i)=tmp(2);
options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accMMSJ(i),tmp]=MatchingTest(dis,dim,tran,tesn,reps,options);
powerMMSJ(i)=tmp(2);
options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accISO(i),tmp]=MatchingTest(dis,dim,tran,tesn,reps,options);
powerISO(i)=tmp(2);
options = struct('nonlinear',2,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accLLE(i),tmp]=MatchingTest(dis,dim,tran,tesn,reps,options);
powerLLE(i)=tmp(2);
options = struct('nonlinear',3,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accLTSA(i),tmp]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
powerLTSA(i)=tmp(2);
end

save(strcat(rootDir,'Data/Results/Brain',num2str(matchingMethod),'Results.mat'),'numRange','K','matchingMethod','numData','iter','dim','reps','accMDS','powerMDS','accMMSJ','powerMMSJ','accISO','powerISO','accLLE','powerLLE','accLTSA','powerLTSA');
    
% %three data
% clear;
% load Wiki_Data.mat
% dis=[TE TF GE];numData=3;
% ss=size(dis,2)/numData;
% TEEuc=SMDS(dis(1:ss, 1:ss),m,0);TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);GEEuc=SMDS(dis(1:ss, 2*ss+1:3*ss),m,0);
% disEuc=[TEEuc TFEuc GEEuc];
% %title('Wikipedia TE and TF and GE Matching Test Using Original Distance')
% %four data
% clear;
% load Wiki_Data.mat
% dis=[TE TF GE GF];numData=4;
% ss=size(dis,2)/numData;
% TEEuc=SMDS(dis(1:ss, 1:ss),m,0);TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);GEEuc=SMDS(dis(1:ss, 2*ss+1:3*ss),m,0);GFEuc=SMDS(dis(1:ss, 3*ss+1:4*ss),m,0);
% disEuc=[TEEuc TFEuc GEEuc GFEuc];
