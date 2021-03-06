function run_Wiki_Data(matchingMethod,reps)
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


load('Wikipedia.mat');

tran=500;dim=10;tesn=100;K=20;numData=2;reps=100;iter=-1;scale=1;m=50;matchingMethod=2;
dis=[TE TF];
ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
disEuc=[TEEuc TFEuc];

% load('BrainCP.mat')
% dis=[distC, distP];
% tran=36;dim=3;tesn=2;K=6;numData=2;reps=100;iter=-1;scale=1;m=5;matchingMethod=1;n=42;
% per=zeros(reps,n);

load('Semipar.mat')
% dis=[distM2g(ind,ind), squareform(pdist(cci(ind)))];
dis=[distM2g(ind,ind), distMigrain(ind,ind)];
tran=79;dim=2;tesn=10;K=7;numData=2;reps=100;iter=-1;scale=1;m=dim;matchingMethod=1;n=109;
per=zeros(reps,n);

ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
disEuc=[TEEuc TFEuc];
for i=1:reps;
per(i,:)=randperm(n);
end
options = struct('nonlinear',0,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accMDS,powerMDS]=MatchingTest(dis,dim,tran,tesn,reps,options);
options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accMMSJ,powerMMSJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accISO,powerISO]=MatchingTest(dis,dim,tran,tesn,reps,options);
options = struct('nonlinear',2,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accLLE,powerLLE]=MatchingTest(dis,dim,tran,tesn,reps,options);
options = struct('nonlinear',3,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
[~,accLTSA,powerLTSA]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);


% ss=size(dis,2)/numData;
% TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
% TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
% disEuc=[TEEuc TFEuc];
% for i=1:reps;
%     per(i,:)=randperm(n);
% end


% 
% load('BrainHippoShape.mat')
% dis=[LMLS, LMRS];
% tran=99;dim=3;tesn=5;K=10;numData=2;reps=100;iter=-1;scale=1;m=10;matchingMethod=1;
% per=zeros(reps,n);
% 
% ss=size(dis,2)/numData;
% TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
% TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
% disEuc=[TEEuc TFEuc];
% for i=1:reps;
%     per(i,:)=randperm(n);
% end
% % y=squareform(pdist(Label));
% % %y=(y>0)+1;
% % y=y+1;
% % for i=1:n
% %     y(i,i)=0;
% % end
% % dis=[LMLS,y];
% %
% 
% %
% %
% % political graphs
% load('PoliticalNetwork.mat')
% tran=22;dim=5;tesn=4;K=6;numData=2;reps=100;iter=-1;scale=1;m=dim;matchingMethod=1;n=34;
% dis=[squareform(pdist(Adj')) squareform(pdist(Label))];
% per=zeros(reps,n);
% 
% ss=size(dis,2)/numData;
% TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
% TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
% disEuc=[TEEuc TFEuc];
% for i=1:reps;
%     per(i,:)=randperm(n);
% end
% TEEuc=SMDS(Adj,m,0);
% TFEuc=[Label,zeros(n,1)]';
% disEuc=[TEEuc TFEuc];
% dis=[squareform(pdist(TEEuc')) squareform(pdist(TFEuc'))];


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
