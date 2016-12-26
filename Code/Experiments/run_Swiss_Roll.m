function run_Swiss_Roll(matchingMethod,tran,reps)
if nargin<1
    matchingMethod=2;
end
if nargin<2
    tran=1000;
end
if nargin<3
    reps=100;
end

fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
addpath(genpath(strcat(rootDir,'Data/')));
load('SwissRoll.mat');

numData=2;tesn=100;K=10;iter=-1;scale=1;dim=2;
Y_data=[Y_data' zeros(1, N)']';
disEuc=[X_data Y_data];
dis=squareform(pdist(disEuc'));
ss=size(dis,2)/2;
dis=[dis(1:ss, 1:ss) dis(ss+1:end, ss+1:end)];
numRange=50:50:tran;
accMDS=zeros(length(numRange),1);powerMDS=zeros(length(numRange),1);
accMMSJ=zeros(length(numRange),1);powerMMSJ=zeros(length(numRange),1);
accISO=zeros(length(numRange),1);powerISO=zeros(length(numRange),1);
accLLE=zeros(length(numRange),1);powerLLE=zeros(length(numRange),1);
accLTSA=zeros(length(numRange),1);powerLTSA=zeros(length(numRange),1);
for i=1:length(numRange)
    tran=numRange(i);
    options = struct('nonlinear',0,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
    [~,accMDS(i),powerMDS(i)]=MatchingTest(dis,dim,tran,tesn,reps,options);
    options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
    [~,accMMSJ(i),powerMMSJ(i)]=MatchingTest(dis,dim,tran,tesn,reps,options);
    options = struct('nonlinear',1,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
    [~,accISO(i),powerISO(i)]=MatchingTest(dis,dim,tran,tesn,reps,options);
    options = struct('nonlinear',2,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
    [~,accLLE(i),powerLLE(i)]=MatchingTest(dis,dim,tran,tesn,reps,options);
    options = struct('nonlinear',3,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
    [~,accLTSA(i),powerLTSA(i)]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
%     options = struct('nonlinear',4,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
%     [solLAP,accLAP,powerLAP]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
%     options = struct('nonlinear',5,'match',matchingMethod,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn);
%     [solHLLE,accHLLE,powerHLLE]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
end
%Seeded trial

save(strcat(rootDir,'Data/Results/SwissRoll',num2str(matchingMethod),'Results.mat'),'numRange','K','matchingMethod','numData','iter','tesn','dim','reps','solMDS','accMDS','powerMDS','solMMSJ','accMMSJ','powerMMSJ','solISO','accISO','powerISO','solLLE','accLLE','powerLLE','solLTSA','accLTSA','powerLTSA');