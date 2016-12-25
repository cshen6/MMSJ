function run_Swiss_Roll_Noise()
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

% dis=[X1P Y1P];%nonlinear vs linear
Y_data=[Y_data' zeros(1, N)']';
%%linear original vs LLE embedded
%dis=[Y1P Y3P];
%disEuc=[Y_data Y_data_LLE];
%%end
% load Sim_Swiss.mat
% load Sim_SwissBroken.mat
% load Sim_Twinpeaks.mat
% load Sim_Inter.mat
%Z_data=[cos(Y_data(1,:)).*sin(Y_data(2,:));sin(Y_data(1,:)).*sin(Y_data(2,:));cos(Y_data(1,:))];
disEuc=[X_data Y_data];
%disEuc=[X_data Z_data];
dis=squareform(pdist(disEuc'));
ss=size(dis,2)/2;
dis=[dis(1:ss, 1:ss) dis(ss+1:end, ss+1:end)];
%pre-process end
tran=1000;numData=2;tesn=100;K=10;iter=-1;reps=100;scale=1;dim=2;
options = struct('nonlinear',0,'match',2,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[sol,power,dCorr]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Original Distance')
options = struct('nonlinear',1,'match',2,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solIsoJ,powerIsoJ,dCorrIsoJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Joint Isomap')
options = struct('nonlinear',1,'match',2,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solIsoS,powerIsoS, dCorrIsoS]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Separate Isomap')
options = struct('nonlinear',2,'match',2,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLLEJ,powerLLEJ, dCorrLLEJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Joint LLE')
options = struct('nonlinear',2,'match',2,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLLES,powerLLES, dCorrLLES]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test Using Separate LLE')
options = struct('nonlinear',3,'match',2,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLTSA,powerLTSA]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
%title('Matching Test Using Joint LLE')
options = struct('nonlinear',4,'match',2,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solLAP,powerLAP]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
options = struct('nonlinear',5,'match',2,'neighborSize',K,'jointSelection',0,'CIsomap',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solHLLE,powerHLLE]=MatchingTestEuc(disEuc,dim,tran,tesn,reps,options);
%Seeded trial
[accIsoJ,acc,accIsoS,accLLES,accLTSA]=ResultsAcc(solIsoJ,sol,solIsoS,solLLES,solLTSA,numData,tran,tesn,reps);

save(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'dim',num2str(dim),'noise',num2str(noise),'.mat'),'tA','test','testN','type','n','dim','noise','rep','powerMLocal','neighbor','pMLocal','pMGC','optimalInd','k','l','C','D','x','y','A','B','mantelH','mcorrH','A_MGC','B_MGC','C_MGC','RC','RD');