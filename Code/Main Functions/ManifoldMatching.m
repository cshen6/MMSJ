function [sol, dCorr]=ManifoldMatching(dis,dim,options)
% This is the main manifold matching function used to match and embed data
% by nonlinear embedding algorithm followed by proper matching.
% It supports joint/separate Isomap/LLE algorithm with joint
% MDS/Procrustes/CCA matching, for any number of given data sets.
% Note that this is the main version used in our paper, and a
% Euclidean-based version is called ManifoldMatchingEuc.m, which is used
% for benchmark purpose in our paper.
%
% For the input, dis should be a distance matrix of size n*(numData*n),
% where numData is the number of data that are matched.
% dim is the embedding and matching dimension.
% Information such as which nonlinear algorithm to use, which matching
% technique to apply, neighborhood size, out-of-sample size,
% joint selection or not, etc., are specified in options (see the code part).
%
% For the output, sol contains the matched data of size dim*(numData*n).
% dCorr returns the distance correlation for in-sample matched data, for
% the original distance, or the shortest-path distance, or LLE embedded
% distance, based on the specified nonlinear algorithm.

%    Notice & Acknowledgement
%
%    This file is part of the Manifold Matching Code, for the following
%    paper:
%
%    Cencheng Shen and Carey E. Priebe, "Manifold Matching using Shortest-Path Distance and Joint Neighborhood Selection",
%    submitted, on arxiv, http://arxiv.org/abs/1412.4098, 2015
%
%    Feel free to modify, extend or distribute this code for non-commercial purposes,
%    as long as proper acknowledgement and citation are made. Please email
%    to cshen6@jhu.edu (prefered email) or cep@jhu.edu for comments and bug reports.
%
%    (C) Cencheng Shen, 2015
%    Johns Hopkins University, Baltimore, Maryland, United States

%Specify options
if nargin < 2
    error('Too few input arguments');
elseif nargin < 3
    options = struct('nonlinear',0,'match',1,'neighborSize',max(10,dim+5),'jointSelection',1,'weight',1,'scaling',1,'numData',2,'maxIter',-1,'oos',0);
end

if ~isfield(options,'weight')
    options.weight = 1; %The weight to be used for raw-stress MDS. By default CMDS is used and weight is not needed.
end
if ~isfield(options,'nonlinear')
    options.nonlinear = 0; %Specify embedding algorithm. 0 for MDS/PCA, 1 for Isomap, 2 for LLE.
end
if ~isfield(options,'match')
    options.match = 1; %Specify matching technique. 0 for joint MDS, 1 for Procrustes, 2 for CCA.
end
if ~isfield(options,'jointSelection')
    options.jointSelection = 1; %Specify joint selection of neighborhood or not in nonlinear embedding. 0 for separate neighbor, otherwise do joint.
end
if ~isfield(options,'neighborSize')
    options.neighborSize = max(10,dim+5); %Set the neighbordhood size for nonlinear embedding.
end
if ~isfield(options,'scaling')
    options.scaling = 1; %Specify scaling or not. 0 for no-scaling, otherwise do scaling based on the Frobenius norm of the first data set.
end
if ~isfield(options,'numData')
    options.numData = 2; %Specify the number of dataset to match.
end
if ~isfield(options,'maxIter')
    options.maxIter = -1; %Specify the iteration maximum. It is used in raw-stress MDS and also Generalized CCA matching. In the default option, CMDS is used in case of MDS; and 30 iterations are used for GCCA.
end
if ~isfield(options,'oos')
    options.oos = 0; %Specify the number of points embedded by out-of-sample (OOS) technique. Note that the last points in each data set will be OOS embedded.
end

%Initialization
match=options.match;
numData=options.numData;
n=size(dis,2)/numData;
tesn=options.oos;
tran=n-tesn;
if options.oos ~= 0 %check if there is any OOS embedding
    %fprintf('The number of OOS-embedded points is %d',options.oos);
    if options.oos == n
        error('The number of OOS points should be smaller than the sample size');
    end
end
solIn=zeros(dim,numData*tran);
sol=zeros(dim,numData*n);
disO=zeros(numData*n,numData*n);
dCorr=0;

%Pre-scale the data according to the Frobenius norm of the first data if needed.
if options.scaling~=0
    norm1=norm(dis(:,1:n),'fro');
    for i=1:numData;
        tmp1=n*(i-1)+1:n*i;
        dis(:,tmp1)=dis(:,tmp1)/norm(dis(:,tmp1),'fro')*norm1;
    end
end

%Nonlinear transformation into shortest-path distance or LLE embedded distance
if options.nonlinear==1
    disIso=NeighborhoodIsomap(dis,numData,options);
    disTmp=zeros(n,numData*n);
    for j=1:numData
        tmp1=n*(j-1)+1:n*j;
        disTmp(:,tmp1)=disIso(tmp1,tmp1);
    end
    dis=disTmp;
end
if options.nonlinear==2
    eucli=NeighborhoodLLE(dis,numData, dim, options);
    disTmp=zeros(n,numData*n);
    for j=1:numData
        tmp1=n*(j-1)+1:n*j;
        disTmp(:,tmp1)=squareform(pdist(eucli(:,tmp1)'));
    end
    dis=disTmp;
end

%Post-scale the distance if needed; this step is not necessary if
%CCA matching is used, but may be useful for other matching with Isomap.
if options.scaling~=0
    for j=1:numData;
        tmp1=n*(j-1)+1:n*j;
        if j==1
            normScale=norm(dis(:,tmp1),'fro');
        else
            dis(:,tmp1)=dis(:,tmp1)/norm(dis(:,tmp1),'fro')*normScale;
        end
    end
end

%Concatenate the Omnibus matrix disO
for i=1:numData;
    tmp1=n*(i-1)+1:n*i;
    disO(tmp1,tmp1)=dis(:,tmp1);
    for j=i+1:numData
        dCorr=dCorr+distCorr(dis(1:tran,n*(i-1)+1:n*(i-1)+tran), dis(1:tran,n*(j-1)+1:n*(j-1)+tran)); %update distance correlation for the in-sample data
    end
end

%In case of one dataset, do the usual MDS/PCA embedding without matching
%and return.
if numData==1;
    sol=SMDS(disO(1:n,1:n),dim,0,options);
    if tesn~=0
        sol=OOSMDS(disO(1:tran+tesn,1:tran+tesn),tran,dim,sol); %OOS the last tesn points if necessary
    end
    return;
end

%For more than one dataset, do Omnibus matching by applying weighted MDS
% or CMDS on the imputed Omnibus Matrix.
if match==0
    weight=options.weight*ones(numData*tran,numData*tran);
    for i=1:numData
        tmp1=tran*(i-1)+1:tran*i;
        for j=i+1:numData
            tmp2=tran*(j-1)+1:tran*j;
            weight(tmp1,tmp2)=eye(tran);
            weight(tmp2,tmp1)=weight(tmp1,tmp2)';
        end
    end
    %Impute the off-diagonal entries of the Omnibus matrix.
    disOIn=zeros(tran*numData,tran*numData); %In-sample omnibus matrix
    for i=1:numData
        tmp1=n*(i-1)+1:n*(i-1)+tran;
        tmp2=tran*(i-1)+1:tran*i;
        disOIn(tmp2,tmp2)=disO(tmp1,tmp1);
    end
    for i=1:numData
        tmp1=tran*(i-1)+1:tran*i;
        for j=i+1:numData
            tmp2=tran*(j-1)+1:tran*j;
            disOIn(tmp1,tmp2)=(disOIn(tmp1,tmp1)+disOIn(tmp2,tmp2))/2;
            disOIn(tmp2,tmp1)=disOIn(tmp1,tmp2)';
        end
    end
    solIn=SMDS(disOIn,dim,weight,options);
    %The equivalent but slower Matlab default MDS is:
    %solIn=mdscale(disOIn,dim,'Criterion','metricstress','Weights',weight,'Start','random')';
    if tesn~=0 %Do OOS for each data set if necessary
        for i=1:numData
            solOOS=OOSMDS(disO(n*(i-1)+1:n*i,n*(i-1)+1:n*i),tran, dim,solIn(:,tran*(i-1)+1:tran*i));
            sol(:, n*(i-1)+1:n*i)=solOOS;
        end;
    else
        sol=solIn;
    end
end

%Do Procrustes matching
if match==1
    for i=1:numData
        tmp1=tran*(i-1)+1:tran*i;
        tmp2=n*(i-1)+1:n*(i-1)+tran;
        solIn(:,tmp1)=SMDS(disO(tmp2,tmp2),dim,0,options);
    end
    [solInP,T]=GeneralProcrustes(solIn,numData); %Generalized Procrustes matching using greedy approximation
    if tesn~=0
        for i=1:numData
            solOOS=OOSMDS(disO(n*(i-1)+1:n*i,n*(i-1)+1:n*i),tran, dim,solIn(:,tran*(i-1)+1:tran*i));
            sol(:, n*(i-1)+1:n*(i-1)+tran) = solInP(:, tran*(i-1)+1:tran*i);
            sol(:, n*(i-1)+tran+1:n*i)=T(:,:,i)*solOOS(:,tran+1:n);
        end
    else
        sol=solInP;
    end
end

%Do CCA matching
if match==2
    for i=1:numData
        tmp1=tran*(i-1)+1:tran*i;
        tmp2=n*(i-1)+1:n*(i-1)+tran;
        solIn(:,tmp1)=SMDS(disO(tmp2,tmp2),dim,0,options);
    end
    [solInP,T]=GeneralCCA(solIn,numData,dim,options); %Generalized CCA matching
    if tesn~=0
        for i=1:numData
            solOOS=OOSMDS(disO(n*(i-1)+1:n*i,n*(i-1)+1:n*i),tran, dim,solIn(:,tran*(i-1)+1:tran*i));
            sol(:, n*(i-1)+1:n*(i-1)+tran) = solInP(:, tran*(i-1)+1:tran*i);
            sol(:, n*(i-1)+tran+1:n*i)=T(:,dim*(i-1)+1:dim*i)'*solOOS(:,tran+1:n);
        end
    else
        sol=solInP;
    end
end

%The auxiliary function to calculate the sample distance correlation.
function corr = distCorr(X,Y)
n=size(X,1);
H=eye(n)-(1/n)*ones(n,n);
X=H*X*H;
Y=H*Y*H;

corr=sum(sum(X.*(Y)));
corr=sqrt(corr/(norm(X,'fro')*norm(Y,'fro')));