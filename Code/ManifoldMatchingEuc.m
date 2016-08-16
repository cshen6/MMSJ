function sol=ManifoldMatchingEuc(dis,dim,options)
% This is the main manifold matching function used to match and embed data
% by nonlinear embedding algorithm followed by proper matching.
% It supports PCA/Isomap/LLE/LTSA/Laplacian eigenmaps/Hessian LLE algorithm
% with joint MDS/Procrustes/CCA matching, for any number of given data sets.
%
% Note that this version is for Euclidean data with separate neighborhood only,
% and the main distance version is in ManifoldMatching.m. Also note that
% all nonlinear algorithms use the Matlab dimension reduction toolbox by L.
% Maaten from http://lvdmaaten.github.io/drtoolbox/, and their
% implementations may be different from us in case of Isomap/LLE/MDS in the distance version.
%
% For the input, dis should be a data of size d*(numData*n), where
% numData is the number of data that are matched, d is the ambient dimension.
% dim is the embedding and matching dimension.
% Information such as which nonlinear algorithm to use, which matching
% technique to apply, neighborhood size, out-of-sample size,
% etc., are specified in options (see the code part).
%
% For the output, sol contains the matched data of size dim*(numData*n).

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
    options = struct('nonlinear',0,'match',1,'neighborSize',max(10,dim+5),'jointSelection',0,'weight',1,'scaling',1,'numData',2,'maxIter',-1,'oos',0);
end

if ~isfield(options,'weight')
    options.weight = 1; %The weight to be used for raw-stress MDS. By default CMDS is used and weight is not needed.
end
if ~isfield(options,'nonlinear')
    options.nonlinear = 0; %Specify embedding algorithm. 0 for MDS/PCA, 1 for Isomap, 2 for LLE, 3 for Laplacian eigenmaps, 4 for LTSA, 5 for Hessian LLE.
end
if ~isfield(options,'match')
    options.match = 1; %Specify matching technique. 0 for MDS/PCA, 1 for Isomap, 2 for LLE.
end
if ~isfield(options,'jointSelection')
    options.jointSelection = 0; %Specify joint selection of neighborhood or not in nonlinear embedding. 0 for separate neighbor, otherwise do joint.
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
disO=zeros(dim,numData*n);

%Pre-scale the data according to the Frobenius norm of the first data if needed
if options.scaling~=0
    norm1=norm(dis(:,1:n),'fro');
    for i=1:numData;
        tmp1=n*(i-1)+1:n*i;
        dis(:,tmp1)=dis(:,tmp1)/norm(dis(:,tmp1),'fro')*norm1;
    end
end

%Nonlinear transformation
nonlinearMode='Laplacian';
switch (options.nonlinear)
    case 0
        nonlinearMode='PCA';
    case 1
        nonlinearMode='Isomap';
    case 2
        nonlinearMode='LLE';
    case 3
        nonlinearMode='LTSA';
    case 4
        nonlinearMode='Laplacian';
    case 5
        nonlinearMode='HLLE';
end

%Embed the data
for j=1:numData
    eucliX=dis(:, n*(j-1)+1:n*j);
    [solX, mapping]=compute_mapping(eucliX(:,1:tran)', nonlinearMode, dim,options.neighborSize); %From Matlab dimension reduction toolbox
    if options.oos~=0
        try
            solXOOS=out_of_sample(eucliX(:,tran+1:n)', mapping); %From Matlab dimension reduction toolbox
        catch
            solXOOS=out_of_sample_est(eucliX(:,tran+1:n)', eucliX(:,1:tran)', solX); %From Matlab dimension reduction toolbox
        end
        solX=[solX' solXOOS']';
    end
    disO(:,n*(j-1)+1:n*j)=solX';
end

%Post-scale the data if needed; note that this step is actually redundant
%for LTSA, Laplacian eigenmaps, HLLE, etc.
if options.scaling~=0
    for j=1:numData;
        tmp1=n*(j-1)+1:n*j;
        if j==1
            normScale=norm(disO(:,tmp1),'fro');
        else
            disO(:,tmp1)=disO(:,tmp1)/norm(disO(:,tmp1),'fro')*normScale;
        end
    end
end

%In case of one dataset, do the usual MDS/PCA embedding without matching.
%and return.
if numData==1;
    sol=compute_mapping(disO', 'PCA', dim)';
    return;
end

%For more than one dataset, do Omnibus matching by applying weighted MDS
% or CMDS on the imputed Omnibus Matrix.
if match==0
    disO=squareform(pdist(disO'));
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
    %solIn=mdscale(disR,dim,'Criterion','metricstress','Weights',weight,'Start','random')';
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
        solIn(:,tran*(i-1)+1:tran*i)=disO(:,n*(i-1)+1:n*(i-1)+tran);
    end
    [solInP,T]=GeneralProcrustes(solIn,numData); %Generalized Procrustes matching using greedy approximation
    if tesn~=0
        for i=1:numData
            solOOS=disO(:,n*(i-1)+tran+1:n*i);
            sol(:, n*(i-1)+1:n*(i-1)+tran) = solInP(:, tran*(i-1)+1:tran*i);
            sol(:, n*(i-1)+tran+1:n*i)=T(:,:,i)*solOOS;
        end
    else
        sol=solInP;
    end
end

%Do CCA matching
if match==2
    for i=1:numData
        solIn(:,tran*(i-1)+1:tran*i)=disO(:,n*(i-1)+1:n*(i-1)+tran);
    end
    [solInP,T]=GeneralCCA(solIn,numData,dim,options); %Generalized CCA matching
    if tesn~=0
        for i=1:numData
            solOOS=disO(:,n*(i-1)+tran+1:n*i);
            sol(:, n*(i-1)+1:n*(i-1)+tran) = solInP(:, tran*(i-1)+1:tran*i);
            sol(:, n*(i-1)+tran+1:n*i)=T(:,dim*(i-1)+1:dim*i)'*solOOS;
        end
    else
        sol=solInP;
    end
end