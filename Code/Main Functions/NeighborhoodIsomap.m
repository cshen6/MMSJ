function D = NeighborhoodIsomap(dis, numData, options)
% This function takes in the distance matrix, and do separate/joint Isomap
% possibly with out-of-sample (OOS) technique and conformal-Isomap (CIso).
%
% For the input, the distance can be of size n*(numData*n) of (numData*n)*(numData*n)
% dim is the embedding dimension,
% and neighborhood size and oos and joint selection can be
% specified in the option.
%
% Output D is the shortest-path distance matrix of size
% (numData*n)*(numData*n).

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

%Initialization and Options
n=size(dis,2)/numData;
if floor(n)~=n
    error('The size of dissimilarity matrix conflcits the number of data specified');
end
if nargin < 2
    error('Too few input arguments');
elseif nargin < 3
    options = struct('jointSelection',1,'CIsomap',0,'neighborSize',max(10,dim+5),'oos',0);
end
if ~isfield(options,'jointSelection')
    options.jointSelection = 1; %0 for separate neighborhood selection, otherwise do joint neighborhood selection
end
if ~isfield(options,'CIsomap')
    options.CIsomap = 0; %0 for usual Isomap, otherwise do CIso
end
if ~isfield(options,'neighborSize')
    options.neighborSize = max(10,dim+5); %default neighorbood size. Note that it should be larger than dim.
end
if ~isfield(options,'oos')
    options.oos = 0;  %The number of oos embedded points. Note that it is used for last points in each dataset.
end
INF = 1000*max(max(dis))*n;
K=options.neighborSize;
tesn=options.oos;
tran=n-tesn;
indNum=zeros(tran,tran,numData); %store the neighborhood graph info for all datasets

%Dis conversion if dis is not yet of size (numData*n)*(numData*n)
D=dis;
if size(dis,1)~=n*numData;
    D=zeros(numData*n,numData*n);
    for j=1:numData;
        tmp1=n*(j-1)+1:n*j;
        D(tmp1,tmp1)=dis(:,tmp1);
    end
end

%Construct neighborhood graph
if options.jointSelection==0;
    %Separate neighborhood selection constructs a graph for each data
    %disp('Separately constructing neighborhood graph');
    for i=1:numData
        tmp1=n*(i-1)+1:n*(i-1)+tran;
        [~,ind]=sort(D(tmp1,tmp1));
        indNum(:,:,i)=ind;
    end
else
    %Joint neighborhood selection constructs the same graph for all data
    %disp('Jointly constructing neighborhood graph');
    tmpSum=0;
    for i=1:numData
        tmp1=n*(i-1)+1:n*(i-1)+tran;
        tmpSum=tmpSum+D(tmp1,tmp1);
    end
    [~,ind]=sort(tmpSum);
    indNum=repmat(ind,[1,1,numData]);
end

%Compute the shortest-path distance for each dataset by Floyd's algorithm
meanDis=ones(tran,numData);
for i=1:numData
    tmp1=n*(i-1)+1:n*(i-1)+tran;
    D1=D(tmp1,tmp1);
    ind=indNum(:,:,i);
    for j=1:tran
        D1(j,ind((2+K):end,j)) = INF;
    end
    
    D1 = min(D1,D1');    % Make sure distance matrix is symmetric
    %disp('Computing shortest paths');
    for j=1:tran
        D1 = min(D1,repmat(D1(:,j),[1 tran])+repmat(D1(j,:),[tran 1]));
    end
    
    %Do C-Isomap if necessary
    meanDis1=ones(tran,1);
    for j=1:tran
        meanDis1(j)=mean(D1(ind(2:K+1,j),j));
    end
    meanDis(:,i)=meanDis1;
    if options.CIsomap~=0
        %disp('Scaling the shortest paths by C-Isomap...');
        D1=D1./(meanDis1*meanDis1').^0.5;
    end
    
    %Do OOS for the last tesn points if necessary
    if tesn~=0
        D(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = NeighborhoodIsomapOOS(D1(1:tran,1:tran),dis(1:tran,n*(i-1)+tran+1:n*i),meanDis1,options);
    else
        D(n*(i-1)+1:n*i,n*(i-1)+1:n*i)=D1;
    end
end

%This auxiliary function calculates the OOS shortest-path distance by ISOMAP
function [D] = NeighborhoodIsomapOOS(DIN, DOOS, meanDis,options)
K=options.neighborSize;
tesn=size(DOOS,2);
tran = size(DIN,1);
n=tran+tesn;
D=zeros(n,n);
D(1:tran,1:tran)=DIN;

INF =  1000*max(max(DIN))*n;
[~, indOOS] = sort(DOOS);

for i=1:tesn
    DOOS(indOOS((2+K):end,i),i) = INF;
end

for i=1:tesn
    for j=1:tran
        DOOS(:,i) = min(DOOS(:,i),DOOS(j,i)+DIN(j,:)');
    end
end

%C-Isomap OOS if necessary
if options.CIsomap~=0
    meanDisOOS=ones(tesn,1);
    for i=1:tesn
        meanDisOOS(i)=mean(DOOS(indOOS(2:K+1,i),i));
    end
    DOOS=DOOS./(meanDis*meanDisOOS').^0.5;
end
D(1:tran,tran+1:n)=DOOS;
D(tran+1:n,1:tran)=DOOS';
