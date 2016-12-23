function [eucli] = NeighborhoodLLE(dis,numData,dim,options)
% This function takes in the distance matrix, and do separate/joint LLE
% possibly with out-of-sample (OOS) technique.
%
% Note that we use the Gram Matrix LLE version for distance input, see:
% S. Roweis and L. Saul, "Think Globally, Fit Locally: Unsupervised Learning of Low Dimensional
% Manifolds", Journal of Machine Learning Research, vol. 4, pp. 119-155, 2003.
%
% For the input, the distance can be of size n*(numData*n) of (numData*n)*(numData*n)
% dim is the embedding dimension,
% and neighborhood size and oos and joint selection can be
% specified in the option.
%
% Output eucli is size dim*(numData*n) containing the embedded data.

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
if nargin < 3
    error('Too few input arguments');
elseif nargin < 4
    options = struct('jointSelection',1,'neighborSize',max(10,dim+5), 'oos',0);
end
if ~isfield(options,'jointSelection')
    options.jointSelection = 1; %0 for separate neighborhood selection, otherwise do joint neighborhood selection
end
if ~isfield(options,'neighborSize')
    options.neighborSize = max(10,dim+5); %default neighorbood size. Note that it should be larger than dim.
end
if ~isfield(options,'oos')
    options.oos = 0; %The number of oos embedded points. Note that it is used for last points in each dataset.
end
K=options.neighborSize;
tesn=options.oos;
tran=n-tesn;

indNum=zeros(K,tran,numData); %store the neighborhood graph info for all datasets
eucli=zeros(dim,numData*n);
tol=1e-3; %always add a regularization term. Note that original LLE adds the term only when K>d, d the original dimension of the data.

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
    %disp('Separately finding K nearest neighbors');
    for i=1:numData
        tmp1=n*(i-1)+1:n*(i-1)+tran;
        [~,ind]=sort(D(tmp1,tmp1));
        indNum(:,:,i)=ind(2:(1+K),:);
    end
else
    %Joint neighborhood selection constructs the same graph for all data
    %disp('Jointly finding K nearest neighbors');
    tmpSum=0;
    for i=1:numData
        tmp1=n*(i-1)+1:n*(i-1)+tran;
        tmpSum=tmpSum+D(tmp1,tmp1);
    end
    [~,ind]=sort(tmpSum);
    ind=ind(2:(1+K),:);
    indNum=repmat(ind,[1,1,numData]);
end

%Solve the reconstruction weights for each data set
%Gram Matrix Version, see in
%disp('Solving for reconstruction weights');
for i=1:numData
    tmp1=n*(i-1)+1:n*(i-1)+tran;
    D1=D(tmp1,tmp1);
    W1=zeros(K,n);
    neighborhood1=indNum(:,:,i);
    for ii=1:tran
        tempSum1=zeros(K+1,1);
        tempSum1(1)=sum(D1(ii,neighborhood1(:,ii)));
        for j=1:K
            tempSum1(j+1)=sum(D1(neighborhood1(j,ii),neighborhood1(:,ii)));
        end
        tempSumT1=sum(tempSum1);
        P1=(repmat(tempSum1',K+1,1)+repmat(tempSum1,1,K+1))/(K+1)-tempSumT1/(K+1)^2-D1([ii neighborhood1(:,ii)'],[ii neighborhood1(:,ii)']);
        P1=P1/2; %The Dot Product Matrix
        G1=-repmat(P1(1,2:K+1),K,1)-repmat(P1(1,2:K+1)',1,K)+P1(2:K+1,2:K+1)+P1(1,1); %The Gram Matrix
        
        G1=G1+eye(K,K)*tol*trace(G1);
        W1(:,ii) = G1\ones(K,1);
        W1(:,ii) = W1(:,ii)/sum(W1(:,ii));
    end
    
    %Compute embedding from eigenvectors of cost matrix M1
    M1 = sparse(1:tran,1:tran,ones(1,tran),tran,tran,4*K*tran);
    for ii=1:tran
        w = W1(:,ii);
        jj = neighborhood1(:,ii);
        M1(ii,jj) = M1(ii,jj) - w';
        M1(jj,ii) = M1(jj,ii) - w;
        M1(jj,jj) = M1(jj,jj) + w*w';
    end;
    
    %Eigenvector calculation
    [Y1,~] = eigs(M1,dim+1,0);
    Y1 = Y1(:,dim:-1:1)'*sqrt(tran); % bottom evect is [1,1,1,1...] with eval 0
    eucli(:,tmp1)=Y1;
    
    %Do OOS for the last tesn points if necessary
    if tesn~=0
        eucli(:, n*(i-1)+tran+1:n*i)=NeighborhoodLLEOOS(D1(1:tran,1:tran),dis(1:tran,n*(i-1)+tran+1:n*i),Y1, options);
    end
end

%LLE OOS implementation
function [eucliOOS] = NeighborhoodLLEOOS(DIN, DOOS, eucli, options)
tran = size(DIN,1);
tesn=size(DOOS,2);
K=options.neighborSize;
dim=size(eucli,1);
tol=1e-3;
eucliOOS=zeros(dim,tesn);
for ii=1:tesn
    tesM=DOOS(:,ii);
    [tmp, indM] = sort(tesM);
    neighborhood1 = indM(1:K);
    
    W1 = zeros(K,1);
    %Gram Matrix Version
    D1=[DIN tesM]';
    D1=[D1 [tesM;0]];
    tempSum1=zeros(K+1,1);
    tempSum1(1)=sum(D1(ii,neighborhood1));
    for j=1:K
        tempSum1(j+1)=sum(D1(neighborhood1(j),neighborhood1));
    end
    tempSumT1=sum(tempSum1);
    P1=(repmat(tempSum1',K+1,1)+repmat(tempSum1,1,K+1))/(K+1)-tempSumT1/(K+1)^2-D1([tran+1 neighborhood1'],[tran+1 neighborhood1']);
    P1=P1/2; %The Dot Product Matrix
    G1=-repmat(P1(1,2:K+1),K,1)-repmat(P1(1,2:K+1)',1,K)+P1(2:K+1,2:K+1)+P1(1,1); %The Gram Matrix
    
    G1=G1+1*eye(K,K)*tol*trace(G1);
    W1 = G1\ones(K,1);
    W1 = W1/sum(W1);
    eucliOOS(:,ii)=eucli(:,neighborhood1)*W1;
end
