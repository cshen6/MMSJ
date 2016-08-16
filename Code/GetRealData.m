function [disW,Label,ptt]=GetRealData(dis,Label,trn,tesn,options) %Data Generator into training and testing pairs
% This file randomizes and re-organizes the given distance, so that the output
% contains the training matched data, testing matched data and testing unmatched datain order.
% It is merely used for easier testing like follows:
% [disW,~,~]=GetRealData(dis,0,trn,tesn)
%
% The input dis should be a size n*(numData*n) distance matrix,
% And label can be used too, or 0 if no label is available,
% trn is the number of training matched pairs,
% tesn specifies the number of testing matched & unmatched pairs.
% Note that trn+(numData+1)*tesn should be no larger than n.
% You can also specify numData and permutation in the options, see below.
%
% The output disW is size (trn+2*tesn) * (numData*(trn+2*tesn)) that can be used in ManifoldMatching.m,
% Label contains the permuted label,
% And ptt returns the permutation.
% Note that the version for Euclidean data is called GetRealDataEuc

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

%Options
n=size(dis,1);
numData=floor(size(dis,2)/n);
if numData<1
    error('The size of dissimilarity matrix is insufficient to do joint embedding');
end
if trn+(numData+1)*tesn > n
    error('Insufficient data for splitting into training and testing');
end

if nargin < 4
    error('Too few input arguments');
elseif nargin < 5
    options = struct('numData',numData,'permutation',0,'scaling',1);
end
if ~isfield(options,'permutation')
    options.permutation = 0; %0 for random permutation, -1 for identity permutation, or any n*1 vector for a user-specified permutation.
end
if ~isfield(options,'scaling')
    options.scaling = 1; %Specify scaling or not. 0 for no-scaling, otherwise do scaling based on the Frobenius norm of the first data set.
end
if ~isfield(options,'numData')
    options.numData = numData;
end
if numData<options.numData %check if the user-specified numData conflicts the user-input data
    error('The size of dissimilarity matrix conflcits the number of data specified');
end
ptt=options.permutation;

%Permute and scale the data if needed
if(ptt==0) %use a random permutation
    ptt=randperm(n);
end
if ptt==-1 %use the identity permutation
    ptt=1:n;
end
norm1=norm(dis(:,1:n),'fro');
for i=1:numData
    tmp1=n*(i-1)+1:n*i;
    dis(:,tmp1)=dis(ptt, ptt+n*(i-1)); %permute the dissimilarity matrix for each dataset
    if options.scaling~=0  %scale all datasets according to the norm of the first data
        dis(:,tmp1)=dis(:,tmp1)/norm(dis(:,tmp1),'fro')*norm1;
    end
end
if (length(Label)==length(ptt))
    Label=Label(ptt); %permute the Label if the label has the same size as the permutation
end

%Create the Omnibus matrix for testing
% disO=zeros(numData*trn,numData*trn);%Training data
% disL=zeros(trn,numData*trn);%Flattened
disW=zeros(trn+2*tesn,numData*(trn+2*tesn)); %This one is used for manifold matching purpose
% disM=zeros(numData*trn,numData*tesn);%Matched
% disU=zeros(numData*trn,numData*tesn);%UnMatched
for i=1:numData;
    %     tmp1=n*(i-1)+1:n*(i-1)+trn;
    %     tmp2=trn*(i-1)+1:trn*i;
    %     tmp3=tesn*(i-1)+1:tesn*i;
    %     disO(tmp2,tmp2)=dis(1:trn,tmp1);
    %     disL(:,tmp2)=dis(1:trn,tmp1);
    %     disM(tmp2,tmp3)=dis(1:trn,n*(i-1)+trn+1:n*(i-1)+trn+tesn);
    %     disU(tmp2,tmp3)=dis(1:trn,tmp3+trn+n*(i-1)+tesn);
    
    indRow=[1: trn+tesn trn+i*tesn+1:trn+(i+1)*tesn ];
    indCol=n*(i-1) + indRow;
    disW(:, (i-1)*(trn+2*tesn)+1: i*(trn+2*tesn))=dis(indRow, indCol);
    %     for j=i+1:numData
    %         tmp2=n*(j-1)+1:n*j;
    %         disO(tmp1,tmp2)=(dis(:,tmp1)+dis(:,tmp2))/2;
    %         disO(tmp2,tmp1)=disO(tmp1,tmp2)';
    %     end
end
