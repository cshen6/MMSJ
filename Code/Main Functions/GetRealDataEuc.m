function [disO,Label,ptt]=GetRealDataEuc(dis,Label,trn,tesn,options) %Data Generator into training and testing pairs
% This file randomizes and re-organizes the given data, so that the output
% contains the training matched data, testing matched data and testing unmatched datain order.
% It is merely used for easier testing like follows:
% [disW,~,~]=GetRealDataEuc(dis,0,trn,tesn)
%
% The input dis should be a size d*(numData*n) data,
% And label can be used too, or 0 if no label is available,
% trn is the number of training matched pairs,
% tesn specifies the number of testing matched & unmatched pairs.
% Note that trn+(numData+1)*tesn should be no larger than n.
% You can also specify numData and permutation in the options, see below.
%
% The output disW is size d * (numData*(trn+2*tesn)) that can be used in ManifoldMatchingEuc.m,
% Label contains the permuted label,
% And ptt returns the permutation.
% Note that the version for distance data is called GetRealData

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
if nargin < 4
    error('Too few input arguments');
elseif nargin < 5
    options = struct('numData',2,'permutation',0,'scaling',1);
end
if ~isfield(options,'permutation')
    options.permutation = 0; %0 for random permutation, -1 for identity permutation, or any n*1 vector for a user-specified permutation.
end
if ~isfield(options,'scaling')
    options.scaling = 1; %Specify scaling or not. 0 for no-scaling, otherwise do scaling based on the Frobenius norm of the first data set.
end
if ~isfield(options,'numData')
    options.numData = 2;
end
numData=options.numData;
if numData*(trn+2*tesn)>size(dis,2)
    error('The size of dissimilarity matrix is insufficient to do joint embedding');
end
ptt=options.permutation;
n=size(dis,2)/numData;
d=size(dis,1);

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
    dis(:,tmp1)=dis(:, ptt+n*(i-1)); %permute the dissimilarity matrix for each dataset
    if options.scaling~=0  %scale all datasets according to the norm of the first data
        dis(:,tmp1)=dis(:,tmp1)/norm(dis(:,tmp1),'fro')*norm1;
    end
end
if (length(Label)==length(ptt))
    Label=Label(ptt); %permute the Label if the label has the same size as the permutation
end

%Create the Omnibus matrix for testing
disO=zeros(d,numData*(trn+2*tesn));
for i=1:numData;
    indRow=[1: trn+tesn trn+i*tesn+1:trn+(i+1)*tesn ];
    indCol=n*(i-1) + indRow;
    disO(:, (i-1)*(trn+2*tesn)+1: i*(trn+2*tesn))=dis(:, indCol);
end