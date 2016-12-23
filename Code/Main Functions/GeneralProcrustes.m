function [sol1,Q]=GeneralProcrustes(sol, numData)
% General Procrustes Runs the Procrustes matching for the
% input, and minimize the Procrustes error in a greedy way.
%
% For example, in case of three data, we try to minimize
% \|Q1X1-X3\|^2+  \|Q2X2-X3\|^2+ \|Q1X1-Q2X2\|^2,
% for which we find Q1 first by ignoring all terms involving Q2, and find Q2 next by
% fixing Q1.
%
% For the input, sol should be size dim* (numData*n), where dim is the dimension, numData is
% the numData for matching, n is the size of each data. Namely the input
% should contain all data in a single input sol.
%
% The output sol1 is of the same size, and contain all Procrustes-transformed data sets in order.
% Q contains the Procrustes transformation in order, which is of size
% dim*dim*numData, with the last tranformation being identity.

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

%Initialization
if numData<2
    error('Less than two datasets, no need for Procrustes');
end
dim=size(sol,1);
n=size(sol,2)/numData;
Q=zeros(dim,dim,numData);
Q(:,:,numData)=eye(dim); %The last transformation is always identity.
sol1=sol;

for i=1:numData-1
    tmp1=n*(i-1)+1:n*i;
    sum=sol(:,n*(numData-1)+1:n*numData)*sol(:,tmp1)';
    for j=1:i-1
        sum=sum+Q(:,:,j)*sol(:, n*(j-1)+1:n*j)*sol(:,tmp1)';
    end
    [U,~,V]=svd(sum);
    Q1=U*V'; %Solving the Procrustes problem by SVD
    sol1(:,tmp1)=Q1*sol(:,tmp1);
    Q(:,:,i)=Q1;
end