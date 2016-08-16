function [sol1, T] = GeneralCCA(sol, numData, dim,options)
% General CCA Runs the Generalized Canonical Correlation Analysis for the
% input, maximize the sum of pairwise correlation iteratively,
% and returns the transformed data as well as the GCCA transformation.
% 
% The algorithm is implemented based on the following paper:
% A. Tenenhaus and M. Tenenhaus, "Regularized Generalized Canonical
% Correlation Analysis", Psychometrika, vol. 76, no. 2, pp. 257-284, 2013.
%
% For the input, sol should be size d* (numData*n), where d is the dimension, numData is
% the numData for matching, n is the size of each data. Namely the input
% should contain all data in a single input sol.
% dim is the GCCA dimension, and should be no larger than d.
%
% The output sol1 is size dim* (numData*n), and contain all matched data sets in order.
% T contains the GCCA transformation in order, which is of size
% d*(dim*numData) and of norm 1 for each vector.

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
if nargin < 3
    error('Too few input arguments');
elseif nargin < 4
    options = struct('maxIter',30);
end
if ~isfield(options,'maxIter') || options.maxIter==-1;
    options.maxIter = 30; % Maximum iteration used.
end

%Initialization
if numData<2
    error('Less than two datasets, no need for CCA');
end
iter=options.maxIter;
d=size(sol,1);
n=size(sol,2)/numData;
Sigma=zeros(d,d,numData);
SigmaInv=zeros(d,d,numData);
T=zeros(d,dim*numData);
sol1=zeros(dim,n*numData);
H=eye(n)-(1/n)*ones(n,n);

%The initial guess of canonical vectors
for i=1:numData
    tmp1=n*(i-1)+1:n*i;
    tmp2=dim*(i-1)+1:dim*i;
    tmp=sol(:,tmp1)*H*H*sol(:,tmp1)'/(n-1);
    Sigma(:,:,i)=tmp;
    SigmaInv(:,:,i)=tmp^-1;
    if i==1
        continue;
    end
    [T1,T2]=canoncorr(sol(:,1:n)',sol(:,tmp1)');
    if i==2
        T(:,1:dim)=T1(:,1:dim);
    end
    T(:,tmp2)=T2(:,1:dim);
end

%Iteratively find the GCCA vectors if more than two datasets
if numData>2
    %Find the canonical coefficients for each dimension
    for dd=1:dim;
        %The starting coefficients
        a=zeros(d,numData);
        for i=1:numData
            tmpA=T(:,dim*(i-1)+dd);
            tmpA=(tmpA'*SigmaInv(:,:,i)*tmpA)^-0.5*SigmaInv(:,:,i)*tmpA;
            a(:,i)=orth(tmpA,T(:,dim*(i-1)+1:dim*(i-1)+dd-1),Sigma(:,:,i),dd-1);
        end
        %Iteration
        for s=1:iter
            for i=1:numData
                tmp1=n*(i-1)+1:n*i;
                v1=0;
                for j=1:numData
                    tmp2=n*(j-1)+1:n*j;
                    if j==i
                        continue;
                    end
                    v1=v1+a(:,j)'*sol(:,tmp2);
                end
                tmpA=sol(:,tmp1)*v1';
                tmpA=(tmpA'*SigmaInv(:,:,i)*tmpA)^-0.5*SigmaInv(:,:,i)*tmpA;
                a(:,i)=orth(tmpA,T(:,dim*(i-1)+1:dim*(i-1)+dd-1),Sigma(:,:,i),dd-1);
            end
        end
        %Normalize all canonical vectors to norm 1
        for i=1:numData
            T(:,dim*(i-1)+dd)=a(:,i)/norm(a(:,i));
        end
    end
end

%Calculate the GCCA-transformed data
for i=1:numData
    tmp1=n*(i-1)+1:n*i;
    tmp2=dim*(i-1)+1:dim*i;
    tmpSol=sol(:,tmp1)';
    tmpSol=tmpSol - repmat(mean(tmpSol),n,1);
    sol1(:,tmp1)=T(:,tmp2)'*tmpSol';
end

%The Gram-Schmidt Orthogonalization
function a=orth(a1,T1,Sigma1,d)
m=size(a1,1);
if (d==0)
    a=a1;
else
    a=zeros(m,1);
    T1=Sigma1*T1;
    for i=1:d;
        b=T1(:,i);
        a=a+a1'*b/(b'*b)*b;
    end
    a=a1-a;
end