function [p,empi,empi2]=plotPower(sol,numData,tesn,div)
% The plotPower function is an auxiliary function to calculate the matching
% test power based on sol.
%
% The input sol should be of size dim*(numData*(training+tesn+tesn)), and
% the first tesn points are matched pairs, while the second tesn points are
% unmatched pairs.
% div is the scale parameter for the type 1 error level, i.e., 1/div is the
% scale. And we usually use 20 in our experiment so that p(2) is the power
% at 0.05 level.
%
% The output p is a vector of size div+1 containing the power at respective levels,
% empi is the empirical distribution of the testing matched pair,
% empi2 is the empirical distribution of the testing unmatched pair.

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

n=size(sol,2)/numData;
tran=n-2*tesn;
p=zeros(div+1,1);
empi=zeros(tesn,1);
empi2=zeros(tesn,1);
ca=zeros(div+1,1);
%Calculate the null distribution
for i=1:tesn
    for j=1:numData
        for k=j+1:numData
            empi(i)=empi(i)+sum((sol(:,(j-1)*n+tran+i)-sol(:,(k-1)*n+tran+i)).^2);
        end
    end
end
%order empirical distribution of the null
empi=sort(empi, 'descend');
for l=1:div
    ca(l)=empi(ceil(tesn/div*l));
end
ca(div+1)=ca(div)*0.8;
%Calculate the alternative distribution
for i=1:tesn
    for j=1:numData
        for k=j+1:numData
            empi2(i)=empi2(i)+sum((sol(:,(j-1)*n+tran+tesn+i)-sol(:,(k-1)*n+tran+tesn+i)).^2);
        end
    end
end
%Output power
for l=1:div+1
    p(l)=mean(empi2>ca(l));
end