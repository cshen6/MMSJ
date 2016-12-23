function [solOOS]=OOSMDS(dis, trn, dim,sol)
% This is the out-of-sample (OOS) embedding function for a distance matrix.
% It only implements OOS for classical multidimensional scaling, but is sufficient to deal with Isomap and LLE
% oos by combining it with the respective functions in NeighborhoodIsomap and
% Neighbohood LLE.
%
% For the input, dis should be a size n*n distance matrix,
% trn is the number of in-sample embedded point,
% dim is the embedding dimension,
% sol is the embedded in-sample data, which should be size dim*trn
%
% The output solOOS is size dim*n, which keeps the first trn points as sol,
% and concatenate the last n-trn points by OOS embedding.

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

n=size(dis,2);
dis=dis.^2; %Square the distance.

%If sol is not provided, do CMDS for in-sample points.
if nargin < 4 || size(sol,1)~=dim
    sol=SMDS(dis(1:trn,1:trn), dim, 0);
end

tv=[ones(1,trn) zeros(1,n-trn)];
H=eye(n)-(1/trn)*ones(n,1)*tv;
dis=-H*dis*H'/2;
dis=dis(1:trn,trn+1:n);
solOOS=(sol*sol')^(-1)*sol*dis;

solOOS = [sol solOOS];