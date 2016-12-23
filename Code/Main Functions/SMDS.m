function [sol1,percent]=SMDS(dis,dim,weight,options)
% This is the MDS function for a distance matrix, which supports both
% classical MDS and raw-stress MDS (SMACOV).
%
% The raw-sress MDS is implemented according to the following paper:
% J. Leeuw and P. Mair, "Multidimensional Scaling Using Majorization:
% SMACOF in R", Journal of Statistical Software, vol. 31, no. 3, pp. 1-30,
% 2009.
%
% For the input, dis should be a size n*n distance matrix,
% weight can be an n*n weight matrix for raw-stress MDS, or 0 for CMDS.
% dim is the embedding dimension,
%
% The output sol1 is size dim*n,
% and the percentage of remaining stress un-explained by raw-stress MDS.

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

% Specify options for raw-stress MDS.
% maxIter limits the maximum iteration, epsilon specifies the stress limit
% to stop iteration. If maxIter is -1, CMDS is used; if maxIter is 0, Matlab
% default stress MDS is used.
if nargin < 3
    error('Too few input arguments');
elseif nargin < 4
    options = struct('maxIter',200,'epsilon',0.1);
end
if ~isfield(options,'maxIter')
    options.maxIter = 200;
end
if ~isfield(options,'epsilon')
    options.epsilon = 0.1;
end

%Initialization
percent=0;
n=size(dis,1);
H=eye(n)-(1/n)*ones(n,n);
iter=options.maxIter;
%sim=0;
dis=dis.^2; %Square the distance

try
    %     if sim==0
    [U,S,~]=svds(-H*dis*H/2,dim);
    %     else
    %         [U,S,~]=svds(dis,dim);
    %     end
    sol1=S(1:dim,1:dim).^(0.5)*(U(:,1:dim))'; %use the CMDS start if possible
    if weight ==0 | iter==-1 %output the CMDS result if weight is 0 or iter is -1
        return;
    end
catch
    disp('Use a randomized start in SMARCOV');
    sol1=rand([dim,n]);%if CMDS fails, use a randomized start
end

if size(weight)~=size(dis)
    weight=ones(size(dis,1),size(dis,2)); %If weight is not corrected specified, use uniform weight
end
if iter==0
    sol1=mdscale(dis,dim,'Criterion','metricstress','Weights',weight,'Start','random')'; %Use the Matlab weighted MDS algorithm if iteration is 0
    return;
end

epsilon=options.epsilon;
V=-weight+diag(sum(weight));
%VI=inv(V+ones(n,n)/n)-ones(n,n)/n;
VI=pinv(V);

DD=squareform(pdist(sol1'));
stress1=sum(sum(weight.*(dis-DD).^2));
percent=stress1;

%Iterative Minimization
for i=1:iter
    B=-weight.*dis;
    for j=1:n
        for k=j:n
            if (DD(j,k)==0)
                B(j,k)=0;
            else
                B(j,k)=B(j,k)/DD(j,k);
            end
            B(k,j)=B(j,k);
        end
    end
    B=B+diag(sum(-B));
    sol2=(VI*B*sol1')';
    DD=squareform(pdist(sol2'));
    stress2=sum(sum(weight.*(dis-DD).^2));
    ee=stress1-stress2;
    if(ee < epsilon)
        %disp('Minimization finished');
        break;
    end
    sol1=sol2;
    stress1=stress2;
end

if i==iter
    %disp('Iteration limits reached.');
end
percent=stress1/percent;