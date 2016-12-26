function [sol,acc,power]=MatchingTestEuc(dis,dim,tran,tesn,reps,options)
%OOS inference test on real data, using both global and local
%dissimilarities (transformed by ISOMAP)
%GE and GF are the dissimilarities.
%m is the estimated true dimension used for regularization, n is the number of points in each dataset.
%dim is the upper limit projecting dimension, tesn is the number of test point.
%div specify the interval between critical level, reps is the number of repeatition. K is the
%choice of neighborhood.
%jointN=0 for separate neighborhood, 1 for joint neighborhood, 2 for nonlinear vs linear
%per is the specified permutation to use for real data inference test, insert 0 for random run.

%Options
%Specify options
if nargin < 4
    error('Too few input arguments');
elseif nargin < 5
    options = struct('nonlinear',0,'match',2,'neighborSize',max(10,dim+5),'jointSelection',0,'CIsomap',0,'weight',1,'scaling',1,'numData',2,'maxIter',-1,'oos',0,'permutation',0);
end

if ~isfield(options,'weight')
    options.weight = 1;
end
if ~isfield(options,'nonlinear')
    options.nonlinear = 0;
end
if ~isfield(options,'match')
    options.match = 2;
end
if ~isfield(options,'jointSelection')
    options.jointSelection = 0;
end
if ~isfield(options,'CIsomap')
    options.CIsomap = 0;
end
if ~isfield(options,'neighborSize')
    options.neighborSize = max(10,dim+5);
end
if ~isfield(options,'scaling')
    options.scaling = 1;
end
if ~isfield(options,'numData')
    options.numData = 2;
end
if ~isfield(options,'maxIter')
    options.maxIter = -1;
end
if ~isfield(options,'oos')
    options.oos = 2*tesn;
end
if ~isfield(options,'permutation')
    options.permutation = 0;
end
%Check number of datasets
if options.numData<2
    error('The size of dissimilarity matrix is insufficient to do matching');
end
per=options.permutation;
numData=options.numData;
n=tran+2*tesn;

div=20;
%Initialization
% x=0:1/div:1;
seqDiv=div+1;
power=zeros(seqDiv,1);
%sol=zeros(dim,n*numData,3,reps);
%dCorr=zeros(reps,1);
%Start testing loop
index=[];
acc=0;
for r=1:reps
    %r
    if per==0
        options.permutation=0;
    else
        options.permutation=per(r,:);
    end
    [disO,~,~]=GetRealDataEuc(dis,0,tran,tesn,options);
    try
        sol=ManifoldMatchingEuc(disO,dim,options);
        power=power+plotPower(sol,numData,tesn,div)/reps;
        acc=acc+ResultsAcc(sol,numData,tran,tesn)/reps;
        power=power(2);
    catch
        index = [index r];
        disp('error, continue to next');
    end
end
%End testing

% %Start plot
% % index
% p1=sum(power(:,:,1),2)/(reps-sum(index~=0));
% p2=sum(power(:,:,2),2)/(reps-sum(index~=0));
% p3=sum(power(:,:,3),2)/(reps-sum(index~=0));
% 
% % p1=mean(power(:,:,1),2);
% % p2=mean(power(:,:,2),2);
% % p3=mean(power(:,:,3),2);
% x=0:0.05:1;
% figure
% hold on
% plot(x,p1,'k-',x,p2,'r-', x,p3,'b-','LineWidth',2);
% legend('Joint MDS','P \circ M','CCA', 'Location','SouthEast');
% xlabel('Type 1 Error Level')
% ylabel('Matching Power')
% axis([0 1 0 1]);
% title('Matching Test using Joint Isomap')
% hold off