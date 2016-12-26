map1=zeros(5,3);
cm = 3;

switch cm
    case 1
        cmap = brewermap(8,'Dark2');
    case 2
        cmap = brewermap(8,'Set2');
    case 3
        gr = [0,1,0];
        ma = [1,0,1];
        cy = [0,1,1];
        cmap(1,:) = gr;
        cmap(2,:) = ma;
        cmap(3,:) = cy;
        cmap(4,:) = [1,0,0];
        cmap(5,:) = [0 0 0];
        %cmap(6,:) = [0.5 0.5 0];
    case 4
        cmap(1,:) = [166,206,227]/255;
        cmap(2,:) = [31,120,180]/255;
        cmap(3,:) = [178,223,138]/255;
    case 5
        cmap(1,:) = [102,194,165]/255;
        cmap(2,:) = [ 252,141,98]/255;
        cmap(3,:) = [141,160,203]/255;
end

map1=cmap;

%figure1-4
set(groot,'defaultAxesColorOrder',map1);

%Matching Power New!
%Generate new data
n=5000;
[X_data, labels, Y_data] = generate_data('brokenswiss', n, 0.05); %1
[X_data, labels, Y_data] = generate_data('twinpeaks', n, 0.05); %2 ok?!
Y_data=[Y_data zeros(1, n)']';

[X_data, labels, Y_data] = generate_data('helix', n, 0.05); %3: working with Lap and Iso only?
[X_data, labels, Y_data] = generate_data('intersect', n, 0.05); %4
Y_data=[Y_data zeros(2, n)']';

[X_data, labels, Y_data] = generate_data('difficult', n, 0.05); %4
Y_data=[Y_data zeros(5, 3125)']';

X_data=X_data';
disEuc=[X_data Y_data];
save SwissBroken.mat;

%Swiss roll
%Debug
clear;
load SwissRoll.mat
X_data=X_data(:,1:5000);
Y_data=Y_data(:,1:5000);
Y_data=[Y_data' zeros(1, 5000)']';
Z_data=[cos(Y_data(1,:)).*sin(Y_data(2,:));sin(Y_data(1,:)).*sin(Y_data(2,:));cos(Y_data(1,:))];
n=1000;numData=2;dim=2;tesn=100;K=10;iter=-1;tran=n;reps=1;
cc=color(1:n+3*tesn);%cc=ones(1,3*n);
%cc=ones(n+3*tesn,1);
disEuc=[X_data(:,1:n+2*tesn) [Y_data(:,1:n+tesn) Y_data(:,n+2*tesn+1:n+3*tesn)]];%nonlinear vs linear
dis=squareform(pdist(disEuc'));
ss=size(dis,2)/2;
dis=[dis(1:ss, 1:ss) dis(ss+1:end, ss+1:end)];
% options = struct('nonlinear',0,'match',2,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',0,'numData',numData,'permutation',per,'oos',0,'maxIter',iter);
% sol=ManifoldMatchingEuc(disEuc,dim,options);
options = struct('nonlinear',1,'match',2,'neighborSize',K,'jointSelection',1,'CIsomap',0,'weight',1,'scaling',0,'numData',numData,'permutation',per,'oos',2*tesn,'maxIter',iter);
sol=ManifoldMatching(dis,dim,options);
power=plotPower(sol,numData,tesn,20);
%check k
dim=3; lim=1; rep1=1; rep2=1000;nn=100;
CorrPermTest([X_data(:,1:nn)' Y_data(:,1:nn)'],nn,dim,lim,rep1,rep2);
%check first data
[p,e1,e2]=plotVelocity([sol(:,1:n) sol(:,n+2*tesn+1:2*n+2*tesn)],options.numData);%matched power
hold on
scatter(sol(1,1:n),sol(2,1:n),20,cc(1:n),'o'); %training
scatter(sol(1,n+2*tesn+1:2*n+2*tesn),sol(2,n+2*tesn+1:2*n+2*tesn),20,cc(1:n),'+');
hold off
[p,e1,e2]=plotVelocity([sol(:,n+1:n+tesn) sol(:,2*n+2*tesn+1:2*n+3*tesn)],options.numData);%matched power
hold on
scatter(sol(1,n+1:n+tesn),sol(2,n+1:n+tesn),30,cc(n+1:n+tesn),'o'); %matched
scatter(sol(1,2*n+2*tesn+1:2*n+3*tesn),sol(2,2*n+2*tesn+1:2*n+3*tesn),30,cc(n+1:n+tesn),'+');
hold off
[p,e1,e2]=plotVelocity([sol(:,n+tesn+1:n+2*tesn) sol(:,2*n+3*tesn+1:2*n+4*tesn)],options.numData);%unmatched power
hold on
scatter(sol(1,n+tesn+1:n+2*tesn),sol(2,n+tesn+1:n+2*tesn),30,cc(n+tesn+1:n+2*tesn),'o'); %unmatched
scatter(sol(1,2*n+3*tesn+1:2*n+4*tesn),sol(2,2*n+3*tesn+1:2*n+4*tesn),30,cc(n+2*tesn+1:n+3*tesn),'+');
hold off
%%%

sol=solIsoJ;
[p0,e0,~]=plotVelocity([sol(:,1:tran) sol(:,tran+2*tesn+1:2*tran+2*tesn)],options.numData);%power
[p1,e1,~]=plotVelocity([sol(:,tran+1:tran+tesn) sol(:,2*tran+2*tesn+1:2*tran+3*tesn)],options.numData);%power
[p2,e2,~]=plotVelocity([sol(:,tran+tesn+1:tran+2*tesn) sol(:,2*tran+3*tesn+1:2*tran+4*tesn)],options.numData);%power
%classfication
LDAModel = fitcdiscr(sol(:,1:tran)',Label(1:tran),'DiscrimType','linear');%'quadratic'    
testLabel=predict(LDAModel,sol(:,tran+1:tran+tesn)');
errorLDA=mean(Label~=testLabel);
%YouTube
%embedding check
ind=1;
[p,e1,e2]=plotVelocity([sol(:,tran+1:tran+tesn,ind) sol(:,2*tran+2*tesn+1:2*tran+3*tesn,ind)],2);%matched testing
[p,e1,e2]=plotVelocity([sol(:,tran+tesn+1:tran+2*tesn,ind) sol(:,2*tran+3*tesn+1:2*tran+4*tesn,ind)],2);%unmatched testing
%
%

%3D curve plotting w.r.t. K and dim
load Wiki_Data.mat;
start=9;en=30;ind=2;dimMax=30;
tran=500;tesn=100;numData=2;reps=100;iter=-1;scale=1;
p1=zeros(dimMax,en-start); p2=zeros(dimMax,en-start); p3=zeros(dimMax,en-start);
p4=zeros(dimMax,en-start); p5=zeros(dimMax,en-start); p6=zeros(dimMax,en-start);
%load Wiki_R.mat;
dis=[TE GE];
for k=10:20;
     k
     for dim=1:dimMax;
         dim
        options = struct('nonlinear',1,'match',2,'neighborSize',k,'jointSelection',1,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
         [sol,power]=MatchingTest(dis,dim,tran,tesn,reps,options);
         p1(dim, k-start) = mean(power(ind,:,1));
         p2(dim, k-start) = mean(power(ind,:,2));
         p3(dim, k-start) = mean(power(ind,:,3));
         options = struct('nonlinear',1,'match',2,'neighborSize',k,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
         [sol,power]=MatchingTest(dis,dim,tran,tesn,reps,options);
         p4(dim, k-start) = mean(power(ind,:,1));
         p5(dim, k-start) = mean(power(ind,:,2));
         p6(dim, k-start) = mean(power(ind,:,3));
         save Wiki_R.mat;
     end
end
k=start+1:en;dd=2:dimMax;
imagesc(k,dd,p2(dd,:));
colorbar();
caxis([0.35 0.55]);
xlabel('Neighborhood Choice')
%xlim([start+1,en]);
%ylim([2,dimMax]);
ylabel('Embedding Dimension')
%zlabel('Inference Power')
title('MMSJ Matching Power Heat-map')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wiki Data
clear;
load Wiki_Data.mat
dim1=10;dim2=10;Kmin=10;Kmax=10;filename='WikiTETF';
dis=[TE TF]; %only
Label=Label+1;
clusterM='Laplacian';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='Kmeans';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='Hclust';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='GMM';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='Classify';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
%
clear;
load BrainCP.mat
dim1=3;dim2=3;Kmin=3;Kmax=10;n=42;filename='CxP';
dis=[distC distP]; %only
y=mvnrnd(zeros(n,dim),eye(dim),n);
distPInd=squareform(pdist(y));
dis=[distPInd distP]; %only
Label=3;
clusterM='Laplacian';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='Kmeans';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='Hclust';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='GMM';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
clusterM='Classify';
EmbedAndCluster(dis,Label,dim1,dim2,Kmin,Kmax,clusterM,filename)
%
clear
load('BrainHippoShape.mat')
dim1=3;dim2=5;Kmin=5;Kmax=20;filename='BrainShape';
dis=[LMLS, LMRS];
y=squareform(pdist(Label));
%y=(y>0)+1;
y=y+1;
for i=1:n
    y(i,i)=0;
end
dis=[LMLS,y];
%Brain Data
clear
load('BrainCP.mat')
dis=[distC, distP];
% m=5;
% x=SMDS(distC, m,0);
% y=SMDS(distP, m,0);
% disEuc=[x y];
dim=5;numData=2;iter=-1;scale=1;
KMAX=10;groups=10;s1=3;
%Cluster compare with truth
% %classfication
% LDAModel = fitcdiscr(sol1(:,1:tran)',Label(1:tran),'DiscrimType','linear');%'quadratic'    
% testLabel=predict(LDAModel,sol1(:,tran+1:tran+tesn)');
% errorLDA1=mean(Label(tran+1:tran+tesn)~=testLabel);
% %
% LDAModel = fitcdiscr(sol2(:,1:tran)',Label(1:tran),'DiscrimType','linear');%'quadratic'    
% testLabel=predict(LDAModel,sol2(:,tran+1:tran+tesn)');
% errorLDA2=mean(Label(tran+1:tran+tesn)~=testLabel);
% %
% LDAModel = fitcdiscr(solIsoJ(:,1:tran)',Label(1:tran),'DiscrimType','linear');%'quadratic'    
% testLabel=predict(LDAModel,solIsoJ(:,tran+1:tran+tesn)');
% errorLDA3=mean(Label(tran+1:tran+tesn)~=testLabel);

%%%%Deep JOFC
load Wiki_Data.mat
tran=500;tesn=100;K=20;numData=2;reps=100;iter=-1;scale=1;m=50;
dim3=[500; 100; 10];
dim2=[200; 10];
dim1=10;
dis=[TE TF];
ss=size(dis,2)/numData;
TEEuc=SMDS(dis(1:ss, 1:ss),m,0);
TFEuc=SMDS(dis(1:ss, ss+1:2*ss),m,0);
disEuc=[TEEuc TFEuc];
dis2=[squareform(pdist(TEEuc')) squareform(pdist(TFEuc'))];
options = struct('nonlinear',0,'match',2,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[sol,power1, dCorr1]=MatchingTest(dis,dim1,tran,tesn,reps,options);
options = struct('nonlinear',0,'match',2,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[sol,power2, dCorr2]=MatchingTest(dis2,dim1,tran,tesn,reps,options);
%options = struct('nonlinear',0,'match',3,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
%[sol,power3, dCorr3]=MatchingTest(dis,dim3,tran,tesn,reps,options);
save('ResultMatchingWikiTETFDeep', 'power1','dCorr1',  'power2','dCorr2','tran','m','dim1','dim2','tesn','K','iter','reps','numData')
%Deep CCA Matching Plot
x=0:0.05:1;
indM=3;
figure
plot(x,mean(power1(:,:,indM),2),'bo-', x,mean(power2(:,:,indM),2),'ro-');
legend('Original Distance  One Layer', 'Original Distance Two Layer','Location','SouthEast');
xlabel('Type 1 Error Level')
ylabel('Matching Power')
title('CCA Matching Test')
ylim([0 1])
%

%title('Matching Test using Original Distance')
options = struct('nonlinear',1,'match',3,'neighborSize',K,'jointSelection',1,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solIsoJ,powerIsoJ, dCorrIsoJ]=MatchingTest(dis,dim,tran,tesn,reps,options);
%title('Matching Test using Joint Isomap')
options = struct('nonlinear',1,'match',3,'neighborSize',K,'jointSelection',0,'weight',1,'scaling',scale,'numData',numData,'maxIter',iter,'permutation',per,'oos',2*tesn); 
[solIsoS,powerIsoS, dCorrIsoS]=MatchingTest(dis,dim,tran,tesn,reps,options);


%%%Classification
clear
load('Wiki_Data.mat')
n=150;
dis=[GE(1:n,1:n), TE(1:n,1:n)];dim=10;trn=100;reps=100;
sol=SMDS(GE,dim,0);
warning('off','all')
[error1,error2]=EmbedAndClassify(sol', Label, dim, trn,reps);
warning('on','all')
match=1;K=5;numData=2;iter=-1;dim=3;nonlinear=1;
options = struct('nonlinear',nonlinear,'match',match,'neighborSize',K,'jointSelection',0,'weight',1,'numData',numData,'maxIter',iter,'oos',0);
[sol, ~]=ManifoldMatching(dis,dim,options);
[error3,error4]=EmbedAndClassify(sol(:,1:n)', Label, dim, trn,reps);
match=0;K=5;numData=1;iter=-1;dim=3;nonlinear=1;
options = struct('nonlinear',nonlinear,'match',match,'neighborSize',K,'jointSelection',0,'weight',1,'numData',numData,'maxIter',iter,'oos',0);
[sol, ~]=ManifoldMatching(dis(1:n,1:n),dim,options);
[error5,error6]=EmbedAndClassify(sol(:,1:n)', Label, dim, trn,reps);

clear
load('BrainHippoShape.mat')
dis=[LMLS, (squareform(pdist(Label))>0)];
%dis=[LMLS, LMRS];
dim=3;trn=100;reps=100;K=4;
sol=SMDS(LMLS,dim,0);
sol2=zeros(dim,n);
for i=1:3
    sol2(:,Label==i)=mean(sol(Label==i),1);
end
dis=[LMLS, (squareform(pdist(sol2')))];
warning('off','all')
[error1,error2]=EmbedAndClassify(sol', Label, dim,trn,reps,1);
warning('on','all')
match=2;numData=2;iter=-1;nonlinear=0;
options = struct('nonlinear',nonlinear,'match',match,'neighborSize',K,'jointSelection',0,'weight',1,'numData',numData,'maxIter',iter,'oos',0);
[sol1, ~]=ManifoldMatching(dis,dim,options);
[error3,error4]=EmbedAndClassify(sol1(:,1:n)', Label, dim, trn,reps,1);
match=2;numData=1;iter=-1;nonlinear=1;
options = struct('nonlinear',nonlinear,'match',match,'neighborSize',K,'jointSelection',0,'weight',1,'numData',numData,'maxIter',iter,'oos',0);
[sol2, ~]=ManifoldMatching(dis(1:n,1:n),dim,options);
[error5,error6]=EmbedAndClassify(sol2(:,1:n)', Label, dim, trn,reps,1);