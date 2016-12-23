%The auxiliary function to calculate the sample distance correlation.
function corr = distCorr(X,Y)
X=squareform(pdist(X'));
Y=squareform(pdist(Y'));
n=size(X,1);
H=eye(n)-(1/n)*ones(n,n);
X=H*X*H;
Y=H*Y*H;

corr=sum(sum(X.*(Y)));
corr=sqrt(corr/(norm(X,'fro')*norm(Y,'fro')));