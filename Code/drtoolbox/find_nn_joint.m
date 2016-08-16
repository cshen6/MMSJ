function [D,ind] = find_nn_joint(X, Y, k)

D1=squareform(pdist(X));
D2=squareform(pdist(Y));
%D2=D2/norm(D2,'fro')*norm(D1,'fro');
D=D1+D2;

[~,ind]=sort(D);

ind=ind(:,1:k);
