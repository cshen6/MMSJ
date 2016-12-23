function plotVelocity(sol,numData)
% The plotVelocity function is an auxiliary function to plot pairs of dataset,
% in a 3D if dim is larger than 2 or 2D plots if dim is 2.
%
% The input sol should be of size dim*(numData*n).
% The first data is plotted by red circle, the second is plotted by blue square,
% and the third is plotted by green plus. The matched pairs are connected by black lines.

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

[dim,n]=size(sol);
n=n/numData;
X=zeros(dim,n,numData);
for i=1:numData
    tmp1=n*(i-1)+1:n*i;
    X(:,:,i)=sol(:,tmp1);
end

if(dim > 2)
    figure
    scatter3(X(1,:,1),X(2,:,1),X(3,:,1),20,'b','s');
    hold on
    if numData>1
        scatter3(X(1,:,2),X(2,:,2),X(3,:,2),20,'r','o');
        for i=1:n
            plot3([X(1,i,1),X(1,i,2)],[X(2,i,1),X(2,i,2)],[X(3,i,1),X(3,i,2)],'k')
        end
    end
    if numData>2
        scatter3(X(1,:,3),X(2,:,3),X(3,:,3),20,'g','+');
        for i=1:n
            plot3([X(1,i,2),X(1,i,3)],[X(2,i,2),X(2,i,3)],[X(3,i,2),X(3,i,3)],'k')
        end
    end
    %quiver3(X(1,:),X(2,:),X(3,:),differ(1,:),differ(2,:),differ(3,:)); %calculate the vector
    %xlabel('The first embedding dimension');
    %ylabel('The second embedding dimension');
    %zlabel('The third embedding dimension');
    hold off
else
    if(dim==2)
        figure
        scatter(X(1,:,1),X(2,:,1),20,'b','s');
        hold on
        if numData>1
            scatter(X(1,:,2),X(2,:,2),20,'r','o');
            for i=1:n
                plot([X(1,i,1),X(1,i,2)],[X(2,i,1),X(2,i,2)],'k')
            end
        end
        if numData>2
            scatter(X(1,:,3),X(2,:,3),20,'g','+');
            for i=1:n
                plot([X(1,i,2),X(1,i,3)],[X(2,i,2),X(2,i,3)],'k')
            end
        end
        %quiver(X(1,:),X(2,:),differ(1,:),differ(2,:));
        %xlabel('The first embedding dimension');
        %ylabel('The second embedding dimension');
        hold off
    end
end