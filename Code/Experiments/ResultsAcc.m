function acc=ResultsAcc(solT,numData,tran,tesn)
acc=0;
sol=zeros(size(solT,1),tesn,numData);
for j=1:numData;
    sol(:,:,j)=solT(:,j*tran+(2*j-2)*tesn+1:j*tran+(2*j-1)*tesn);
end
for i=1:tesn
    for j=1:numData;
        flagg=true;
        for tt=1:numData
            if tt~=j
                dd=sol(:,:,tt)-repmat(sol(:,i,j),1,tesn);
                dd=sum(dd.^2);
                [~,ind]=sort(dd,'ascend');
                if (ind(1)~=i)
                    flagg=false;
                    break;
                end
            end
        end
        if (flagg)
            acc=acc+1/tesn/numData;
        end
    end
end