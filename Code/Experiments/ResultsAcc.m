function [accIsoJ,acc,accIsoS,accLLES,accLTSA]=ResultsAcc(solIsoJ,sol,solIsoS,solLLES,solLTSA,numData,tran,tesn,reps)

% [accIsoJ,acc,accIsoS,accLLES,accLTSA]=ResultsAcc(solIsoJ,sol,solIsoS,solLLES,solLTSA,numData,tran,tesn,reps);

accIsoJ=accFind(solIsoJ,numData,tran,tesn,reps);
acc=accFind(sol,numData,tran,tesn,reps);
accIsoS=accFind(solIsoS,numData,tran,tesn,reps);
accLLES=accFind(solLLES,numData,tran,tesn,reps);
accLTSA=accFind(solLTSA,numData,tran,tesn,reps);

function acc=accFind(solT,numData,tran,tesn,reps)
acc=0;
sol=zeros(size(solT,1),tesn,numData);
for rr=1:reps
    for j=1:numData;
        sol(:,:,j)=solT(:,j*tran+(2*j-2)*tesn+1:j*tran+(2*j-1)*tesn,2,rr);
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
                acc=acc+1/tesn/reps/numData;
            end           
        end
    end
end