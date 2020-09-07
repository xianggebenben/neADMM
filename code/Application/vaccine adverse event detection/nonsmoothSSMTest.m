function [CM,AUC,ROCX,ROCY,AUPR,PRX,PRY]=nonsmoothSSMTest(beta,test,testInstanceIndex,testLabel)
l = length(testLabel);
k = size(test,2);
test =[ones(size(test,1),1),test];
S=test*beta';
maxS=zeros(l,1);
maxInstance=zeros(l,k+1);
for i=1:l
    index = find(testInstanceIndex==i);
    if ~isempty(index)
        [maxS(i),pos]=max(S(index));
        maxInstance(i,:) =test(index(pos),:);
    else
        maxS(i)=beta(1);
        maxInstance(i,:) =[1,zeros(1,k)];
    end
end
CM=confusionmat(testLabel,(maxS>=0)*1);
[ROCX,ROCY,~,AUC]=perfcurve(testLabel,1./(1+exp(-maxS)),1);
[PRX,PRY,~,AUPR]=perfcurve(testLabel,1./(1+exp(-maxS)),1,'xCrit', 'reca', 'yCrit', 'prec');
end