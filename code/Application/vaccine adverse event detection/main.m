clear;
load data.mat;
rho =1;
lambda = 1;
nfold=5;
nu=0;
rng('default');
indices = crossvalind('Kfold', length(label), nfold);
beta=zeros(1,size(data,2)+1);
m=0;
    testIdx = find(indices == 1);
    trainIdx = find(indices ~= 1);
    train=[];
    trainInstanceIndex=[];
    for j=1:length(trainIdx)
        index=find(InstanceIndex==trainIdx(j));
        train=[train;data(index,:)];
        trainInstanceIndex=[trainInstanceIndex,j*ones(1,length(index))];
    end
    test=[];
    testInstanceIndex=[];
    for j=1:length(testIdx)
        index=find(InstanceIndex==testIdx(j));
        test=[test;data(index,:)];
        testInstanceIndex=[testInstanceIndex,j*ones(1,length(index))];
    end
    n=length(trainIdx);
[r_history,s_history,beta1,obj_history1,~] =nonsmoothSSM(train,[],trainInstanceIndex,[],label(trainIdx),rho,nu,lambda,n,m,beta);
[beta2,obj_history2]=selectorVariable(train,trainInstanceIndex,label(trainIdx),lambda,n,beta);
[CM1,AUC1,ROCX1,ROCY1,AUPR1,PRX1,PRY1]=nonsmoothSSMTest(beta1,test,testInstanceIndex,label(testIdx));
acc1=(CM1(1,1)+CM1(2,2))/sum(sum(CM1));
pr1=CM1(2,2)/sum(CM1(:,2));
re1=CM1(2,2)/sum(CM1(2,:));
fs1=2*pr1*re1/(pr1+re1);
[CM2,AUC2,ROCX2,ROCY2,AUPR2,PRX2,PRY2]=nonsmoothSSMTest(beta2,test,testInstanceIndex,label(testIdx));
acc2=(CM2(1,1)+CM2(2,2))/sum(sum(CM2));
pr2=CM2(2,2)/sum(CM2(:,2));
re2=CM2(2,2)/sum(CM2(2,:));
fs2=2*pr2*re2/(pr2+re2);
save('result.mat');