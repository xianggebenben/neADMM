function [r_history,s_history,beta,obj_history,rho] = nonsmoothSSM(xlInst,xuInst,xlIdx,xuIdx,yl,rho,nu,lambda,n,m,beta0)
%  threeStageAdmm Solve L1 regularized multi-instance and semi-supervised logistic regression via ADMM and FISTA
%
% % minimize  obj =sum(log(1+exp(epsilon)))-sum(yl.*epsilon(1:n))-nu*sum(max(epsilon((n+1):end),0))+lambda*norm(beta,1);
% where epsilon=max(beta*xlInst[xlIdx==i,:]') for user i,
% xlInst  and xuInst are all instance vectors for labeled and unlabeled users
% respectively.
% xlIdx and xuIdx are mapping from users to instances  for labeled
% and unlabeled users respectively.
% yl is the label set for labeled users.
% rho is the augmented Lagrangian parameter.
% nu is a tuning parameter
% lambda is a tuning parameter.
% n is the number of labeled users.
% m is the number of unlabeled users.
% beta0 is an initialization for beta(it could be omitted).
% return:
% r_history is a history record for L2 norm of primal residual r.
% s_history is a history record for L2 norm of dual residual s.
% beta is the objective parameter.
% obj is the corresponding objective function values for optimal beta.
% rho is the final augmented Lagrangian parameter.

%global constants and defaults
k = size(xlInst,2);
MAX_ITER = 100;
MAX_COUNT = 50;
TOLERANCE1   = 10e-3;
TOLERANCE2   = 10e-3;
%threeStageAdmm solver

%preprocessing instance sets; adding new instances with all zeros to users that
%do not have any messages.
newIndex =setdiff(1:n,unique(xlIdx));
xlInst =[xlInst;zeros(length(newIndex),k)];
xlIdx =[xlIdx,newIndex];
newIndex =setdiff(1:m,unique(xuIdx));
xuInst =[xuInst;zeros(length(newIndex),k)];
xuIdx =[xuIdx,newIndex];
xlInst =[ones(size(xlInst,1),1),xlInst];
xuInst =[ones(size(xuInst,1),1),xuInst];
idx =[xlIdx,xuIdx+n];
inst =[xlInst;xuInst];
[idx,index]=sort(idx,2,'ascend');
inst =inst(index,:);
clear newIndex xlInst xlIdx xuInst xuIdx index;
%initialization
if exist('beta0','var')
    beta = beta0;
else
    rng(100);
    beta = 2*rand(1,k+1)-1;
end
S =beta*inst';
epsilon=max_S(idx,S,n+m);
maxS =epsilon;
len =size(inst,1);
u1 = zeros(1,n+m);
u2 = zeros(1,len);
r = 0;
s = 0;
r_history =[];
s_history =[];
mu =4;
rMin =9999;
sMin=9999;
count =0;
user =n+m;
for iter = 1:MAX_ITER
        %rho-update
%              if(r>mu*s)
%                  rho =rho*2;
%              else if (mu*r<s)
%                      rho = rho/2;
%                  end
%              end
    %rho =rho +2/MAX_ITER;
    % epsilon-update
    epsilon = update_epsilon(yl,maxS,u1,rho,epsilon,nu);
    % beta-update
    [beta,FitInfo] = lasso(inst(:,2:end),S+u2, 'Alpha',1,'Lambda',lambda/(rho*length(S)));
    beta=[FitInfo.Intercept,beta'];
    %beta=update_beta(S,u2,inst,beta,lambda,rho);
    % S-update
    S_old = S;
    S=update_S(inst,idx,epsilon,beta,u1,u2,user);
    maxS_old = maxS;
    maxS=max_S(idx,S,n+m);
    % compute prime and dual residuals
    r1 = epsilon-maxS;
    r2 = S -beta*inst';
    s1 =rho*(maxS_old-maxS);
    s2 =rho*(S-S_old)*inst;
    u1 = u1 + r1;
    u2 = u2 + r2;
    r =norm([r1,r2]);
    r_history(iter)=r;
    s =norm([s1,s2]);
    s_history(iter)=s;
    temp_S=beta*inst';
    temp_epsilon=max_S(idx,temp_S,n+m);
    obj =objective(temp_epsilon,yl,lambda,beta,n,nu);
    obj_history(iter)=obj;
    % record the smallest r and s
        if ( r<rMin)
            rMin=r;
            count =0;
        else if (s<sMin)
                sMin =s;
                count =0;
            else
                count = count +1;
            end
        end
    %  termination checks
    if (r<=TOLERANCE1 &&  s<=TOLERANCE2) || count ==MAX_COUNT
        break;
    end
end
end

function obj = objective(epsilon,yl,lambda,beta,n,nu)
obj =sum(log(1+exp(epsilon)))-epsilon(1:n)*yl-nu*sum(max(epsilon((n+1):end),0))+lambda*norm(beta(2:end),1);
end
function epsilon = update_epsilon(yl,maxS,u1,rho,epsilon0,nu)
% solve the epsilon update
% via FISTA
% global constants and defaults
n = length(yl);
yepsilon = 10e10;
MAX_ITER = 500;
TOLERANCE =10e-5;
epsilon = epsilon0;
y = @(epsilon) (sum(log(1 + exp(epsilon)))-epsilon(1:n)*yl-nu*sum(max(epsilon(n+1:end),0)) + (rho/2)*norm(epsilon - maxS + u1,2)^2);
lambda =1;
zeta =epsilon;
eta = 0.5;
% FISTA
for iter = 1:MAX_ITER
    yepsilon_old =yepsilon;
    yepsilon = y(epsilon);
    if(abs(yepsilon-yepsilon_old)<TOLERANCE)
        break;
    end
    lambda_old =lambda;
    lambda =(1+sqrt(1+4*lambda^2))/2;
    gamma =(1-lambda_old)/lambda;
    gradient=exp(epsilon)./(exp(epsilon)+1);
    gradient(1:n)=gradient(1:n)-yl';
    gradient(n+1:end)=gradient(n+1:end)+rho*(epsilon(n+1:end)-maxS(n+1:end)+u1(n+1:end));
    zeta_old =zeta;
    zeta(1:n) = (rho*(maxS(1:n)-u1(1:n))+(epsilon(1:n)-eta*gradient(1:n))/eta)/(rho+1/eta);
    % unlabeled users
    for i =(n+1):length(epsilon0)
        v = @(zeta)(-nu*max(zeta,0)+1/(2*eta)*(zeta-epsilon(i)+eta*gradient(i))^2);
        candidate =[0,epsilon(i)-eta*gradient(i),epsilon(i)-eta*gradient(i)+nu*eta];
        [~,pos]=min([v(candidate(1)),v(candidate(2)),v(candidate(3))]);
        zeta(i)=candidate(pos);
    end
    epsilon =(1-gamma)*zeta+gamma*zeta_old;
end
end

function S=update_S(inst,idx,epsilon,beta,u1,u2,user)
% solve S-update
S=beta*inst'-u2;
for i=1:user
    index =find(idx==i);
    l=length(index);
    phi =S(index);
    temp =epsilon(i)+u1(i);
    [phi_sort,seq] =sort(phi,'descend');
    avg =zeros(1,l);
    avg(1)=(phi_sort(1)+temp)/2;
    for j=2:l
        avg(j) =(avg(j-1)*j+phi_sort(j))/(j+1);
    end
    % find the largest suitable j*
    pos =find(avg(1:end-1)>phi_sort(2:end), 1, 'first' );
    if isempty(pos)
        pos=l;
    end
    S(index(seq(1:pos(1))))=avg(pos(1));
end
end
function maxS =max_S(idx,S,n)
% obtain the maximal instance for each user
maxS =zeros(1,n);
for i=1:n
    index = idx==i;
    maxS(i) = max(S(index));
end
end