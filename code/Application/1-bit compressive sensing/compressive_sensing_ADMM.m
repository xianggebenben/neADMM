function [r_history,s_history,x_hat,x,obj]=compressive_sensing_ADMM(N,M,K,lambda)
x=zeros(N,1);
rng(100);
%% nonzero index
index =randi(N,1,K);
rng(100);
R =normrnd(0,1,1,K);
x(index)=R;
rng(100);
phi =normrnd(0,1,M,N);
y=phi*x;
Y=diag(y);
 t=Y*phi;
 y1=0;
 y2=zeros(M,1);
 rng(0);
 x_hat=normrnd(0,1,N,1);
 x_hat=x_hat/norm(x_hat);
 w=x_hat;
 y3=zeros(N,1);
 rho=1;
 iter_max=100;
 tolerance=10e-5;
for i=1:iter_max
    i
    w_old=w;
    % update x_hat
    x_hat=update_x(rho,w,y1,y3);
    % update z
    z=update_z(lambda,y2,rho,w,t);
    % update w
    w=update_w(w,x_hat,rho,y3,t,z,y2);
    r1=norm(x_hat)^2-1;
    r2=t*w-z;
    r3=w-x_hat;
    r=norm([r1;r2;r3])
    s1=rho*(w_old-w);
    s2=rho*t*(w_old-w);
    s=norm([s1;s2])
    y1=y1+rho*r1;
    y2=y2+rho*r2;
    y3=y3+rho*r3;
    r_history(i)=r;
    s_history(i)=s;
    if(r<tolerance && s<tolerance)
        break;
    end
end
x_hat=x_hat/norm(x_hat);
    obj=object(lambda,t,x_hat);
end
function obj=object(lambda,t,x)
 obj=sum(abs(x))+lambda/2*sum(min(t*x,0).^2);
end

function z=update_z(lambda,y2,rho,w,t)
z=(rho*t*w+y2)/rho;
temp=(rho*t*w+y2)/(rho+lambda);
z(z<0)=temp(z<0);
end
function w=update_w(w0,x,rho,y3,t,z,y2)
% solve the w update
% via FISTA
% global constants and defaults
yw = 10e10;
MAX_ITER = 1000;
TOLERANCE =10e-6;
w = w0;
y = @(w) rho/2*norm(w-x+y3/rho)^2+rho/2*norm(t*w-z+y2/rho)^2+norm(w,1);
lambda =1;
zeta =w;
eta = 10e-10;
% FISTA
for iter = 1:MAX_ITER
    yw_old =yw;
    yw = y(w);
    if(abs(yw-yw_old)<TOLERANCE)
        break;
    end
    lambda_old =lambda;
    lambda =(1+sqrt(1+4*lambda^2))/2;
    gamma =(1-lambda_old)/lambda;
    gradient=rho*t'*(t*w-z+y2/rho)+rho*(w-x+y3/rho);
    zeta_old =zeta;
    zeta0=zeta-eta*gradient;
    zeta = max(abs(zeta0)-eta,0).*sign(zeta0);
    w =(1-gamma)*zeta+gamma*zeta_old;
end
end
function x=update_x(rho,w,y1,y3)
% u=\Vert x\Vert
% if u^2+1+y_1/rho>0
    p=1+y1/rho;
    q=-norm(2*(w+y3/rho),2);
    delta =(p/3)^3+(q/2)^2;
    if delta>0
    A=nthroot(-q/2+sqrt(delta),3);
    B=nthroot(-q/2-sqrt(delta),3);
        u=A+B;
    else
    A=(-q/2+delta^(1/2))^(1/3);
    B=(-q/2-delta^(1/2))^(1/3);
        % pay attention -sqrt(3)/2i=i*sqrt(3)/2;
        omega=-1/2-sqrt(3)/2i;
        u1 = A + B;
        u2 = real(omega*A+omega^2*B);
        u3 = real(omega^2*A+omega*B);
        u=[u1,u2,u3];
    end
    u=u(u.^2+1+y1/rho>0);
    x1=[];
    if ~isempty(u)
      for i=1:length(u)
    x1(:,i)=2*(w+y3/rho)/(u^2+1+y1/rho);
      end
    end
 % else u^2+1+y_1/rho<0
    p=1+y1/rho;
    q=norm(2*(w+y3/rho),2);
    delta =(p/3)^3+(q/2)^2;
    if delta>0
    A=nthroot(-q/2+sqrt(delta),3);
    B=nthroot(-q/2-sqrt(delta),3);
        u=A+B;
    else
    A=(-q/2+delta^(1/2))^(1/3);
    B=(-q/2-delta^(1/2))^(1/3);
        % pay attention -sqrt(3)/2i=i*sqrt(3)/2;
        omega=-1/2-sqrt(3)/2i;
        u1 = A + B;
        u2 = real(omega*A+omega^2*B);
        u3 = real(omega^2*A+omega*B);
        u=[u1,u2,u3];
    end
    u=u(u.^2+1+y1/rho<0);
    x2=[];
        if(~isempty(u))
    for i=1:length(u)        
    x2(:,i)=2*(w+y3/rho)./(u(i)^2+1+y1/rho);
    end
        end
    x=[x1,x2];
    f=@(x) (norm(x,2)^2-1+y1/rho)^2+norm(w-x+y3/rho,2)^2;
    for i=1:size(x,2)
        fval(i)=f(x(:,i));
    end
    [~,idx]=min(fval);
    x=x(:,idx);
end