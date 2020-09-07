clear;
MAX_ITER = 1000;
x=0;
y=0;
z=0;
r=0;
s=0;
rho =1;
tolerance=10e-10;
mu=2;
for i = 1:MAX_ITER
    %rho-update
%         if(abs(r)>mu*abs(s))
%             rho=rho*2;
%         else if (abs(s)>mu*abs(r))
%                 rho=rho/2;
%             end
%         end
    rho =rho +0.01;
    %X-update
    p=z^2-1+y/rho;
    q=1/(2*rho);
    delta =(p/3)^3+(q/2)^2;
    if delta>0
    A=nthroot(-q/2+sqrt(delta),3);
    B=nthroot(-q/2-sqrt(delta),3);
        x=A+B;
    else
    A=(-q/2+delta^(1/2))^(1/3);
    B=(-q/2-delta^(1/2))^(1/3);
        % pay attention -sqrt(3)/2i=i*sqrt(3)/2;
        omega=-1/2-sqrt(3)/2i;
        x1 = A + B;
        x2 = real(omega*A+omega^2*B);
        x3 = real(omega^2*A+omega*B);
        candi=[x1,x2,x3];
        value =candi+y*candi.^2+rho/2*(candi.^2+z^2-1).^2;
        idx =find(value==min(value),1,'first');
        x=candi(idx);
    end
    %Z-update
    z_old =z;
    p=x^2-1+y/rho;
    q=1/(2*rho);
    delta =(p/3)^3+(q/2)^2;
    if delta>0
    A=nthroot(-q/2+sqrt(delta),3);
    B=nthroot(-q/2-sqrt(delta),3);
        z=A+B;
    else
    A=(-q/2+delta^(1/2))^(1/3);
    B=(-q/2-delta^(1/2))^(1/3);
    % pay attention -sqrt(3)/2i=i*sqrt(3)/2;
     omega=-1/2-sqrt(3)/2i;
        z1 = A + B;
        z2 = real(omega*A+omega^2*B);
        z3 = real(omega^2*A+omega*B);
        candi=[z1,z2,z3];
        value =candi+y*candi.^2+rho/2*(candi.^2+x^2-1).^2;
        idx =find(value==min(value),1,'first');
        z=candi(idx);
    end
    r=x^2+z^2-1;
    s=1/2*rho*(x)^(-2)*(z^2-z_old^2);
    p=x+z;
    y=y+rho*r;
    r_history(i)=abs(r);
    s_history(i)=abs(s);
    p_history(i)=p;
end
figure;
plot(1:i,r_history,'r');
hold on;
plot(1:i,s_history,'b*');
hold off;
legend('r','s');
xlabel('Iteration');
ylabel('Residual');
figure;
plot(1:i,p_history,'b');
xlabel('Iteration');
ylabel('Objective Value');