clear;
MAX_ITER = 30;
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
    t1=-(2*(sqrt(z)-1)+y)/(2*(rho+1));
    if(t1>=0)
        x=t1^2;
    else
        x=0;
    end
    %Z-update
    z_old =z;
    t2=-(2*(sqrt(x)-1)+y)/(2*(rho+1));
    if(t2>=0)
        z=t2^2;
    else
        z=0;
    end
    r=sqrt(x)+sqrt(z)-1;
    s=1/2*rho*(x)^(-1/2)*(sqrt(z)-sqrt(z_old));
    p=x+z;
    y=y+rho*r;
    r_history(i)=abs(r);
    s_history(i)=abs(s);
    p_history(i)=p;
end
figure(1);
plot(1:i,r_history,'r');
hold on;
plot(1:i,s_history,'b*');
hold off;
legend('r','s');
xlabel('Iteration');
ylabel('Residual');
figure(2);
plot(1:i,p_history);
xlabel('Iteration');
ylabel('Objective Value');