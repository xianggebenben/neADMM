function f=objective(x,t)
lambda=0.01;
 f=sum(abs(x))+lambda/2*sum(min(t*x,0).^2);
