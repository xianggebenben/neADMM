clear;
obj=[];
for m=1:4
for n=1:10
    m
    n
[~,~,~,~,obj(m,n)]=compressive_sensing_ADMM(512,50*n,16*m,0.01);
end
end