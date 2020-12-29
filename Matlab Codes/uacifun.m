function uaci=uacifun(x,y)
[m,n,p]=size(x);
z=abs(double(x)-double(y));
uaci=sum(sum(sum(z)))/(255*m*n*p)*100;
end