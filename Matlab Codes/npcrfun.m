function npcr1=npcrfun(x,y,m,n)
z=x==y;
cap=double(z(z==0));
npcr1=size(cap,1)/(n*m)*100;
end