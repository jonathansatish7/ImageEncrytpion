function corrc(a,Z)
%%%%Correlation coeff Original image horizontal direction
xi = a(:,1:end-1);  % original image
yi = Z(:,2:end);  % encrypted image
randIndex = randperm(numel(xi));                                   
randIndex = randIndex(1:2000);   
xRandi = xi(randIndex);            
yRandi = yi(randIndex); 
r_xyi = corr2(xRandi(:),yRandi(:));
disp(r_xyi);
% R=corrcoef(xRandi(:),yRandi(:));
% R_soh=R(2)^2;
figure;
scatter(xRandi(:),yRandi(:),'.');
%Correlation coeff encrypted image horizontal direction
x = Z(:,1:end-1);  
y = Z(:,2:end);  
randIndex = randperm(numel(x));                               
randIndex = randIndex(1:2000);   
xRand = x(randIndex);  
yRand = y(randIndex);
% R=corrcoef(xRand(:),yRand(:));
% R_seh=R(2)^2;
% disp(R_seh);
r_xy = corr2(xRand(:),yRand(:));
disp(r_xy);
figure;
scatter(xRand(:),yRand(:),'.');
%Correlation coeff Original image vertical direction
xiv = a(1:end-1,1:end-1,1);  
yiv = a(2:end,2:end,1); 
randIndex = randperm(numel(xiv));                                   
randIndex = randIndex(1:2000);   
xRandiv = xiv(randIndex);            
yRandiv = yiv(randIndex); 
% R=corrcoef(xRandiv(:),yRandiv(:));
% R_sov=R(2)^2;
% disp(R_sov);
r_xyiv =corr2(xRandiv(:),yRandiv(:));
disp(r_xyiv);
figure;
scatter(xRandiv(:),yRandiv(:),'.');
%Correlation coeff encrypted image vert direction
xv = Z(1:end-1,1:end-1,1);  
yv = Z(2:end,2:end,1); 
randIndex = randperm(numel(xv));                                   
randIndex = randIndex(1:2000);   
xRandv = xv(randIndex);            
yRandv = yv(randIndex); 
% R=corrcoef(xRandv(:),yRandv(:));
% R_sev=R(2)^2;
% disp(R_sev);
r_xyv = corr2(xRandv(:),yRandv(:));
disp(r_xyv);
figure;
scatter(xRandv(:),yRandv(:),'.');
%Correlation coeff original image diagonal direction
xid = a(2:end,1:end-1,1);  
yid = a(1:end-1,2:end,1);
randIndex = randperm(numel(xid));                                   
randIndex = randIndex(1:2000);   
xRandid = xid(randIndex);            
yRandid = yid(randIndex); 
% R=corrcoef(xRandid(:),yRandid(:));
% R_sod=R(2)^2;
% disp(R_sod);
r_xyid =corr2(xRandid(:),yRandid(:));
disp(r_xyid);
figure;
scatter(xRandid(:),yRandid(:),'.');

%Correlation coeff encrypted image diagonal direction
xd = Z(2:end,1:end-1,1);  
yd = Z(1:end-1,2:end,1); 
randIndex = randperm(numel(xd));                                   
randIndex = randIndex(1:2000);   
xRandd = xd(randIndex);            
yRandd = yd(randIndex); 
%  R=corrcoef(xRandd(:),yRandd(:));
% R_sed=R(2)^2;
% disp(R_sed);
r_xyd = corr2(xRandd(:),yRandd(:));
disp(r_xyd);
figure;
scatter(xRandd(:),yRandd(:),'.');