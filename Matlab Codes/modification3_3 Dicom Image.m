clc 
clear all
close all
DICOM_IMAGE=dicomread('image_OT.dcm');%Loading of the DICOM image
DICOM_IMAGE=imresize(DICOM_IMAGE,[256 256]);
figure(1)
imshow(DICOM_IMAGE,[]);
title('DICOM image');
[row,col]=size(DICOM_IMAGE);
u=3.989;n=(row*col)-1;           %Logistic Map initial conditions       
x0=0.499;
x(1)=x0;        %Chaotic sequence generarion through Logistic Map
for i=1:n
    x(i+1)=u*x(i)*(1-x(i));
end
a = reshape(x,[row,col]);
h =[mod(a* 10^17,256)+1];

%%%%%%%%% converting to binary 8 bit %%%%%%%%%

for i=1:1:row
    for j=1:1:col
      finalOutput{i,j} = dec2bin(h(i,j),8);
    end
end 

%%%%%%%%%%converting acc to rule 3 %%%%%%%%%

codebook1 = containers.Map({'00','11','01','10'},{'A','T','C','G'}); %// Lookup
C1 = cellfun(@(x) values(codebook1, {x(1:2),x(3:4),x(5:6),x(7:8)}), ...
               finalOutput, 'uni', 0);
A = cellfun(@cell2mat, C1, 'uni', 0);

%%%%%%%%% converting to binary 8 bit %%%%%%%%%
for i=1:1:row
    for j=1:1:col
      finalOutput1{i,j} = dec2bin(DICOM_IMAGE(i,j),8);
    end
end 

%%%%%%%%%%converting acc to rule 3 %%%%%%%%%

codebook2 = containers.Map({'00','11','01','10'},{'A','T','C','G'}); %// Lookup
C2 = cellfun(@(x) values(codebook2, {x(1:2),x(3:4),x(5:6),x(7:8)}), ...
               finalOutput1, 'uni', 0);
B = cellfun(@cell2mat, C2, 'uni', 0);



%%%%%%% DNA xor operation %%%%%%%%%%%%

for i=1:row%%%%%% row Size
    for j=1:col%%%%%%%%%%% column size
        
         img_a = A{i,j};
         key_b= B{i,j};
         
        for i__0 = 1:4     %%%%%%% no.of encoded bits per pixel since its medical image we have 16 bits which are encoded to 8 bits
            
    if (img_a(i__0)=='A')&&(key_b(i__0)=='A')
        dna_a(i__0)='A';
    else if (img_a(i__0)=='A')&&(key_b(i__0)=='C')
            dna_a(i__0)='C';
        else if (img_a(i__0)=='A')&&(key_b(i__0)=='G')
                dna_a(i__0)='G';
            else if (img_a(i__0)=='A')&&(key_b(i__0)=='T')
                    dna_a(i__0)='T';
                    
                else if (img_a(i__0)=='C')&&(key_b(i__0)=='A')
                        dna_a(i__0)='C';
                        else if (img_a(i__0)=='C')&&(key_b(i__0)=='C')
                                dna_a(i__0)='A';
                            else if (img_a(i__0)=='C')&&(key_b(i__0)=='G')
                                    dna_a(i__0)='T';
                                else if (img_a(i__0)=='C')&&(key_b(i__0)=='T')
                                        dna_a(i__0)='G';
                                        
                                        else if (img_a(i__0)=='G')&&(key_b(i__0)=='A')
                                                dna_a(i__0)='G';
                                                else if (img_a(i__0)=='G')&&(key_b(i__0)=='C')
                                                        dna_a(i__0)='T';
                                                        else if (img_a(i__0)=='G')&&(key_b(i__0)=='G')
                                                                dna_a(i__0)='A';
                                                                else if (img_a(i__0)=='G')&&(key_b(i__0)=='T')
                                                                        dna_a(i__0)='C';
                                                                        
                                                                        else if (img_a(i__0)=='T')&&(key_b(i__0)=='A')
                                                                                dna_a(i__0)='T';
                                                                                else if (img_a(i__0)=='T')&&(key_b(i__0)=='C')
                                                                                        dna_a(i__0)='G';
                                                                                        else if (img_a(i__0)=='T')&&(key_b(i__0)=='G')
                                                                                                dna_a(i__0)='C';
                                                                                                else  (img_a(i__0)=='T')&&(key_b(i__0)=='T')
                                                                                                        dna_a(i__0)='A';
                                                                                            end
                                                                                    end
                                                                            end
                                                                    end
                                                            end
                                                    end
                                            end
                                    end
                                end
                            end
                    end
                end
            end
        end
    end
        end
        dna_add{i,j}=dna_a;   %%%%%%%%%% Forming the matrix again
    end

end

%%%%%%%%%% converting to binary and then to decimal to disp image rule 8 %%%%%%%%

codebook3= containers.Map({'A','T','C','G'},{'11','00','10','01'}); %// Lookup
C3 = cellfun(@(x) values(codebook3, {x(1),x(2),x(3),x(4)}), ...
             dna_add, 'uni', 0);
B11= cellfun(@cell2mat, C3, 'uni', 0);
% disp(B);

for i=1:1:row
    for j=1:1:col
        encrypted(i,j,1) = bin2dec(B11(i,j));
    end
end
figure(2);
subplot(1,1,1);imshow(encrypted,[]);
%  imwrite(encrypted,'D:\dna_fingerprint_enc.bmp')
x1(1)=0.567;
r11=3.987;
for k=1:n
    x1(k+1)=r11*x1(k)*(1-x1(k));
end
h1 =[mod(x1* 10^17,256)+1];
a_1= reshape(h1,[row,col]);
a_1=uint8(a_1);

for i=1:1:row
    for j=1:1:col
        shuffled_im(i,j) = bitxor(encrypted(i,j),a_1(i,j));
    end
end
figure(3)
imshow(shuffled_im);
title('Shuffled image after logistic map');
shuffled_im_1D=shuffled_im(:);
yy=zeros(1,row*col);
yy(1)=shuffled_im_1D(1);
for i=2:row*col
    yy(i)=bitxor(shuffled_im_1D(i-1),shuffled_im_1D(i));  %XOR encryption
end
shuffled_im_2D=reshape(yy,row,col);    %Wavelet decomposition of the image 
[row_1,col_1]=size(shuffled_im_2D);
a1=double(shuffled_im_2D);
lshaar = liftwave('haar','int2int');
els = {'p',[-0.125 0.125],0};
lsnew = addlift(lshaar,els);
[cAint,cHint,cVint,cDint] = lwt2(a1,lsnew);
figure(4)
subplot(2,2,1)
imshow(cAint,[]);
title('LL part');
subplot(2,2,2)
imshow(cHint,[]);
title('LH part');
subplot(2,2,3)
imshow(cVint,[]);
title('HL part');
subplot(2,2,4)
imshow(cDint,[]);
title('HH part');
[cA1int,cH1int,cV1int,cD1int]=lwt2(cAint,lsnew);  %Second level wavelet decompostion 
[cA2int,cH2int,cV2int,cD2int]=lwt2(cHint,lsnew);
[cA3int,cH3int,cV3int,cD3int]=lwt2(cVint,lsnew);
[cA4int,cH4int,cV4int,cD4int]=lwt2(cDint,lsnew);
[row_2,col_2]=size(cA1int);
u=3.925;n=row_2*col_2;                      %Scrambling process done to the approximation values
x0_1=0.1;
x_1(1)=u*x0_1*(1-x0_1);
for i=2:n
    x_1(i)=u*x_1(i-1)*(1-x_1(i-1));
end
k=1;
asc_x1=sort(x_1);
N1_1=length(x_1);
for i=1:N1_1
    for j=1:N1_1
        if(asc_x1(i)==x_1(j))
            new_mat_1(1,k)=j;
            k=k+1;
        end
    end
end
cA1int_1D=(cA1int(:))';
cA2int_1D=(cA2int(:))';
cA3int_1D=(cA3int(:))';
cA4int_1D=(cA4int(:))';
for i=1:N1_1
    cA1int_shuffle(i)=cA1int_1D(new_mat_1(i));
    cA2int_shuffle(i)=cA2int_1D(new_mat_1(i));
    cA3int_shuffle(i)=cA3int_1D(new_mat_1(i));
    cA4int_shuffle(i)=cA4int_1D(new_mat_1(i));
end
len=1;
for i=1:row_2
    for j=1:col_2
        encrypted_matrix1(i,j) = cA1int_shuffle(1,len);
        encrypted_matrix2(i,j) = cA2int_shuffle(1,len);
        encrypted_matrix3(i,j) = cA3int_shuffle(1,len);
        encrypted_matrix4(i,j) = cA4int_shuffle(1,len);
        len=len+1;
    end
end
figure(5)
subplot(2,2,1)
imshow(encrypted_matrix1,[]);
title('shuffled LL1');
subplot(2,2,2)
imshow(encrypted_matrix2,[]);
title('shuffled LL2');
subplot(2,2,3)
imshow(encrypted_matrix3,[]);
title('shuffled LL3');
subplot(2,2,4)
imshow(encrypted_matrix4,[]);
title('shuffled LL4');
encrypted_matrix1=double(encrypted_matrix1);
encrypted_matrix2=double(encrypted_matrix2);
encrypted_matrix3=double(encrypted_matrix3);
encrypted_matrix4=double(encrypted_matrix4);
inverse1=ilwt2(encrypted_matrix1,cH1int,cV1int,cD1int,lsnew);
inverse2=ilwt2(encrypted_matrix2,cH2int,cV2int,cD2int,lsnew);
inverse3=ilwt2(encrypted_matrix3,cH3int,cV3int,cD3int,lsnew);
inverse4=ilwt2(encrypted_matrix4,cH4int,cV4int,cD4int,lsnew);
inverse_im=ilwt2(inverse1,inverse2,inverse3,inverse4,lsnew);   %The encrypted image 
figure(6)
imshow(inverse_im,[]);
title('Encrypted Image');
figure(7)
histogram(inverse_im);
title('Histogram of encrypted image');
%Decryption
[inverse1_1,inverse2_1,inverse3_1,inverse4_1]=lwt2(inverse_im,lsnew); %Wavelet decomposition during Decryption
[encrypted_matrix1_1,cH1int_1,cV1int_1,cD1int_1]=lwt2(inverse1_1,lsnew);
[encrypted_matrix2_1,cH2int_1,cV2int_1,cD2int_1]=lwt2(inverse2_1,lsnew);
[encrypted_matrix3_1,cH3int_1,cV3int_1,cD3int_1]=lwt2(inverse3_1,lsnew);
[encrypted_matrix4_1,cH4int_1,cV4int_1,cD4int_1]=lwt2(inverse4_1,lsnew);
[row_3,col_3]=size(encrypted_matrix1_1);
len=1;
for i=1:row_3
    for j=1:col_3
        encrypted_matrix1_2(1,len)=encrypted_matrix1_1(i,j);  
        encrypted_matrix2_2(1,len)=encrypted_matrix2_1(i,j);
        encrypted_matrix3_2(1,len)=encrypted_matrix3_1(i,j);
        encrypted_matrix4_2(1,len)=encrypted_matrix4_1(i,j);
        len=len+1;
    end
end
for i=1:N1_1
    encrypted_matrix1_3(new_mat_1(i))=encrypted_matrix1_2(i);
    encrypted_matrix2_3(new_mat_1(i))=encrypted_matrix2_2(i);
    encrypted_matrix3_3(new_mat_1(i))=encrypted_matrix3_2(i);
    encrypted_matrix4_3(new_mat_1(i))=encrypted_matrix4_2(i);
end
cA1int_3=reshape(encrypted_matrix1_3,row_3,col_3);  
cA2int_3=reshape(encrypted_matrix2_3,row_3,col_3);
cA3int_3=reshape(encrypted_matrix3_3,row_3,col_3);
cA4int_3=reshape(encrypted_matrix4_3,row_3,col_3);
cAint_1=ilwt2(cA1int_3,cH1int,cV1int,cD1int,lsnew);
cHint_1=ilwt2(cA2int_3,cH2int,cV2int,cD2int,lsnew);
cVint_1=ilwt2(cA3int_3,cH3int,cV3int,cD3int,lsnew);
cDint_1=ilwt2(cA4int_3,cH4int,cV4int,cD4int,lsnew);
decrypted_im=ilwt2(cAint_1,cHint_1,cVint_1,cDint_1,lsnew);
decrypted_im=uint8(decrypted_im);
decrypted_im=double(decrypted_im);
for i=1:row
    for j=1:col
        c1(i,j)=a1(i,j)-decrypted_im(i,j);
    end
end
for i=1:row
    for j=1:col
        decrypted_im1(i,j)=decrypted_im(i,j)+c1(i,j);
    end
end
decrypted_im_1=decrypted_im1(:);
zz=zeros(1,row*col);
zz(1)=decrypted_im_1(1);
for i=2:row*col
    zz(i)=bitxor(zz(i-1),decrypted_im_1(i));
end
z=reshape(zz,row,col);

for i=1:1:row
    for j=1:1:col
        encrypted1(i,j) = bitxor(z(i,j),a_1(i,j));
    end
end

%%%%%%%%% converting to binary 8 bit %%%%%%%%%
for i=1:1:row
    for j=1:1:col
      encr_bin{i,j} = dec2bin(encrypted1(i,j),8);
    end
end 

%%%%%%%%%%converting acc to rule 1 %%%%%%%%%

codebook4 = containers.Map({'11','00','10','01'},{'A','T','C','G'}); %// Lookup
C_2 = cellfun(@(x) values(codebook4, {x(1:2),x(3:4),x(5:6),x(7:8)}), ...
              encr_bin, 'uni', 0);
B_1 = cellfun(@cell2mat, C_2, 'uni', 0);

n=(row*col)-1;
x=zeros(1,n);
r1=3.989;
x0=0.499;
x(1)=x0;
for k=1:n
    x(k+1)=r1*x(k)*(1-x(k));
end
a = reshape(x,[row,col]);
h =[mod(a* 10^17,256)+1];

%%%%%%%% converting to binary 8 bit %%%%%%%%%

for i=1:1:row
    for j=1:1:col
      finalOutput_2{i,j} = dec2bin(h(i,j),8);
    end
end 


%%%%%%%%%converting acc to rule 8 %%%%%%%%%

codebook5 = containers.Map({'00','11','01','10'},{'A','T','C','G'}); %// Lookup
C_1 = cellfun(@(x) values(codebook5, {x(1:2),x(3:4),x(5:6),x(7:8)}), ...
               finalOutput_2, 'uni', 0);
A_1 = cellfun(@cell2mat, C_1, 'uni', 0);

for i=1:row%%%%%% row Size
    for j=1:col%%%%%%%%%%% column size
        
         img_a = B_1{i,j};
         key_b= A_1{i,j};
         
        for i__0 = 1:4     %%%%%%% no.of encoded bits per pixel since its medical image we have 16 bits which are encoded to 8 bits
       
            
            
   
    if (img_a(i__0)=='A')&&(key_b(i__0)=='A')
        dna_a(i__0)='A';
    else if (img_a(i__0)=='A')&&(key_b(i__0)=='C')
            dna_a(i__0)='C';
        else if (img_a(i__0)=='A')&&(key_b(i__0)=='G')
                dna_a(i__0)='G';
            else if (img_a(i__0)=='A')&&(key_b(i__0)=='T')
                    dna_a(i__0)='T';
                    
                else if (img_a(i__0)=='C')&&(key_b(i__0)=='A')
                        dna_a(i__0)='C';
                        else if (img_a(i__0)=='C')&&(key_b(i__0)=='C')
                                dna_a(i__0)='A';
                            else if (img_a(i__0)=='C')&&(key_b(i__0)=='G')
                                    dna_a(i__0)='T';
                                else if (img_a(i__0)=='C')&&(key_b(i__0)=='T')
                                        dna_a(i__0)='G';
                                        
                                        else if (img_a(i__0)=='G')&&(key_b(i__0)=='A')
                                                dna_a(i__0)='G';
                                                else if (img_a(i__0)=='G')&&(key_b(i__0)=='C')
                                                        dna_a(i__0)='T';
                                                        else if (img_a(i__0)=='G')&&(key_b(i__0)=='G')
                                                                dna_a(i__0)='A';
                                                                else if (img_a(i__0)=='G')&&(key_b(i__0)=='T')
                                                                        dna_a(i__0)='C';
                                                                        
                                                                        else if (img_a(i__0)=='T')&&(key_b(i__0)=='A')
                                                                                dna_a(i__0)='T';
                                                                                else if (img_a(i__0)=='T')&&(key_b(i__0)=='C')
                                                                                        dna_a(i__0)='G';
                                                                                        else if (img_a(i__0)=='T')&&(key_b(i__0)=='G')
                                                                                                dna_a(i__0)='C';
                                                                                                else  (img_a(i__0)=='T')&&(key_b(i__0)=='T')
                                                                                                        dna_a(i__0)='A';
                                                                                            end
                                                                                    end
                                                                            end
                                                                    end
                                                            end
                                                    end
                                            end
                                    end
                                end
                            end
                    end
                end
            end
        end
    end
        end
        dna_add_1{i,j}=dna_a;   %%%%%%%%%% Forming the matrix again
    end

end

%%%%%%%%%% converting to binary and then to decimal to disp image%%%%%%%%

codebook6= containers.Map({'A','T','C','G'},{'00','11','01','10'}); %// Lookup
C_3 = cellfun(@(x) values(codebook6, {x(1),x(2),x(3),x(4)}), ...
             dna_add_1, 'uni', 0);
B_11= cellfun(@cell2mat, C_3, 'uni', 0);
% disp(B);

for i=1:1:row
    for j=1:1:col
        decrypted(i,j,1) = bin2dec(B_11(i,j));
        
    end
end
figure(8)
imshow(decrypted,[]);
title('Decrypted image');    
%%% Performance Analysis
uaci=uacifun(inverse_im,DICOM_IMAGE)
npcr=npcrfun(DICOM_IMAGE,inverse_im,256,256)
corrc(inverse_im,DICOM_IMAGE)
%PSNR Peak Signal Noise Ratio
sum=0;
for i=1:row
    for j=1:col
        k=(double(DICOM_IMAGE(i,j))-double(inverse_im(i,j)));
        sum=sum+k*k;
    end
end
MSE=((sum/(256*256)));

PSNR=20*log10(255)-10*log10(MSE)
