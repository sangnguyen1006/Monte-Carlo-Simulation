% Chuong trinh MCML bang Matlab
clear all
tic
% Buoc song
lambda=780;

% Nguong dung photon
Wthr=0.001;
% ban kinh chum tia
spot=1.5;%{mm}
% nguong xo so
m=10; %{nm}
% Mo phong don gian cho 1 lop vat lieu ----------------------------

N=10000;% So photon mo phong
c1=500;c2=500;dr=0.05;dz=0.05; % so o phan chia theo r va z (1 hat chiem 0.05 o)
ma=0.14; % he so hap thu 
ms=10;% he so tan xa
g=0.9;% he so bat dang huong
%//////////////////////////////////////////////////////////////////////////////////////////////////////
k=0;
n1=1;%{chiet suat khong khi}
n2=1.33; %{chiet suat mo}
tetac=asin(n1/n2); % goc phan xa toan phan (dung trong zl<0)
Rsp=((n1-n2)/(n1+n2))^2;% he so phan xa cua song toi tinh W
%///////////////////////////////////////////////////////////////////////////////////////////////////////
Rdif=Rsp;survive=0;
zam=0;totalref=0;

Q(c1,c2)=0;% Khoi dau ma tran bang 0    
% Bat dau m o phong 1 photon
%//////////////////////////////////////////////////////////////////////////////////////////////// 
for nn=1:N
    W=1-Rsp; %trong luong photon (1-Rsp:he so truyen qua)
    %khoi tao vi tri va huong di ban dau
    x1=spot*sqrt(-log(rand)/2); %random vi tri ban dau theo spot
    y1=0;
    z1=0;
    mx=0;my=0;mz=1; %d(0,0,1) vec to don vi di thang goc
%////////////////////////////////////////////////////////////////////////////////////////////////    
    while W ~= 0 
        k=k+1;
        step=-log(rand)/(ma+ms);% mo phong quang duong tu do S
        x1=x1+mx*step; %x1,y1,z1: vi tri moi cua photon (r=r+d*s),d la vec to huong di moi
        y1=y1+my*step; %mx,my,mz: huong di moi cua photon
        z1=z1+mz*step; %vi tri moi = vi tri cu + huong*quang duong tu do
%///////////////////////////////////////////////////////////////////////////////////////////////
        if z1<=0 % di nguoc ra ngoai
            z1=-z1;
            zam=zam+1;
            tetai=acos(abs(mz));
            if tetai<tetac
            tetar=asin(n2*sin(tetai)/n1);%{Snell's law}
            del=tetai-tetar;
            sum=tetai+tetar;
            Reflectance=0.5*(sin(del)/sin(sum))^2+0.5*(tan(del)/tan(sum))^2;
            Rdif=Rdif+(1-Reflectance)*W;
            else
            Reflectance=1;
            totalref=totalref+1;
            end
            W=W*Reflectance;
        end;
%///////////////////////////////////////////////////////////////////////////////////////////////////
        %Ghi
        %vi tri trong ma tran Q(C1,C2)
        r=sqrt(x1*x1+y1*y1);
        
        i=round(r/dr+0.5); %ban kinh
        j=round(z1/dz+0.5); %do sau
        
        dQ=W*ma/(ma+ms); %nang luong da duoc tich luy khi photon den 1 vi tri 
        W=W*ms/(ma+ms); %trong luong cua photon tiep tuc di chuyen
        if (i<=c1)&(j<=c2) 
            Q(i,j)= Q(i,j)+dQ; %tich luy nang luong tai 1 vi tri (x1,y1,z1)
        end
%////////////////////////////////////////////////////////////////////////////////////////////////////        
        %tinh goc tan xa va goc toi
        teta=acos((1+g^2-((1-g^2)/(1-g+2*g*rand))^2)/2/g); %goc tan xa
        fi=2*pi*rand; %goc toi
        %tinh toan vecto moi (huong di chuyen moi cua photon) sau khi tich
        %luy nang luong dQ va di tiep voi nang luong W 
        if abs(mz)>0.9999 %truong hop dat biet
            mx=sin(teta)*cos(fi);
            my=sin(teta)*sin(fi);
            mz=mz*cos(teta)/abs(mz);
        else
            mx1=mx;
            my1=my;
            mz1=mz;
            %truong hop tong quat
            mx=sin(teta)*(mx1*mz1*cos(fi)-my1*sin(fi))/sqrt(1-mz1*mz1)+mx1*cos(teta);
            my=sin(teta)*(my1*mz1*cos(fi)+mx1*sin(fi))/sqrt(1-mz1*mz1)+my1*cos(teta);
            mz=-sin(teta)*cos(fi)*sqrt(1-mz1*mz1)+mz1*cos(teta);
        end;
%////////////////////////////////////////////////////////////////////////////////////////////////////        
        %cham dut photon 
        if W<Wthr 
            if rand<=1/m %1 c? h?i song sot trong m co hoi duoc trao
                W=m*W;
                survive=survive+1;
            else
                W=0;
            end
        end
    end
end

%CHUYEN NANG LUONG
%Sau khi thuc hien chuong trinh C:\MATLAB\bin\monte.m voi N photon, 
%trong so duoc ghi trong phan tu luoi Q(i,j).
%chuong trinh nay chuyen thanh mat do nang luong voi don vi J/mm3 
%voi gia thiet nang luong cua N photon la 1J
%the tich cua phan thu luoi
V=zeros;
for i=1:500
   V(i)=(2*i+1)*pi*dr^2*dz;	%(mm3)
   V=V'
end
H=zeros(500,500);
for i=1:500
   H(i,:)=Q(i,:)/N/V(i);
end
toc
% plot
figure(1)
A=[0.2 0.15 0.1 0.08 0.06 0.04 0.02 0.008 0.004 0.002 ];
T=contour(H(1:100,1:150),A);
clabel(T,'manual');
title('su phan bo mat do nang luong, J/mm3');
xlabel(['do sau x ' num2str(dz) 'mm' ]);
ylabel(['ban kinh x ' num2str(dr) 'mm']);
gtext(['buoc song ' num2str(lambda) 'nm']);
gtext(['so photon  N=' num2str(N)]);
gtext('chum Gauss, ban kinh vet 1.5mm');

P=5e-3;
F=H/ma*P*100;%(fluence rate W/cm^2)
figure(2)
A=[0.2 0.15 0.1 0.08 0.06 0.04 0.02 0.008 0.004 0.002 0.0008 ];
FF=contour(F(1:100,1:150));
clabel(T,'manual');
title('su phan bo fluence rate, W/mm2');
xlabel(['do sau x ' num2str(dz) 'mm' ]);
ylabel(['ban kinh x ' num2str(dr) 'mm']);
gtext(['buoc song ' num2str(lambda) 'nm']);
gtext(['so photon  N=' num2str(N)]);
gtext('chum Gauss, ban kinh vet 1.5mm,cong suat 5mW');

