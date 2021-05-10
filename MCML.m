function MCML()
global dz dr
global H1 F 
global lamda power
%dat bien toan cuc
%Mo phong cong suat den
p=0.0002; %Cong suat den LED (watt)
%w_l=[700,720];%pho lamda (nm)
%delta_lamda= 1; %do chia
t=2; %thoi gian chieu sang (s)
ttype='musc_interp.mat'; %loai mo 

led_interp();%chay chuong trinh noi suy led_interp.m
%load('led_interp.mat');
p=p*t; %nang luong den chieu trong khoang thoi gian t (J)
d=length(lamda);
e_n(2,d)=0;% khoi tao ma tran nang luong va so hat
N=0;
load(ttype);%lay du lieu tu file '.mat'

Wthr=0.001;%Nguong dung photon
spot=1.5; %ban kinh chum tia (nm)
m=10;%nguong xo so (nm)
c1=500;c2=500;dr=0.05;dz=0.05; % so o phan chia theo r va z
k=0;
n1=1;%{chiet suat khong khi}
n2=1.38; %{chiet suat mo}
tetac=asin(n1/n2); % goc phan xa toan phan
Rsp=((n1-n2)/(n1+n2))^2;% he so phan xa cua song toi
Rdif=Rsp;survive=0;
zam=0;totalref=0;
Q1(c1,c2)=0;% Khoi dau ma tran bang 0 
Q2(c1,c2)=0;

%tinh nang luong va so hat
tic
for i=1:d
    e_n(1,i) =((6.626*10.^(-34))*3*10.^8)/(lamda(i)*10.^(-9));%nang luong cua tung buoc song
    e_n(2,i)=round((p*power(i))/e_n(1,i)); %tinh so hat cho tung buoc song
    e_n(2,i)=round(e_n(2,i)/(10.^(10))); %giam so hat 
    N=N+e_n(2,i);
end

for a=1:d
    n=e_n(2,a);
    local=lamda(a)- musc_interp(1,1)+1;
    ma=musc_interp(local,2);
    ms=musc_interp(local,3);
    g=musc_interp(local,4);
    
    %mo phong 
    for nn=1:n
        W=1-Rsp;
        x1=spot*sqrt(-log(rand)/2);
        y1=0;
        z1=0;
        mx=0;my=0;mz=1;
        
        while W ~= 0 
              k=k+1;
              step=-log(rand)/(ma+ms);% mo phong quang duong tu do
              x1=x1+mx*step;
              y1=y1+my*step;
              z1=z1+mz*step;
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
              end
              %Ghi
              r=sqrt(x1*x1+y1*y1);
              i=round(r/dr+0.5);
              j=round(z1/dz+0.5);
              dQ=W*ma/(ma+ms);
              W=W*ms/(ma+ms);
              if (i<=c1)&&(j<=c2) 
                  Q1(i,j)= Q1(i,j)+dQ;
                  Q2(i,j)= Q2(i,j)+dQ/ma;   
              end
              teta=acos((1+g^2-((1-g^2)/(1-g+2*g*rand))^2)/2/g);
              fi=2*pi*rand;      
              if abs(mz)>0.9999
                 mx=sin(teta)*cos(fi);
                 my=sin(teta)*sin(fi);
                 mz=mz*cos(teta)/abs(mz);
              else
                 mx1=mx;
                 my1=my;
                 mz1=mz;
                 mx=sin(teta)*(mx1*mz1*cos(fi)-my1*sin(fi))/sqrt(1-mz1*mz1)+mx1*cos(teta);
                 my=sin(teta)*(my1*mz1*cos(fi)+mx1*sin(fi))/sqrt(1-mz1*mz1)+my1*cos(teta);
                 mz=-sin(teta)*cos(fi)*sqrt(1-mz1*mz1)+mz1*cos(teta);
             end
             if W<Wthr
                if rand<=1/m
                   W=m*W;
                   survive=survive+1;
                else
                   W=0;
                end
             end
        end
    end
     
end
toc

%chuyen doi nang luong
V=zeros;
for i=1:500
   V(i)=(2*i+1)*pi*dr^2*dz;	%(mm3)
   V=V';
end
H1=zeros(500,500);
H2=zeros(500,500);

for i=1:500
   H1(i,:)=(Q1(i,:)/N/V(i))*10.^(10);
   H2(i,:)=(Q2(i,:)/N/V(i))*10.^(10); 
end

P=5e-3;
F=H2*P*100;%(fluence rate W/cm^2)

look();%chay file look.m de hien thi

end



    