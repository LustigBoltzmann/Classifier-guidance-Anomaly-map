clear all 
alpha0=0*pi/180;  %%%空载接触角
Nb=14;             %%%滚珠数量
rb=4*10^(-3);              %%%滚珠半径
rd=22.5*10^(-3);          %%%节圆半径
%rL=12*10^(-6);    %%%径向间隙
rL=3*10^(-6);
%A0=0*10^(-6);           %%%空载时内外滚道曲率中心的相对距离
h=100*10^(-6);            %%%深度h

%h = 50*10^(-6);       %%% 故障轴承缺陷深度


phif=270*pi/180;    %%%phif =270+12.85°故障位置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phifd=phif*180/pi;    %%%phif =270+12.85°c
%Ro=22.4*10^(-3);                 %%%Ro外滚道半径
Ro=rb+rd; 
%k=1.45*10^10 ;  %%%赫兹接触的挠度系s数
K=5453*10^6;   %%%赫兹刚度
mi=2;
mo=2;
kox=212.2*10^6;
koy=212.2*10^6;
cox=1648.2;
c=200;
coy=1648.2;
fs=33.33;
fc=fs*(1-rb*cos(alpha0)/rd)/2;
fo=Nb*(1-rb*cos(alpha0)/rd)*fs/2;      %外圈故障特征频率
To=1/fo;                          %外圈故障特征周期
fi=Nb*(1+rb*cos(alpha0)/rd)*fs/2;      %内圈故障特征频率
fbpo=40.3;      %%%固有频率
Fx=0;                  %%%静力载荷
Fy=-100;
Fz=0;
Mx=0;                  %%%%力矩
My=0;
% steps=0.05;
% tmm=1/(360*fc)*steps;
% Fs=1/tmm;
Fs=65536;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%采样频率
Fs=8192;
Fs = 20000;
%Fs=8192;
tmm=1/Fs;
steps=tmm*(360*fc);
amount=10240*2;   
amount = 10240*2+2048%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%仿真点数
%****初值
%x1(1)=0; x2(1)=0; x3(1)=0; x4(1)=0; x5(1)=0; x6(1)=0; x7(1)=0; x8(1)=0;
x1=0; x2=0; x3=-17.5*10^-6; x4=0; x5=0; x6=0; x7=0; x8=0;
%***解微分方程程序

for i=1:amount
    Dphif(i)=(1)*pi/180;      %%%此时的Δ?f=5°或者40°或者70   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%故障宽度
        %%%%%%%%%%%%%%%%% Eq.(4)
    phicd(i)=steps*(i-1);
    tmms(i)=steps*(i-1)/(360*fc);
    phic(i)=(steps*(i-1))*pi/180;

    phi1(i)=phic(i)+2*pi*(1-1)/Nb;
    phi2(i)=phic(i)+2*pi*(2-1)/Nb;
    phi3(i)=phic(i)+2*pi*(3-1)/Nb;
    phi4(i)=phic(i)+2*pi*(4-1)/Nb;
    phi5(i)=phic(i)+2*pi*(5-1)/Nb;
    phi6(i)=phic(i)+2*pi*(6-1)/Nb;
    phi7(i)=phic(i)+2*pi*(7-1)/Nb;
    phi8(i)=phic(i)+2*pi*(8-1)/Nb;
    phi9(i)=phic(i)+2*pi*(9-1)/Nb;
    phi10(i)=phic(i)+2*pi*(10-1)/Nb;
    phi11(i)=phic(i)+2*pi*(11-1)/Nb;
    phi12(i)=phic(i)+2*pi*(12-1)/Nb;
    phi13(i)=phic(i)+2*pi*(13-1)/Nb;
    phi14(i)=phic(i)+2*pi*(14-1)/Nb;

    %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Eq.(6)
if  mod(phi1(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi1(i),2*pi)<=phif;
    a1(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi1(i),2*pi)-phif+0.5*Dphif(i))^2);
    d1(i)=min(a1(i),h); 
    %d1(i)=h;
    elseif mod(phi1(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi1(i),2*pi)>phif;
        b1(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi1(i),2*pi)+0.5*Dphif(i))^2);
        d1(i)=min(b1(i),h);
        %d1(i)=h;
    else 
        d1(i)=0;
end;   

if  mod(phi2(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi2(i),2*pi)<=phif;
    a2(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi2(i),2*pi)-phif+0.5*Dphif(i))^2);
    d2(i)=min(a2(i),h); 
    %d2(i)=h;
    elseif mod(phi2(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi2(i),2*pi)>phif;
        b2(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi2(i),2*pi)+0.5*Dphif(i))^2);
        d2(i)=min(b2(i),h);
        %d2(i)=h;
    else 
        d2(i)=0;
end;   
if  mod(phi3(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi3(i),2*pi)<=phif;
    a3(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi3(i),2*pi)-phif+0.5*Dphif(i))^2);
    d3(i)=min(a3(i),h); 
    %d3(i)=h;
    elseif mod(phi3(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi3(i),2*pi)>phif;
        b3(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi3(i),2*pi)+0.5*Dphif(i))^2);
        d3(i)=min(b3(i),h);
        %d3(i)=h;
    else 
        d3(i)=0;
end;   
if  mod(phi4(i),2*pi)>(phif-0.5*Dphif(i)) && mod(phi4(i),2*pi)<=phif;
    a4(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi4(i),2*pi)-phif+0.5*Dphif(i))^2);
    d4(i)=min(a4(i),h); 
    %d4(i)=h;
    elseif mod(phi4(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi4(i),2*pi)>phif;
        b4(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi4(i),2*pi)+0.5*Dphif(i))^2);
        d4(i)=min(b4(i),h);
        %d4(i)=h;
    else 
        d4(i)=0;
end;   
if  mod(phi5(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi5(i),2*pi)<=phif;
    a5(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi5(i),2*pi)-phif+0.5*Dphif(i))^2);
    d5(i)=min(a5(i),h); 
    %d5(i)=h;
    elseif mod(phi5(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi5(i),2*pi)>phif;
        b5(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi5(i),2*pi)+0.5*Dphif(i))^2);
        d5(i)=min(b5(i),h);
        %d5(i)=h;
    else 
        d5(i)=0;
end;   
if  mod(phi6(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi6(i),2*pi)<=phif;
    a6(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi6(i),2*pi)-phif+0.5*Dphif(i))^2);
    d6(i)=min(a6(i),h); 
    %d6(i)=h;
    elseif mod(phi6(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi6(i),2*pi)>phif;
        b6(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi6(i),2*pi)+0.5*Dphif(i))^2);
        d6(i)=min(b6(i),h);
        %d6(i)=h;
    else 
        d6(i)=0;
end;   
if  mod(phi7(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi7(i),2*pi)<=phif;
    a7(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi7(i),2*pi)-phif+0.5*Dphif(i))^2);
    d7(i)=min(a7(i),h); 
    %d7(i)=h;
    elseif mod(phi7(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi7(i),2*pi)>phif;
        b7(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi7(i),2*pi)+0.5*Dphif(i))^2);
        d7(i)=min(b7(i),h);
        %d7(i)=h;
    else 
        d7(i)=0;
end;   
if  mod(phi8(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi8(i),2*pi)<=phif;
    a8(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi8(i),2*pi)-phif+0.5*Dphif(i))^2);
    d8(i)=min(a8(i),h); 
    %d8(i)=h;
    elseif mod(phi8(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi8(i),2*pi)>phif;
        b8(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi8(i),2*pi)+0.5*Dphif(i))^2);
        d8(i)=min(b8(i),h);
        %d8(i)=h;
    else 
        d8(i)=0;
end;   
if  mod(phi9(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi9(i),2*pi)<=phif;
    a9(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi9(i),2*pi)-phif+0.5*Dphif(i))^2);
    d9(i)=min(a9(i),h); 
    %d9(i)=h;
    elseif mod(phi9(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi9(i),2*pi)>phif;
        b9(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi9(i),2*pi)+0.5*Dphif(i))^2);
        d9(i)=min(b9(i),h);
        %d9(i)=h;
    else 
        d9(i)=0;
end;   
if  mod(phi10(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi10(i),2*pi)<=phif;
    a10(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi10(i),2*pi)-phif+0.5*Dphif(i))^2);
    d10(i)=min(a10(i),h); 
    %d10(i)=h;
    elseif mod(phi10(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi10(i),2*pi)>phif;
        b10(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi10(i),2*pi)+0.5*Dphif(i))^2);
        d10(i)=min(b10(i),h);
        %d10(i)=h;
    else 
        d10(i)=0;
end;   
if  mod(phi11(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi11(i),2*pi)<=phif;
    a11(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi11(i),2*pi)-phif+0.5*Dphif(i))^2);
    d11(i)=min(a11(i),h); 
    %d11(i)=h;
    elseif mod(phi11(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi11(i),2*pi)>phif;
        b11(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi11(i),2*pi)+0.5*Dphif(i))^2);
        d11(i)=min(b11(i),h);
        %d11(i)=h;
    else 
        d11(i)=0;
end;   
if  mod(phi12(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi12(i),2*pi)<=phif;
    a12(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi12(i),2*pi)-phif+0.5*Dphif(i))^2);
    d12(i)=min(a12(i),h); 
    %d12(i)=h;
    elseif mod(phi12(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi12(i),2*pi)>phif;
        b12(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi12(i),2*pi)+0.5*Dphif(i))^2);
        d12(i)=min(b12(i),h);
        %d12(i)=h;
    else 
        d12(i)=0;
end;   
if  mod(phi13(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi13(i),2*pi)<=phif;
    a13(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi13(i),2*pi)-phif+0.5*Dphif(i))^2);
    d13(i)=min(a13(i),h); 
    %d13(i)=h;
    elseif mod(phi13(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi13(i),2*pi)>phif;
        b13(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi13(i),2*pi)+0.5*Dphif(i))^2);
        d13(i)=min(b13(i),h);
        %d13(i)=h;
    else 
        d13(i)=0;
end;   
if  mod(phi14(i),2*pi)>(phif-0.5*Dphif(i)) & mod(phi14(i),2*pi)<=phif;
    a14(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(mod(phi14(i),2*pi)-phif+0.5*Dphif(i))^2);
    d14(i)=min(a14(i),h); 
    %d14(i)=h;
    elseif mod(phi14(i),2*pi)<(phif+0.5*Dphif(i))& mod(phi14(i),2*pi)>phif;
        b14(i)=rb-sqrt((rb)^2-0.25*(Ro)^2*(phif-mod(phi14(i),2*pi)+0.5*Dphif(i))^2);
        d14(i)=min(b14(i),h);
        %d14(i)=h;
    else 
        d14(i)=0;
end;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Eq.

if (x1-x5)*cos(phi1(i))+(x3-x7)*sin(phi1(i))-rL-d1(i)>0;
    delta1(i)=(x1-x5)*cos(phi1(i))+(x3-x7)*sin(phi1(i))-rL-d1(i);
    gama1=1;
else 
    delta1(i)=0;
    gama1=0;
end


if (x1-x5)*cos(phi2(i))+(x3-x7)*sin(phi2(i))-rL-d2(i)>0;
    delta2(i)=(x1-x5)*cos(phi2(i))+(x3-x7)*sin(phi2(i))-rL-d2(i);
    gama2=1;
else 
    delta2(i)=0;
    gama2=0;
end

if (x1-x5)*cos(phi3(i))+(x3-x7)*sin(phi3(i))-rL-d3(i)>0;
    delta3(i)=(x1-x5)*cos(phi3(i))+(x3-x7)*sin(phi3(i))-rL-d3(i);
    gama3=1;
else 
    delta3(i)=0;
    gama3=0;
end


if (x1-x5)*cos(phi4(i))+(x3-x7)*sin(phi4(i))-rL-d4(i)>0;
    delta4(i)=(x1-x5)*cos(phi4(i))+(x3-x7)*sin(phi4(i))-rL-d4(i);
    gama4=1;
else 
    delta4(i)=0;
    gama4=0;
end


if (x1-x5)*cos(phi5(i))+(x3-x7)*sin(phi5(i))-rL-d5(i)>0;
    delta5(i)=(x1-x5)*cos(phi5(i))+(x3-x7)*sin(phi5(i))-rL-d5(i);
    gama5=1;
else 
    delta5(i)=0;
    gama5=0;
end

if (x1-x5)*cos(phi6(i))+(x3-x7)*sin(phi6(i))-rL-d6(i)>0;
    delta6(i)=(x1-x5)*cos(phi6(i))+(x3-x7)*sin(phi6(i))-rL-d6(i);
    gama6=1;
else 
    delta6(i)=0;
    gama6=0;
end

if (x1-x5)*cos(phi7(i))+(x3-x7)*sin(phi7(i))-rL-d7(i)>0;
    delta7(i)=(x1-x5)*cos(phi7(i))+(x3-x7)*sin(phi7(i))-rL-d7(i);
    gama7=1;
else 
    delta7(i)=0;
    gama7=0;
end

if (x1-x5)*cos(phi8(i))+(x3-x7)*sin(phi8(i))-rL-d8(i)>0;
    delta8(i)=(x1-x5)*cos(phi8(i))+(x3-x7)*sin(phi8(i))-rL-d8(i);
    gama8=1;
else 
    delta8(i)=0;
    gama8=0;
end

if (x1-x5)*cos(phi9(i))+(x3-x7)*sin(phi9(i))-rL-d9(i)>0;
    delta9(i)=(x1-x5)*cos(phi9(i))+(x3-x7)*sin(phi9(i))-rL-d9(i);
    gama9=1;
else 
    delta9(i)=0;
    gama9=0;
end

if (x1-x5)*cos(phi10(i))+(x3-x7)*sin(phi10(i))-rL-d10(i)>0;
    delta10(i)=(x1-x5)*cos(phi10(i))+(x3-x7)*sin(phi10(i))-rL-d10(i);
    gama10=1;
else 
    delta10(i)=0;
    gama10=0;
end

if (x1-x5)*cos(phi11(i))+(x3-x7)*sin(phi11(i))-rL-d11(i)>0;
    delta11(i)=(x1-x5)*cos(phi11(i))+(x3-x7)*sin(phi11(i))-rL-d11(i);
    gama11=1;
else 
    delta11(i)=0;
    gama11=0;
end

if (x1-x5)*cos(phi12(i))+(x3-x7)*sin(phi12(i))-rL-d12(i)>0;
    delta12(i)=(x1-x5)*cos(phi12(i))+(x3-x7)*sin(phi12(i))-rL-d12(i);
    gama12=1;
else 
    delta12(i)=0;
    gama12=0;
end

if (x1-x5)*cos(phi13(i))+(x3-x7)*sin(phi13(i))-rL-d13(i)>0;
    delta13(i)=(x1-x5)*cos(phi13(i))+(x3-x7)*sin(phi13(i))-rL-d13(i);
    gama13=1;
else 
    delta13(i)=0;
    gama13=0;
end

if (x1-x5)*cos(phi14(i))+(x3-x7)*sin(phi14(i))-rL-d14(i)>0;
    delta14(i)=(x1-x5)*cos(phi14(i))+(x3-x7)*sin(phi14(i))-rL-d14(i);
    gama14=1;
else 
    delta14(i)=0;
    gama14=0;
end
%sdelta1(i)=(x(2)-x(6))*cos(phi1(i))+(x(4)-x(8))*sin(phi1(i));
%sdelta2(i)=(x(2)-x(6))*cos(phi2(i))+(x(4)-x(8))*sin(phi2(i));
%sdelta3(i)=(x(2)-x(6))*cos(phi3(i))+(x(4)-x(8))*sin(phi3(i));
%sdelta4(i)=(x(2)-x(6))*cos(phi4(i))+(x(4)-x(8))*sin(phi4(i));
%sdelta5(i)=(x(2)-x(6))*cos(phi5(i))+(x(4)-x(8))*sin(phi5(i));
%sdelta6(i)=(x(2)-x(6))*cos(phi6(i))+(x(4)-x(8))*sin(phi6(i));
%sdelta7(i)=(x(2)-x(6))*cos(phi7(i))+(x(4)-x(8))*sin(phi7(i));
%sdelta8(i)=(x(2)-x(6))*cos(phi8(i))+(x(4)-x(8))*sin(phi8(i));
%sdelta9(i)=(x(2)-x(6))*cos(phi9(i))+(x(4)-x(8))*sin(phi9(i));
%sdelta10(i)=(x(2)-x(6))*cos(phi10(i))+(x(4)-x(8))*sin(phi10(i));
%sdelta11(i)=(x(2)-x(6))*cos(phi11(i))+(x(4)-x(8))*sin(phi11(i));
%sdelta12(i)=(x(2)-x(6))*cos(phi12(i))+(x(4)-x(8))*sin(phi12(i));
%sdelta13(i)=(x(2)-x(6))*cos(phi13(i))+(x(4)-x(8))*sin(phi13(i));
%sdelta14(i)=(x(2)-x(6))*cos(phi14(i))+(x(4)-x(8))*sin(phi14(i));
%cos(phi1(i))*(K*((x(1)-x(5))*cos(phi1(i))-rL-d1(i))+c*((x(2)-x(6))*cos(phi1(i))-rL-d1(i)))+cos(phi2(i))*(K*((x(1)-x(5))*cos(phi2(i))-rL-d2(i))+c*((x(2)-x(6))*cos(phi2(i))-rL-d2(i)))+cos(phi3(i))*(K*((x(1)-x(5))*cos(phi3(i))-rL-d3(i))+c*((x(2)-x(6))*cos(phi3(i))-rL-d3(i)))+cos(phi4(i))*(K*((x(1)-x(5))*cos(phi4(i))+-rL-d4(i))+c*((x(2)-x(6))*cos(phi4(i))-rL-d4(i)))+cos(phi5(i))*(K*((x(1)-x(5))*cos(phi5(i))-rL-d5(i))+c*((x(2)-x(6))*cos(phi5(i))-rL-d5(i)))+cos(phi6(i))*(K*((x(1)-x(5))*cos(phi6(i))-rL-d6(i))+c*((x(2)-x(6))*cos(phi6(i))-rL-d6(i)))+cos(phi7(i))*(K*((x(1)-x(5))*cos(phi7(i))-rL-d7(i))+c*((x(2)-x(6))*cos(phi7(i))-rL-d7(i)))+cos(phi8(i))*(K*((x(1)-x(5))*cos(phi8(i))-rL-d8(i))+c*((x(2)-x(6))*cos(phi8(i))-rL-d8(i)))+cos(phi9(i))*(K*((x(1)-x(5))*cos(phi9(i))-rL-d9(i))+c*((x(2)-x(6))*cos(phi9(i))-rL-d9(i)))+cos(phi10(i))*(K*((x(1)-x(5))*cos(phi10(i))-rL-d10(i))+c*((x(2)-x(6))*cos(phi10(i))-rL-d10(i)))+cos(phi11(i))*(K*((x(1)-x(5))*cos(phi11(i))-rL-d11(i))+c*((x(2)-x(6))*cos(phi11(i))-rL-d11(i)))+cos(phi12(i))*(K*((x(1)-x(5))*cos(phi12(i))-rL-d12(i))+c*((x(2)-x(6))*cos(phi12(i))-rL-d12(i)))+cos(phi13(i))*(K*((x(1)-x(5))*cos(phi13(i))-rL-d13(i))+c*((x(2)-x(6))*cos(phi13(i))-rL-d13(i)))+cos(phi14(i))*(K*((x(1)-x(5))*cos(phi14(i))-rL-d14(i))+c*((x(2)-x(6))*cos(phi14(i))-rL-d14(i)))

%cos(phi1(i))*(K*((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i))+c*((x(2)-x(6))*cos(phi1(i))+(0-x(8))*sin(phi1(i))-rL-d1(i)))+cos(phi2(i))*(K*((x(1)-x(5))*cos(phi2(i))+(x(3)-x(7))*sin(phi2(i))-rL-d2(i))+c*((x(2)-x(6))*cos(phi2(i))+(0-x(8))*sin(phi2(i))-rL-d2(i)))+cos(phi3(i))*(K*((x(1)-x(5))*cos(phi3(i))+(x(3)-x(7))*sin(phi3(i))-rL-d3(i))+c*((x(2)-x(6))*cos(phi3(i))+(0-x(8))*sin(phi3(i))-rL-d3(i)))+cos(phi4(i))*(K*((x(1)-x(5))*cos(phi4(i))+(x(3)-x(7))*sin(phi4(i))-rL-d4(i))+c*((x(2)-x(6))*cos(phi4(i))+(0-x(8))*sin(phi4(i))-rL-d4(i)))+cos(phi5(i))*(K*((x(1)-x(5))*cos(phi5(i))+(x(3)-x(7))*sin(phi5(i))-rL-d5(i))+c*((x(2)-x(6))*cos(phi5(i))+(0-x(8))*sin(phi5(i))-rL-d5(i)))+cos(phi6(i))*(K*((x(1)-x(5))*cos(phi6(i))+(x(3)-x(7))*sin(phi6(i))-rL-d6(i))+c*((x(2)-x(6))*cos(phi6(i))+(0-x(8))*sin(phi6(i))-rL-d6(i)))+cos(phi7(i))*(K*((x(1)-x(5))*cos(phi7(i))+(x(3)-x(7))*sin(phi7(i))-rL-d7(i))+c*((x(2)-x(6))*cos(phi7(i))+(0-x(8))*sin(phi7(i))-rL-d7(i)))+cos(phi8(i))*(K*((x(1)-x(5))*cos(phi8(i))+(x(3)-x(7))*sin(phi8(i))-rL-d8(i))+c*((x(2)-x(6))*cos(phi8(i))+(0-x(8))*sin(phi8(i))-rL-d8(i)))+cos(phi9(i))*(K*((x(1)-x(5))*cos(phi9(i))+(x(3)-x(7))*sin(phi9(i))-rL-d9(i))+c*((x(2)-x(6))*cos(phi9(i))+(0-x(8))*sin(phi9(i))-rL-d9(i)))+cos(phi10(i))*(K*((x(1)-x(5))*cos(phi10(i))+(x(3)-x(7))*sin(phi10(i))-rL-d10(i))+c*((x(2)-x(6))*cos(phi10(i))+(0-x(8))*sin(phi10(i))-rL-d10(i)))+cos(phi11(i))*(K*((x(1)-x(5))*cos(phi11(i))+(x(3)-x(7))*sin(phi11(i))-rL-d11(i))+c*((x(2)-x(6))*cos(phi11(i))+(0-x(8))*sin(phi11(i))-rL-d11(i)))+cos(phi12(i))*(K*((x(1)-x(5))*cos(phi12(i))+(x(3)-x(7))*sin(phi12(i))-rL-d12(i))+c*((x(2)-x(6))*cos(phi12(i))+(0-x(8))*sin(phi12(i))-rL-d12(i)))+cos(phi13(i))*(K*((x(1)-x(5))*cos(phi13(i))+(x(3)-x(7))*sin(phi13(i))-rL-d13(i))+c*((x(2)-x(6))*cos(phi13(i))+(0-x(8))*sin(phi13(i))-rL-d13(i)))+cos(phi14(i))*(K*((x(1)-x(5))*cos(phi14(i))+(x(3)-x(7))*sin(phi14(i))-rL-d14(i))+c*((x(2)-x(6))*cos(phi14(i))+(0-x(8))*sin(phi14(i))-rL-d14(i)))

%sin(phi1(i))*(K*((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i))+c*((x(2)-x(6))*cos(phi1(i))+(0-x(8))*sin(phi1(i))-rL-d1(i)))+sin(phi2(i))*(K*((x(1)-x(5))*cos(phi2(i))+(x(3)-x(7))*sin(phi2(i))-rL-d2(i))+c*((x(2)-x(6))*cos(phi2(i))+(0-x(8))*sin(phi2(i))-rL-d2(i)))+sin(phi3(i))*(K*((x(1)-x(5))*cos(phi3(i))+(x(3)-x(7))*sin(phi3(i))-rL-d3(i))+c*((x(2)-x(6))*cos(phi3(i))+(0-x(8))*sin(phi3(i))-rL-d3(i)))+sin(phi4(i))*(K*((x(1)-x(5))*cos(phi4(i))+(x(3)-x(7))*sin(phi4(i))-rL-d4(i))+c*((x(2)-x(6))*cos(phi4(i))+(0-x(8))*sin(phi4(i))-rL-d4(i)))+sin(phi5(i))*(K*((x(1)-x(5))*cos(phi5(i))+(x(3)-x(7))*sin(phi5(i))-rL-d5(i))+c*((x(2)-x(6))*cos(phi5(i))+(0-x(8))*sin(phi5(i))-rL-d5(i)))+sin(phi6(i))*(K*((x(1)-x(5))*cos(phi6(i))+(x(3)-x(7))*sin(phi6(i))-rL-d6(i))+c*((x(2)-x(6))*cos(phi6(i))+(0-x(8))*sin(phi6(i))-rL-d6(i)))+sin(phi7(i))*(K*((x(1)-x(5))*cos(phi7(i))+(x(3)-x(7))*sin(phi7(i))-rL-d7(i))+c*((x(2)-x(6))*cos(phi7(i))+(0-x(8))*sin(phi7(i))-rL-d7(i)))+sin(phi8(i))*(K*((x(1)-x(5))*cos(phi8(i))+(x(3)-x(7))*sin(phi8(i))-rL-d8(i))+c*((x(2)-x(6))*cos(phi8(i))+(0-x(8))*sin(phi8(i))-rL-d8(i)))+sin(phi9(i))*(K*((x(1)-x(5))*cos(phi9(i))+(x(3)-x(7))*sin(phi9(i))-rL-d9(i))+c*((x(2)-x(6))*cos(phi9(i))+(0-x(8))*sin(phi9(i))-rL-d9(i)))+sin(phi10(i))*(K*((x(1)-x(5))*cos(phi10(i))+(x(3)-x(7))*sin(phi10(i))-rL-d10(i))+c*((x(2)-x(6))*cos(phi10(i))+(0-x(8))*sin(phi10(i))-rL-d10(i)))+sin(phi11(i))*(K*((x(1)-x(5))*cos(phi11(i))+(x(3)-x(7))*sin(phi11(i))-rL-d11(i))+c*((x(2)-x(6))*cos(phi11(i))+(0-x(8))*sin(phi11(i))-rL-d11(i)))+sin(phi12(i))*(K*((x(1)-x(5))*cos(phi12(i))+(x(3)-x(7))*sin(phi12(i))-rL-d12(i))+c*((x(2)-x(6))*cos(phi12(i))+(0-x(8))*sin(phi12(i))-rL-d12(i)))+sin(phi13(i))*(K*((x(1)-x(5))*cos(phi13(i))+(x(3)-x(7))*sin(phi13(i))-rL-d13(i))+c*((x(2)-x(6))*cos(phi13(i))+(0-x(8))*sin(phi13(i))-rL-d13(i)))+sin(phi14(i))*(K*((x(1)-x(5))*cos(phi14(i))+(x(3)-x(7))*sin(phi14(i))-rL-d14(i))+c*((x(2)-x(6))*cos(phi14(i))+(0-x(8))*sin(phi14(i))-rL-d14(i)))

%kox*cos(phi1(i))*((x1-x5)*cos(phi1(i))+(x3-x7)*sin(phi1(i))-rL-d1(i))+cox*cos(phi1(i))*((x2-x6)*cos(phi1(i))+(x4-x8)*sin(phi1(i))-rL-d1(i));
%koy*sin(phi1(i))*((x1-x5)*cos(phi1(i))+(x3-x7)*sin(phi1(i))-rL-d1(i))+coy*sin(phi1(i))*((x2-x6)*cos(phi1(i))+(x4-x8)*sin(phi1(i))-rL-d1(i));
%kox*cos(phi13(i))*((x1-x5)*cos(phi13(i))+(x3-x7)*sin(phi13(i))-rL-d13(i))+cox*cos(phi13(i))*((x2-x6)*cos(phi13(i))+(x4-x8)*sin(phi13(i))-rL-d13(i))+kox*cos(phi14(i))*((x1-x5)*cos(phi14(i))+(x3-x7)*sin(phi14(i))-rL-d14(i))+cox*cos(phi14(i))*((x2-x6)*cos(phi14(i))+(x4-x8)*sin(phi14(i))-rL-d14(i))+kox*cos(phi1(i))*((x1-x5)*cos(phi1(i))+(x3-x7)*sin(phi1(i))-rL-d1(i))+cox*cos(phi1(i))*((x2-x6)*cos(phi1(i))+(x4-x8)*sin(phi1(i))-rL-d1(i))+kox*cos(phi2(i))*((x1-x5)*cos(phi2(i))+(x3-x7)*sin(phi2(i))-rL-d2(i))+cox*cos(phi2(i))*((x2-x6)*cos(phi2(i))+(x4-x8)*sin(phi2(i))-rL-d2(i))+kox*cos(phi3(i))*((x1-x5)*cos(phi3(i))+(x3-x7)*sin(phi3(i))-rL-d3(i))+cox*cos(phi3(i))*((x2-x6)*cos(phi3(i))+(x4-x8)*sin(phi3(i))-rL-d3(i))+kox*cos(phi4(i))*((x1-x5)*cos(phi4(i))+(x3-x7)*sin(phi4(i))-rL-d4(i))+cox*cos(phi4(i))*((x2-x6)*cos(phi4(i))+(x4-x8)*sin(phi4(i))-rL-d4(i))+kox*cos(phi5(i))*((x1-x5)*cos(phi5(i))+(x3-x7)*sin(phi5(i))-rL-d5(i))+cox*cos(phi5(i))*((x2-x6)*cos(phi5(i))+(x4-x8)*sin(phi5(i))-rL-d5(i))+kox*cos(phi6(i))*((x1-x5)*cos(phi6(i))+(x3-x7)*sin(phi6(i))-rL-d6(i))+cox*cos(phi6(i))*((x2-x6)*cos(phi6(i))+(x4-x8)*sin(phi6(i))-rL-d6(i))+kox*cos(phi7(i))*((x1-x5)*cos(phi7(i))+(x3-x7)*sin(phi7(i))-rL-d7(i))+cox*cos(phi7(i))*((x2-x6)*cos(phi7(i))+(x4-x8)*sin(phi7(i))-rL-d7(i))+kox*cos(phi8(i))*((x1-x5)*cos(phi8(i))+(x3-x7)*sin(phi8(i))-rL-d8(i))+cox*cos(phi8(i))*((x2-x6)*cos(phi8(i))+(x4-x8)*sin(phi8(i))-rL-d8(i))+kox*cos(phi9(i))*((x1-x5)*cos(phi9(i))+(x3-x7)*sin(phi9(i))-rL-d9(i))+cox*cos(phi9(i))*((x2-x6)*cos(phi9(i))+(x4-x8)*sin(phi9(i))-rL-d9(i))+kox*cos(phi10(i))*((x1-x5)*cos(phi10(i))+(x3-x7)*sin(phi10(i))-rL-d10(i))+cox*cos(phi10(i))*((x2-x6)*cos(phi10(i))+(x4-x8)*sin(phi10(i))-rL-d10(i))+kox*cos(phi11(i))*((x1-x5)*cos(phi11(i))+(x3-x7)*sin(phi11(i))-rL-d11(i))+cox*cos(phi11(i))*((x2-x6)*cos(phi11(i))+(x4-x8)*sin(phi11(i))-rL-d11(i))+kox*cos(phi12(i))*((x1-x5)*cos(phi12(i))+(x3-x7)*sin(phi12(i))-rL-d12(i))+cox*cos(phi12(i))*((x2-x6)*cos(phi12(i))+(x4-x8)*sin(phi12(i))-rL-d12(i))
%koy*sin(phi1(i))*((x1-x5)*cos(phi1(i))+(x3-x7)*sin(phi1(i))-rL-d1(i))+coy*sin(phi1(i))*((x2-x6)*cos(phi1(i))+(x4-x8)*sin(phi1(i))-rL-d1(i))+koy*sin(phi2(i))*((x1-x5)*cos(phi2(i))+(x3-x7)*sin(phi2(i))-rL-d1(i))+coy*sin(phi2(i))*((x2-x6)*cos(phi2(i))+(x4-x8)*sin(phi2(i))-rL-d2(i))+koy*sin(phi3(i))*((x1-x5)*cos(phi3(i))+(x3-x7)*sin(phi3(i))-rL-d1(i))+coy*sin(phi3(i))*((x2-x6)*cos(phi3(i))+(x4-x8)*sin(phi3(i))-rL-d3(i))+koy*sin(phi4(i))*((x1-x5)*cos(phi4(i))+(x3-x7)*sin(phi4(i))-rL-d1(i))+coy*sin(phi4(i))*((x2-x6)*cos(phi4(i))+(x4-x8)*sin(phi4(i))-rL-d4(i))+koy*sin(phi5(i))*((x1-x5)*cos(phi5(i))+(x3-x7)*sin(phi5(i))-rL-d1(i))+coy*sin(phi5(i))*((x2-x6)*cos(phi5(i))+(x4-x8)*sin(phi5(i))-rL-d5(i))+koy*sin(phi6(i))*((x1-x5)*cos(phi6(i))+(x3-x7)*sin(phi6(i))-rL-d1(i))+coy*sin(phi6(i))*((x2-x6)*cos(phi6(i))+(x4-x8)*sin(phi6(i))-rL-d6(i))+koy*sin(phi7(i))*((x1-x5)*cos(phi7(i))+(x3-x7)*sin(phi7(i))-rL-d1(i))+coy*sin(phi7(i))*((x2-x6)*cos(phi7(i))+(x4-x8)*sin(phi7(i))-rL-d7(i))+koy*sin(phi8(i))*((x1-x5)*cos(phi8(i))+(x3-x7)*sin(phi8(i))-rL-d1(i))+coy*sin(phi8(i))*((x2-x6)*cos(phi8(i))+(x4-x8)*sin(phi8(i))-rL-d8(i))+koy*sin(phi9(i))*((x1-x5)*cos(phi9(i))+(x3-x7)*sin(phi9(i))-rL-d1(i))+coy*sin(phi9(i))*((x2-x6)*cos(phi9(i))+(x4-x8)*sin(phi9(i))-rL-d9(i))+koy*sin(phi10(i))*((x1-x5)*cos(phi10(i))+(x3-x7)*sin(phi10(i))-rL-d1(i))+coy*sin(phi10(i))*((x2-x6)*cos(phi10(i))+(x4-x8)*sin(phi10(i))-rL-d10(i))+koy*sin(phi11(i))*((x1-x5)*cos(phi11(i))+(x3-x7)*sin(phi11(i))-rL-d1(i))+coy*sin(phi11(i))*((x2-x6)*cos(phi11(i))+(x4-x8)*sin(phi11(i))-rL-d11(i))+koy*sin(phi12(i))*((x1-x5)*cos(phi12(i))+(x3-x7)*sin(phi12(i))-rL-d1(i))+coy*sin(phi12(i))*((x2-x6)*cos(phi12(i))+(x4-x8)*sin(phi12(i))-rL-d12(i))+koy*sin(phi13(i))*((x1-x5)*cos(phi13(i))+(x3-x7)*sin(phi13(i))-rL-d1(i))+coy*sin(phi13(i))*((x2-x6)*cos(phi13(i))+(x4-x8)*sin(phi13(i))-rL-d13(i))+koy*sin(phi14(i))*((x1-x5)*cos(phi14(i))+(x3-x7)*sin(phi14(i))-rL-d1(i))+coy*sin(phi14(i))*((x2-x6)*cos(phi14(i))+(x4-x8)*sin(phi14(i))-rL-d14(i))

odefun=@(t,x)[x(2);
            Fx/mi-(gama1*(K*(((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i)))^1.5+c*((x(2)-x(6))*cos(phi1(i))+(x(4)-x(8))*sin(phi1(i))))*cos(phi1(i))+gama2*(K*(((x(1)-x(5))*cos(phi2(i))+(x(3)-x(7))*sin(phi2(i))-rL-d2(i)))^1.5+c*((x(2)-x(6))*cos(phi2(i))+(x(4)-x(8))*sin(phi2(i))))*cos(phi2(i))+gama3*(K*(((x(1)-x(5))*cos(phi3(i))+(x(3)-x(7))*sin(phi3(i))-rL-d3(i)))^1.5+c*((x(2)-x(6))*cos(phi3(i))+(x(4)-x(8))*sin(phi3(i))))*cos(phi3(i))+gama4*(K*(((x(1)-x(5))*cos(phi4(i))+(x(3)-x(7))*sin(phi4(i))-rL-d4(i)))^1.5+c*((x(2)-x(6))*cos(phi4(i))+(x(4)-x(8))*sin(phi4(i))))*cos(phi4(i))+gama5*(K*(((x(1)-x(5))*cos(phi5(i))+(x(3)-x(7))*sin(phi5(i))-rL-d5(i)))^1.5+c*((x(2)-x(6))*cos(phi5(i))+(x(4)-x(8))*sin(phi5(i))))*cos(phi5(i))+gama6*(K*(((x(1)-x(5))*cos(phi6(i))+(x(3)-x(7))*sin(phi6(i))-rL-d6(i)))^1.5+c*((x(2)-x(6))*cos(phi6(i))+(x(4)-x(8))*sin(phi6(i))))*cos(phi6(i))+gama7*(K*(((x(1)-x(5))*cos(phi7(i))+(x(3)-x(7))*sin(phi7(i))-rL-d7(i)))^1.5+c*((x(2)-x(6))*cos(phi7(i))+(x(4)-x(8))*sin(phi7(i))))*cos(phi7(i))+gama8*(K*(((x(1)-x(5))*cos(phi8(i))+(x(3)-x(7))*sin(phi8(i))-rL-d8(i)))^1.5+c*((x(2)-x(6))*cos(phi8(i))+(x(4)-x(8))*sin(phi8(i))))*cos(phi8(i))+gama9*(K*(((x(1)-x(5))*cos(phi9(i))+(x(3)-x(7))*sin(phi9(i))-rL-d9(i)))^1.5+c*((x(2)-x(6))*cos(phi9(i))+(x(4)-x(8))*sin(phi9(i))))*cos(phi9(i))+gama10*(K*(((x(1)-x(5))*cos(phi10(i))+(x(3)-x(7))*sin(phi10(i))-rL-d10(i)))^1.5+c*((x(2)-x(6))*cos(phi10(i))+(x(4)-x(8))*sin(phi10(i))))*cos(phi10(i))+gama11*(K*(((x(1)-x(5))*cos(phi11(i))+(x(3)-x(7))*sin(phi11(i))-rL-d11(i)))^1.5+c*((x(2)-x(6))*cos(phi11(i))+(x(4)-x(8))*sin(phi11(i))))*cos(phi11(i))+gama12*(K*(((x(1)-x(5))*cos(phi12(i))+(x(3)-x(7))*sin(phi12(i))-rL-d12(i)))^1.5+c*((x(2)-x(6))*cos(phi12(i))+(x(4)-x(8))*sin(phi12(i))))*cos(phi12(i))+gama13*(K*(((x(1)-x(5))*cos(phi13(i))+(x(3)-x(7))*sin(phi13(i))-rL-d13(i)))^1.5+c*((x(2)-x(6))*cos(phi13(i))+(x(4)-x(8))*sin(phi13(i))))*cos(phi13(i))+gama14*(K*(((x(1)-x(5))*cos(phi14(i))+(x(3)-x(7))*sin(phi14(i))-rL-d14(i)))^1.5+c*((x(2)-x(6))*cos(phi14(i))+(x(4)-x(8))*sin(phi14(i))))*cos(phi14(i)))/mi;
              x(4);
            Fy/mi-(gama1*(K*(((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i)))^1.5+c*((x(2)-x(6))*cos(phi1(i))+(x(4)-x(8))*sin(phi1(i))))*sin(phi1(i))+gama2*(K*(((x(1)-x(5))*cos(phi2(i))+(x(3)-x(7))*sin(phi2(i))-rL-d2(i)))^1.5+c*((x(2)-x(6))*cos(phi2(i))+(x(4)-x(8))*sin(phi2(i))))*sin(phi2(i))+gama3*(K*(((x(1)-x(5))*cos(phi3(i))+(x(3)-x(7))*sin(phi3(i))-rL-d3(i)))^1.5+c*((x(2)-x(6))*cos(phi3(i))+(x(4)-x(8))*sin(phi3(i))))*sin(phi3(i))+gama4*(K*(((x(1)-x(5))*cos(phi4(i))+(x(3)-x(7))*sin(phi4(i))-rL-d4(i)))^1.5+c*((x(2)-x(6))*cos(phi4(i))+(x(4)-x(8))*sin(phi4(i))))*sin(phi4(i))+gama5*(K*(((x(1)-x(5))*cos(phi5(i))+(x(3)-x(7))*sin(phi5(i))-rL-d5(i)))^1.5+c*((x(2)-x(6))*cos(phi5(i))+(x(4)-x(8))*sin(phi5(i))))*sin(phi5(i))+gama6*(K*(((x(1)-x(5))*cos(phi6(i))+(x(3)-x(7))*sin(phi6(i))-rL-d6(i)))^1.5+c*((x(2)-x(6))*cos(phi6(i))+(x(4)-x(8))*sin(phi6(i))))*sin(phi6(i))+gama7*(K*(((x(1)-x(5))*cos(phi7(i))+(x(3)-x(7))*sin(phi7(i))-rL-d7(i)))^1.5+c*((x(2)-x(6))*cos(phi7(i))+(x(4)-x(8))*sin(phi7(i))))*sin(phi7(i))+gama8*(K*(((x(1)-x(5))*cos(phi8(i))+(x(3)-x(7))*sin(phi8(i))-rL-d8(i)))^1.5+c*((x(2)-x(6))*cos(phi8(i))+(x(4)-x(8))*sin(phi8(i))))*sin(phi8(i))+gama9*(K*(((x(1)-x(5))*cos(phi9(i))+(x(3)-x(7))*sin(phi9(i))-rL-d9(i)))^1.5+c*((x(2)-x(6))*cos(phi9(i))+(x(4)-x(8))*sin(phi9(i))))*sin(phi9(i))+gama10*(K*(((x(1)-x(5))*cos(phi10(i))+(x(3)-x(7))*sin(phi10(i))-rL-d10(i)))^1.5+c*((x(2)-x(6))*cos(phi10(i))+(x(4)-x(8))*sin(phi10(i))))*sin(phi10(i))+gama11*(K*(((x(1)-x(5))*cos(phi11(i))+(x(3)-x(7))*sin(phi11(i))-rL-d11(i)))^1.5+c*((x(2)-x(6))*cos(phi11(i))+(x(4)-x(8))*sin(phi11(i))))*sin(phi11(i))+gama12*(K*(((x(1)-x(5))*cos(phi12(i))+(x(3)-x(7))*sin(phi12(i))-rL-d12(i)))^1.5+c*((x(2)-x(6))*cos(phi12(i))+(x(4)-x(8))*sin(phi12(i))))*sin(phi12(i))+gama13*(K*(((x(1)-x(5))*cos(phi13(i))+(x(3)-x(7))*sin(phi13(i))-rL-d13(i)))^1.5+c*((x(2)-x(6))*cos(phi13(i))+(x(4)-x(8))*sin(phi13(i))))*sin(phi13(i))+gama14*(K*(((x(1)-x(5))*cos(phi14(i))+(x(3)-x(7))*sin(phi14(i))-rL-d14(i)))^1.5+c*((x(2)-x(6))*cos(phi14(i))+(x(4)-x(8))*sin(phi14(i))))*sin(phi14(i)))/mi;
              x(6);
            -cox*x(6)/mo-kox*x(5)/mo+(gama1*(K*(((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i)))^1.5+c*((x(2)-x(6))*cos(phi1(i))+(x(4)-x(8))*sin(phi1(i))))*cos(phi1(i))+gama2*(K*(((x(1)-x(5))*cos(phi2(i))+(x(3)-x(7))*sin(phi2(i))-rL-d2(i)))^1.5+c*((x(2)-x(6))*cos(phi2(i))+(x(4)-x(8))*sin(phi2(i))))*cos(phi2(i))+gama3*(K*(((x(1)-x(5))*cos(phi3(i))+(x(3)-x(7))*sin(phi3(i))-rL-d3(i)))^1.5+c*((x(2)-x(6))*cos(phi3(i))+(x(4)-x(8))*sin(phi3(i))))*cos(phi3(i))+gama4*(K*(((x(1)-x(5))*cos(phi4(i))+(x(3)-x(7))*sin(phi4(i))-rL-d4(i)))^1.5+c*((x(2)-x(6))*cos(phi4(i))+(x(4)-x(8))*sin(phi4(i))))*cos(phi4(i))+gama5*(K*(((x(1)-x(5))*cos(phi5(i))+(x(3)-x(7))*sin(phi5(i))-rL-d5(i)))^1.5+c*((x(2)-x(6))*cos(phi5(i))+(x(4)-x(8))*sin(phi5(i))))*cos(phi5(i))+gama6*(K*(((x(1)-x(5))*cos(phi6(i))+(x(3)-x(7))*sin(phi6(i))-rL-d6(i)))^1.5+c*((x(2)-x(6))*cos(phi6(i))+(x(4)-x(8))*sin(phi6(i))))*cos(phi6(i))+gama7*(K*(((x(1)-x(5))*cos(phi7(i))+(x(3)-x(7))*sin(phi7(i))-rL-d7(i)))^1.5+c*((x(2)-x(6))*cos(phi7(i))+(x(4)-x(8))*sin(phi7(i))))*cos(phi7(i))+gama8*(K*(((x(1)-x(5))*cos(phi8(i))+(x(3)-x(7))*sin(phi8(i))-rL-d8(i)))^1.5+c*((x(2)-x(6))*cos(phi8(i))+(x(4)-x(8))*sin(phi8(i))))*cos(phi8(i))+gama9*(K*(((x(1)-x(5))*cos(phi9(i))+(x(3)-x(7))*sin(phi9(i))-rL-d9(i)))^1.5+c*((x(2)-x(6))*cos(phi9(i))+(x(4)-x(8))*sin(phi9(i))))*cos(phi9(i))+gama10*(K*(((x(1)-x(5))*cos(phi10(i))+(x(3)-x(7))*sin(phi10(i))-rL-d10(i)))^1.5+c*((x(2)-x(6))*cos(phi10(i))+(x(4)-x(8))*sin(phi10(i))))*cos(phi10(i))+gama11*(K*(((x(1)-x(5))*cos(phi11(i))+(x(3)-x(7))*sin(phi11(i))-rL-d11(i)))^1.5+c*((x(2)-x(6))*cos(phi11(i))+(x(4)-x(8))*sin(phi11(i))))*cos(phi11(i))+gama12*(K*(((x(1)-x(5))*cos(phi12(i))+(x(3)-x(7))*sin(phi12(i))-rL-d12(i)))^1.5+c*((x(2)-x(6))*cos(phi12(i))+(x(4)-x(8))*sin(phi12(i))))*cos(phi12(i))+gama13*(K*(((x(1)-x(5))*cos(phi13(i))+(x(3)-x(7))*sin(phi13(i))-rL-d13(i)))^1.5+c*((x(2)-x(6))*cos(phi13(i))+(x(4)-x(8))*sin(phi13(i))))*cos(phi13(i))+gama14*(K*(((x(1)-x(5))*cos(phi14(i))+(x(3)-x(7))*sin(phi14(i))-rL-d14(i)))^1.5+c*((x(2)-x(6))*cos(phi14(i))+(x(4)-x(8))*sin(phi14(i))))*cos(phi14(i)))/mo;
              x(8);
            -coy*x(8)/mo-koy*x(7)/mo+(gama1*(K*(((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i)))^1.5+c*((x(2)-x(6))*cos(phi1(i))+(x(4)-x(8))*sin(phi1(i))))*sin(phi1(i))+gama2*(K*(((x(1)-x(5))*cos(phi2(i))+(x(3)-x(7))*sin(phi2(i))-rL-d2(i)))^1.5+c*((x(2)-x(6))*cos(phi2(i))+(x(4)-x(8))*sin(phi2(i))))*sin(phi2(i))+gama3*(K*(((x(1)-x(5))*cos(phi3(i))+(x(3)-x(7))*sin(phi3(i))-rL-d3(i)))^1.5+c*((x(2)-x(6))*cos(phi3(i))+(x(4)-x(8))*sin(phi3(i))))*sin(phi3(i))+gama4*(K*(((x(1)-x(5))*cos(phi4(i))+(x(3)-x(7))*sin(phi4(i))-rL-d4(i)))^1.5+c*((x(2)-x(6))*cos(phi4(i))+(x(4)-x(8))*sin(phi4(i))))*sin(phi4(i))+gama5*(K*(((x(1)-x(5))*cos(phi5(i))+(x(3)-x(7))*sin(phi5(i))-rL-d5(i)))^1.5+c*((x(2)-x(6))*cos(phi5(i))+(x(4)-x(8))*sin(phi5(i))))*sin(phi5(i))+gama6*(K*(((x(1)-x(5))*cos(phi6(i))+(x(3)-x(7))*sin(phi6(i))-rL-d6(i)))^1.5+c*((x(2)-x(6))*cos(phi6(i))+(x(4)-x(8))*sin(phi6(i))))*sin(phi6(i))+gama7*(K*(((x(1)-x(5))*cos(phi7(i))+(x(3)-x(7))*sin(phi7(i))-rL-d7(i)))^1.5+c*((x(2)-x(6))*cos(phi7(i))+(x(4)-x(8))*sin(phi7(i))))*sin(phi7(i))+gama8*(K*(((x(1)-x(5))*cos(phi8(i))+(x(3)-x(7))*sin(phi8(i))-rL-d8(i)))^1.5+c*((x(2)-x(6))*cos(phi8(i))+(x(4)-x(8))*sin(phi8(i))))*sin(phi8(i))+gama9*(K*(((x(1)-x(5))*cos(phi9(i))+(x(3)-x(7))*sin(phi9(i))-rL-d9(i)))^1.5+c*((x(2)-x(6))*cos(phi9(i))+(x(4)-x(8))*sin(phi9(i))))*sin(phi9(i))+gama10*(K*(((x(1)-x(5))*cos(phi10(i))+(x(3)-x(7))*sin(phi10(i))-rL-d10(i)))^1.5+c*((x(2)-x(6))*cos(phi10(i))+(x(4)-x(8))*sin(phi10(i))))*sin(phi10(i))+gama11*(K*(((x(1)-x(5))*cos(phi11(i))+(x(3)-x(7))*sin(phi11(i))-rL-d11(i)))^1.5+c*((x(2)-x(6))*cos(phi11(i))+(x(4)-x(8))*sin(phi11(i))))*sin(phi11(i))+gama12*(K*(((x(1)-x(5))*cos(phi12(i))+(x(3)-x(7))*sin(phi12(i))-rL-d12(i)))^1.5+c*((x(2)-x(6))*cos(phi12(i))+(x(4)-x(8))*sin(phi12(i))))*sin(phi12(i))+gama13*(K*(((x(1)-x(5))*cos(phi13(i))+(x(3)-x(7))*sin(phi13(i))-rL-d13(i)))^1.5+c*((x(2)-x(6))*cos(phi13(i))+(x(4)-x(8))*sin(phi13(i))))*sin(phi13(i))+gama14*(K*(((x(1)-x(5))*cos(phi14(i))+(x(3)-x(7))*sin(phi14(i))-rL-d14(i)))^1.5+c*((x(2)-x(6))*cos(phi14(i))+(x(4)-x(8))*sin(phi14(i))))*sin(phi14(i)))/mo];

%(cos(phi1(i))*(K*gama11((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i))^1.5)+cos(phi1(i))*(K*gama11((x(1)-x(5))*cos(phi1(i))+(x(3)-x(7))*sin(phi1(i))-rL-d1(i))^1.5));

x0=[x1 x2 x3 x4 x5 x6 x7 x8 ]; %初值
tspan=(i-1)*tmm:tmm:i*tmm;    %求解区间
%options=odeset('reltol',1e-10,'abstol',1e-10); %精度
options=odeset('relTol',1e-3,'Maxstep',0.1);

sol=ode15s(odefun,tspan,x0,options);
xint=i*tmm;                      %取t=i*tmm时各状态变量的值
sxint= deval(sol,xint);

%以下代入式，求xo的二阶导
axo(i)=-cox*sxint(6)/mo-kox*sxint(5)/mo+(gama1*(K*(((sxint(1)-sxint(5))*cos(phi1(i))+(sxint(3)-sxint(7))*sin(phi1(i))-rL-d1(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi1(i))+(sxint(4)-sxint(8))*sin(phi1(i))))*cos(phi1(i))+gama2*(K*(((sxint(1)-sxint(5))*cos(phi2(i))+(sxint(3)-sxint(7))*sin(phi2(i))-rL-d2(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi2(i))+(sxint(4)-sxint(8))*sin(phi2(i))))*cos(phi2(i))+gama3*(K*(((sxint(1)-sxint(5))*cos(phi3(i))+(sxint(3)-sxint(7))*sin(phi3(i))-rL-d3(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi3(i))+(sxint(4)-sxint(8))*sin(phi3(i))))*cos(phi3(i))+gama4*(K*(((sxint(1)-sxint(5))*cos(phi4(i))+(sxint(3)-sxint(7))*sin(phi4(i))-rL-d4(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi4(i))+(sxint(4)-sxint(8))*sin(phi4(i))))*cos(phi4(i))+gama5*(K*(((sxint(1)-sxint(5))*cos(phi5(i))+(sxint(3)-sxint(7))*sin(phi5(i))-rL-d5(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi5(i))+(sxint(4)-sxint(8))*sin(phi5(i))))*cos(phi5(i))+gama6*(K*(((sxint(1)-sxint(5))*cos(phi6(i))+(sxint(3)-sxint(7))*sin(phi6(i))-rL-d6(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi6(i))+(sxint(4)-sxint(8))*sin(phi6(i))))*cos(phi6(i))+gama7*(K*(((sxint(1)-sxint(5))*cos(phi7(i))+(sxint(3)-sxint(7))*sin(phi7(i))-rL-d7(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi7(i))+(sxint(4)-sxint(8))*sin(phi7(i))))*cos(phi7(i))+gama8*(K*(((sxint(1)-sxint(5))*cos(phi8(i))+(sxint(3)-sxint(7))*sin(phi8(i))-rL-d8(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi8(i))+(sxint(4)-sxint(8))*sin(phi8(i))))*cos(phi8(i))+gama9*(K*(((sxint(1)-sxint(5))*cos(phi9(i))+(sxint(3)-sxint(7))*sin(phi9(i))-rL-d9(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi9(i))+(sxint(4)-sxint(8))*sin(phi9(i))))*cos(phi9(i))+gama10*(K*(((sxint(1)-sxint(5))*cos(phi10(i))+(sxint(3)-sxint(7))*sin(phi10(i))-rL-d10(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi10(i))+(sxint(4)-sxint(8))*sin(phi10(i))))*cos(phi10(i))+gama11*(K*(((sxint(1)-sxint(5))*cos(phi11(i))+(sxint(3)-sxint(7))*sin(phi11(i))-rL-d11(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi11(i))+(sxint(4)-sxint(8))*sin(phi11(i))))*cos(phi11(i))+gama12*(K*(((sxint(1)-sxint(5))*cos(phi12(i))+(sxint(3)-sxint(7))*sin(phi12(i))-rL-d12(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi12(i))+(sxint(4)-sxint(8))*sin(phi12(i))))*cos(phi12(i))+gama13*(K*(((sxint(1)-sxint(5))*cos(phi13(i))+(sxint(3)-sxint(7))*sin(phi13(i))-rL-d13(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi13(i))+(sxint(4)-sxint(8))*sin(phi13(i))))*cos(phi13(i))+gama14*(K*(((sxint(1)-sxint(5))*cos(phi14(i))+(sxint(3)-sxint(7))*sin(phi14(i))-rL-d14(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi14(i))+(sxint(4)-sxint(8))*sin(phi14(i))))*cos(phi14(i)))/mo;
%axo(i)= -cox*sxint(5)/mo-kox*sxint(5)/mo+(cos(phi1(i))*(K*((sxint(1)-sxint(5))*cos(phi1(i))+(sxint(3)-sxint(7))*sin(phi1(i))-rL-d1(i))+c*((sxint(2)-sxint(6))*cos(phi1(i))+(sxint(4)-sxint(8))*sin(phi1(i))-rL-d1(i))^1.5)+cos(phi2(i))*(K*((sxint(1)-sxint(5))*cos(phi2(i))+(sxint(3)-sxint(7))*sin(phi2(i))-rL-d2(i))+c*((sxint(2)-sxint(6))*cos(phi2(i))+(sxint(4)-sxint(8))*sin(phi2(i))-rL-d2(i))^1.5)+cos(phi3(i))*(K*((sxint(1)-sxint(5))*cos(phi3(i))+(sxint(3)-sxint(7))*sin(phi3(i))-rL-d3(i))+c*((sxint(2)-sxint(6))*cos(phi3(i))+(sxint(4)-sxint(8))*sin(phi3(i))-rL-d3(i))^1.5)+cos(phi4(i))*(K*((sxint(1)-sxint(5))*cos(phi4(i))+(sxint(3)-sxint(7))*sin(phi4(i))-rL-d4(i))+c*((sxint(2)-sxint(6))*cos(phi4(i))+(sxint(4)-sxint(8))*sin(phi4(i))-rL-d4(i))^1.5)+cos(phi5(i))*(K*((sxint(1)-sxint(5))*cos(phi5(i))+(sxint(3)-sxint(7))*sin(phi5(i))-rL-d5(i))+c*((sxint(2)-sxint(6))*cos(phi5(i))+(sxint(4)-sxint(8))*sin(phi5(i))-rL-d5(i))^1.5)+cos(phi6(i))*(K*((sxint(1)-sxint(5))*cos(phi6(i))+(sxint(3)-sxint(7))*sin(phi6(i))-rL-d6(i))+c*((sxint(2)-sxint(6))*cos(phi6(i))+(sxint(4)-sxint(8))*sin(phi6(i))-rL-d6(i))^1.5)+cos(phi7(i))*(K*((sxint(1)-sxint(5))*cos(phi7(i))+(sxint(3)-sxint(7))*sin(phi7(i))-rL-d7(i))+c*((sxint(2)-sxint(6))*cos(phi7(i))+(sxint(4)-sxint(8))*sin(phi7(i))-rL-d7(i))^1.5)+cos(phi8(i))*(K*((sxint(1)-sxint(5))*cos(phi8(i))+(sxint(3)-sxint(7))*sin(phi8(i))-rL-d8(i))+c*((sxint(2)-sxint(6))*cos(phi8(i))+(sxint(4)-sxint(8))*sin(phi8(i))-rL-d8(i))^1.5)+cos(phi9(i))*(K*((sxint(1)-sxint(5))*cos(phi9(i))+(sxint(3)-sxint(7))*sin(phi9(i))-rL-d9(i))+c*((sxint(2)-sxint(6))*cos(phi9(i))+(sxint(4)-sxint(8))*sin(phi9(i))-rL-d9(i))^1.5)+cos(phi10(i))*(K*((sxint(1)-sxint(5))*cos(phi10(i))+(sxint(3)-sxint(7))*sin(phi10(i))-rL-d10(i))+c*((sxint(2)-sxint(6))*cos(phi10(i))+(sxint(4)-sxint(8))*sin(phi10(i))-rL-d10(i))^1.5)+cos(phi11(i))*(K*((sxint(1)-sxint(5))*cos(phi11(i))+(sxint(3)-sxint(7))*sin(phi11(i))-rL-d11(i))+c*((sxint(2)-sxint(6))*cos(phi11(i))+(sxint(4)-sxint(8))*sin(phi11(i))-rL-d11(i))^1.5)+cos(phi12(i))*(K*((sxint(1)-sxint(5))*cos(phi12(i))+(sxint(3)-sxint(7))*sin(phi12(i))-rL-d12(i))+c*((sxint(2)-sxint(6))*cos(phi12(i))+(sxint(4)-sxint(8))*sin(phi12(i))-rL-d12(i))^1.5)+cos(phi13(i))*(K*((sxint(1)-sxint(5))*cos(phi13(i))+(sxint(3)-sxint(7))*sin(phi13(i))-rL-d13(i))+c*((sxint(2)-sxint(6))*cos(phi13(i))+(sxint(4)-sxint(8))*sin(phi13(i))-rL-d13(i))^1.5)+cos(phi14(i))*(K*((sxint(1)-sxint(5))*cos(phi14(i))+(sxint(3)-sxint(7))*sin(phi14(i))-rL-d14(i))+c*((sxint(2)-sxint(6))*cos(phi14(i))+(sxint(4)-sxint(8))*sin(phi14(i))-rL-d14(i))^1.5))/mo;
%以下代入式4.4，求yo的二阶导
ayo(i)=-coy*sxint(8)/mo-koy*sxint(7)/mo+(gama1*(K*(((sxint(1)-sxint(5))*cos(phi1(i))+(sxint(3)-sxint(7))*sin(phi1(i))-rL-d1(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi1(i))+(sxint(4)-sxint(8))*sin(phi1(i))))*sin(phi1(i))+gama2*(K*(((sxint(1)-sxint(5))*cos(phi2(i))+(sxint(3)-sxint(7))*sin(phi2(i))-rL-d2(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi2(i))+(sxint(4)-sxint(8))*sin(phi2(i))))*sin(phi2(i))+gama3*(K*(((sxint(1)-sxint(5))*cos(phi3(i))+(sxint(3)-sxint(7))*sin(phi3(i))-rL-d3(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi3(i))+(sxint(4)-sxint(8))*sin(phi3(i))))*sin(phi3(i))+gama4*(K*(((sxint(1)-sxint(5))*cos(phi4(i))+(sxint(3)-sxint(7))*sin(phi4(i))-rL-d4(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi4(i))+(sxint(4)-sxint(8))*sin(phi4(i))))*sin(phi4(i))+gama5*(K*(((sxint(1)-sxint(5))*cos(phi5(i))+(sxint(3)-sxint(7))*sin(phi5(i))-rL-d5(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi5(i))+(sxint(4)-sxint(8))*sin(phi5(i))))*sin(phi5(i))+gama6*(K*(((sxint(1)-sxint(5))*cos(phi6(i))+(sxint(3)-sxint(7))*sin(phi6(i))-rL-d6(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi6(i))+(sxint(4)-sxint(8))*sin(phi6(i))))*sin(phi6(i))+gama7*(K*(((sxint(1)-sxint(5))*cos(phi7(i))+(sxint(3)-sxint(7))*sin(phi7(i))-rL-d7(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi7(i))+(sxint(4)-sxint(8))*sin(phi7(i))))*sin(phi7(i))+gama8*(K*(((sxint(1)-sxint(5))*cos(phi8(i))+(sxint(3)-sxint(7))*sin(phi8(i))-rL-d8(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi8(i))+(sxint(4)-sxint(8))*sin(phi8(i))))*sin(phi8(i))+gama9*(K*(((sxint(1)-sxint(5))*cos(phi9(i))+(sxint(3)-sxint(7))*sin(phi9(i))-rL-d9(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi9(i))+(sxint(4)-sxint(8))*sin(phi9(i))))*sin(phi9(i))+gama10*(K*(((sxint(1)-sxint(5))*cos(phi10(i))+(sxint(3)-sxint(7))*sin(phi10(i))-rL-d10(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi10(i))+(sxint(4)-sxint(8))*sin(phi10(i))))*sin(phi10(i))+gama11*(K*(((sxint(1)-sxint(5))*cos(phi11(i))+(sxint(3)-sxint(7))*sin(phi11(i))-rL-d11(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi11(i))+(sxint(4)-sxint(8))*sin(phi11(i))))*sin(phi11(i))+gama12*(K*(((sxint(1)-sxint(5))*cos(phi12(i))+(sxint(3)-sxint(7))*sin(phi12(i))-rL-d12(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi12(i))+(sxint(4)-sxint(8))*sin(phi12(i))))*sin(phi12(i))+gama13*(K*(((sxint(1)-sxint(5))*cos(phi13(i))+(sxint(3)-sxint(7))*sin(phi13(i))-rL-d13(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi13(i))+(sxint(4)-sxint(8))*sin(phi13(i))))*sin(phi13(i))+gama14*(K*(((sxint(1)-sxint(5))*cos(phi14(i))+(sxint(3)-sxint(7))*sin(phi14(i))-rL-d14(i)))^1.5+c*((sxint(2)-sxint(6))*cos(phi14(i))+(sxint(4)-sxint(8))*sin(phi14(i))))*sin(phi14(i)))/mo;
yi(i)=sxint(3);
xi(i)=sxint(1);
yo(i)=sxint(7);
xo(i)=sxint(5);
x1=sxint(1); x2=sxint(2); x3=sxint(3); x4=sxint(4); x5=sxint(5); x6=sxint(6);
x7=sxint(7); x8=sxint(8);

%循环进度显示
  jindu=i/amount*100;
  %if mod(jindu,5)<1;
  disp([num2str(jindu) '%']);
 % end  
  end
% ac=ayo-axo;     
% %Xt=ayo;
% Xt=axo(4097:6144);
% dt=1/Fs;
% N=length(Xt);
% f=[0:N/2-1]/(N*dt);
% t=[4097-1:6144-1]/Fs;       %%%%%%%%%%%%与Xt的取值关联
% Xt1=hamming(N)'.*Xt(1:N);     
% Xf=fft(Xt1);
% Xf1=abs(Xf)*2/N;

%load('F:\python file\rotor_imputation\simu_dataset\simu_axy.mat')
d1 = real(ayo(2000:end));d1_ = ones(3000,512);
for i =1:4000
    d1_(i,:) = d1(5*(i-1)+1:512+5*(i-1));
end
% save(['D:\Work in EPFL\diffusion_model\AE_Diffusion\AE_data_ims\data_bo.mat','d1_')
%% =======================================完整信号作图=============================================
set(0,'defaultfigurecolor','w');%设置背景板的颜色为白色
ax_=axo(10241:end);
ay_=ayo(10241:end);
N=length(ax_);
dt=1/Fs;
t=[0:N-1]*dt;

figure(1)
set(gca,'fontname','times new roman','fontsize',12)
subplot(2,1,1);
ylabel('\itacceleration\rm(m/s^2)','fontname','times new roman','fontsize',15,'FontWeight','bold');
xlabel('\ittime\rm(s)','fontname','times new roman','fontsize',15,'FontWeight','bold');
title('\itoriginal signal','fontname','times new roman','fontsize',15,'FontWeight','bold');
hold on
plot(t,ax_,'m-');
subplot(2,1,2);
ylabel('\itacceleration\rm(m/s^2)','fontname','times new roman','fontsize',15,'FontWeight','bold');
xlabel('\ittime\rm(s)','fontname','times new roman','fontsize',15,'FontWeight','bold');
title('\itoriginal signal','fontname','times new roman','fontsize',15,'FontWeight','bold');
hold on
plot(t,ay_,'c-');

y_fft=fft(ax_);
p2_y_fft=abs(y_fft/N);
p1_y_fft=p2_y_fft(1:N/2+1);
p1_y_fft(2:end-1)=2*p2_y_fft(2:N/2);
f=(Fs/N)*(0:N/2);
y_fft_=fft(ay_);
p2_y_fft_=abs(y_fft_/N);
p1_y_fft_=p2_y_fft_(1:N/2+1);
p1_y_fft_(2:end-1)=2*p2_y_fft_(2:N/2);
%f=(Fs/N)*(0:N/2);

figure(2)
set(gca,'fontname','times new roman','fontsize',12)
subplot(2,1,1)
ylabel('\itAmplitude\rm(m/Hz*s^2)','fontname','times new roman','fontsize',15);
xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',15);
title('\itx- signal','fontname','times new roman','fontsize',15);
hold on
plot(f,p1_y_fft,'r')

subplot(2,1,2)
ylabel('\itAmplitude\rm(m/Hz*s^2)','fontname','times new roman','fontsize',15);
xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',15);
title('\ity- signal','fontname','times new roman','fontsize',15);
hold on
plot(f,p1_y_fft_,'r')


figure(3)
hax_=abs(hilbert(ax_));
hay_=abs(hilbert(ay_));
y_fft=fft(hax_);
p2_y_fft=abs(y_fft/N);
p1_y_fft=p2_y_fft(1:N/2+1);
p1_y_fft(2:end-1)=2*p2_y_fft(2:N/2);
f=(Fs/N)*(0:N/2);
y_fft_=fft(hay_);
p2_y_fft_=abs(y_fft_/N);
p1_y_fft_=p2_y_fft_(1:N/2+1);
p1_y_fft_(2:end-1)=2*p2_y_fft_(2:N/2);
%f=(Fs/N)*(0:N/2);
% % set(gca,'fontname','times new roman','fontsize',12)
subplot(2,1,1)
ylabel('\itAmplitude\rm(m/(s^2\cdotHz))','fontname','times new roman','fontsize',12,'FontWeight','bold');
xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',12,'FontWeight','bold');
%xticks([0,1.4],'fontsize',10);yticks([-0.05,0,0.05],'fontsize',10);
title('\itx- signal','fontname','times new roman','fontsize',15,'FontWeight','bold');
hold on
plot(f,p1_y_fft,'r');
hold on
axis([0,4500,0,0.020]);
set(gca,'fontname','times new roman','fontsize',12);
subplot(2,1,2)
ylabel('\itAmplitude\rm(m/(s^2\cdotHz))','fontname','times new roman','fontsize',12,'FontWeight','bold');
xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',12,'FontWeight','bold');
title('\ity- signal','fontname','times new roman','fontsize',15,'FontWeight','bold');
hold on
plot(f,p1_y_fft_,'r');
hold on
axis([0,4500,0,1.50]);
set(gca,'fontname','times new roman','fontsize',12);

%% =======================================不完整信号作图=============================================
% %ax_=axo(1967:end);
% %ay_=ayo(1967:end);
% ax_m=[ax_(1:2048),zeros(1,8192-2048),ax_(8193:end)];%intact part-x
% ax_m_1=ax_(1:2048);
% ax_m_2=ax_(8193:end);
% ax_e=ax_(2049:8192);%missing part-x
% 
% ay_m=[ay_(1:2048),zeros(1,8192-2048),ay_(8193:end)];%intact part-x
% ay_m_1=ay_(1:2048);
% ay_m_2=ay_(8193:end);
% ay_e=ay_(2049:8192);%missing part -y
% N=length(ax_);
% N_=length(ay_);
% dt=1/Fs;
% tx=[0:N-1]*dt;
% ty=[0:N_-1]*dt;
% 
% figure(4)
% 
% subplot(2,1,1);
% ylabel('\itacceleration\rm(m/s^2)','fontname','times new roman','fontsize',19,'FontWeight','bold');
% xlabel('\ittime\rm(s)','fontname','times new roman','fontsize',19,'FontWeight','bold');
% title('\ithorizontal signal(x)','fontname','times new roman','fontsize',19,'FontWeight','bold');
% hold on
% %plot(tx,ax_m,'m-');
% plot(tx(1:2048),ax_m_1,'m');
% hold on
% plot(tx(2049:8192),ax_e,'b');
% hold on
% plot(tx(8193:end),ax_m_2,'m');
% hold on
% axis([0,1.4,-0.1,0.1]);
% set(gca,'fontname','times new roman','fontsize',15);
% subplot(2,1,2);
% ylabel('\itacceleration\rm(m/s^2)','fontname','times new roman','fontsize',19,'FontWeight','bold');
% xlabel('\ittime\rm(s)','fontname','times new roman','fontsize',19,'FontWeight','bold');
% title('\itvertical signal(y)','fontname','times new roman','fontsize',19,'FontWeight','bold');
% hold on
% plot(ty(1:2048),ay_m_1,'c','linewidth',0.7);
% hold on
% plot(ty(2049:8192),ay_e,'b','linewidth',0.7);
% hold on
% plot(ty(8193:end),ay_m_2,'c','linewidth',0.7);
% hold on
% axis([0,1.4,-20,20]);
% set(gca,'fontname','times new roman','fontsize',15);
% y_fft=fft(ax_m);
% p2_y_fft=abs(y_fft/N);
% p1_y_fft=p2_y_fft(1:N/2+1);
% p1_y_fft(2:end-1)=2*p2_y_fft(2:N/2);
% fx=(Fs/N)*(0:N/2);
% fy=(Fs/N_)*(0:N_/2); 
% y_fft_=fft(ay_m);
% p2_y_fft_=abs(y_fft_/N_);
% p1_y_fft_=p2_y_fft_(1:N_/2+1);
% p1_y_fft_(2:end-1)=2*p2_y_fft_(2:N_/2);
% %f=(Fs/N)*(0:N/2);
% 
% figure(5)
% set(gca,'fontname','times new roman','fontsize',12)
% subplot(2,1,1)
% ylabel('\itAmplitude','fontname','times new roman','fontsize',15);
% xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',15);
% title('\itx- signal','fontname','times new roman','fontsize',15);
% hold on
% plot(fx,p1_y_fft,'r')
% subplot(2,1,2)
% ylabel('\itAmplitude','fontname','times new roman','fontsize',15);
% xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',15);
% title('\ity- signal','fontname','times new roman','fontsize',15);
% hold on
% plot(fy,p1_y_fft_,'r')
% 
% figure(6)
% hax_m=abs(hilbert(ax_m));
% hay_=abs(hilbert(ay_m));
% y_fft=fft(hax_m);
% p2_y_fft=abs(y_fft/N);
% p1_y_fft=p2_y_fft(1:N/2+1);
% p1_y_fft(2:end-1)=2*p2_y_fft(2:N/2);
% y_fft_=fft(hay_);
% p2_y_fft_=abs(y_fft_/N_);
% p1_y_fft_=p2_y_fft_(1:N_/2+1);
% p1_y_fft_(2:end-1)=2*p2_y_fft_(2:N_/2);
% %f=(Fs/N)*(0:N/2);
% % set(gca,'fontname','times new roman','fontsize',12)
% subplot(2,1,1)
% ylabel('\itAmplitude\rm(m/(s^2\cdotHz))','fontname','times new roman','fontsize',12,'FontWeight','bold');
% xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',12,'FontWeight','bold');
% title('\itx- signal','fontname','times new roman','fontsize',15,'FontWeight','bold');
% hold on
% plot(fx,p1_y_fft,'r');
% hold on
% axis([0,4500,0,0.015]);
% set(gca,'fontname','times new roman','fontsize',12);
% subplot(2,1,2)
% ylabel('\itAmplitude\rm(m/(s^2\cdotHz))','fontname','times new roman','fontsize',12,'FontWeight','bold');
% xlabel('\itfrequency\rm(Hz)','fontname','times new roman','fontsize',12,'FontWeight','bold');
% title('\ity- signal','fontname','times new roman','fontsize',15,'FontWeight','bold');
% hold on
% plot(fy,p1_y_fft_,'r');
% hold on
% axis([0,4500,0,1.2]);
% set(gca,'fontname','times new roman','fontsize',12);
% %xticks('fontsize',10);yticks('fontsize',10)
%%
% figure(1)
% patch([0.3653 0.3653 0.4369 0.4369],[-54.8 -5.4 -5.4 -54.8],'r');
% patch([0.5736 0.5736 0.4369 0.4369],[-54.8 -5.4 -5.4 -54.8],'r');
% patch([1.68 1.68 1.719 1.719],[-54.8 -5.4 -5.4 -54.8],'r');
% patch([1.719 1.719 1.764 1.764],[-54.8 -5.4 -5.4 -54.8],'r');
% %patch([0.533 0.533 0.616 0.616],[-49.8 -0.4 -0.4 -49.8],'r');
% %patch([0.616 0.616 0.6786 0.6786],[-49.8 -0.4 -0.4 -49.8],'r');
% %patch([1.743 1.743 1.791 1.791],[-49.8 -0.4 -0.4 -49.8],'r');
% %patch([1.791 1.791 1.839 1.839],[-49.8 -0.4 -0.4 -49.8],'r');
% hold on;
% plot(f/1000,10*log10(Xf1(1:N/2)));
% ylabel('Acceleration(dB)');
% xlabel('Frequency(kHz)');  
% axis([0.3 2.1 -50 0]);
% set(gcf,'Units','centimeters','Position',[10 10 9.2 8.5]);%设置图片大小为7cm×7cm
% set(gcf,'color','white');%设定边框背景为白色
% set(gca,'looseInset',[0 0 0 0.05]);%设定边框距离为0
% annotation('textbox',...
%     [0 1 0 0],...
%     'String',{'(b)'},...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% hold off;
% figure(2)
% %subplot(211);
% plot(t,Xt);
% xlabel('时间t/s'); ylabel('振幅');title('时域波形图');
% %axis([0 0.035 -2 2]);
% %grid on;
% %subplot(212);
% figure(3)
% plot(f,Xf1(1:N/2));
% xlabel('频率f/Hz'); ylabel('幅值');title('频谱图');
% %axis([0 8000 0 120]);
% grid on;
% 
% %希尔伯特包络解调
% figure(4)
% y3=abs(hilbert(Xt));
% y4=abs(fft(y3-mean(y3)))*2/N;
% plot([0:N-1]/(N*dt),y4);
% title('Hilbert demodulation spectrum');xlabel('Frequency/Hz'); ylabel('Amplitude/(m/s^2)');
% %axis([0 500 0 8]);
%      
% % figure (5);
% % %plot(tmms(1/(14*fc*tmm)+1:3/(14*fc*tmm)+1),axo(1/(14*fc*tmm)+1:3/(14*fc*tmm)+1),tmms(1/(14*fc*tmm)+1:3/(14*fc*tmm)+1)+2/(14*fc),axo(1/(14*fc*tmm)+1:3/(14*fc*tmm)+1));
% % plot(tmms(4.5/(14*fc*tmm)+1:6144*2),ayo(4.5/(14*fc*tmm)+1:6144*2));
% % % plot(phic*180/pi, kxx);
% % ylabel('加速度(m/s^2)');
% % xlabel('时间(s)');  
% % axis([0.0248*4.5 0.187 -3.5 6]);
% % set(gcf,'Units','centimeters','Position',[10 10 9.2*2 8.5]);%设置图片大小为7cm×7cm
% % set(gcf,'color','white');%设定边框背景为白色
% % set(gca,'looseInset',[0 0 0 0.05]);%设定边框距离为0
% % annotation('textbox',...
% %     [0 1 0 0],...
% %     'String',{'(a)'},...
% %     'FitBoxToText','off',...
% %     'EdgeColor',[1 1 1]);