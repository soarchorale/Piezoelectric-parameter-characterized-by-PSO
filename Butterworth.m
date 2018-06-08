%实现了对频率为20和200Hz单频叠加cos信号的低通滤波，使输出仅含有20Hz分量

clear all


fs=1200;        %采样频率
N=120;     % N/fs 秒数据
n=0:N-1;
t=n/fs;     %时间
fL=20; 
fH=200;    

% f_m=xlsread('pzt_round20_1M.xlsx','A2:A249');
% z_m=xlsread('pzt_round20_1M.xlsx','F2:F249');
% phase_m=xlsread('pzt_round20_1M.xlsx','G2:G249');

% f_m=load('pin_pmn_30pt011_f.txt');
% z_m=load('pin_pmn_30pt011_z.txt');
% phase_m=load('pin_pmn_30pt011_p.txt');

% f_m=load('Mn_pinpmnpt_big_f.txt');
% z_m=load('Mn_pinpmnpt_big_z.txt');
% phase_m=load('Mn_pinpmnpt_big_p.txt');

f_m=load('pzt_roundplate_10k_3M_f.txt');
z_m=load('pzt_roundplate_10k_3M_z.txt');
phase_m=load('pzt_roundplate_10k_3M_p.txt');

t=f_m;
s=z_m;

% sL=cos(2*pi*fL*t);    
% sH=cos(2*pi*fH*t);   
% s=sL+sH;

% s_in_f=fft(s);
% i=1:250;
% 
% 
% plot(i,s_in_f(i),'r');
% 
figure(1);
plot(t,log(s),'b');
title('输入信号'); 
xlabel('t/s');
ylabel('幅度');

%设计低通滤波器：
Wp = 65/fs; Ws = 500/fs;    %截止频率50Hz,阻带截止频率100Hz,采样频率200Hz
[n,Wn] = buttord(Wp,Ws,1,50);   %阻带衰减大于50db,通带纹波小于1db
%估算得到Butterworth低通滤波器的最小阶数N和3dB截止频率Wn
[a,b]=butter(n,Wn);     %设计Butterworth低通滤波器
% [h,f]=freqz(a,b,'whole',fs);        %求数字低通滤波器的频率响应
% f=(0:length(f)-1)'*fs/length(f);    %进行对应的频率转换
% figure(2);
% plot(f,abs(h));     %绘制Butterworth低通滤波器的幅频响应图
% title('巴氏低通滤波器'); 
% grid;
sF=filter(a,b,s);       %叠加函数s经过低通滤波器以后的新函数
figure(3);
plot(t,log(sF));     %绘制叠加函数S经过低通滤波器以后的时域图形
xlabel('t/s');
ylabel('幅度');
% 
% SF=fft(sF,N);       %对叠加函数S经过低通滤波器以后的新函数进行256点的基—2快速傅立叶变换
% mag=abs(SF);        %求幅值
% f=(0:length(SF)-1)'*fs/length(SF);    %进行对应的频率转换
% figure(4);
% plot(f,mag);        %绘制叠加函数S经过低通滤波器以后的频谱图
% title('低通滤波后的频谱图');

% figure(5);
% plot(t,sL);        
% title('低频信号');
