%ʵ���˶�Ƶ��Ϊ20��200Hz��Ƶ����cos�źŵĵ�ͨ�˲���ʹ���������20Hz����

clear all


fs=1200;        %����Ƶ��
N=120;     % N/fs ������
n=0:N-1;
t=n/fs;     %ʱ��
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
title('�����ź�'); 
xlabel('t/s');
ylabel('����');

%��Ƶ�ͨ�˲�����
Wp = 65/fs; Ws = 500/fs;    %��ֹƵ��50Hz,�����ֹƵ��100Hz,����Ƶ��200Hz
[n,Wn] = buttord(Wp,Ws,1,50);   %���˥������50db,ͨ���Ʋ�С��1db
%����õ�Butterworth��ͨ�˲�������С����N��3dB��ֹƵ��Wn
[a,b]=butter(n,Wn);     %���Butterworth��ͨ�˲���
% [h,f]=freqz(a,b,'whole',fs);        %�����ֵ�ͨ�˲�����Ƶ����Ӧ
% f=(0:length(f)-1)'*fs/length(f);    %���ж�Ӧ��Ƶ��ת��
% figure(2);
% plot(f,abs(h));     %����Butterworth��ͨ�˲����ķ�Ƶ��Ӧͼ
% title('���ϵ�ͨ�˲���'); 
% grid;
sF=filter(a,b,s);       %���Ӻ���s������ͨ�˲����Ժ���º���
figure(3);
plot(t,log(sF));     %���Ƶ��Ӻ���S������ͨ�˲����Ժ��ʱ��ͼ��
xlabel('t/s');
ylabel('����');
% 
% SF=fft(sF,N);       %�Ե��Ӻ���S������ͨ�˲����Ժ���º�������256��Ļ���2���ٸ���Ҷ�任
% mag=abs(SF);        %���ֵ
% f=(0:length(SF)-1)'*fs/length(SF);    %���ж�Ӧ��Ƶ��ת��
% figure(4);
% plot(f,mag);        %���Ƶ��Ӻ���S������ͨ�˲����Ժ��Ƶ��ͼ
% title('��ͨ�˲����Ƶ��ͼ');

% figure(5);
% plot(t,sL);        
% title('��Ƶ�ź�');
