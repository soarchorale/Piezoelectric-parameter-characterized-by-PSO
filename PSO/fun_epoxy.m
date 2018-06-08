%��֤PSO�㷨��SA �㷨�Ľ����ȡ������ע��aa֮��Ĵ��� ��ע��AA
%���㸴�ϲ��ϵ���Ͻ����ȡ������ע��A֮��Ĵ��� ��ע�� aa

function [R_f,X_f,Z_f,phase_f,ff]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c )

% clear all
% epl33s=5.556e-9;
% tan_ct_e=0.08;
% k_t=0.602;
% tan_ct_k= - 0.0129;
% c33_D=6.926e+10;
% tan_ct_c=0.022;

%****************************  a
%��֤PSO��SA�㷨ʱȡ������ע��
t=1.6e-3;
rou=5395;
tan_ct_m=tan_ct_c;
A=pi*(17.2e-3/2)^2;%SA��ԭʼ����
%**************************** a

%****************************  A
%���ѹ�縴�ϲ���ʱȡ������ע��
% t=3.264e-3;
% rou=3.0715e+03;
% tan_ct_m=tan_ct_c;
% A=(12.75e-3)^2*pi;%���ϲ��ϲ���
%************************** A

%**************************** A
%���ϲ����迹�ײ�������
% f_m=load('2f.txt')/1000;
% z_m=load('2z.txt');
% phase_m=load('2p.txt');%���ϲ���
%**************************** A

%**************************** A
%����������ʼ�����ֹ��
% start_point=240;
% check_point=93;
% ff=zeros(1,check_point);
% zz_m=zeros(1,check_point);
% pphase_m=zeros(1,check_point);
% R_f=zeros(1,check_point);
% X_f=zeros(1,check_point);
% Z_f=zeros(1,check_point);
% phase_f=zeros(1,check_point);%���ϲ���
%**************************** A

%**************************** a
for i=1:410
    f0=0.8e+6;
    ff(i)=f0/1000 + i ;   % 1M Hz
    f=f0 + i*1000;%SA��ԭʼ����
%**************************** a

%**************************** A
% for i= 1:check_point
%     ff(i)=f_m(start_point+i-1);
%     zz_m(i)=z_m(start_point+i-1);
%     pphase_m(i)=phase_m(start_point+i-1);
%     f=f_m(start_point+i-1)*1e+3;%���ϲ���
%**************************** A

M=pi*f*t * sqrt(rou/c33_D) *(1- 3* tan_ct_m^2 /8);

N=  (pi*f*t *tan_ct_m/2)* sqrt(rou/c33_D) *(1- 5* tan_ct_m^2 /16);

T= (tanh(N)* sec(M)^2*(1-tan_ct_k^2) - 2* tan_ct_k* tan(M) *sech(N)^2) *(M - N* tan_ct_e);

U= ( tan(M)* sech(N)^2*(1-tan_ct_k^2) + 2* tan_ct_k*tanh(N)* sec(M)^2) *(N + M* tan_ct_e);

F=( tanh(N)* sec(M)^2*(1-tan_ct_k^2) - 2* tan_ct_k*tan(M)* sech(N)^2) *(N + M* tan_ct_e);

G=(tan(M)* sech(N)^2*(1-tan_ct_k^2) + 2* tan_ct_k* tanh(N) *sec(M)^2) *(M - N* tan_ct_e);

S= (M^2 + N^2)* (1+ tan_ct_e^2) * ( 1+ tan(M)^2 * tanh(N)^2);

R_f(i)= t*tan_ct_e /(2*pi*f * epl33s *A *(1+ tan_ct_e^2)) + t * k_t^2 *(T-U)/(2*pi*f * epl33s *A * S);

X_f(i)= - t /(2*pi*f * epl33s *A *(1+ tan_ct_e^2)) + t * k_t^2 *(F+G)/(2*pi*f * epl33s *A * S);

Z_f(i)= sqrt (R_f(i)^2 + X_f(i)^2 );

phase_f(i)= atan (X_f(i) / R_f(i)) * 180/pi;

end
  
% figure(1);
%plot(ff,log(Z_f),'b',ff,log(zz_m),'r');legend('���','����');
% plot(ff,log(R_f),'g',ff,log(X_f),'r');
% figure(2);
% plot(ff,phase_f,'b',ff,pphase_m,'r');legend('���','����');
