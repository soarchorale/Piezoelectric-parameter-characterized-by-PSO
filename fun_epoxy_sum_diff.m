%��֤PSO�㷨��SA �㷨�Ľ����ȡ������ע��aa֮��Ĵ��� ��ע��AA
%���㸴�ϲ��ϵ���Ͻ����ȡ������ע��A֮��Ĵ��� ��ע�� aa
function [sum_diff_Z]=  fun_epoxy_sum_diff(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c )

 [R_cal,X_cal,Z_cal,phase_cal,ff]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );%����ֵ

%%%����ֵ��ʵ��[Z_m,phase_m,ff]�ı�����ֵΪ��������ֵ
%**************************** a
epl33s=5.556e-9;
tan_ct_e=0.08;
k_t=0.602;
tan_ct_k= - 0.0129;
c33_D=6.926e+10;
tan_ct_c=0.022;%SA��ԭʼ����
%**************************** a

[R_m,X_m,Z_m,phase_m,ff]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );

%**************************** A
% f=load('2f.txt')/1000;
% z=load('2z.txt');
% r=load('2r.txt');
% x=load('2x.txt');
% phase=load('2p.txt');%���ϲ���
%**************************** A

%**************************** A
% start_point=240;

% for i=1:93
%     ff(i)=f(start_point+i-1);
%     Z_m(i)=z(start_point+i-1);
%     X_m(i)=x(start_point+i-1);
%     R_m(i)=r(start_point+i-1);
%     phase_m(i)=phase(start_point+i-1);
% end%���ϲ���
%**************************** A
sum_diff_Z=0;
ff_num=length(ff);
for j=1:ff_num
    sum_diff_Z = sum_diff_Z+(R_cal(j) - R_m(j) )^2+(X_cal(j) - X_m(j) )^2;% + (Z_cal(j) - Z_m(j) )^2
end
end
