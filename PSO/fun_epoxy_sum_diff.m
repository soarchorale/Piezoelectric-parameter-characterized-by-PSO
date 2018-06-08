%验证PSO算法和SA 算法的结果请取消所有注释aa之间的代码 并注释AA
%计算复合材料的拟合结果请取消所有注释A之间的代码 并注释 aa
function [sum_diff_Z]=  fun_epoxy_sum_diff(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c )

 [R_cal,X_cal,Z_cal,phase_cal,ff]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );%迭代值

%%%测量值：实际[Z_m,phase_m,ff]的变量赋值为仪器测量值
%**************************** a
epl33s=5.556e-9;
tan_ct_e=0.08;
k_t=0.602;
tan_ct_k= - 0.0129;
c33_D=6.926e+10;
tan_ct_c=0.022;%SA中原始数据
%**************************** a

[R_m,X_m,Z_m,phase_m,ff]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );

%**************************** A
% f=load('2f.txt')/1000;
% z=load('2z.txt');
% r=load('2r.txt');
% x=load('2x.txt');
% phase=load('2p.txt');%复合材料
%**************************** A

%**************************** A
% start_point=240;

% for i=1:93
%     ff(i)=f(start_point+i-1);
%     Z_m(i)=z(start_point+i-1);
%     X_m(i)=x(start_point+i-1);
%     R_m(i)=r(start_point+i-1);
%     phase_m(i)=phase(start_point+i-1);
% end%复合材料
%**************************** A
sum_diff_Z=0;
ff_num=length(ff);
for j=1:ff_num
    sum_diff_Z = sum_diff_Z+(R_cal(j) - R_m(j) )^2+(X_cal(j) - X_m(j) )^2;% + (Z_cal(j) - Z_m(j) )^2
end
end
