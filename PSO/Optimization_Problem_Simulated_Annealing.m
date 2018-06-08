%验证PSO算法和SA 算法的结果请取消所有注释aa之间的代码 并注释AA
%计算复合材料的拟合结果请取消所有注释A之间的代码 并注释 aa
clc;
clear;
%%初始赋值
T0=1e+7;%0.1
T_min=1e-2;
r=0.98;%退火速率
q=0;
tl=100;     %在每个温度下的最大迭代次数
Iter=450;

 % Bound=[1.6668 	7.2228;0.0240 	0.104;  0.1806 	0.7826; -0.0039 	-0.01677;  2.0778 	9.0038  ; 0.0066 	0.0286]	%粒子的坐标范围	
% %粒子的坐标范围:数1   数2         数3     数4      数5     数6
%Bound=[1.6668 7.2228;0.0240 	0.104;  0.1806 	0.7826; -0.0039 	-0.01677;  2.0778 9.0038  ; 0.0066 	0.0286];	%30%
  Bound=[5.0004	6.1116;0.072	0.088;0.5418	0.6622;-0.01161 -0.01419;6.2334	7.6186;0.0198	0.0242];%10%
  %Bound=[4.4448	6.6672;0.064	0.096;0.4816	0.7224;-0.01032 -0.01548;5.5408	8.3112;0.0176	0.0264];%20%
%Bound=[1  8; 0  0.1;    0.1  1;  -0.1  -0.005;  2  10; 0.01 0.1]%原始区间
% Bound=[6 8;0 0.1;0.4 0.6;  -0.1 0.1; 3  4; 0 0.05];%凑区间
% Bound=[5.44488	5.66712;0.0784	0.0816;0.58996	0.61404;-0.012642	-0.013158;6.78748	7.06452;0.02156	0.02244];
%定义粒子上下边界
lowerbound=Bound(:,1);
upperbound=Bound(:,2);

	%初始化各粒子初始位置，在有效范围内随机选数	
    population=rand*(Bound(:,2)-Bound(:,1)) + Bound(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%population(1,:)=round(population(1,;)*1000)/1000;
%population(2,;)=round(population(2,;)*1000)/1000;

epl33s  =  population(1) *1e-9;
tan_ct_e  =  population(2);
k_t  =  population(3);
tan_ct_k  =  population(4);
c33_D  =  population(5)*1e+10;
tan_ct_c  =  population(6);

%*************************
%此段可验证文中的初值选择的拟合结果
% epl33s  =  3*1e-9;
% tan_ct_e  =  -0.5;
% k_t  =  0.5;
% tan_ct_k  = 0.5;
% c33_D  =  5*1e+10;
% tan_ct_c  =  0.9;
%****************************

%Var0=rand(1,6)*(-6)+[1 3]; %初始变量赋值，随机在区间[-3,3]中取值
Z0= fun_epoxy_sum_diff(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );
Var_best=population;
Z_best=Z0;
T=T0;
Time=ceil(double(solve(['0.1*(0.98)^x=',num2str(T_min)])));
obj=zeros(Time,1);
i=0;
M=zeros(6,tl);
N=zeros(1,tl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%读取样本数据
%*************************** a
epl33s=5.556e-9;
tan_ct_e=0.08;
k_t=0.602;
tan_ct_k= - 0.0129;
c33_D=6.926e+10;
tan_ct_c=0.022;
[rr_m,xx_m,zz_m,phase_m,ff_m]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );%原始数据
%*************************** a

%*************************** A
% f_m=load('2f.txt')/1000;
% z_m=load('2z.txt');
% phase_m=load('2p.txt');
% x_m=load('2x.txt');
% r_m=load('2r.txt');%复合材料

% ff_m=f_m(240:333);
% zz_m=z_m(240:333);
% xx_m=x_m(240:333);
% rr_m=r_m(240:333);%复合材料
%*************************** A

% d0=zeros(1,Iter);
% f0=zeros(6,Iter);
% Best=zeros(1,Iter);
while T>T_min && i<=Iter && Z_best>100
    i=i+1;
    u(i)=i;
%     if i>300
        T=T*r;
%     elseif i<=300
%         q=(T-T_min)/(T_min*T*tl);
%         T=T/(1+T*q);
%     end

    for k=1:tl
        Var_new=get_newanswer(Var_best);
        Z_new=fun_epoxy_sum_diff(Var_new(1)*1e-9,Var_new(2),Var_new(3),Var_new(4),Var_new(5)*1e+10,Var_new(6));
        d=Z_new-Z_best;
        if d<0
           Z_best=Z_new;
           Var_best=Var_new;
        else
            if exp(-d/T)>rand
                Z_best=Z_new;
                Var_best=Var_new; 
            end
        end
         M(:,k)=Var_best;
         N(k)=Z_best;
    end
    [d0(i),n]=min(N);      %找出当前温度下最优解
    f0(:,i)=M(:,n);
    if i==1||d0(i)<obj(i-1)             
       obj(i)=d0(i);          %如果当前温度下最优解小于上一温度，则记录当前解
    else
       obj(i)=obj(i-1);    %如果当前温度下最优解大于上一温度，则记录上一温度解
    end 
    [Best(i),q]=min(d0);
    Var_best=f0(:,q);
%     Var_best=M(:,n);
%     Best(i)=Z_best;

%      plot(obj,'ro-');
    figure(1);
%   subplot(2,1,1);
    plot(u,d0,'bo-',u,Best,'ro-');legend('current solutions','best solutions','Location','NorthWest');
    xlabel({'Generation','(a)'},'FontSize',14, 'FontWeight','bold','Color','k');
    ylabel('Curent Solutions');
    title({'Optimization process'},'FontSize',14, 'FontWeight','bold','Color','k');
    text(0.7,0.95,['d0 = ', num2str(d0(i))],'Units','normalized');
    text(0.7,0.85,['Best = ', num2str(Best(i))],'Units','normalized');
    set(gca,'linewidth',1.5);
    set(gca,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold')
%     plot(Best,'ro-');
%     xlabel('迭代次数');
%       ylabel('最优解');
%     title('优化过程');
%     subplot(3,1,2);
%     plot(Best,'ro-');
%     text(0.5,0.95,['Best = ', num2str(Best(i))],'Units','normalized');
         
    [R_c,X_c,Z_c,phase_c,ff]=  fun_epoxy(Var_best(1)*1e-9,Var_best(2),Var_best(3),Var_best(4),Var_best(5)*1e+10,Var_best(6));
%   subplot(2,1,2);
    figure(2)
    plot(ff,log(Z_c),'mo-',ff_m,log(zz_m),'r');
    legend('Z_f Fitting Curve','Z_f Data Curve','Location','NorthWest');
    xlabel({'Frequency (kHz)','(b)'},'FontSize',14, 'FontWeight','bold','Color','k');
    ylabel({'Resistance (Ω)'},'FontSize',14);
    title({'Fitting Result'},'FontSize',14, 'FontWeight','bold','Color','k');
    set(gca,'linewidth',1.5,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold');
    figure(3)
    plot(ff,phase_c,'bo-',ff_m,phase_m,'g');
%   figure(3)
%   semilogy(ff,R_c,'bo-',ff_m,rr_m,'g',ff,X_c,'mo-',ff_m,xx_m,'r');legend('R_f Fitting Curve','R_f Data Curve','X_f Fitting Curve','X_f Data Curve','Location','NorthWest');xlabel({'Frequency (kHz)','(c)'},'FontSize',14, 'FontWeight','bold','Color','k');ylabel({'Resistance (Ω)'},'FontSize',14, 'FontWeight','bold','Color','k');title({'Fitting Result'},'FontSize',14, 'FontWeight','bold','Color','k');
   
%     subplot(2,1,1);
%     plot(ff,log(R_c),'ro-',ff_m,log(rr_m),'b');legend('fitting curve','data curve','Location','NorthWest');xlabel('Frequency(kHz)');ylabel('R(f)');title('R_f Fitting Result');title({'R_f fitting result'},'FontSize',14, 'FontWeight','bold','Color','k');
%     set(gca,'linewidth',1.5);
%     set(gca,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold');
%     subplot(2,1,2);
%     hold on;
%     plot(ff,log(X_c),'ro-',ff_m,log(xx_m),'b');legend('fitting curve','data curve','Location','SouthWest');xlabel({'Frequency (kHz)','(b)'},'FontSize',14, 'FontWeight','bold','Color','k');ylabel('X(f)');title('X_f Fitting Result');title({'X_f fitting result'},'FontSize',14, 'FontWeight','bold','Color','k');
%     set(gca,'linewidth',1.5);
%     set(gca,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold');
%     hold off;
    drawnow;
end

i
Var_best
min(Best(i))







            
        


