%验证PSO算法和SA 算法的结果请取消所有注释aa之间的代码 并注释AA
%计算复合材料的拟合结果请取消所有注释A之间的代码 并注释 aa
clear all;
clc;

flag=0;		    %停止程序标志
oldbestval=0;	%记录旧的适应度值	
samecounter=0; 	%记录得到相同适应度值的迭代次数
iteration = 0;	%迭代次数
MaxIter =128;  	%最大迭代次数
PopSize=96;   	%种群大小 （加大）
c1 = 0.85;      %学习因子原0.65
c2 = 0.85;      %学习因子    加
w=0.65; 		%惯性因子原0.75	降
%k=0.729;       %收敛因子
%粒子的坐标范围:数1   数2         数3     数4      数5     数6
%Bound=[1.6668 7.2228;0.0240 	0.104;  0.1806 	0.7826; -0.0039 	-0.01677;  2.0778 9.0038  ; 0.0066 	0.0286]	%30%
%Bound=[5.0004	6.1116;0.072	0.088;0.5418	0.6622;-0.01161 -0.01419;6.2334	7.6186;0.0198	0.0242]%10%
Bound=[4.4448	6.6672;0.064	0.096;0.4816	0.7224;-0.01032 -0.01548;5.5408	8.3112;0.0176	0.0264]%20%
%Bound=[1  8; 0  0.1;    0.1  1;  -0.1  -0.005;  2  10; 0.01 0.1]%原始区间
%Bound=[0 10;0 1;0.2 0.9;  -0.1 0.5; 0  20; 0 0.1];%最佳结果
%Bound=[6.3414	7.7506;0.0567	0.0693;0.45657	0.55803; -0.0396	-0.0484; 3.15288	3.85352;0.03519	0.04301];%复合材料10%
% % 当PopSize=6   	%种群大小，可产生population 6维数，每列为一组数，从左到右共6组，
% % population = 
% %     3.5600    5.2230    1.5830    9.8950    6.2450    4.8110
% %     0.1520    0.1330    0.1430    0.1230    0.1580    0.1760
% %    25.2982   26.4053   22.0907   23.7982   27.8333   26.8085
% %     0.2461    0.2568    0.2794    0.2059    0.2603    0.2050
% %    34.1537   33.0500   38.7437   30.1501   37.6795   39.7084
% %     0.3990    0.3789    0.3439    0.3498    0.3214    0.3643
% % 对应 fvalue为按行输入，共6个数

Ndim = length(Bound);  %空间维数

for i=1:PopSize			 %定义粒子上下边界
    lowerbound(:,i)=Bound(:,1);
    upperbound(:,i)=Bound(:,2);
end

for i=1:Ndim		%初始化各粒子初始位置，在有效范围内随机选数	
    population(i,:)=rand(1, PopSize)*(Bound(i,2)-Bound(i,1)) + Bound(i,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population(1,:)=round(population(1,:)*1000)/1000;
% population(2,:)=round(population(2,:)*1000)/1000;

for i=1:Ndim		%初始化各粒子最大速度，使粒子不能越出边界
    vmax(i,:)=(Bound(i,2)-Bound(i,1))/2;
end

 for i=1:Ndim		%初始化各粒子速度
    velocity(i, :) =vmax(i,1)*rand(1, PopSize);
 end
velocity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fvalue=input('初始附值:下行自动产生Fvalue')

for i = 1:PopSize		%计算各粒子的适应度值
%         fvalue(i) = population(1,i)^2+population(2,i)^2+population(3,i)^2;  
%%%%%%后加计算适应度程序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epl33s  =  population(1,i) *1e-9;%e-9
tan_ct_e   =  population(2,i);
k_t  =  population(3,i);
tan_ct_k  =  population(4,i);
c33_D  =  population(5,i)   *1e+10;%e+10
tan_ct_c  =  population(6,i);

[sum_diff_Z]=  fun_epoxy_sum_diff(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );

fvalue(i) =sum_diff_Z;

end
fvalue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbest = population;   		%记录各粒子的个体极值点位置
fpbest = fvalue;      		%记录最佳适应度值
[fbestval,index] = min(fvalue);    	% 找出全局极值和相应的序号   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%读取样品数据
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
% x_m=load('2x.txt');
% r_m=load('2r.txt');
% phase_m=load('2p.txt');%复合材料

% ff_m=f_m(240:333);
% zz_m=z_m(240:333);
% xx_m=x_m(240:333);
% rr_m=r_m(240:333);
% phase_m=phase_m(240:333);%复合材料
%*************************** A

trace=zeros(iteration,2);   
Best=zeros(1,iteration);
while(flag == 0) && (iteration < MaxIter)	%寻找最优主程序
    
    iteration = iteration +1;		%迭代次数递增
    %w=(1+1/(0.5+log(iteration)))/2;  %惯性因子的一种随迭代改进法
  
    for i=1:PopSize			%更新全局极值点位置
        gbest(:,i) = population(:,index);
    end
     
    R1 = rand(Ndim, PopSize);		%重新设置随机性
    R2 = rand(Ndim, PopSize);
%   R3 = rand(Ndim, PopSize);
%   velocity = k*(velocity + c1*R1.*(pbest-population) +
%   c2*R2.*(gbest-population)); %加入收敛因子k
%   population = population + velocity;
    velocity = w*velocity + c1*R1.*(pbest-population) + c2*R2.*(gbest-population);% 更新各粒子速度    				
    population = population + velocity;	% 更新各粒子位置
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    population(1,:)=round(population(1,:)*1000)/1000;
    population(2,:)=round(population(2,:)*1000)/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OutFlag = population<=lowerbound | population>=upperbound; 	% 逸出标志
    population = population - OutFlag.*velocity;		% 阻止逸出
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      population
%      fvalue=input('初始附值:下行自动产生Fvalue')
     
    for i = 1:PopSize			% 更新各粒子适应度值
%           fvalue(i) = population(1,i)^2+population(2,i)^2+population(3,i)^2;
%%%%%后加计算适应度程序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epl33s  =  population(1,i) *1e-9;
tan_ct_e   =  population(2,i);
k_t  =  population(3,i);
tan_ct_k  =  population(4,i);
c33_D  =  population(5,i)   *1e+10;
tan_ct_c  =  population(6,i);

[sum_diff_Z]=  fun_epoxy_sum_diff(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );

fvalue(i) =sum_diff_Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    end
%     fvalue
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    changeColumns = fvalue <fpbest;	% 更新后的适应度值优于更新前的，记录序号
    pbest(:, find(changeColumns)) = population(:, find(changeColumns));    % 更新个体极值点位置
    fpbest = fpbest.*( ~changeColumns) + fvalue.*changeColumns ;           %更新个体极值
    [fbestval, index] = min(fvalue); 	%更新全局极值和相应的序号
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   trace(iteration,1)=min(fvalue);
   trace(iteration,2)=sum(fvalue)/length(fvalue);
   population(:,index);
   size(trace);
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
    if floor(fbestval*1e30)==oldbestval	 %比较更新前和更新后的放大的适应度值;
    samecounter=samecounter+1;		 %相等时记录加一;
    else
        oldbestval=floor(fbestval*1e30);	%不相等时更新放大的适应度值，并记录清零;
        samecounter=0;
    end 

    if samecounter >= 20		%多次迭代的适应度值相近时程序停止
        flag=1;
    end

     Best(iteration) =fbestval;		% 输出及描出找到的全局极值
	 
     figure(1)
     plot(Best,'ro-');xlabel({'generation','(a)'},'FontSize',14, 'FontWeight','bold','Color','k'); ylabel('f(x)');
     text(0.5,0.95,['Best = ', num2str(Best(iteration))],'Units','normalized');   
     set(gca,'linewidth',1.5);
     set(gca,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold');
     solutions=population(:,index);
     [R_c,X_c,Z_c,phase_c,ff]=  fun_epoxy(solutions(1)*1e-9,solutions(2),solutions(3),solutions(4),solutions(5)*1e+10,solutions(6));
         
     figure(2)
     plot(ff,log(Z_c),'mo-',ff_m,log(zz_m),'r');
     legend('Z_f Fitting Curve','Z_f Data Curve','Phase_f Fitting Curve','phase_f Data Curve','Location','NorthWest');
     xlabel({'Frequency (kHz)','(b)'},'FontSize',14, 'FontWeight','bold','Color','k');
     ylabel({'Resistance (Ω)'},'FontSize',14);
     title({'Fitting Result'},'FontSize',14, 'FontWeight','bold','Color','k');
     set(gca,'linewidth',1.5,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold');
     figure(3)
     plot(ff,phase_c,'bo-',ff_m,phase_m,'g');
     legend('Phase_f Fitting Curve','phase_f Data Curve');
     xlabel({'Frequency (kHz)','(b)'},'FontSize',14, 'FontWeight','bold','Color','k');
     ylabel({'degree'},'FontSize',14);
     title({'Fitting Result'},'FontSize',14, 'FontWeight','bold','Color','k');
     set(gca,'linewidth',1.5,'FontName','Times New Roman','FontSize',13, 'FontWeight','bold');
     drawnow; 
end		

solutions
fbestval

