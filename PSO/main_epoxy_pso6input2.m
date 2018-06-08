%��֤PSO�㷨��SA �㷨�Ľ����ȡ������ע��aa֮��Ĵ��� ��ע��AA
%���㸴�ϲ��ϵ���Ͻ����ȡ������ע��A֮��Ĵ��� ��ע�� aa
clear all;
clc;

flag=0;		    %ֹͣ�����־
oldbestval=0;	%��¼�ɵ���Ӧ��ֵ	
samecounter=0; 	%��¼�õ���ͬ��Ӧ��ֵ�ĵ�������
iteration = 0;	%��������
MaxIter =128;  	%����������
PopSize=96;   	%��Ⱥ��С ���Ӵ�
c1 = 0.85;      %ѧϰ����ԭ0.65
c2 = 0.85;      %ѧϰ����    ��
w=0.65; 		%��������ԭ0.75	��
%k=0.729;       %��������
%���ӵ����귶Χ:��1   ��2         ��3     ��4      ��5     ��6
%Bound=[1.6668 7.2228;0.0240 	0.104;  0.1806 	0.7826; -0.0039 	-0.01677;  2.0778 9.0038  ; 0.0066 	0.0286]	%30%
%Bound=[5.0004	6.1116;0.072	0.088;0.5418	0.6622;-0.01161 -0.01419;6.2334	7.6186;0.0198	0.0242]%10%
Bound=[4.4448	6.6672;0.064	0.096;0.4816	0.7224;-0.01032 -0.01548;5.5408	8.3112;0.0176	0.0264]%20%
%Bound=[1  8; 0  0.1;    0.1  1;  -0.1  -0.005;  2  10; 0.01 0.1]%ԭʼ����
%Bound=[0 10;0 1;0.2 0.9;  -0.1 0.5; 0  20; 0 0.1];%��ѽ��
%Bound=[6.3414	7.7506;0.0567	0.0693;0.45657	0.55803; -0.0396	-0.0484; 3.15288	3.85352;0.03519	0.04301];%���ϲ���10%
% % ��PopSize=6   	%��Ⱥ��С���ɲ���population 6ά����ÿ��Ϊһ�����������ҹ�6�飬
% % population = 
% %     3.5600    5.2230    1.5830    9.8950    6.2450    4.8110
% %     0.1520    0.1330    0.1430    0.1230    0.1580    0.1760
% %    25.2982   26.4053   22.0907   23.7982   27.8333   26.8085
% %     0.2461    0.2568    0.2794    0.2059    0.2603    0.2050
% %    34.1537   33.0500   38.7437   30.1501   37.6795   39.7084
% %     0.3990    0.3789    0.3439    0.3498    0.3214    0.3643
% % ��Ӧ fvalueΪ�������룬��6����

Ndim = length(Bound);  %�ռ�ά��

for i=1:PopSize			 %�����������±߽�
    lowerbound(:,i)=Bound(:,1);
    upperbound(:,i)=Bound(:,2);
end

for i=1:Ndim		%��ʼ�������ӳ�ʼλ�ã�����Ч��Χ�����ѡ��	
    population(i,:)=rand(1, PopSize)*(Bound(i,2)-Bound(i,1)) + Bound(i,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population(1,:)=round(population(1,:)*1000)/1000;
% population(2,:)=round(population(2,:)*1000)/1000;

for i=1:Ndim		%��ʼ������������ٶȣ�ʹ���Ӳ���Խ���߽�
    vmax(i,:)=(Bound(i,2)-Bound(i,1))/2;
end

 for i=1:Ndim		%��ʼ���������ٶ�
    velocity(i, :) =vmax(i,1)*rand(1, PopSize);
 end
velocity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fvalue=input('��ʼ��ֵ:�����Զ�����Fvalue')

for i = 1:PopSize		%��������ӵ���Ӧ��ֵ
%         fvalue(i) = population(1,i)^2+population(2,i)^2+population(3,i)^2;  
%%%%%%��Ӽ�����Ӧ�ȳ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
pbest = population;   		%��¼�����ӵĸ��弫ֵ��λ��
fpbest = fvalue;      		%��¼�����Ӧ��ֵ
[fbestval,index] = min(fvalue);    	% �ҳ�ȫ�ּ�ֵ����Ӧ�����   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ȡ��Ʒ����
%*************************** a
epl33s=5.556e-9;
tan_ct_e=0.08;
k_t=0.602;
tan_ct_k= - 0.0129;
c33_D=6.926e+10;
tan_ct_c=0.022;
[rr_m,xx_m,zz_m,phase_m,ff_m]=  fun_epoxy(epl33s, tan_ct_e, k_t,tan_ct_k,c33_D, tan_ct_c );%ԭʼ����
%*************************** a

%*************************** A
% f_m=load('2f.txt')/1000;
% z_m=load('2z.txt');
% x_m=load('2x.txt');
% r_m=load('2r.txt');
% phase_m=load('2p.txt');%���ϲ���

% ff_m=f_m(240:333);
% zz_m=z_m(240:333);
% xx_m=x_m(240:333);
% rr_m=r_m(240:333);
% phase_m=phase_m(240:333);%���ϲ���
%*************************** A

trace=zeros(iteration,2);   
Best=zeros(1,iteration);
while(flag == 0) && (iteration < MaxIter)	%Ѱ������������
    
    iteration = iteration +1;		%������������
    %w=(1+1/(0.5+log(iteration)))/2;  %�������ӵ�һ��������Ľ���
  
    for i=1:PopSize			%����ȫ�ּ�ֵ��λ��
        gbest(:,i) = population(:,index);
    end
     
    R1 = rand(Ndim, PopSize);		%�������������
    R2 = rand(Ndim, PopSize);
%   R3 = rand(Ndim, PopSize);
%   velocity = k*(velocity + c1*R1.*(pbest-population) +
%   c2*R2.*(gbest-population)); %������������k
%   population = population + velocity;
    velocity = w*velocity + c1*R1.*(pbest-population) + c2*R2.*(gbest-population);% ���¸������ٶ�    				
    population = population + velocity;	% ���¸�����λ��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    population(1,:)=round(population(1,:)*1000)/1000;
    population(2,:)=round(population(2,:)*1000)/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OutFlag = population<=lowerbound | population>=upperbound; 	% �ݳ���־
    population = population - OutFlag.*velocity;		% ��ֹ�ݳ�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      population
%      fvalue=input('��ʼ��ֵ:�����Զ�����Fvalue')
     
    for i = 1:PopSize			% ���¸�������Ӧ��ֵ
%           fvalue(i) = population(1,i)^2+population(2,i)^2+population(3,i)^2;
%%%%%��Ӽ�����Ӧ�ȳ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 
    changeColumns = fvalue <fpbest;	% ���º����Ӧ��ֵ���ڸ���ǰ�ģ���¼���
    pbest(:, find(changeColumns)) = population(:, find(changeColumns));    % ���¸��弫ֵ��λ��
    fpbest = fpbest.*( ~changeColumns) + fvalue.*changeColumns ;           %���¸��弫ֵ
    [fbestval, index] = min(fvalue); 	%����ȫ�ּ�ֵ����Ӧ�����
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   trace(iteration,1)=min(fvalue);
   trace(iteration,2)=sum(fvalue)/length(fvalue);
   population(:,index);
   size(trace);
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
    if floor(fbestval*1e30)==oldbestval	 %�Ƚϸ���ǰ�͸��º�ķŴ����Ӧ��ֵ;
    samecounter=samecounter+1;		 %���ʱ��¼��һ;
    else
        oldbestval=floor(fbestval*1e30);	%�����ʱ���·Ŵ����Ӧ��ֵ������¼����;
        samecounter=0;
    end 

    if samecounter >= 20		%��ε�������Ӧ��ֵ���ʱ����ֹͣ
        flag=1;
    end

     Best(iteration) =fbestval;		% ���������ҵ���ȫ�ּ�ֵ
	 
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
     ylabel({'Resistance (��)'},'FontSize',14);
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

