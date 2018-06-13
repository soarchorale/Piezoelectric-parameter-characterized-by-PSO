%此处bound的界限要和主程序中一致
function Var_new=get_newanswer(Var_best)
Var_new=Var_best+0.2*(2*rand(6,1)-[1;1;1;1;1;1]);%以0.02为精度随机搜索邻域
 Bound=[5.0004	6.1116;0.072	0.088; 0.5418	0.6622;-0.01161	-0.01419; 6.2334	7.6186;0.0198	0.0242];%10%
%Bound=[1.6668 	7.2228;0.0240 	0.104;  0.1806 	0.7826; -0.0039 	-0.01677;  2.0778 	9.0038  ; 0.0066 	0.0286];	%粒子的坐标范围	
%Bound=[1.6668 7.2228;0.0240 	0.104;  0.1806 	0.7826; -0.0039 	-0.01677;  2.0778 9.0038  ; 0.0066 	0.0286];	%30%
% Bound=[4.4448	6.6672;0.064	0.096;0.4816	0.7224;-0.01032 -0.01548;5.5408	8.3112;0.0176	0.0264];%20%
%Bound=[1  8; 0  0.1;    0.1  1;  -0.1  -0.005;  2  10; 0.01 0.1];%原始区间
% Bound=[5.44488	5.66712;0.0784	0.0816;0.58996	0.61404;-0.012642	-0.013158;6.78748	7.06452;0.02156	0.02244];
%Bound=[6 8;0 0.1;0.4 0.6;  -0.1 0.1; 3  4; 0 0.05];;%凑区间
lowerbound=Bound(:,1);%定义粒子上下边界
upperbound=Bound(:,2);
outupper=Var_new>=upperbound;
outlower=Var_new<=lowerbound;
while max(outupper)||max(outlower) % 逸出标志并阻止溢出
    if max(outupper)==1
    Var_new=Var_new-outupper.*(0.2*rand(6,1));
    outupper=Var_new>=upperbound;
    else
        if max(outlower)==1
    Var_new=Var_new+outlower.*(0.2*rand(6,1));
    outlower=Var_new<=lowerbound;
        end
    end
end
end