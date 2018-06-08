% clear all
% fa=1.15e+6;
% fr=1.04e+6;
% fp=fa;
% rou=7500;
% c=3.01e-9;
% l=10e-3;
% w=10e-3;
% t=1.945e-3;
% kt33=0.5*pi*fr/fa*atan(0.5*pi*(fa-fr)/fa)
% e33=t*c/(l*w)
% c33=4*rou*fp^2*t^2

X=Z_c';
Y=zz_m';

B=regress(X,Y)