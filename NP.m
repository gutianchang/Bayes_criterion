%奈曼-皮尔逊准则（Neyman-Pearson criterion）恒虚警检测
%H0: x =  n
%H1: x = 1 + n 
%n是均值为0，方差为1的高斯观测噪声
%H0和H1的先验概率未知
%代价因子未知
%虚警概率为0.1
%观测次数为1
a=0.1;

syms r;%检测门限
eqn=(1-normcdf(r,0,1))==a;
solx = solve(eqn,r);
fprintf('奈曼-皮尔逊准则的检测门限为\n%f\n',solx);

PH1_H1=1-normcdf(solx,1,1);
fprintf('正确预警的概率P(H1|H1)为\n%f\n',PH1_H1);