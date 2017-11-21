%最小平均错误概率准则（minimum mean probability of error criterion）
%假设二元信号模型为：
%H0: x = -1 + n
%H1: x = 1 + n 
%n是均值为0，方差为0.5的高斯观测噪声
%H0和H1的先验概率分别是0.3和0.7
%代价因子:c00 = 0;c10 = 1; c11 = 0; c01 = 1
%观测次数为1

clear;close all;
syms x;%观测量
mu0 = -1;
mu1 = 1;
sigma = (0.5)^(0.5);
PH0 = 0.3;
PH1 = 0.7;
c00 = 0;
c10 = 1; 
c11 = 0; 
c01 = 1;
%计算似然比函数
lamdax = exp(4*x);
lx = x;%两边取对数

%计算似然比门限cita
cita = (PH0*(c10-c00))/(PH1*(c01-c11));
gama = 1/4*log(cita);%两边取对数，gama即为贝叶斯检测中的最佳门限
fprintf('最小错误概率判决的检测门限为\n%f\n',gama);

%计算四个条件概率和检测门限的关系
PH0_H1 = zeros(1,601);
PH0_H0 = zeros(1,601);
PH1_H1 = zeros(1,601);
PH1_H0 = zeros(1,601);
C = zeros(1,601);   %平均代价，即平均错误概率
for n=1:601
    PH0_H0(n) = normcdf(0.01*(n-301),mu0,sigma);
    PH1_H0(n) = 1-normcdf(0.01*(n-301),mu0,sigma);
    PH0_H1(n) = normcdf(0.01*(n-301),mu1,sigma);
    PH1_H1(n) = 1-normcdf(0.01*(n-301),mu1,sigma);
    C(n) = PH0*c10 + PH1*c11 + PH1*(c01-c11)*PH0_H1(n) - PH0*(c10-c00)*PH0_H0(n);
end
n=1:601;
figure(1)
plot(0.01*(n-301),PH0_H0);hold on;
plot(0.01*(n-301),PH1_H0);hold on;
plot(0.01*(n-301),PH0_H1);hold on;
plot(0.01*(n-301),PH1_H1);hold on;
title("条件概率和检测门限的关系");legend("P(H0|H0)","P(H1|H0)","P(H0|H1)","P(H1|H1)");
xlabel("检测门限");ylabel("概率");
figure(2)
plot(0.01*(n-301),C);hold on;
[~,t] = min(C);%寻找最小平均代价
plot(0.01*(n(t)-301),C(t),'rs','MarkerSize',6);%在图上标出该点 
str = ['P(' num2str(0.01*(n(t)-301)) ',' num2str(C(t)) ')'];
text(0.01*(n(t)-301),C(t),str) %在图上标出平均代价和其对应的门限
title("平均错误概率和检测门限的关系");xlabel("检测门限");ylabel("平均错误概率");
