%极小化极大准则（minimax criterion）
%假设二元信号模型为：
%H0: x = -1 + n
%H1: x = 1 + n 
%n是均值为0，方差为0.5的高斯观测噪声
%H0和H1的先验概率未知
%代价因子:c00 = 1;c10 = 4; c11 = 2; c01 = 8
%观测次数为1

clear;close all;
syms x;%观测量
mu0 = -1;
mu1 = 1;
sigma = (0.5)^(0.5);
c00 = 1;
c10 = 4; 
c11 = 2; 
c01 = 8;

PH0_H1 = zeros(1,601);
PH0_H0 = zeros(1,601);
PH1_H1 = zeros(1,601);
PH1_H0 = zeros(1,601);
C_MAX = zeros(1,601);   %极大平均代价
for n=1:601
    PH0_H0(n) = normcdf(0.01*(n-301),mu0,sigma);
    PH1_H0(n) = 1-normcdf(0.01*(n-301),mu0,sigma);
    PH0_H1(n) = normcdf(0.01*(n-301),mu1,sigma);
    PH1_H1(n) = 1-normcdf(0.01*(n-301),mu1,sigma);
    %对所有可能的先验概率的值，求平均代价
    C=zeros(1,1001);
    for k=1:1:1001
        PH1=0.001*(k-1);
        PH0=1-PH1;
        C(k) = PH0*c10 + PH1*c11 + PH1*(c01-c11)*PH0_H1(n) - PH0*(c10-c00)*PH0_H0(n);
    end
    %求每个先验概率对应的极大平均代价
    [C_MAX(n),i]=max(C(:));
end
[C_MIN_MAX,j]=min(C_MAX(:));%找到极小的极大平均代价
%求判决门限
x=0.01*(j-301);
fprintf('极大化极小的检测门限为\n%f\n',x);
n=1:601;
figure(1)
plot(0.01*(n-301),PH0_H0);hold on;
plot(0.01*(n-301),PH1_H0);hold on;
plot(0.01*(n-301),PH0_H1);hold on;
plot(0.01*(n-301),PH1_H1);hold on;
title("条件概率和检测门限的关系");legend("P(H0|H0)","P(H1|H0)","P(H0|H1)","P(H1|H1)");
xlabel("检测门限");ylabel("概率");
figure(2)
plot(0.01*(n-301),C_MAX);hold on;
[~,t] = min(C_MAX);%寻找极小的极大
plot(0.01*(n(t)-301),C_MAX(t),'rs','MarkerSize',6);%在图上标出该点 
str = ['P(' num2str(0.01*(n(t)-301)) ',' num2str(C_MAX(t)) ')'];
text(0.01*(n(t)-301),C_MAX(t),str) %在图上标出平均代价和其对应的门限
title("极大平均代价和检测门限的关系");xlabel("检测门限");ylabel("极大平均代价");





