%��С������׼��minimax criterion��
%�����Ԫ�ź�ģ��Ϊ��
%H0: x = -1 + n
%H1: x = 1 + n 
%n�Ǿ�ֵΪ0������Ϊ0.5�ĸ�˹�۲�����
%H0��H1���������δ֪
%��������:c00 = 1;c10 = 4; c11 = 2; c01 = 8
%�۲����Ϊ1

clear;close all;
syms x;%�۲���
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
C_MAX = zeros(1,601);   %����ƽ������
for n=1:601
    PH0_H0(n) = normcdf(0.01*(n-301),mu0,sigma);
    PH1_H0(n) = 1-normcdf(0.01*(n-301),mu0,sigma);
    PH0_H1(n) = normcdf(0.01*(n-301),mu1,sigma);
    PH1_H1(n) = 1-normcdf(0.01*(n-301),mu1,sigma);
    %�����п��ܵ�������ʵ�ֵ����ƽ������
    C=zeros(1,1001);
    for k=1:1:1001
        PH1=0.001*(k-1);
        PH0=1-PH1;
        C(k) = PH0*c10 + PH1*c11 + PH1*(c01-c11)*PH0_H1(n) - PH0*(c10-c00)*PH0_H0(n);
    end
    %��ÿ��������ʶ�Ӧ�ļ���ƽ������
    [C_MAX(n),i]=max(C(:));
end
[C_MIN_MAX,j]=min(C_MAX(:));%�ҵ���С�ļ���ƽ������
%���о�����
x=0.01*(j-301);
fprintf('���󻯼�С�ļ������Ϊ\n%f\n',x);
n=1:601;
figure(1)
plot(0.01*(n-301),PH0_H0);hold on;
plot(0.01*(n-301),PH1_H0);hold on;
plot(0.01*(n-301),PH0_H1);hold on;
plot(0.01*(n-301),PH1_H1);hold on;
title("�������ʺͼ�����޵Ĺ�ϵ");legend("P(H0|H0)","P(H1|H0)","P(H0|H1)","P(H1|H1)");
xlabel("�������");ylabel("����");
figure(2)
plot(0.01*(n-301),C_MAX);hold on;
[~,t] = min(C_MAX);%Ѱ�Ҽ�С�ļ���
plot(0.01*(n(t)-301),C_MAX(t),'rs','MarkerSize',6);%��ͼ�ϱ���õ� 
str = ['P(' num2str(0.01*(n(t)-301)) ',' num2str(C_MAX(t)) ')'];
text(0.01*(n(t)-301),C_MAX(t),str) %��ͼ�ϱ��ƽ�����ۺ����Ӧ������
title("����ƽ�����ۺͼ�����޵Ĺ�ϵ");xlabel("�������");ylabel("����ƽ������");





