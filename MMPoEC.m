%��Сƽ���������׼��minimum mean probability of error criterion��
%�����Ԫ�ź�ģ��Ϊ��
%H0: x = -1 + n
%H1: x = 1 + n 
%n�Ǿ�ֵΪ0������Ϊ0.5�ĸ�˹�۲�����
%H0��H1��������ʷֱ���0.3��0.7
%��������:c00 = 0;c10 = 1; c11 = 0; c01 = 1
%�۲����Ϊ1

clear;close all;
syms x;%�۲���
mu0 = -1;
mu1 = 1;
sigma = (0.5)^(0.5);
PH0 = 0.3;
PH1 = 0.7;
c00 = 0;
c10 = 1; 
c11 = 0; 
c01 = 1;
%������Ȼ�Ⱥ���
lamdax = exp(4*x);
lx = x;%����ȡ����

%������Ȼ������cita
cita = (PH0*(c10-c00))/(PH1*(c01-c11));
gama = 1/4*log(cita);%����ȡ������gama��Ϊ��Ҷ˹����е��������
fprintf('��С��������о��ļ������Ϊ\n%f\n',gama);

%�����ĸ��������ʺͼ�����޵Ĺ�ϵ
PH0_H1 = zeros(1,601);
PH0_H0 = zeros(1,601);
PH1_H1 = zeros(1,601);
PH1_H0 = zeros(1,601);
C = zeros(1,601);   %ƽ�����ۣ���ƽ���������
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
title("�������ʺͼ�����޵Ĺ�ϵ");legend("P(H0|H0)","P(H1|H0)","P(H0|H1)","P(H1|H1)");
xlabel("�������");ylabel("����");
figure(2)
plot(0.01*(n-301),C);hold on;
[~,t] = min(C);%Ѱ����Сƽ������
plot(0.01*(n(t)-301),C(t),'rs','MarkerSize',6);%��ͼ�ϱ���õ� 
str = ['P(' num2str(0.01*(n(t)-301)) ',' num2str(C(t)) ')'];
text(0.01*(n(t)-301),C(t),str) %��ͼ�ϱ��ƽ�����ۺ����Ӧ������
title("ƽ��������ʺͼ�����޵Ĺ�ϵ");xlabel("�������");ylabel("ƽ���������");
