%����-Ƥ��ѷ׼��Neyman-Pearson criterion�����龯���
%H0: x =  n
%H1: x = 1 + n 
%n�Ǿ�ֵΪ0������Ϊ1�ĸ�˹�۲�����
%H0��H1���������δ֪
%��������δ֪
%�龯����Ϊ0.1
%�۲����Ϊ1
a=0.1;

syms r;%�������
eqn=(1-normcdf(r,0,1))==a;
solx = solve(eqn,r);
fprintf('����-Ƥ��ѷ׼��ļ������Ϊ\n%f\n',solx);

PH1_H1=1-normcdf(solx,1,1);
fprintf('��ȷԤ���ĸ���P(H1|H1)Ϊ\n%f\n',PH1_H1);