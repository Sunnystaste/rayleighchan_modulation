function[t,gt]=noise_nb(fc,B,N0,t,T)
%noise_nb.m
%խ����˹����
% INPUTS:
% fc:�ز�Ƶ��
% B:խ����˹�źŵĴ���
% N0:�������߹������ܶ�
% t:ʱ�����ֵ
% OUTPUTS:
% t:ʱ�����ֵ
%T=5
dt=t(2)-t(1);
t=[0:dt:T];             %ʱ������
%��������ΪN0*B�ĸ�˹������
P=N0*B;        %�����dB,10log10(Ac^2/2/(N0*B)
st=sqrt(P)*randn(1,length(t));     %�����˹�ź�
[f,sf]=T2F(t,st);                  %�����˹�źŵĸ���Ҷ�任
%������������ͨ��խ����ͨϵͳ
[t,gt]=bpf(f,sf,fc-B/2,fc+B/2);   %�����˹�źž���խ����ͨϵͳ
gt=real(gt);
