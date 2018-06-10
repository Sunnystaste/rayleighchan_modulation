function[t,gt]=noise_nb(fc,B,N0,t,T)
%noise_nb.m
%窄带高斯过程
% INPUTS:
% fc:载波频率
% B:窄带高斯信号的带宽
% N0:噪声单边功率谱密度
% t:时间采样值
% OUTPUTS:
% t:时间采样值
%T=5
dt=t(2)-t(1);
t=[0:dt:T];             %时间向量
%产生功率为N0*B的高斯白噪声
P=N0*B;        %信噪比dB,10log10(Ac^2/2/(N0*B)
st=sqrt(P)*randn(1,length(t));     %宽带高斯信号
[f,sf]=T2F(t,st);                  %宽带高斯信号的傅立叶变换
%将上述白噪声通过窄带带通系统
[t,gt]=bpf(f,sf,fc-B/2,fc+B/2);   %宽带高斯信号经过窄带带通系统
gt=real(gt);
