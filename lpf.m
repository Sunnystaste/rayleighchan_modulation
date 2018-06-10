function [t,st]=lpf(f,sf,B)
% This function filter an input data using a
% lowpass filter
% Inputs:
%   f: input  frequency samples
%   sf: input data spectrum samples
%   B: lowpass's bandwidth with a rectangle lowpass
% Outputs:
%   t: time samples
%   st: output data 's time samples
df=f(2)-f(1);
T=1/df;
hf=zeros(1,length(f));
bf=[-floor(B/df):floor(B/df)]+floor(length(f)/2); %��ͨ�˲�����Ƶ�ʷ�Χ
hf(bf)=1;
%yf=hf.*sf.*exp(-j*2*pi*f*0.5*T);
yf=hf.*sf;
[t,st]=F2T(f,yf);