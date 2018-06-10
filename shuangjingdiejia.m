clear all
close all
Rs=9600;             %符号频率
ts=1/Rs;            %符号间隔
k=100;               %10个码元
A=1;                %bpsk信号幅值
fc=9600;              %载波速率
Fs=960000;            %采样频率
Ts=1/Fs;            %采样间隔
tau=5*Ts;
theta=tau*2*pi*Rs;

Fd=100;              %Doppler频偏，以Hz为单位
Npc = 1/Rs*Fs;                                                
l = 0;
for i=1:k,
   for j = l:l+Npc-1 
         s3(1,j+1) = A*cos(2*pi*(fc+2*Fd)*j/Fs+2*theta); 
         s2(1,j+1) = A*cos(2*pi*(fc+Fd)*j/Fs+theta);
         s1(1,j+1) = A*cos(2*pi*fc*j/Fs);  
   end  
   l = l+Npc;  
end
s=s1+s2+s3;
figure;
x1=1:length(s);
h1=stem(x1,s);