clear all
close all
Rs=9600;             %����Ƶ��
ts=1/Rs;            %���ż��
k=100;               %10����Ԫ
A=1;                %bpsk�źŷ�ֵ
fc=9600;              %�ز�����
Fs=960000;            %����Ƶ��
Ts=1/Fs;            %�������
tau=5*Ts;
theta=tau*2*pi*Rs;

Fd=100;              %DopplerƵƫ����HzΪ��λ
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