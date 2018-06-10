clear all
close all
echo off
clc
%-----------------------------升余弦滤波器
Rs=100;             %符号频率
ts=1/Rs;            %符号间隔
k=1000;               %k个码元
A=1;                %bpsk信号幅值
fc=100;              %载波速率
Fs=10000;            %采样频率
Ts=1/Fs;            %采样间隔

snr_in_dB=-20:1:40;
for i=1:length(snr_in_dB),
    SNR=exp(snr_in_dB(i)*log(10)/10);            % 信噪比
    theo(i)=0.5*(1-1/sqrt(1+1/SNR));   % RayLeigh 信道误比特率理论值
end;

for i=1:length(snr_in_dB)
    snr=10^(snr_in_dB(i)/10);                          % signal to noise ratio
    sgma=sqrt(A/snr);                             % noise variance
    
    Fd=100;              %Doppler频偏，以Hz为单位
    tau=[0,5*Ts,10*Ts];          %多径延时，以s为单位
    pdb=[0,-5,-8];          %各径功率，以dB位单位
    h=rayleighchan(Ts,Fd,tau,pdb);


    Fd2=50;              %Doppler频偏，以Hz为单位
    tau2=[0,5*Ts];          %多径延时，以s为单位
    pdf2=[0,-5];          %各径功率，以dB位单位
    h0=rayleighchan(Ts,Fd2,tau2,pdf2);

    %-------------------------------通过信道

    code = randint(1,k); 
    N = k/Rs*Fs;              
    Npc = 1/Rs*Fs;                                                
    l = 0;  
    bpsk = zeros(1,N);  
    ct=zeros(1,N);
    n=sgma*randn(1,N);    % N normal distributed r.v with 0, variance sgma
    for j=1:k,
       for m = l:l+Npc-1  
           if code(1,j)==1  
             ct(1,m+1) = A*cos(2*pi*fc*m/Fs);
             bpsk(1,m+1) = A*cos(2*pi*fc*m/Fs);  
           elseif code(1,j)==0 
             ct(1,m+1) = A*cos(2*pi*fc*m/Fs);
             bpsk(1,m+1) = A*cos(2*pi*fc*m/Fs + pi);  
           end  
       end  
       l = l+Npc;  
    end
    
    bpsk=bpsk+n;
    code;


    y1=filter(h,bpsk);
    y1_fft=fft(y1);
    y1_abs=abs(y1_fft);

    y2=filter(h0,bpsk);
    y2_fft=fft(y2);
    y2_abs=abs(y2_fft);
    
    avr=mean(abs(y1).^2);
    db=-10*log10(avr);
    
    h1=zeros(1,length(y1));
    h2=zeros(1,length(y1));
    for j=1:length(y1),
            if y1(1,j)/bpsk(1,j)>0  
               h1(1,j)=1;
            else
               h1(1,j)= -1;
            end   
            if y2(1,j)/bpsk(1,j)>0  
               h2(1,j)= 1;
            else
               h2(1,j)= -1;
            end
    end

    out1=bpsk.*ct;
    Wc=2*100/Fs;                                          %截止频率 10000Hz
    [b,a]=butter(5,Wc);
    Signal_Filter=filter(b,a,out1);
    out2=y1.*ct.*h1;
    Signal_Filter2=filter(b,a,out2);
    out3=y2.*ct.*h2;
    Signal_Filter3=filter(b,a,out3);

    %{
    
    figure;
    subplot(3,1,1);
    x1=1:length(code);
    h1=stem(x1,code-1);
    title('信源信号s[n]');

    subplot(3,1,2);
    x2=1:length(bpsk);
    h2=stem(x2,bpsk);
    title('bpsk调制之后的信号s(t)');

    bpsk_fft=fft(bpsk);
    bpsk_abs=abs(bpsk_fft);  
    subplot(3,1,3);
    x3=1:length(bpsk_abs);
    h3=stem(x3,bpsk_abs);
    title('bpsk调制之后的信号s(t)的频谱');
    
    figure;
    subplot(2,2,1);
    x21=1:length(y1);
    h=stem(x21,y1);
    title('过瑞利信道1之后的信号r1(t)');

    subplot(2,2,2);
    x22=1:length(y1_abs);
    h22=stem(x22,y1_abs);
    title('过瑞利信道1之后的信号r1(t)的频谱');

    subplot(2,2,3);
    x23=1:length(y2);
    h=stem(x23,y2);
    title('过瑞利信道2之后的信号r2(t)');

    subplot(2,2,4);
    x24=1:length(y2_abs);
    h24=stem(x24,y2_abs);
    title('过瑞利信道2之后的信号r2(t)的频谱');

    figure;
    subplot(3,2,1);
    x31=1:length(Signal_Filter);
    h=stem(x31,Signal_Filter);
    title('out(t)过低通滤波器');

    subplot(3,2,2);
    x32=1:length(out1);
    h=stem(x32,out1);
    title('bpsk解调之后的信号out(t)');

    subplot(3,2,3);
    x33=1:length(Signal_Filter2);
    h=stem(x33,Signal_Filter2);
    title('out1(t)过低通滤波器');

    subplot(3,2,4);
    x34=1:length(out2);
    h=stem(x34,out2);
    title('r1(t)按bpsk解调之后的信号out1(t)');

    subplot(3,2,5);
    x35=1:length(Signal_Filter3);
    h=stem(x35,Signal_Filter3);
    title('out2(t)过低通滤波器');

    subplot(3,2,6);
    x36=1:length(out3);
    h=stem(x36,out3);
    title('r2(t)按bpsk解调之后的信号out2(t)');
    %}


    for j=1:k,
        if(Signal_Filter2(1,0.8*Fs/Rs+Fs/Rs*(j-1))>0)
            decode(1,j)=1;
        else
            decode(1,j)=0;
        end
    end
    decode;
    er=decode-code;
    number=0;
    for j=1:k,
            if(er(1,j)~=0)
               number=number+1;
            end
    end
    ber(1,i)=number/k;
end
ber
figure;
semilogy(snr_in_dB,theo,'r-');
hold on;
semilogy(snr_in_dB,ber);
hold on;
xlabel('Eb/No');
ylabel('Pe')
title('误码率');
legend('BPSK瑞利理论误码率','BPSK瑞利误码率仿真(无分集)');

figure;
plot(snr_in_dB,ber);
hold on;
plot(snr_in_dB,theo,'r-');
hold on;
xlabel('Eb/No');
ylabel('Pe')
title('误码率');
legend('BPSK瑞利理论误码率','BPSK瑞利误码率仿真(无分集)');
