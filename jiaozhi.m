clear all;
close all;
echo off;
fs=1e6;                 %抽样频率1Mhz
fm=5e3;                   %基带频率\Rs\Rb=5Khz
n0=1e2;                     %码元数1k个
code = randi([0,1],1,n0); %随机二进制数据
code_fenji=code;%时间分集(两次)
code_hebing=[code,code_fenji];%交织之后的信号
n=2*n0;
Npc = 1/fm*fs;           %单个码元对应的采样点数   
final=Npc*n;             %采样点数
fc=1e5;                   %载波频率100Khz
t=0:1/fs:n/fm-1/fs;           %t时间
A=1;                     %信号幅值
l = 0;                   %计数mark  
bpsk = zeros(1,final);%bpsk信号   
ct=zeros(1,final);           %载波信号
code0=zeros(1,final);           %基带方波信号
decode=zeros(1,n);           %普通解码
decode5=zeros(1,n);           %普通解码
decode0=zeros(1,n);          %分集解码
decode1=zeros(1,n);          %分集解码
decode00=zeros(1,n);          %分集解码
decode000=zeros(1,n);          %分集解码
%生成BPSK序列
%生成载波序列
%生成基带方波序列
for j=1:n,
       for m = l:l+Npc-1  
           if code_hebing(1,j)==1
             code0(1,m+1)=1;
             ct(1,m+1) = A*cos(2*pi*fc*m/fs);
             bpsk(1,m+1) = A*cos(2*pi*fc*m/fs);  
           elseif code_hebing(1,j)==0 
             code0(1,m+1)=-1;
             ct(1,m+1) = A*cos(2*pi*fc*m/fs);
             bpsk(1,m+1) = A*cos(2*pi*fc*m/fs + pi);  
           end  
       end  
       l = l+Npc;  
end


figure
subplot(2,2,1);
plot(t,code0);
axis([0 n0/fm -2 2]);
title('等概二进制信源');
xlabel('t');
grid on;

subplot(2,2,2);
plot(abs(fft(code0)));%产生2psk信号的频谱
grid on;
title('矩形波频谱');xlabel('f');

subplot(2,2,3);
plot(t,bpsk);
title('bpsk信号');
axis([0 n/fm -2 2]);
xlabel('t');
grid on;


subplot(2,2,4);
plot(abs(fft(bpsk)));%产生2psk信号的频谱
grid on;
title('2psk信号频谱');xlabel('f');
%}

dpsk=bpsk.*ct;        %相干解调
[f,sf]=T2F(t,dpsk);
[t,dpsk0]=lpf(f,sf,2*fm);

Fd=55.5;              %Doppler频偏，以Hz为单位
tau=[0,5/fc,10/fc];          %多径延时，以s为单位
pdb=[0,-5,-8];          %各径功率，以dB位单位
h=rayleighchan(1/fs,Fd,tau,pdb);
y1=filter(h,bpsk);
h1=zeros(1,length(y1));


%生成估计序列
for j=1:length(y1),
    %估计序列
    if y1(1,j)/bpsk(1,j)>0  
               h1(1,j)=1;
    else
               h1(1,j)= -1;  
    end  
end
y11=y1.*ct.*h1;%瑞利信号解调1
%过LPF带宽为5Rs
[f0,sf0]=T2F(t,y11);
[t,ruili]=lpf(f0,sf0,2*fm);

r1=ruili(1:length(ruili)/2);
t0=0:1/fs:n0/fm-1/fs;           %分集时间
r2=ruili(length(ruili)/2+1:length(ruili));



subplot(2,2,1);
plot(t0,dpsk0(1:length(dpsk0)/2));
xlabel('t');
title('2psk解调信号之后过LPF（解分集）');

subplot(2,2,2);
plot(t,dpsk0);
xlabel('t');
title('2psk解调信号之后过LPF(未解分集)');

subplot(2,2,3);
plot(t0,r1);
xlabel('t');
title('1号时间段上的分集信号r1(t)');


subplot(2,2,4);
plot(t0,r2);
xlabel('t');
title('2号时间段上的分集信号r2(t)');


figure;
subplot(2,2,1)
plot(t0,r1);
xlabel('t');
title('r1(t)');

subplot(2,2,2)
plot(t0,r2);
xlabel('t');
title('r2(t)');

ruili0=zeros(1,length(r1));
for i=1:length(r1)%选择合并
    if(real(r1(i))>real(r2(i)))
        ruili0(i)=r1(i);
    else
        ruili0(i)=r2(i);
    end
end

subplot(2,2,3)
plot(t0,ruili0);
xlabel('t');
title('选择合并信号r0(t)');

ruili000=zeros(1,length(r1));
for i=1:length(r1)%最大比合并
        sum=abs(r1(i))+abs(r2(i));
        ruili000(i)=r1(i)*r1(i)/sum+r2(i)*r2(i)/sum;
end
 
subplot(2,2,4)
plot(t0,ruili000);
xlabel('t');
title('最大比合并信号r000(t)');
%{
snr_in_dB=-20:2:10;
theo=zeros(1,length(snr_in_dB));
theo_gause=zeros(1,length(snr_in_dB));
%snr=10.^(snr_in_dB/10);                          % snr数值
%sgma=sqrt(A./snr);                             % 噪声幅值
%N0=sgma'*randn(1,final);    %生成高斯噪声矩阵

for i=1:length(snr_in_dB),
    SNR=exp(snr_in_dB(i)*log(10)/10);            % 信噪比
    theo(i)=0.5*(1-1/sqrt(1+1/SNR));   % RayLeigh 信道误比特率理论值
end;

for i=1:length(snr_in_dB),
    SNR=exp(snr_in_dB(i)*log(10)/10);            % 信噪比
    theo_gause(i)=0.5*erfc(sqrt(SNR));       % AWGN 信道误比特率理论值
end;
%误码率序列初始化
ber=zeros(1,length(snr_in_dB));
ber0=zeros(1,length(snr_in_dB));
ber5=zeros(1,length(snr_in_dB));
ber1=zeros(1,length(snr_in_dB));
ber00=zeros(1,length(snr_in_dB));
ber000=zeros(1,length(snr_in_dB));

for k=1:length(snr_in_dB)
    %{
    out=y1.*ct.*h1+N0(k,:);
    out2=y2.*ct.*h2+N0(k,:);
    out3=y3.*ct.*h3+N0(k,:);
    out4=y4.*ct.*h4+N0(k,:);
    %}
    %瑞利分析
    %同上列分析
    %信号加噪
    n1 = wgn(1,final,-snr_in_dB(k));
    n2 = wgn(1,final,-snr_in_dB(k));
    n3 = wgn(1,final,-snr_in_dB(k));
    n4 = wgn(1,final,-snr_in_dB(k));
    n5 = wgn(1,final,-snr_in_dB(k));
    %{
    y11=awgn(y1,snr_in_dB(k));
    y22=awgn(y2,snr_in_dB(k));
    y33=awgn(y3,snr_in_dB(k));
    y44=awgn(y4,snr_in_dB(k));
    %}
    y11=y1+n1;
    y22=y2+n2;
    y33=y3+n3;
    y44=y4+n4;
    out=y11.*ct.*h1;%瑞利信号解调1
    out2=y22.*ct.*h2;%瑞利信号解调2
    out3=y33.*ct.*h3;%瑞利信号解调3
    out4=y44.*ct.*h4;%瑞利信号解调3
    
    %过低通滤波器
    [f0,sf0]=T2F(t,out);
    [t,ruili]=lpf(f0,sf0,10*fm);

    [f1,sf1]=T2F(t,out2);
    [t,ruili1]=lpf(f1,sf1,10*fm);
    
    [f3,sf3]=T2F(t,out3);
    [t,ruili3]=lpf(f3,sf3,10*fm);

    [f4,sf4]=T2F(t,out4);
    [t,ruili4]=lpf(f4,sf4,10*fm);
    
    %高斯分析------------------------------  
    %bpsk0=awgn(bpsk,snr_in_dB(k));        %加噪声
    bpsk0=bpsk+n5;
    dpsk=bpsk0.*ct;        %相干解调
    [f,sf]=T2F(t,dpsk);    %过低通滤波器
    [t,dpsk0]=lpf(f,sf,10*fm);
    
 
    %采样判决
    for j=1:n,
        if(dpsk0(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode5(1,j)=1;
        else
                decode5(1,j)=0;
        end
    end
    %计算误码率
    er5=decode5-code;
    number5=0;
    for j=1:n,
        if(er5(1,j)~=0)
            number5=number5+1;
        end
    end
    ber5(1,k)=number5/n;
    
    
    
    
    %瑞利分析------------------------------
    %无分集分析
    %采样判决
    for j=1:n,
        if(ruili(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode(1,j)=1;
        else
                decode(1,j)=0;
        end
    end
    %计算误码率
    er=decode-code;
    number=0;
    for j=1:n,
        if(er(1,j)~=0)
            number=number+1;
        end
    end
    ber(1,k)=number/n;
    
    
    %分集情况分析1
    ruili0=zeros(1,length(ruili1));
    for i=1:length(ruili)%选择合并4
        if(abs(ruili(i))>abs(ruili1(i))&&abs(ruili(i))>abs(ruili3(i))&&abs(ruili(i))>abs(ruili4(i)))
            ruili0(i)=ruili(i);
        elseif(abs(ruili1(i))>abs(ruili(i))&&abs(ruili1(i))>abs(ruili3(i))&&abs(ruili1(i))>abs(ruili4(i)))
            ruili0(i)=ruili1(i);
        elseif(abs(ruili3(i))>abs(ruili(i))&&abs(ruili3(i))>abs(ruili1(i))&&abs(ruili3(i))>abs(ruili4(i)))
            ruili0(i)=ruili3(i);
        elseif(abs(ruili4(i))>abs(ruili(i))&&abs(ruili4(i))>abs(ruili1(i))&&abs(ruili4(i))>abs(ruili3(i)))
            ruili0(i)=ruili4(i);
        else
            ruili0(i)=ruili(i);
        end
    end
    %采样判决
    for j=1:n,
        if(ruili0(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode0(1,j)=1;
        else
                decode0(1,j)=0;
        end
    end
    %计算误码率
    er0=decode0-code;
    number0=0;
    for j=1:n,
        if(er0(1,j)~=0)
            number0=number0+1;
        end
    end
    ber0(1,k)=number0/n;
    
    
    %分集情况分析2
    ruili00=zeros(1,length(ruili1));
    for i=1:length(ruili)%选择合并2
        if(abs(ruili(i))>abs(ruili1(i)))
            ruili00(i)=ruili(i);
        else
            ruili00(i)=ruili1(i);
        end
    end
    
    %采样判决
    for j=1:n,
        if(ruili00(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode1(1,j)=1;
        else
                decode1(1,j)=0;
        end
    end
    %计算误码率
    er1=decode1-code;
    number1=0;
    for j=1:n,
        if(er1(1,j)~=0)
            number1=number1+1;
        end
    end
    ber1(1,k)=number1/n;
    
    %分集情况分析3
    ruili000=zeros(1,length(ruili1));
    for i=1:length(ruili)%最大比合并
        sum=abs(ruili(i))^2+abs(ruili1(i))^2+abs(ruili3(i))^2+abs(ruili4(i))^2;
        ruili000(i)=ruili(i)*abs(ruili(i))^2/sum+ruili1(i)*abs(ruili1(i))^2/sum+ruili3(i)*abs(ruili3(i))^2/sum+ruili4(i)*abs(ruili4(i))^2/sum;
    end
    %采样判决
    for j=1:n,
        if(ruili000(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode00(1,j)=1;
        else
                decode00(1,j)=0;
        end
    end
    %计算误码率
    er00=decode00-code;
    number00=0;
    for j=1:n,
        if(er00(1,j)~=0)
            number00=number00+1;
        end
    end
    ber00(1,k)=number00/n;
    
    %分集情况分析4
    ruili0000=zeros(1,length(ruili1));
    for i=1:length(ruili)%等增益合并
        ruili0000(i)=(ruili(i)+ruili1(i)+ruili3(i)+ruili4(i))/4;
    end
    %采样判决
    for j=1:n,
        if(ruili0000(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode000(1,j)=1;
        else
                decode000(1,j)=0;
        end
    end
    %计算误码率
    er000=decode000-code;
    number000=0;
    for j=1:n,
        if(er000(1,j)~=0)
            number000=number000+1;
        end
    end
    ber000(1,k)=number000/n;
end
figure;
plot(snr_in_dB,theo_gause,'--b',snr_in_dB,ber5,'pr',snr_in_dB,ber,'*r',snr_in_dB,theo,'-b');
legend('BPSK高斯理论误码率','BPSK高斯误码率仿真','BPSK瑞利理论误码率','BPSK瑞利误码率仿真(无分集)');
figure;
plot(snr_in_dB,ber000,'dr',snr_in_dB,ber00,'^k',snr_in_dB,ber1,'sb',snr_in_dB,ber0,'og',snr_in_dB,ber,'*r',snr_in_dB,theo,'-b');
legend('BPSK瑞利误码率仿真(等增益合并)','BPSK瑞利误码率仿真(最大比合并)','BPSK瑞利误码率仿真(选择合并2)','BPSK瑞利误码率仿真(选择合并4)','BPSK瑞利误码率仿真(无分集)','BPSK瑞利理论误码率');
%}