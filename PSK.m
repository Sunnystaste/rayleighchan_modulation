clear all;
close all;
echo off;
fs=1e6;                 %����Ƶ��1Mhz
fm=5e3;                   %����Ƶ��\Rs\Rb=5Khz
n=1e2;                     %��Ԫ��1k��
code = randi([0,1],1,n); %�������������
Npc = 1/fm*fs;           %������Ԫ��Ӧ�Ĳ�������   
final=Npc*n;             %��������
fc=1e5;                   %�ز�Ƶ��100Khz
t=0:1/fs:n/fm-1/fs;           %tʱ��
A=1;                     %�źŷ�ֵ
l = 0;                   %����mark  
bpsk = zeros(1,final);%bpsk�ź�   
ct=zeros(1,final);           %�ز��ź�
code0=zeros(1,final);           %���������ź�
decode=zeros(1,n);           %��ͨ����
decode5=zeros(1,n);           %��ͨ����
decode0=zeros(1,n);          %�ּ�����
decode1=zeros(1,n);          %�ּ�����
decode00=zeros(1,n);          %�ּ�����
decode000=zeros(1,n);          %�ּ�����
%����BPSK����
%�����ز�����
%���ɻ�����������
for j=1:n,
       for m = l:l+Npc-1  
           if code(1,j)==1
             code0(1,m+1)=1;
             ct(1,m+1) = A*cos(2*pi*fc*m/fs);
             bpsk(1,m+1) = A*cos(2*pi*fc*m/fs);  
           elseif code(1,j)==0 
             code0(1,m+1)=-1;
             ct(1,m+1) = A*cos(2*pi*fc*m/fs);
             bpsk(1,m+1) = A*cos(2*pi*fc*m/fs + pi);  
           end  
       end  
       l = l+Npc;  
end


%{
figure
subplot(2,2,1);
plot(t,code0);
axis([0 n/fm -2 2]);
title('�ȸŶ�������Դ');
xlabel('t');
grid on;

subplot(2,2,2);
plot(abs(fft(code0)));%����2psk�źŵ�Ƶ��
grid on;
title('���β�Ƶ��');xlabel('f');

subplot(2,2,3);
plot(t,bpsk);
title('bpsk�ź�');
axis([0 n/fm -2 2]);
xlabel('t');
grid on;


subplot(2,2,4);
plot(abs(fft(bpsk)));%����2psk�źŵ�Ƶ��
grid on;
title('2psk�ź�Ƶ��');xlabel('f');
%}

dpsk=bpsk.*ct;        %��ɽ��
[f,sf]=T2F(t,dpsk);
[t,dpsk0]=lpf(f,sf,2*fm);

%�ּ�1
Fd=55.5;              %DopplerƵƫ����HzΪ��λ
tau=[0,5/fc,10/fc];          %�ྶ��ʱ����sΪ��λ
pdb=[0,-5,-8];          %�������ʣ���dBλ��λ
h=rayleighchan(1/fs,Fd,tau,pdb);
y1=filter(h,bpsk);
h1=zeros(1,length(y1));

%�ּ�2
Fd0=100.3;              %DopplerƵƫ����HzΪ��λ
tau0=[0,1/fc,2/fc];          %�ྶ��ʱ����sΪ��λ
pdb0=[0,-4,-6];          %�������ʣ���dBλ��λ
h0=rayleighchan(1/fs,Fd,tau0,pdb0);
y2=filter(h0,bpsk);
h2=zeros(1,length(y2));

%�ּ�3
Fd3=1.11;              %DopplerƵƫ����HzΪ��λ
tau3=[0,3/fc,4/fc];          %�ྶ��ʱ����sΪ��λ
pdb3=[0,-8,-9];          %�������ʣ���dBλ��λ
h3=rayleighchan(1/fs,Fd,tau3,pdb3);
y3=filter(h3,bpsk);
h3=zeros(1,length(y3));

%�ּ�4
Fd4=243.5;              %DopplerƵƫ����HzΪ��λ
tau4=[0,5/fc,1/fc];          %�ྶ��ʱ����sΪ��λ
pdb4=[0,-8,-3];          %�������ʣ���dBλ��λ
h4=rayleighchan(1/fs,Fd,tau4,pdb4);
y4=filter(h4,bpsk);
h4=zeros(1,length(y4));

%���ɹ�������
for j=1:length(y1),
    %�ּ�1��������
    if y1(1,j)/bpsk(1,j)>0  
               h1(1,j)=1;
    else
               h1(1,j)= -1;
    end  
    %�ּ�2��������
    if y2(1,j)/bpsk(1,j)>0  
               h2(1,j)=1;
    else
               h2(1,j)= -1;
    end 
    %�ּ�3��������
    if y3(1,j)/bpsk(1,j)>0  
               h3(1,j)=1;
    else
               h3(1,j)= -1;
    end
    %�ּ�4��������
    if y4(1,j)/bpsk(1,j)>0  
               h4(1,j)=1;
    else
               h4(1,j)= -1;
    end  
end
y11=y1.*ct.*h1;%�����źŽ��1
y22=y2.*ct.*h2;%�����źŽ��2
y33=y3.*ct.*h3;%�����źŽ��3
y44=y4.*ct.*h4;%�����źŽ��3
avr=mean(abs(y1).^2);
db=-10*log10(avr);
%��LPF����Ϊ5Rs
[f0,sf0]=T2F(t,y11);
[t,ruili]=lpf(f0,sf0,2*fm);

[f1,sf1]=T2F(t,y22);
[t,ruili1]=lpf(f1,sf1,2*fm);

[f3,sf3]=T2F(t,y33);
[t,ruili3]=lpf(f3,sf3,2*fm);

[f4,sf4]=T2F(t,y44);
[t,ruili4]=lpf(f4,sf4,2*fm);
%{
figure
subplot(3,2,1);
plot(t,dpsk);
xlabel('t');
axis([0 n/fm -2 2]);
title('2psk����ź�');

subplot(3,2,2);
plot(t,dpsk0);
xlabel('t');
title('2psk����ź�֮���LPF');

subplot(3,2,3);
plot(t,y11);
xlabel('t');
title('�������ŵ�1֮����ź�r1(t)*�ز�*�����ź�');


subplot(3,2,4);
plot(t,ruili);
xlabel('t');
title('r1(t)*�ز�*�����ź�֮���LPF');

subplot(3,2,5);
plot(abs(fft(dpsk0)));
grid on;
title('2psk����źŹ�LPF��Ƶ��');xlabel('f');

subplot(3,2,6);
plot(abs(fft(ruili)));
grid on;
title('�������ŵ�1֮����ź�r1(t)*�ز�*�����ź�Ƶ��');xlabel('f');

figure;
subplot(2,2,1)
plot(t,ruili);
xlabel('t');
title('r1(t)*�ز�*�����ź�֮���LPF');

subplot(2,2,2)
plot(t,ruili1);
xlabel('t');
title('r2(t)*�ز�*�����ź�֮���LPF');
%}
ruili0=zeros(1,length(ruili1));
for i=1:length(ruili)%ѡ��ϲ�
    if(real(ruili(i))>real(ruili1(i)))
        ruili0(i)=ruili(i);
    else
        ruili0(i)=ruili1(i);
    end
end
%{
subplot(2,2,3)
plot(t,ruili0);
xlabel('t');
title('ѡ��ϲ��ź�r0(t)');
%}
ruili000=zeros(1,length(ruili1));
for i=1:length(ruili)%���Ⱥϲ�
        sum=abs(ruili(i))+abs(ruili1(i))+abs(ruili3(i))+abs(ruili4(i));
        ruili000(i)=(ruili(i)+ruili1(i)+ruili3(i)+ruili4(i))/sum;
end
 %{
subplot(2,2,4)
plot(t,ruili000);
xlabel('t');
title('���Ⱥϲ��ź�r000(t)');
%}
snr_in_dB=-20:2:10;
theo=zeros(1,length(snr_in_dB));
theo_gause=zeros(1,length(snr_in_dB));
%snr=10.^(snr_in_dB/10);                          % snr��ֵ
%sgma=sqrt(A./snr);                             % ������ֵ
%N0=sgma'*randn(1,final);    %���ɸ�˹��������

for i=1:length(snr_in_dB),
    SNR=exp(snr_in_dB(i)*log(10)/10);            % �����
    theo(i)=0.5*(1-1/sqrt(1+1/SNR));   % RayLeigh �ŵ������������ֵ
end;

for i=1:length(snr_in_dB),
    SNR=exp(snr_in_dB(i)*log(10)/10);            % �����
    theo_gause(i)=0.5*erfc(sqrt(SNR));       % AWGN �ŵ������������ֵ
end;
%���������г�ʼ��
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
    %��������
    %ͬ���з���
    %�źż���
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
    out=y11.*ct.*h1;%�����źŽ��1
    out2=y22.*ct.*h2;%�����źŽ��2
    out3=y33.*ct.*h3;%�����źŽ��3
    out4=y44.*ct.*h4;%�����źŽ��3
    
    %����ͨ�˲���
    [f0,sf0]=T2F(t,out);
    [t,ruili]=lpf(f0,sf0,10*fm);

    [f1,sf1]=T2F(t,out2);
    [t,ruili1]=lpf(f1,sf1,10*fm);
    
    [f3,sf3]=T2F(t,out3);
    [t,ruili3]=lpf(f3,sf3,10*fm);

    [f4,sf4]=T2F(t,out4);
    [t,ruili4]=lpf(f4,sf4,10*fm);
    
    %��˹����------------------------------  
    %bpsk0=awgn(bpsk,snr_in_dB(k));        %������
    bpsk0=bpsk+n5;
    dpsk=bpsk0.*ct;        %��ɽ��
    [f,sf]=T2F(t,dpsk);    %����ͨ�˲���
    [t,dpsk0]=lpf(f,sf,10*fm);
    
 
    %�����о�
    for j=1:n,
        if(dpsk0(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode5(1,j)=1;
        else
                decode5(1,j)=0;
        end
    end
    %����������
    er5=decode5-code;
    number5=0;
    for j=1:n,
        if(er5(1,j)~=0)
            number5=number5+1;
        end
    end
    ber5(1,k)=number5/n;
    
    
    
    
    %��������------------------------------
    %�޷ּ�����
    %�����о�
    for j=1:n,
        if(ruili(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode(1,j)=1;
        else
                decode(1,j)=0;
        end
    end
    %����������
    er=decode-code;
    number=0;
    for j=1:n,
        if(er(1,j)~=0)
            number=number+1;
        end
    end
    ber(1,k)=number/n;
    
    
    %�ּ��������1
    ruili0=zeros(1,length(ruili1));
    for i=1:length(ruili)%ѡ��ϲ�4
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
    %�����о�
    for j=1:n,
        if(ruili0(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode0(1,j)=1;
        else
                decode0(1,j)=0;
        end
    end
    %����������
    er0=decode0-code;
    number0=0;
    for j=1:n,
        if(er0(1,j)~=0)
            number0=number0+1;
        end
    end
    ber0(1,k)=number0/n;
    
    
    %�ּ��������2
    ruili00=zeros(1,length(ruili1));
    for i=1:length(ruili)%ѡ��ϲ�2
        if(abs(ruili(i))>abs(ruili1(i)))
            ruili00(i)=ruili(i);
        else
            ruili00(i)=ruili1(i);
        end
    end
    
    %�����о�
    for j=1:n,
        if(ruili00(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode1(1,j)=1;
        else
                decode1(1,j)=0;
        end
    end
    %����������
    er1=decode1-code;
    number1=0;
    for j=1:n,
        if(er1(1,j)~=0)
            number1=number1+1;
        end
    end
    ber1(1,k)=number1/n;
    
    %�ּ��������3
    ruili000=zeros(1,length(ruili1));
    for i=1:length(ruili)%���Ⱥϲ�
        sum=abs(ruili(i))^2+abs(ruili1(i))^2+abs(ruili3(i))^2+abs(ruili4(i))^2;
        ruili000(i)=ruili(i)*abs(ruili(i))^2/sum+ruili1(i)*abs(ruili1(i))^2/sum+ruili3(i)*abs(ruili3(i))^2/sum+ruili4(i)*abs(ruili4(i))^2/sum;
    end
    %�����о�
    for j=1:n,
        if(ruili000(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode00(1,j)=1;
        else
                decode00(1,j)=0;
        end
    end
    %����������
    er00=decode00-code;
    number00=0;
    for j=1:n,
        if(er00(1,j)~=0)
            number00=number00+1;
        end
    end
    ber00(1,k)=number00/n;
    
    %�ּ��������4
    ruili0000=zeros(1,length(ruili1));
    for i=1:length(ruili)%������ϲ�
        ruili0000(i)=(ruili(i)+ruili1(i)+ruili3(i)+ruili4(i))/4;
    end
    %�����о�
    for j=1:n,
        if(ruili0000(1,0.5*fs/fm+fs/fm*(j-1))>0)
                decode000(1,j)=1;
        else
                decode000(1,j)=0;
        end
    end
    %����������
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
legend('BPSK��˹����������','BPSK��˹�����ʷ���','BPSK��������������','BPSK���������ʷ���(�޷ּ�)');
figure;
plot(snr_in_dB,ber000,'dr',snr_in_dB,ber00,'^k',snr_in_dB,ber1,'sb',snr_in_dB,ber0,'og',snr_in_dB,ber,'*r',snr_in_dB,theo,'-b');
legend('BPSK���������ʷ���(������ϲ�)','BPSK���������ʷ���(���Ⱥϲ�)','BPSK���������ʷ���(ѡ��ϲ�2)','BPSK���������ʷ���(ѡ��ϲ�4)','BPSK���������ʷ���(�޷ּ�)','BPSK��������������');