clc; clear all; close all;
% QPSK�������У�����ͼ
mapper=[1/sqrt(2)+1i/sqrt(2) -1/sqrt(2)+1i/sqrt(2) 
       1/sqrt(2)-1i/sqrt(2) -1/sqrt(2)-1i/sqrt(2)];
N_OFDM_Frame=100;           % OFDM�źŵĸ���
N_Subcarriers=1024;         % ���ز�����
ber_begin=0;ber_end=15;     % ��ʼ����ֹ�����
snr_draw = [0,5,10,15];     % ����0��5��10��15db������µ������Ϣ
title_n = ["0dB","5dB","10dB","15dB"];
tool1 = 1;tool2 = 1;tool3 = 1;
for EsN0=ber_begin:ber_end
    Num_Error_Symbol=0;
    for n=1:N_OFDM_Frame
        %% ��������������У�QPSK��
        InputBits=randi([0 1], 2, N_Subcarriers);
        IndexSymbol = InputBits(1,:)*2+InputBits(2,:)+1;%����ά�ź�ת����һά��������
        %% �����ź�ӳ�� 
        Tx_FreqDomain=mapper(IndexSymbol); % QPSK�ź�
        %% FFT���������ź�ת��Ϊʱ��
        Tx_TimeDomain=sqrt(N_Subcarriers) *ifft(Tx_FreqDomain);
        figure(1)
        if EsN0 == ber_begin
            plot(1:N_Subcarriers,real(Tx_TimeDomain));title("OFDMʱ���źŲ���")
        end
        %% ͨ����˹���԰������ŵ�
        Rx_TimeDomain = awgn(Tx_TimeDomain,EsN0);
        %% �źŽ�� IDFT���ź���ʱ��תΪƵ��
        Rx_FreqDomain=fft(Rx_TimeDomain)/sqrt(N_Subcarriers);
        %% �����Ȼ�������źŵĽ�ӳ��
        Distance = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 4
                % ��������ź���ԭ�źž���
                Distance(nm, ns)=(real(Rx_FreqDomain(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex]=min(Distance);  % �ҳ�ÿһ���е���С���뼰���Ӧ��������
        N_Error_Symbol=length(find(OutputIndex-IndexSymbol));% ÿһ�ε������� ���ҵ���Ϊ0�ĸ���
        Num_Error_Symbol=Num_Error_Symbol+N_Error_Symbol;% ���������� N_Subcarriers�����ִμӺ�
    end
    SER_Ep_QPSK(EsN0+1)=Num_Error_Symbol/ (N_OFDM_Frame*N_Subcarriers); %��ÿһ������ȼ���ʵ�ʷ���������
    SER_Th_QPSK(EsN0+1)=1-(1-0.5*erfc (sqrt (10^ (EsN0/10)/2)))^2; %��ÿһ������ȼ�������������
    N=length(snr_draw);
    figure(2) % �����ź�����ͼ
    for i=1:1:N
        if(EsN0==snr_draw(i))
            subplot(2,2,tool1);
            plot(Rx_FreqDomain,'*r');
            axis([-5, 5, -5, 5]);
            title(title_n(tool1));
            tool1=tool1+1;
        end
    end
    figure(6) %�����źŲ���
    for i=1:1:N   
        if(EsN0==snr_draw(i))
            subplot(2,2,tool3);
            plot(1:1:N_Subcarriers,real(Tx_TimeDomain));
            title(title_n(tool3));
            tool3=tool3+1;
        end
    end
    figure(3) %�����źŲ���
    for i=1:1:N   
        if(EsN0==snr_draw(i))
            subplot(2,2,tool2);
            plot(1:1:N_Subcarriers,real(Rx_TimeDomain));
            title(title_n(tool2));
            tool2=tool2+1;
        end
    end
end
figure(4) % ����������
semilogy(ber_begin:ber_end,SER_Th_QPSK,'-b*' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_QPSK,'-ro' );grid on; hold on;
title("����������ͼ")
%legend('QPSK Theoretical','QPSK Simulation');hold on;
axis([ber_begin ber_end 10^-6 1]);
xlabel('E_s/N_0(dB)');
ylabel( 'SER')
figure(5);
plot(mapper,'*r');
title('�����ź�����ͼ');
axis([-5, 5, -5, 5]);