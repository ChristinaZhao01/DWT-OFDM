clc; clear all; close all;
% QPSK索引序列，星座图
mapper=[1/sqrt(2)+1i/sqrt(2) -1/sqrt(2)+1i/sqrt(2) 
       1/sqrt(2)-1i/sqrt(2) -1/sqrt(2)-1i/sqrt(2)];
N_OFDM_Frame=100;           % OFDM信号的个数
N_Subcarriers=1024;         % 子载波长度
ber_begin=0;ber_end=15;     % 起始、截止信噪比
snr_draw = [0,5,10,15];     % 画出0、5、10、15db信噪比下的相关信息
title_n = ["0dB","5dB","10dB","15dB"];
tool1 = 1;tool2 = 1;tool3 = 1;
for EsN0=ber_begin:ber_end
    Num_Error_Symbol=0;
    for n=1:N_OFDM_Frame
        %% 随机生成索引序列（QPSK）
        InputBits=randi([0 1], 2, N_Subcarriers);
        IndexSymbol = InputBits(1,:)*2+InputBits(2,:)+1;%将二维信号转换成一维索引序列
        %% 进行信号映射 
        Tx_FreqDomain=mapper(IndexSymbol); % QPSK信号
        %% FFT，将发送信号转换为时域
        Tx_TimeDomain=sqrt(N_Subcarriers) *ifft(Tx_FreqDomain);
        figure(1)
        if EsN0 == ber_begin
            plot(1:N_Subcarriers,real(Tx_TimeDomain));title("OFDM时域信号波形")
        end
        %% 通过高斯加性白噪声信道
        Rx_TimeDomain = awgn(Tx_TimeDomain,EsN0);
        %% 信号解调 IDFT将信号由时域转为频域
        Rx_FreqDomain=fft(Rx_TimeDomain)/sqrt(N_Subcarriers);
        %% 最大似然法进行信号的解映射
        Distance = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 4
                % 计算接收信号与原信号距离
                Distance(nm, ns)=(real(Rx_FreqDomain(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex]=min(Distance);  % 找出每一列中的最小距离及其对应的行索引
        N_Error_Symbol=length(find(OutputIndex-IndexSymbol));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol=Num_Error_Symbol+N_Error_Symbol;% 计算误码数 N_Subcarriers的总轮次加和
    end
    SER_Ep_QPSK(EsN0+1)=Num_Error_Symbol/ (N_OFDM_Frame*N_Subcarriers); %对每一个信噪比计算实际仿真误码率
    SER_Th_QPSK(EsN0+1)=1-(1-0.5*erfc (sqrt (10^ (EsN0/10)/2)))^2; %对每一个信噪比计算理论误码率
    N=length(snr_draw);
    figure(2) % 接收信号星座图
    for i=1:1:N
        if(EsN0==snr_draw(i))
            subplot(2,2,tool1);
            plot(Rx_FreqDomain,'*r');
            axis([-5, 5, -5, 5]);
            title(title_n(tool1));
            tool1=tool1+1;
        end
    end
    figure(6) %发送信号波形
    for i=1:1:N   
        if(EsN0==snr_draw(i))
            subplot(2,2,tool3);
            plot(1:1:N_Subcarriers,real(Tx_TimeDomain));
            title(title_n(tool3));
            tool3=tool3+1;
        end
    end
    figure(3) %接收信号波形
    for i=1:1:N   
        if(EsN0==snr_draw(i))
            subplot(2,2,tool2);
            plot(1:1:N_Subcarriers,real(Rx_TimeDomain));
            title(title_n(tool2));
            tool2=tool2+1;
        end
    end
end
figure(4) % 误码率曲线
semilogy(ber_begin:ber_end,SER_Th_QPSK,'-b*' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_QPSK,'-ro' );grid on; hold on;
title("误码率曲线图")
%legend('QPSK Theoretical','QPSK Simulation');hold on;
axis([ber_begin ber_end 10^-6 1]);
xlabel('E_s/N_0(dB)');
ylabel( 'SER')
figure(5);
plot(mapper,'*r');
title('发送信号星座图');
axis([-5, 5, -5, 5]);