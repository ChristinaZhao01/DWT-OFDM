function OFDM_16QAM_AWGN
%16QAM索引序列，星座图
mapper = [1/(3*sqrt(2))+1i/(3*sqrt(2)) 1/sqrt(2)+1i/(3*sqrt(2)) 1/(3*sqrt(2))+1i/sqrt(2) 1/sqrt(2)+1i/sqrt(2)
          -1/(3*sqrt(2))+1i/(3*sqrt(2)) -1/sqrt(2)+1i/(3*sqrt(2)) -1/(3*sqrt(2))+1i/sqrt(2) -1/sqrt(2)+1i/sqrt(2)
          1/(3*sqrt(2))-1i/(3*sqrt(2)) 1/sqrt(2)-1i/(3*sqrt(2)) 1/(3*sqrt(2))-1i/sqrt(2) 1/sqrt(2)-1i/sqrt(2)
          -1/(3*sqrt(2))-1i/(3*sqrt(2)) -1/sqrt(2)-1i/(3*sqrt(2)) -1/(3*sqrt(2))-1i/sqrt(2) -1/sqrt(2)-1i/sqrt(2)]
N_OFDM_Frame=100; % OFDM信号的个数
N_Subcarriers=1024;% 子载波长度
ber_begin=0;ber_end=15;% 起始、截止信噪比
snr_draw = [0,5,10,15];
title_n = ["0dB","5dB","10dB","15dB"];
tool1 = 1;tool2 = 1;
for EsN0=ber_begin:ber_end
    Num_Error_Symbol=0;
    for n=1:N_OFDM_Frame
        %随机生成索引序列（16QAM）
        IndexSymbol=randi([1 16], 1, N_Subcarriers);
          % 进行信号映射 
        Tx_FreqDomain=mapper(IndexSymbol); % 16QAM信号
        Tx_TimeDomain=sqrt(N_Subcarriers) *ifft(Tx_FreqDomain);%FFT，将发送信号转换为时域
        figure(6)
        if EsN0 == ber_begin
            plot(1:N_Subcarriers,real(Tx_TimeDomain));title("OFDM时域信号波形")
        end
        %通过高斯白噪声信道
        Rx_TimeDomain = awgn(Tx_TimeDomain,EsN0);
        Rx_FreqDomain=fft(Rx_TimeDomain)/sqrt(N_Subcarriers);% 信号解调 IDFT将信号由时域转为频域
        Distance = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 16
                %计算接收信号与原信号距离
                Distance(nm, ns)=(real(Rx_FreqDomain(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex]=min(Distance);  % 找出每一列中的最小距离及其对应的行索引
        N_Error_Symbol=length(find(OutputIndex-IndexSymbol));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol=Num_Error_Symbol+N_Error_Symbol;%计算误码数 N_Subcarriers的总轮次加和
    end
    SER_Ep_16QAM(EsN0+1)=Num_Error_Symbol/ (N_OFDM_Frame*N_Subcarriers);  %对每一个信噪比计算实际仿真误码率
    SER_Th_16QAM(EsN0+1)=1-(1-(3/4)*erfc(sqrt((3/2)*10^(EsN0/10)/(16-1)))).^2;%对每一个信噪比计算理论误码率
    N=length(snr_draw);
    figure(7) %接收信号星座图
    for i=1:1:N
        if(EsN0==snr_draw(i))
            subplot(2,2,tool1);
            plot(Rx_FreqDomain,'*r');
            axis([-5, 5, -5, 5]);
            title(title_n(tool1));
            tool1=tool1+1;
        end
    end
    figure(8) %接收信号波形
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
semilogy(ber_begin:ber_end,SER_Th_16QAM,'-k*' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_16QAM,'-mo' );grid on; hold on;
title("误码率曲线图")
%legend('16QAM Theoretical','16QAM Simulation');
axis([ber_begin ber_end 10^-6 1]);
xlabel('E_s/N_0(dB)');
ylabel( 'SER')
figure(9);
plot(mapper,'*r');
title('发送信号星座图');
axis([-5, 5, -5, 5]);