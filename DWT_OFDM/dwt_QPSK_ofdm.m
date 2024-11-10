clc; clear all; close all;
% 发端
% QPSK索引序列，星座图
mapper=[1/sqrt(2)+1i/sqrt(2) -1/sqrt(2)+1i/sqrt(2) 
       1/sqrt(2)-1i/sqrt(2) -1/sqrt(2)-1i/sqrt(2)];
N_OFDM_Frame=100;           % OFDM信号的个数
N_Subcarriers=1024;         % 子载波长度
ber_begin=0;ber_end=15;     % 起始、截止信噪比
snr_draw = [0,5,10,15];     % 画出0、5、10、15db信噪比下的相关信息
title_n = ["0dB","5dB","10dB","15dB"];
tool1 = 1;tool2 = 1;tool3 = 1
for EsN0=ber_begin:ber_end
    Num_Error_Symbol=0;
    for n=1:N_OFDM_Frame
        %%  随机生成索引序列（QPSK）
        InputBits=randi([0 1], 2, N_Subcarriers);
        IndexSymbol = InputBits(1,:)*2+InputBits(2,:)+1;%将二维信号转换成一维索引序列
        %% 进行信号映射 
        Tx_FreqDomain=mapper(IndexSymbol); % QPSK信号
        %Tx_FreqDomain = Tx_FreqDomain.';
        %%  idwt，将发送信号转换为时域
        Tx_TimeDomain =idwt(Tx_FreqDomain,0*Tx_FreqDomain, 'haar'); % 采用haar小波基分解
        len = length(Tx_TimeDomain)
        %%  通过瑞利信道 
        x = randn(1,2*N_Subcarriers);
        y = randn(1,2*N_Subcarriers);
        h = (x+1i.*y)/sqrt(2); 
        for a=1:2*N_Subcarriers
            H(1,a)=sqrt((x(a))^2+(y(a))^2); 
        end
        fadeSig3 = H.*Tx_TimeDomain;
        Rx_TimeDomain = awgn(fadeSig3,EsN0);
        %% 通过高斯加性白噪声信道
        %Rx_TimeDomain = awgn(Tx_TimeDomain,EsN0);
        %%  瑞利信道
        %h = (randn(2*N_Subcarriers,1)+i*randn(2*N_Subcarriers,1))/sqrt(2);
        %Rx_TimeDomain = h.*Tx_TimeDomain;
        %RayleighChan = comm.RayleighChannel();
        %Rx_TimeDomain = RayleighChan(Tx_TimeDomain);
        %h_rayleigh = sqrt(1/2) *(randn(1,1) + 1i*randn(1,1));
        %Rx_TimeDomain = Tx_TimeDomain*h_rayleigh;
        %%  莱斯信道
        %RicianChan = comm.RicianChannel();
        %Rx_TimeDomain = RicianChan(Tx_TimeDomain);
        k = 10000;
        h_rician = sqrt(k/(k+1)) + sqrt(1/(k+1))* H;
        Rx_TimeDomain = Tx_TimeDomain.*h_rician;
        %% 信号解调 dwt将信号由时域转为频域
        [ca1, cd1] = dwt(Rx_TimeDomain, 'haar'); % 采用haar小波基分解
        a1 = upcoef('a', ca1, 'haar', 1, N_Subcarriers); % 从系数得到近似信号
        d1 = upcoef('d', cd1, 'haar', 1, N_Subcarriers); % 从系数得到细节信号

        Rx_FreqDomain = ca1; % 重构信号
        figure(2)
        subplot(2, 2, 1); plot(ca1); title('一层小波分解的低频信息');
        subplot(2, 2, 2); plot(cd1); title('一层小波分解的高频信息');
        subplot(2, 2, 3); plot(Rx_TimeDomain, 'r-'); title('一层小波分解的重构信号');
        subplot(2, 2, 4); plot(Tx_TimeDomain, 'r-'); title('一层小波分解的重构信号');
         %% 最大似然法进行信号的解映射
        Distance = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 4
                % 计算接收信号与原信号距离
                Distance(nm, ns)=(real(Rx_FreqDomain(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex]=min(Distance);  % 找出每一列中的最小距离及其对应的行索引
        %%  计算误比特数
        N_Error_Symbol=length(find(OutputIndex-IndexSymbol));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol=Num_Error_Symbol+N_Error_Symbol;% 计算误码数 N_Subcarriers的总轮次加和
    end
    SER_Ep_QPSK(EsN0+1)=Num_Error_Symbol/ (N_OFDM_Frame*N_Subcarriers); %对每一个信噪比计算实际仿真误码率
    SER_Th_QPSK(EsN0+1)=1-(1-0.5*erfc (sqrt (10^ (EsN0/10)/2)))^2; %对每一个信噪比计算理论误码率
    N=length(snr_draw);
    figure(3) % 接收信号星座图
    for i=1:1:N
        if(EsN0==snr_draw(i))
            subplot(2,2,tool1);
            plot(Rx_FreqDomain,'*r');
            axis([-5, 5, -5, 5]);
            title(title_n(tool1));
            tool1=tool1+1;
        end
    end
    figure(1) %发送信号波形
    for i=1:1:N   
        if(EsN0==snr_draw(i))
            subplot(2,2,tool2);
            plot(1:1:len,real(Tx_TimeDomain));
            title(title_n(tool2));
            tool2=tool2+1;
        end
    end
    figure(4) %接收信号波形
    for i=1:1:N   
        if(EsN0==snr_draw(i))
            subplot(2,2,tool3);
            plot(1:1:len,real(Rx_TimeDomain));
            title(title_n(tool3));
            tool3=tool3+1;
        end
    end
end
figure(5) % 误码率曲线
semilogy(ber_begin:ber_end,SER_Th_QPSK,'-b*' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_QPSK,'-ro' );grid on; hold on;
title("误码率曲线图")
%legend('QPSK Theoretical','QPSK Simulation');hold on;
axis([ber_begin ber_end 10^-6 1]);
xlabel('E_s/N_0(dB)');
ylabel( 'SER')
figure(6);
plot(mapper,'*r');
title('发送信号星座图');
axis([-5, 5, -5, 5]);