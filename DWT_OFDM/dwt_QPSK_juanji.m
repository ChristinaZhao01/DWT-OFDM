clc; clear all; close all;
% 发端
% QPSK索引序列，星座图
mapper=[1/sqrt(2)+1i/sqrt(2) -1/sqrt(2)+1i/sqrt(2) 
       1/sqrt(2)-1i/sqrt(2) -1/sqrt(2)-1i/sqrt(2)];
N_OFDM_Frame=10;            % OFDM信号的个数
N_Subcarriers=1024;           % 子载波长度
M=4;                        %4PSK调制
ber_begin=0;ber_end=25;     % 起始、截止信噪比
snr_draw = [0,10,15,25];     % 画出0、5、10、15db信噪比下的相关信息
title_n = ["0dB","10dB","15dB","25dB"];
tool1 = 1;tool2 = 1;tool3 = 1;tool4 = 1;
L=7;                %卷积码约束长度
tblen=5*L;          %Viterbi译码器回溯深度

for EsN0=ber_begin:ber_end
    Num_Error_Symbol=0;
    Num_Error_Symbol_q = 0;
    Num_Error_Symbol_f = 0;
    Num_Error_Symbol_r = 0;
    Num_Error_Symbol_q_r = 0;
    for n=1:N_OFDM_Frame
        %%  随机生成索引序列（QPSK）
        InputBits=randi([0 1], 1, N_Subcarriers);
        %%  信道编码（卷积码编码）
        trellis = poly2trellis(7,[133 171]);       %(2,1,7)卷积编码
        code_data=convenc(InputBits,trellis);
        %% 串并变换
        data_temp1= reshape(code_data,log2(M),[])';             %以每组2比特进行分组，M=4
        data_temp2= bi2de(data_temp1)+1;                        %二进制转化为十进制
        %% 进行信号映射 
        Tx_FreqDomain=mapper(data_temp2);                       % QPSK信号
        Tx_FreqDomain1 = Tx_FreqDomain;
        %%  fft-ofdm
        Tx_TimeDomain1=sqrt(N_Subcarriers) *ifft(Tx_FreqDomain1);%FFT，将发送信号转换为时域
        Tx_TimeDomain1 = Tx_TimeDomain1.';        %并串转换
        %%  idwt，将发送信号转换为时域
        Tx_TimeDomain =idwt(Tx_FreqDomain,0*Tx_FreqDomain, 'haar'); % 采用db1小波基分解
        len = length(Tx_TimeDomain)
        Tx_TimeDomain = Tx_TimeDomain.';        %并串转换
        %% 信号通过瑞利信道
        x = randn(1,2*N_Subcarriers);
        y = randn(1,2*N_Subcarriers);
        for a=1:2*N_Subcarriers
            H(1,a)=sqrt((x(a))^2+(y(a))^2); 
        end
        fadeSig3 = H.*Tx_TimeDomain;
        Rx_TimeDomain_r = awgn(fadeSig3,EsN0);
        Rx_TimeDomain_r = Rx_TimeDomain_r.';
        %% 通过高斯加性白噪声信道
        Rx_TimeDomain = awgn(Tx_TimeDomain,EsN0);
        Rx_TimeDomain = Rx_TimeDomain.';
        Rx_TimeDomain1 = awgn(Tx_TimeDomain1,EsN0)
        Rx_TimeDomain1 = Rx_TimeDomain1.';
        %%  fft
        Rx_FreqDomain1=fft(Rx_TimeDomain1)/sqrt(N_Subcarriers);% 信号解调 FFT将信号由时域转为频域
        %% 信号解调 dwt将信号由时域转为频域
        [ca1, cd1] = dwt(Rx_TimeDomain, 'haar'); % 采用haar小波基分解
        [ca2, cd2] = dwt(Rx_TimeDomain_r, 'haar'); % 采用haar小波基分解
        %a1 = upcoef('a', ca1, 'haar', 1, N_Subcarriers); % 从系数得到近似信号
        %d1 = upcoef('d', cd1, 'haar', 1, N_Subcarriers); % 从系数得到细节信号
        Rx_FreqDomain = ca1; % 重构信号
        Rx_FreqDomain_r = ca2;
        figure(2)
        subplot(2, 2, 1); plot(ca1); title('一层小波分解的低频信息');
        subplot(2, 2, 2); plot(cd1); title('一层小波分解的高频信息');
        subplot(2, 2, 3); plot(Rx_TimeDomain, 'r-'); title('一层小波分解的重构信号');
        subplot(2, 2, 4); plot(Tx_TimeDomain, 'r-'); title('一层小波分解的重构信号');
         %% 最大似然法进行信号的解映射(通过高斯信道）
        Distance = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 4
                % 计算接收信号与原信号距离
                Distance(nm, ns)=(real(Rx_FreqDomain(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex]=min(Distance);  % 找出每一列中的最小距离及其对应的行索引
         %% 最大似然法进行信号的解映射（通过瑞利信道）
        Distance2 = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 4
                % 计算接收信号与原信号距离
                Distance2(nm, ns)=(real(Rx_FreqDomain_r(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain_r(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex2]=min(Distance2);  % 找出每一列中的最小距离及其对应的行索引
        %%  对fft-ofdm进行解映射
        Distance1 = zeros(length(mapper),length(N_Subcarriers));
        for ns=1:N_Subcarriers
            for nm = 1 : 4
                % 计算接收信号与原信号距离
                Distance1(nm, ns)=(real(Rx_FreqDomain1(ns))-real(mapper(nm)))^2+(imag(Rx_FreqDomain1(ns))-imag(mapper(nm)))^2; 
            end
        end
        [Z,OutputIndex1]=min(Distance1);  % 找出每一列中的最小距离及其对应的行索引
        data1 = reshape(OutputIndex1,[],1);
        data2 = de2bi(data1-1);
        Bit = reshape(data2',1,[]);
        N_Error_Symbolf=length(find(Bit-code_data));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol_f=Num_Error_Symbol_f+N_Error_Symbolf;% 计算误码数 N_Subcarriers的总轮次加和
        %%  维特比译码(通过高斯信道）
        De_data1 = reshape(OutputIndex,[],1);
        De_data2 = de2bi(De_data1-1);
        De_Bit = reshape(De_data2',1,[]);
        rx_c_de = vitdec(De_Bit,trellis,tblen,'trunc','hard');   %硬判决
        %%  计算误比特数（通过高斯信道）
        %   译码前的误比特数
        N_Error_Symbol_q=length(find(De_Bit-code_data));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol_q=Num_Error_Symbol_q+N_Error_Symbol_q;% 计算误码数 N_Subcarriers的总轮次加和
        %   译码后的误比特数
        N_Error_Symbol=length(find(rx_c_de-InputBits));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol=Num_Error_Symbol+N_Error_Symbol;% 计算误码数 N_Subcarriers的总轮次加和
        %%  维特比译码（通过瑞利信道）
        da1 = reshape(OutputIndex2,[],1);
        da2 = de2bi(da1-1);
        db = reshape(da2',1,[]);
        rx_c_de2 = vitdec(db,trellis,tblen,'trunc','hard');   %硬判决
        %%  计算误比特数（通过瑞利信道）
        %   译码前的误比特数
        N_Error_Symbol_q_r=length(find(db-code_data));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol_q_r=Num_Error_Symbol_q_r+N_Error_Symbol_q_r;% 计算误码数 N_Subcarriers的总轮次加和
        %   译码后的误比特数
        N_Error_Symbol_r=length(find(rx_c_de2-InputBits));% 每一次的误码数 即找到不为0的个数
        Num_Error_Symbol_r=Num_Error_Symbol_r+N_Error_Symbol_r;% 计算误码数 N_Subcarriers的总轮次加和
    end
    SER_Ep_QPSK(EsN0+1)=Num_Error_Symbol/ (N_OFDM_Frame*N_Subcarriers); %对每一个信噪比计算实际仿真误码率
    SER_Ep_QPSK_q(EsN0+1)=Num_Error_Symbol_q/ (N_OFDM_Frame*N_Subcarriers); %对每一个信噪比计算实际仿真误码率
    SER_Ep_QPSK_f(EsN0+1)=Num_Error_Symbol_f/ (N_OFDM_Frame*N_Subcarriers);
    SER_Ep_QPSK_r(EsN0+1)=Num_Error_Symbol_r/ (N_OFDM_Frame*N_Subcarriers);
    SER_Ep_QPSK_q_r(EsN0+1)=Num_Error_Symbol_q_r/ (N_OFDM_Frame*N_Subcarriers);
    SER_Th_QPSK(EsN0+1)=1-(1-0.5*erfc (sqrt (10^ (EsN0/10)/2)))^2; %对每一个信噪比计算理论误码率
    N=length(snr_draw);
    figure(3) % 接收信号星座图（通过高斯信道）
    for i=1:1:N
        if(EsN0==snr_draw(i))
            subplot(2,2,tool1);
            plot(Rx_FreqDomain,'*r');
            axis([-5, 5, -5, 5]);
            title(title_n(tool1));
            tool1=tool1+1;
        end
    end
    figure(7) % 接收信号星座图（通过瑞利信道）
    for i=1:1:N
        if(EsN0==snr_draw(i))
            subplot(2,2,tool4);
            plot(Rx_FreqDomain_r,'*r');
            axis([-5, 5, -5, 5]);
            title(title_n(tool4));
            tool4=tool4+1;
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
semilogy(ber_begin:ber_end,SER_Ep_QPSK_q,'-mo' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_QPSK_f,'-g*' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_QPSK_r,'-ko' );grid on; hold on;
semilogy(ber_begin:ber_end,SER_Ep_QPSK_q_r,'-k*' );grid on; hold on;
title("误码率曲线图")
axis([ber_begin,ber_end,10^-7 ,1]);
xlabel('E_s/N_0(dB)');
ylabel( 'SER')
figure(6);
plot(mapper,'*r');
title('发送信号星座图');
axis([-5, 5, -5, 5]);