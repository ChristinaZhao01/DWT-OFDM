# DWT-OFDM Signal Transmission System Simulation and Performance Analysis

## Project contents

- Implemented a DWT-OFDM signal system by improving the conventional FFT-OFDM signal system using Discrete Wavelet Transform (DWT) and simulated and analysed the BER of the two systems using MATLAB.
- Converted the input binary signal into a low speed parallel transmission over multiple subcarriers, modulated using QPSK or 16QAM, converted from frequency to time domain using IFFT and IDWT, transmitted serially over noisy channels and finally the process is reversed at the receiver end to recover the original binary signal.
- Analysed the performance of the DWT-OFDM signal transmission system under Gaussian and Rayleigh channels respectively.
- Compared the BER of the DTW-OFDM signal system with and without convolutional coding of the input signal sequence at the transmitter side and Viterbi decoding at the receiver side to provide a comprehensive evaluation and optimisation of the system performance.

## Programming Platform


 - Matlab

## Simulation result diagram
### Graph of total QPSK BER for different cases
![image](https://github.com/ChristinaZhao01/DWT-OFDM/blob/main/QPSK%E8%AF%AF%E7%A0%81%E7%8E%87%E6%80%BB%E6%9B%B2%E7%BA%BF%E5%9B%BE.png)
### Graph of total 16QAM BER for different cases
![image](https://github.com/ChristinaZhao01/DWT-OFDM/blob/main/16QAM%E8%AF%AF%E7%A0%81%E7%8E%87%E6%80%BB%E6%9B%B2%E7%BA%BF%E5%9B%BE.png)
