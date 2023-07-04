function para = parameter_init()
%% System Parameters
para.M = 200;                                                               % number of measurements
para.K = 500;                                                               % number of total devices
para.Ka = 50;                                                               % number of active devices
para.Nr = 512;                                                              % number of total antennas
para.Num_sub = 2;                                                           % number of subarrays
para.Nr_sub = para.Nr/para.Num_sub;                                         % number of antennas in each subarray
para.Sa = para.Ka / para.K;                                                 % sparsity ratio
para.Nbits = 2;                                                             % number of quantization bits
para.Nmax = 1;                                                              % quantization range
para.intvl = 32;                                                            % interval of high-resolution ADCs
para.fc = 6e9;                                                              % carrier frenquency
para.c = 3e8;                                                               % speed of light
para.Bw = 1e6;                                                              % system bandwidth
para.Ps_dBm = 5;                                                            % maximum transmit power in dBm
para.Psd_n = -174;                                                          % AWGN Psd dBm/Hz
para.Pn = 10^(para.Psd_n/10)*1e-3*para.Bw;                                  % noise power
para.r_max = 100;                                                           % maximum of distance between device and BS
para.r_min = 10;                                                            % minimum of distance between device and BS

%% Algorithm Parameters
para.T1 = 30;                                                               % maximum number of iterations in State 1
para.T2 = 20;                                                               % maximum number of iterations in State 2
para.damping = 0.4;                                                         % damping factor
para.threshold = 0.5;                                                       % threshold
para.niter = para.T1+para.T2;
end
