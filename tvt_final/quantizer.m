function [sig_quantize, delta, codebook] = quantizer(sig, Nbits, sigMax, sigMin)
% uniform scalar quantization for real signals
sig_size = size(sig);
sig = sig(:);
delta = (sigMax-sigMin) / (2^Nbits);                                         % quantization interval 
b = (-2^Nbits/2 + 1) : 2^Nbits/2;
codebook = (-1/2 + b).*delta;                                                % quantization codebook generation 
index = 2^Nbits/2 + sign(sig).*min(2^Nbits/2,ceil(abs(sig)/delta)) + ceil((1-sign(sig))/2); % corresponding index in codebook
sig_quantize = reshape(codebook(index),sig_size); 
end

