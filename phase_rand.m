function [x_t,x_amp,sym_phase] = phase_rand(x,permutation)

% function [truecorr, pvals, nullcorrs] = PHASE_RAND_CORR(x,y,nscram, tail)
%
% This function calculates the correlation between 
%     [1] a phase-scrambled version of each column of 'x'
% and 
%     [2] the intact vector 'y'
% to produce a distribution of null correlations in which we have controlled
% for the power spectrum (and thus temporal autocorrelation) of the input time-series.
%
% INPUT
%   -x = [Nsamp by K] matrix, each of whose K columns will be phase-scrambled
%   -phermutation = 1-perform phase permutation, 0-generate random variable 
%
% OUTPUT
%   -x_t:  phase scrambled x
%   -x_amp: amplitudes of Fourier components of x
%   -sym_phase: symmetrized randomized phases
%
% Author: CJ Honey
% Version: 0.1, April 2010  (phase_scram_corr_Nvs1)
% Version: 0.2, March 2011    -- fixed bug with zero-th phase component;
%                             -- randomizes rather than scrambles phase 
%                             based on feedback from Jochen Weber
% modify Erez 2013
% cleaned up MLN 2016


% Extract number of samples and number of signals
[Nsamp, K] = size(x);  

% Transform the vectors-to-be-scrambled to the frequency domain
Fx = fft(x); 

% Udentify indices of positive and negative frequency components of the fft
% we need to know these so that we can symmetrize phase of neg and pos freq
if mod(Nsamp,2) == 0
    posfreqs = 2:(Nsamp/2);
    negfreqs = Nsamp : -1 : (Nsamp/2)+2;
else
    posfreqs = 2:(Nsamp+1)/2;
    negfreqs = Nsamp : -1 : (Nsamp+1)/2 + 1;
end

x_amp = abs(Fx);  %get the amplitude of the Fourier components
x_phase = atan2(imag(Fx), real(Fx)); %get the phases of the Fourier components [NB: must use 'atan2', not 'atan' to get the sign of the angle right]
J = sqrt(-1);  %define the vertical vector in the complex plane

% will contain symmetrized randomized phases for each bootstrap
sym_phase = x_phase;%zeros(Nsamp,K);



% Phase scramble
if permutation
    [~,rp] = sort(rand(Nsamp,1));
    x_phase=x_phase(rp,:);
    new_phase=x_phase(1:length(posfreqs),:);
else
    new_phase=2*pi*rand(length(posfreqs),K);
end

sym_phase(posfreqs,:) = sym_phase(posfreqs,:)+new_phase;
sym_phase(negfreqs,:) = sym_phase(negfreqs,:)-new_phase;

z = x_amp.*exp(J.*sym_phase); %generate (symmetric)-phase-scrambled Fourier components
x_t = real(ifft(z)); %invert the fft to generate a phase-scrambled version of x







