clear all
close all
clc


%% Script performs Data-driven Kronecker-product wavelet-Fisz Despeckling 




x = 'DATA'; % .mat file first dimension (time), second and third (space)


smin=0; smax=255;

x=(x-min(x(:)))*smax/(max(x(:))-min(x(:)) )+smin;

time_averaged = squeeze(mean(x,1));
Reconstructed = ddhf_est_Iso_Kronecker_soft_unk(x,time_averaged);


