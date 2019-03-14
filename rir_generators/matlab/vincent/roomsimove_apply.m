function x=roomsimove_apply(time,HH,s,fs)

% ROOMSIMOVE_APPLY Filters a signal by precomputed simulated room impulse
% responses, performing further spatial interpolation
%
% x=roomsimove_apply(time,HH,s,fs)
%
% Input:
% time: n_samples x 1 vector listing the sampling instants
% HH: H_length x channels x n_samples matrix containing the associated filters
% s: s_length x 1 vector containing the signal to be filtered (s_length
% must be larger or equal to time(end)*fs)
% fs: sampling frequency in Hz
%
% Output:
% x: s_length x channels matrix containing the filtered signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Roomsimove, http://homepages.loria.fr/evincent/software/Roomsimove.zip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Spatial interpolation %%%
[H_length,channels,n_samples]=size(HH);
% Interpolation factor between the sampled points
interp_factor=3;
itime=zeros((n_samples-1)*interp_factor+2,1); % interpolated filters
iHH=zeros(H_length,channels,(n_samples-1)*interp_factor+1); % change times
itime(1)=time(1); itime(end)=time(end);
iHH(:,:,1)=HH(:,:,1);
for n=1:n_samples-1,
    for i=1:interp_factor,
        itime((n-1)*interp_factor+i+1)=((interp_factor-i+.5)*time(n)+(i-.5)*time(n+1))/interp_factor;
        iHH(:,:,(n-1)*interp_factor+i+1)=((interp_factor-i)*HH(:,:,n)+i*HH(:,:,n+1))/interp_factor;
    end
end

%%% Filtering %%%
s_length=length(s);
x=zeros(s_length+H_length-1,channels);
for n=1:(n_samples-1)*interp_factor+1,
    tmin=max(1,floor(itime(n)*fs)+1);
    tmax=min(s_length,floor(itime(n+1)*fs));
    x(tmin:tmax+H_length-1,:)=x(tmin:tmax+H_length-1,:)+fftfilt(iHH(:,:,n),[s(tmin:tmax); zeros(H_length-1,1)]);
end
x=x(1:s_length,:);

return