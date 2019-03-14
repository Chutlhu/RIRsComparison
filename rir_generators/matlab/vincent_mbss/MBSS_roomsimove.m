function [time,HH]=MBSS_roomsimove(fs,room_size,F_abs,A,sensor_xyz,sensor_off,sensor_type,ptime,source_xyz)

% ROOMSIMOVE Compute shoebox room filters for a moving source described
% by its position at discrete instants
%
% [time,HH]=MBSS_roomsimove(room_sensor_config,source_config,fs)
%
% Inputs :
% fs : sampling frequency (in Hz)
% room_size : 1 x 3 , room dimensions (in meters)
% F_abs : nfreq x 1, frequencies to define absorption coefficients
% A : nfreq x 6, frequency-dependant absorption coefficients for each
% surface of the room
% sensor_xyz : 3 x nchan, cartesian coordinates of the nchan microphones
% sensor_off : 3 x nchan, Sensor directions (azimuth, elevation and roll
% offset in degrees)
% sensor_type : 1 x nchan, sensor type (1 = omnidirectional / 2 = cardioid)
% ptime : 1 x timeTick, time stamps array for corresponding source position
% (in seconds)
% source_xyz : 3 x timeTick, corresponding source position over ptime
%
% Outputs :
% time: n_samples x 1 vector listing the sampling instants
% HH: H_length x channels x n_samples matrix containing the associated filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008 Emmanuel Vincent
% Copyright 2015 Ewen Camberlein and Romain Lebarbenchon
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Roomsimove, http://www.irisa.fr/metiss/members/evincent/Roomsimove.zip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Differences since copyright 2008 (roomsimove.m):
% - input parameters (no usage of .txt file for parameters storage)

%%% Configuration data %%%
% Maximum distance in m between the sampled points
sample_dist=1e-2*16000/fs;
% When the true source moves a certain distance, the variation of distance
% between the image sources and the sensors is at most equal to this
% distance. Hence it is sufficient to check that this distance corresponds
% to a small enough delay to ensure that the movement is accurately sampled
% for all image sources.
% Distance    Delay    L2 error on a pure delay filter, integrated over the path (worst case with delay close to zero)
%   1mm    0.05 sample    33 dB without interpolation
%                         61 dB with maximum interpolation (interp_factor >= 26)
%   4mm    0.2 sample     21 dB without interpolation
%                         37 dB with maximum interpolation (interp_factor >= 7)
%   1cm    0.5 sample     13 dB without interpolation
%                         22 dB with maximum interpolation (interp_factor >= 3)
% Reading all data into variables

n_breaks=length(ptime);

%%% Spatial sampling %%%
% First breakpoint
time=ptime(1);
HH=MBSS_roomsimove_single(fs,room_size,F_abs,A,sensor_xyz(:,:,1),sensor_off,sensor_type,source_xyz(:,1));
% Subsequent breakpoints
for b=1:n_breaks-1
    fprintf('    - %d / %d filter(s) generation\n',b,n_breaks-1);
    
    % Moving source distance
    if(size(source_xyz,2)>1)
        breakdist_source  = norm(source_xyz(:,b+1)-source_xyz(:,b),2);
    else
        breakdist_source = 0;
    end
    
    % Moving array (sensors) distance
    if(size(sensor_xyz,3)>1)
        breakdist_sensors = norm(mean(sensor_xyz(:,:,b+1),2)-mean(sensor_xyz(:,:,b),2),2);
    else
        breakdist_sensors = 0;
    end
    
    % Moving source only
    if breakdist_source > 0 && breakdist_sensors == 0 
        n_samples=ceil(breakdist_source/sample_dist);
        for s=1:n_samples
            time(end+1)=((n_samples-s)*ptime(b)+s*ptime(b+1))/n_samples;
            HH(:,:,end+1)=MBSS_roomsimove_single(fs,room_size,F_abs,A,sensor_xyz,sensor_off,sensor_type,((n_samples-s)*source_xyz(:,b)+s*source_xyz(:,b+1))/n_samples);
        end
    % Moving sensors only
    elseif breakdist_source == 0 && breakdist_sensors > 0 
        n_samples=ceil(breakdist_sensors/sample_dist);
        for s=1:n_samples
            time(end+1)=((n_samples-s)*ptime(b)+s*ptime(b+1))/n_samples;
            HH(:,:,end+1)=MBSS_roomsimove_single(fs,room_size,F_abs,A,((n_samples-s)*sensor_xyz(:,:,b)+s*sensor_xyz(:,:,b+1))/n_samples,sensor_off,sensor_type,source_xyz);
        end
    % Moving sensors and source
    elseif breakdist_source > 0 && breakdist_sensors > 0 
        
        n_samples=ceil(min(breakdist_sensors,breakdist_source)/sample_dist);
        for s=1:n_samples
            time(end+1)=((n_samples-s)*ptime(b)+s*ptime(b+1))/n_samples;
            HH(:,:,end+1)=MBSS_roomsimove_single(fs,room_size,F_abs,A,((n_samples-s)*sensor_xyz(:,:,b)+s*sensor_xyz(:,:,b+1))/n_samples,sensor_off,sensor_type,((n_samples-s)*source_xyz(:,b)+s*source_xyz(:,b+1))/n_samples);
        end
    % Too small (distance) movement => we keep previous filter
    else
        time(end+1)=ptime(b+1);
        HH(:,:,end+1)=HH(:,:,end);
    end
end

return;