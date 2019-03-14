function [time,HH]=roomsimove(room_sensor_config,source_config,fs)

% ROOMSIMOVE Compute shoebox room filters for a moving source described
% by its position at discrete instants obtained from a configuration file
%
% [time,HH]=roomsimove(room_sensor_config,source_config,fs)
%
% Input:
% room_sensor_config: room/sensor configuration file
% source_config: source movements configuration file
% fs: sampling frequency in Hz
%
% Output:
% time: n_samples x 1 vector listing the sampling instants
% HH: H_length x channels x n_samples matrix containing the associated filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008-2016 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Roomsimove, http://homepages.loria.fr/evincent/software/Roomsimove.zip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% Reading source config
[param,value]=textread(source_config,'%9s%[^\n]','commentstyle','matlab');
ptime=eval(['[' value{strcmp(param,'ptime')} ']'])';
px=eval(['[' value{strcmp(param,'px')} ']'])'; py=eval(['[' value{strcmp(param,'py')} ']'])'; pz=eval(['[' value{strcmp(param,'pz')} ']'])';
source_xyz=[px py pz]';
if any(strcmp(param,'pd')),
    pa=eval(['[' value{strcmp(param,'pa')} ']'])'; pe=eval(['[' value{strcmp(param,'pe')} ']'])'; pr=eval(['[' value{strcmp(param,'pr')} ']'])';
    source_off=[pa pe pr]';
    source_dir=eval(value{strcmp(param,'pd')});
end
n_breaks=size(source_xyz,2);

%%% Spatial sampling %%%
% First breakpoint
time=ptime(1);
HH=roomsimove_single(room_sensor_config,source_xyz(:,1),source_off(:,1),source_dir);
% Subsequent breakpoints
for b=1:n_breaks-1,
    breakdist=norm(source_xyz(:,b+1)-source_xyz(:,b),2);
    if breakdist > 0,
        n_samples=ceil(breakdist/sample_dist);
        for s=1:n_samples,
            time(end+1)=((n_samples-s)*ptime(b)+s*ptime(b+1))/n_samples;
            HH(:,:,end+1)=roomsimove_single(room_sensor_config,((n_samples-s)*source_xyz(:,b)+s*source_xyz(:,b+1))/n_samples);
        end
    else
        time(end+1)=ptime(b+1);
        HH(:,:,end+1)=HH(:,:,end);
    end
end

return