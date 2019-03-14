function H=roomsimove_single(room_sensor_config,source_xyz,source_off,source_dir)

% ROOMSIMOVE_SINGLE Compute shoebox room filters for a single source
% position given a room/sensor configuration file
%
% H=roomsimove_single(room_sensor_config,source_xyz)
%
% Input:
% room_sensor_config: room/sensor configuration file
% source_xyz: 3 x 1 vector representing the source position (x, y, z) 
% source_off: 3 x 1 vector representing the source orientation (azimuth,
% elevation, roll)
% source_dir: string corresponding to the source directivity
%
% Output:
% H: H_length x channels matrix containing the filters for all sensors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2003 Douglas R. Campbell
% Copyright 2008-2016 Emmanuel Vincent
% This software is a stripped-down version of the Roomsim toolbox version
% 3.3 by Douglas R. Campbell , which was previously available at
% http://media.paisley.ac.uk/~campbell/Roomsim/
% which has been extended in order to take source directivity into account
% It is distributed under the terms of the GNU Public License version 3
% (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Roomsimove, http://homepages.loria.fr/evincent/software/Roomsimove.zip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Differences with Roomsim version 3.3:
% - three or more sensors allowed
% - source directivity taken into account
% - faster implementation
% - distance and air attenuation applied to all image sources whatever
% their distance and modeled by a variable filter length (usually shorter)
% - no low-pass filter included in the fractional delay filters
% - no check of sensor/source positions (must be in the room and distant
% from at least 5 cm)
% - no check of simulation time or memory overflow
%
% Default parameters:
% humidity=40
% temperature=20 (corresponding to c=343m/s)
% order=-1 (default reflection order)
% H_length=-1 (default filter length)
% air_F=1 (air attenuation)
% dist_F=1 (distance attenuation)
% Fc_HP=20 (high-pass filter above 20Hz)
% smooth_F=1 (fractional delay filters)
% alpha_F=1 (surface opacity = reflectivity)


%%% Configuration data %%%
% Reading source config
if nargin < 3,
    source_off=zeros(size(source_xyz));
    source_dir='omnidirectional';
end
source_dir=load(source_dir);
source_dir=source_dir.S3D;
c_psi=cos(pi/180*source_off(1)); s_psi=sin(pi/180*source_off(1));
c_theta=cos(-pi/180.*source_off(2)); s_theta=sin(-pi/180.*source_off(2));
c_phi=cos(pi/180.*source_off(3)); s_phi=sin(pi/180.*source_off(3));
tm_source=[c_theta.*c_psi c_theta.*s_psi -s_theta;...
    s_phi.*s_theta.*c_psi-c_phi.*s_psi s_phi.*s_theta.*s_psi+c_phi.*c_psi s_phi.*c_theta;...
    c_phi.*s_theta.*c_psi+s_phi.*s_psi c_phi.*s_theta.*s_psi-s_phi.*c_psi c_phi.*c_theta];
% Reading room config
[param,value]=textread(room_sensor_config,'%9s%[^\n]','commentstyle','matlab');
Fs=eval(value{strcmp(param,'Fs')});
room_size=eval(['[' value{strcmp(param,'room_size')} ']'])';
F_abs=eval(['[' value{strcmp(param,'F_abs')} ']'])';
Ax1=eval(['[' value{strcmp(param,'Ax1')} ']'])'; Ax2=eval(['[' value{strcmp(param,'Ax2')} ']'])'; 
Ay1=eval(['[' value{strcmp(param,'Ay1')} ']'])'; Ay2=eval(['[' value{strcmp(param,'Ay2')} ']'])'; 
Az1=eval(['[' value{strcmp(param,'Az1')} ']'])'; Az2=eval(['[' value{strcmp(param,'Az2')} ']'])'; 
% Room absorption from 0 to Nyquist
A=[Ax1 Ax2 Ay1 Ay2 Az1 Az2]; clear Ax1 Ax2 Ay1 Ay2 Az1 Az2;
A=A(F_abs<=Fs/2,:); F_abs=F_abs(F_abs<=Fs/2);
if F_abs(1)~=0,
    A=[A(1,:); A]; F_abs=[0; F_abs];
end
if F_abs(end)~=Fs/2,
    A=[A; A(end,:)]; F_abs=[F_abs; Fs/2];
end
RR=length(F_abs);
B=realsqrt(1-A);
bx1=B(:,1); bx2=B(:,2); by1=B(:,3); by2=B(:,4); bz1=B(:,5); bz2=B(:,6);
% Reading sensor config
channels=0;
while any(strcmp(param,['sp' int2str(channels+1)])),
    channels=channels+1;
end
sensor_xyz=zeros(3,channels);
sensor_off=zeros(3,channels);
tm_sensor=zeros(3,3,channels);
for sensor_No=1:channels,
    sensor_xyz(:,sensor_No)=eval(['[' value{strcmp(param,['sp' int2str(sensor_No)])} ']'])';
    sensor_off(:,sensor_No)=eval(['[' value{strcmp(param,['so' int2str(sensor_No)])} ']'])';
    c_psi=cos(pi/180*sensor_off(1,sensor_No)); s_psi=sin(pi/180*sensor_off(1,sensor_No));
    c_theta=cos(-pi/180.*sensor_off(2,sensor_No)); s_theta=sin(-pi/180.*sensor_off(2,sensor_No));
    c_phi=cos(pi/180.*sensor_off(3,sensor_No)); s_phi=sin(pi/180.*sensor_off(3,sensor_No));
    tm_sensor(:,:,sensor_No)=[c_theta.*c_psi c_theta.*s_psi -s_theta;...
        s_phi.*s_theta.*c_psi-c_phi.*s_psi s_phi.*s_theta.*s_psi+c_phi.*c_psi s_phi.*c_theta;...
        c_phi.*s_theta.*c_psi+s_phi.*s_psi c_phi.*s_theta.*s_psi-s_phi.*c_psi c_phi.*c_theta];
end
% Air absorption and reverberation time
m_air=6.875e-4*(F_abs.'/1000).^(1.7);
atten_air=exp(-0.5*m_air.'); %attenuation factors for one metre travelled in air
Lx=room_size(1); Ly=room_size(2); Lz=room_size(3);
V_room=Lx*Ly*Lz; % Volume of room m^3
Sxz=Lx*Lz; Syz=Ly*Lz; Sxy=Lx*Ly; S=2*(Sxz+Syz+Sxy); % Total area of shoebox room surfaces
Se=Syz*(A(:,1)+A(:,2))+Sxz*(A(:,3)+A(:,4))+Sxy*(A(:,5)+A(:,6)); %Effective absorbing area of room surfaces at each frequency
a_bar=Se./S; % Mean absorption of each room surface
RT60=0.1611*V_room./(4*m_air'*V_room-S*log(1-a_bar)); % Norris-Eyring estimate adjusted for air absorption


%%% Simulated impulse responses %%%
% Constants
Two_pi=2*pi; % Compute here for efficiency
T=1/Fs; % Sampling Period
nyquist=Fs/2; % Half sampling frequency
Fs_c=Fs/343; % Samples per metre
% Reflection order and impulse response length
H_length = fix(max(RT60)*Fs); % H_length = longest reverberation time in samples (rounded down to integer)
range=H_length/Fs_c; % H_length in metres
order_x = ceil(range./(2.*Lx)); %  Number in +x direction
order_y = ceil(range./(2.*Ly)); %  Number in +y direction
order_z = ceil(range./(2.*Lz)); %  Number in +z direction
n_isources = (2.*order_x+1).*(2.*order_y+1).*(2.*order_z+1).*8; %Maximum number of image sources
delay_s=Fs_c*sqrt(sum((source_xyz*ones(1,channels)-sensor_xyz).^2));
H_length=max(H_length,ceil(max(max(delay_s)))+200); % Ensure H_length > 200 points so that a full CIPIC or MIT HRIR can be viewed
% Interpolation filter for fractional delays
N_frac = 32; % Order of FIR fractional delay filter
Tw = N_frac*T; % Window duration (seconds)
Two_pi_Tw=Two_pi/Tw; % Compute here for efficiency
t=(-Tw/2:T:Tw/2)'; % Filter time window NB column vector of length (N_frac+1) symmetrical about t=0
pad_frac = zeros(N_frac,1); % Column vector of zero values for post-padding
% Second order high-pass IIR filter to remove DC buildup (nominal -4dB cut-off at 20 Hz)
w=2*pi*20;
r1=exp(-w*T); r2=r1;
b1=-(1+r2); b2=r2; %Numerator coefficients (fix zeros)
a1=2*r1*cos(w*T); a2=-r1*r1; %Denominator coefficients (fix poles)
HP_gain=(1-b1+b2)/(1+a1-a2); %Normalisation gain
b_HP=[1 b1 b2]/HP_gain;
a_HP=[1 -a1 -a2];
% Further constants
Two_Lx=2*room_size(1); % Twice Length (Depth)
Two_Ly=2*room_size(2); % Twice Width
Two_Lz=2*room_size(3); % Twice Height
isource_ident=[-1 -1 -1; -1 -1 1; -1 1 -1; -1 1 1; 1 -1 -1; 1 -1 1; 1 1 -1; 1 1 1]; %codes the eight permutations of x+/-xp,y+/-yp,z+/-zp (the source to receiver vector components) where [-1 -1 -1] identifies the parent source.
surface_coeff=[0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]; % Includes/excludes bx,by,bz depending on 0/1 state.
qq=surface_coeff(:,1); %  for bx1
jj=surface_coeff(:,2); %  for by1
kk=surface_coeff(:,3); %  for bz1
F_abs_N=F_abs./nyquist; % Normalise the standard absorption frequency range for surfaces, (0 to 1) = (0 to Fs/2)
N_refl=2*round(nyquist/F_abs(2)); % Required length of FIR filter modelling impulse response of surface(+air)
Half_I=N_refl/2; % Half length of FIR filter model
Half_I_plusone = Half_I+1; % Half length shift required for FIR filter model of surface impulse response
window = 0.5.*(1 - cos(2*pi*(0:N_refl).'./N_refl)); % Compute the (N_refl+1) point column vector Hann window
% Image locations and impulse responses
isource_xyz = zeros(3,n_isources); % image source co-ordinates
refl = zeros(RR,n_isources); % surface reflection impulse amplitude (MEMORY CRITICAL POINT)
xx=isource_ident(:,1)*source_xyz(1); % partial x coord of image.
yy=isource_ident(:,2)*source_xyz(2); % partial y coord of image.
zz=isource_ident(:,3)*source_xyz(3); % partial z coord of image. 
xx_yy_zz=[xx yy zz]';
n_images=0; %number of significant images of each parent source
% Frequency dependent surface reflection and coordinates and distance for each image
for n=-order_x:order_x,
    bx2_abs_n=bx2.^abs(n); % Compute here for efficiency
    Two_n_Lx=n*Two_Lx; % Compute here for efficiency
    for l=-order_y:order_y,
        bx2y2_abs_nl=bx2_abs_n.*(by2.^abs(l)); % Compute here for efficiency
        Two_l_Ly=l*Two_Ly; % Compute here for efficiency
        for m=-order_z:order_z,
            bx2y2z2_abs_nlm=bx2y2_abs_nl.*(bz2.^abs(m)); % Compute here for efficiency
            Two_m_Lz=m*Two_Lz; % Compute here for efficiency
            Two_nlm_Lxyz = [Two_n_Lx; Two_l_Ly; Two_m_Lz]; % Concatenate here for efficiency
            for permu=1:8,
                n_images=n_images+1; %Accumulate count of the image sources
                % calculate xyz coordinates of image source n_images
                isource_xyz(:,n_images)=Two_nlm_Lxyz - xx_yy_zz(:,permu);
                % minimal delay to sensors in samples
                delay=min(Fs_c*sqrt(sum((isource_xyz(:,n_images)*ones(1,channels)-sensor_xyz).^2)));
                if delay <= H_length, % compute only for image sources within impulse response length
                    refl(:,n_images)=bx1.^abs(n-qq(permu)).*by1.^abs(l-jj(permu)).*bz1.^abs(m-kk(permu)).*bx2y2z2_abs_nlm;                    
                    if sum(refl(:,n_images)) < 1E-6, % (NB refl always +ve for air to surface, otherwise need abs here)
                        n_images=n_images-1; % Delete image sources with a sum of reflection coeffs below 1*10^-6 i.e. -120dB
                    end;
                else
                    n_images=n_images-1; % Delete image sources with a delay > impulse response length H_length
                end
            end
        end
    end
end
% Complete impulse response for the source
isource_xyz = isource_xyz(:,1:n_images); % Re-Allocate array for image source co-ordinates (discard trailing zero values)
refl = refl(:,1:n_images); % Re-Allocate array for surface reflection impulse amplitude (discard trailing zero values)
H = zeros(H_length,channels);
for sensor_No = 1:channels, % For each sensor
    % Get the sensor direction-dependent impulse responses
    sensor_dir=eval(value{strcmp(param,['sd' int2str(sensor_No)])});
    sensor_dir=load(sensor_dir);
    sensor_dir=sensor_dir.S3D;
    for is = 1:n_images, % for each of the n_images image sources
        b_refl = refl(:,is);
        xyz=isource_xyz(:,is)-sensor_xyz(:,sensor_No); % Position vector from sensor_No to source(is)                      
        dist = norm(xyz,2); % Distance (m) between image source(is) and sensor_No
        b_refl = b_refl./dist; % Include effect of distance (ie. 1/R) attenuation
        b_refl=b_refl.*(atten_air.^dist); % Include the absorption due to air
        % Estimate the values of reflection coefficient at the linear interpolated grid points
        b_refl=interp1q(F_abs_N,b_refl,1/Half_I*(0:Half_I).');
        % Half spectrum of data b_refl is now made conjugate-symmetric about Nyquist frequency,
        % and last data point discarded to make periodic spectrum corresponding to a real data sequence.
        b_refl=b_refl([1:Half_I+1 Half_I:-1:2]);
        % Transform surface data from frequency response to impulse response.
        h_refl = real(ifft(b_refl,N_refl)); % IFFT to calculate impulse response column vector of length N_refl samples
        % Make the impulse realisable (half length shift) and Hann window it 
        h_refl = window.*[h_refl(Half_I_plusone:N_refl); h_refl(1:Half_I_plusone)];
        if (n_images==1) || max(abs(h_refl(1:Half_I_plusone))) >= 1E-5, % For primary sources, and image sources with impulse response peak magnitudes >= -100dB (1/100000)
            % Fractional delay filter
            delay = Fs_c*dist; % delay in samples = (Samples per metre)*Distance
            rdelay=round(delay); % Extract integer delay (concatenated later with impulse response)
            t_Td=t-(delay-rdelay).*T; % Take account of fractional delay  -0.5 < D < +0.5 sample period
            hsf=.5*(1+cos(Two_pi_Tw*t_Td)).*sinc(Fs*t_Td); % Compute delayed filter impulse response for sensor
            h=filter(hsf,1,[h_refl; pad_frac]); % Convolve channel signals
            len_h=length(h); % length of impulse response modelling image source response
            adjust_delay = rdelay - ceil(len_h./2); % Half length shift to remove delay due to impulse response
            % Sensor filter
            xyz_source=tm_sensor(:,:,sensor_No)*xyz;% position vector from each sensor location to each image source in sensor axes system
            hyp=sqrt(xyz_source(1)^2+xyz_source(2)^2); % Distance (m) between sensor_No and proj of image source on xy plane 
            elevation=atan(xyz_source(3)./(hyp+eps)); % Calculate -pi/2 <= elevation <= +pi/2 rads
            azimuth=atan2(xyz_source(2),xyz_source(1)); % Calculate -pi <= azimuth <= +pi rad
            e_index=round(elevation.*180/pi)+91; 
            a_index=round(azimuth.*180/pi)+181;
            sensor_ir=sensor_dir{e_index,a_index};
            h=filter(sensor_ir,1,[h; zeros(length(sensor_ir)-1,1)]);
            % Source filter
            xyz_sensor=-tm_source*xyz;% position vector from each image source location to each sensor in source axes system
            hyp=sqrt(xyz_sensor(1)^2+xyz_sensor(2)^2); % Distance (m) between image source and proj of sensor_No on xy plane
            elevation=atan(xyz_sensor(3)./(hyp+eps)); % Calculate -pi/2 <= elevation <= +pi/2 rads
            azimuth=atan2(xyz_sensor(2),xyz_sensor(1)); % Calculate -pi <= azimuth <= +pi rad
            e_index=round(elevation.*180/pi)+91; 
            a_index=round(azimuth.*180/pi)+181;
            source_ir=source_dir{e_index,a_index};
            h=filter(source_ir,1,[h; zeros(length(source_ir)-1,1)]);
            len_h=length(h);
            %Accumulate the impulse responses from each image source within an array of length H_length
            start_index_Hp = max(adjust_delay+1+(adjust_delay>=0),1);
            stop_index_Hp = min(adjust_delay+1+len_h,H_length);
            start_index_h = max(-adjust_delay,1);
            stop_index_h = start_index_h + (stop_index_Hp - start_index_Hp);
            %Add whole or part of impulse response
            H(start_index_Hp:stop_index_Hp, sensor_No)= H(start_index_Hp:stop_index_Hp, sensor_No) + h(start_index_h:stop_index_h);
        end
    end
    %High-pass filtering
    H(:,sensor_No)=filter(b_HP,a_HP,H(:,sensor_No));
end

return
