function rirs = generate_rirs_with(method, room, source, array)
% generate_rirs_with Description
%    rirs = generate_rirs_with(method, room, source, array)
%
% Long description
%
%
switch method
case 'schimmel'
    rirs = generate_with_schimmel(room, source, array);
case 'vincent_mbss'
    rirs = generate_with_vincent(room, source, array);
case 'habets'
    rirs = generate_with_habets(room, source, array);
otherwise
    error([ method 'is not implemented'])
end
disp(method)
disp(size(rirs))
assert(size(rirs,1) == size(array.pos,1))

% function end: 'generate_rirs_with'
end

function rirs = generate_with_schimmel(in_room, in_source, in_array)
    room.dimension   = [ in_room.size(1),
                         in_room.size(2),
                         in_room.size(3)];   % room dimension (x,y,z)
    room.humidity    = 0.40;         % relative humidity (0,...,1)
    room.temperature = 20;           % room temperature (celsius)
    % surface absorption/diffusion coefficients
    room.surface.frequency  = [ 125, 250, 500, 1000, 2000, 4000, 8000];
    room.surface.absorption = [ repmat(in_room.walls_abs.south,  1,7); % south
                                repmat(in_room.walls_abs.north,  1,7); % north
                                    repmat(in_room.walls_abs.west,   1,7); % west
                                repmat(in_room.walls_abs.east,   1,7); % east
                                repmat(in_room.walls_abs.floor,  1,7); % floor
                                repmat(in_room.walls_abs.ceiling,1,7)];% ceiling
    room.surface.diffusion = [  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    roomsim_params.room = room;
    %    simulation options
    options.fs                  = in_room.Fs;              % sampling frequency in Hz
    options.responseduration    = in_room.max_sim_smpl/in_room.Fs; % duration of impulse response
    options.bandsperoctave      = 1;                    % simulation frequency accuracy (1, 2, 3, or 4 bands/octave)
    options.referencefrequency  = 125;                  % reference frequency for frequency octaves
    options.airabsorption       = true;                 % apply air absorption?
    options.distanceattenuation = true;                 % apply distance attenuation?
    options.subsampleaccuracy   = false;                % apply subsample accuracy?
    options.highpasscutoff      = 0;                    % 3dB frequency of high-pass filter (0=none)
    options.verbose             = false;                 % print status messages?
    % specular reflections simulation options
    options.simulatespecular    = true;                 % simulate specular reflections?
    options.reflectionorder     = [ in_room.k_refl,
                                    in_room.k_refl,
                                    in_room.k_refl]; % maximum specular reflection order (x,y,z)
    % diffuse reflections simulation options
    options.simulatediffuse     = in_room.do_diffusion;    % simulate diffuse reflections?
    options.numberofrays        = 2000;                 % number of rays in simulation (20*K^2)
    options.diffusetimestep     = 0.010;                % time resolution in diffuse energy histogram (seconds)
    options.rayenergyfloordB    = -80;                  % ray energy threshold (dB, with respect to initial energy)
    options.uncorrelatednoise   = true;                 % use uncorrelated poisson arrivals for binaural impulse responses?
    roomsim_params.options = options;
    % source
    [N,D] = size(in_source.pos);
    if D ~= 3
        error('Wrong dimension')
    end
    source(1).location    = in_source.pos;       % location of source (x,y,z; meters)
    source(1).orientation = [0, 0, 0];         % orientation of source (yaw,pitch,roll; degrees)
    source(1).description = 'omnidirectional';       % source type
    roomsim_params.source = source;
    disp(roomsim_params.source)
    % array
    [M,D] = size(in_array.pos);
    if D ~= 3
        error('Wrong dimensions')
    end
    for i = 1:M
        receiver(i).location    = in_array.pos(i,:);         % location of receiver (x,y,z; meters)
        receiver(i).orientation = [ 0 0 0 ];           % orientation of receiver (yaw,pitch,roll; degrees)
        receiver(i).description = 'omnidirectional';  % receiver type
    end
    roomsim_params.receiver = receiver;
    disp(roomsim_params.receiver)

    roomsim load omnidirectional;   % loads the receiver model
    roomsim load omnidirectional;   % loads the source model

    rirs =roomsim(roomsim_params);
    if iscell(rirs)
        rirs = cell2mat(rirs);
        rirs = single(rirs);
    end
    % check
    if any(isnan(rirs(:)))
        error('Wrong Impulse Response');
    end

    rirs = rirs';
end

function rirs = generate_with_vincent(in_room, in_source, in_array)
    fs  = in_room.Fs;
    k_refl = in_room.k_refl;
    room_size = in_room.size;
    F_abs = [ 125, 250, 500, 1000, 2000, 4000, 8000]';
    A = ...
       [repmat(in_room.walls_abs.south,  1,7); % south
        repmat(in_room.walls_abs.north,  1,7); % north
        repmat(in_room.walls_abs.west,   1,7); % west
        repmat(in_room.walls_abs.east,   1,7); % east
        repmat(in_room.walls_abs.floor,  1,7); % floor
        repmat(in_room.walls_abs.ceiling,1,7)]';% ceiling
    N = size(in_array.pos,1);
    sensor_xyz = in_array.pos';
    sensor_off =  repmat(0,3,N);
    sensor_type = repmat(1,1,N);
    M = size(in_source,1);
    source_xyz = in_source.pos';
    rir_length = in_room.max_sim_smpl;

    rirs = MBSS_roomsimove_single(fs,room_size,F_abs,A,sensor_xyz,sensor_off,sensor_type,source_xyz,k_refl,rir_length);
    rirs = rirs';
end

function rirs = generate_with_habets(in_room, in_source, in_array)
    c = 343;                  % Sound velocity (m/s)
    fs = in_room.Fs;          % Sample frequency (samples/s)
    r = in_array.pos;         % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
    s = in_source.pos;        % Source position [x y z] (m)
    L = in_room.size;         % Room dimensions [x y z] (m)
    beta = ...                % absorption coefficient
        [in_room.walls_abs.south,
        in_room.walls_abs.north,
        in_room.walls_abs.west,
        in_room.walls_abs.east,
        in_room.walls_abs.floor,
        in_room.walls_abs.ceiling]';
    n = in_room.max_sim_smpl; % Number of samples
    mtype = 'omnidirectional';% Type of microphone
    order = in_room.k_refl;   % -1 equals maximum reflection order!
    dim = 3;                  % Room dimension
    orientation = 0;          % Microphone orientation (rad)
    hp_filter = 1;            % Enable high-pass filter

    rirs = rir_generator(c, fs, r, s, L, beta, n, mtype, order, dim, orientation, hp_filter);

end
