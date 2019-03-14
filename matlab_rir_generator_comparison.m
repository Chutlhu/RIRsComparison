MATLAB_RIRS_SIMs = {'schimmel', 'habets', 'vincent_mbss'};
n_simulators = length(MATLAB_RIRS_SIMs);

PATH_TO_SIMs = './rir_generators/matlab/';

% all the coordinates are in x, y, z
% all the angle are in r, az, el

%% audio scene
%% signals
room.Fs = 48000;
room.max_sim_smpl = floor(0.5*room.Fs);

% room:
room.size = [3,5,4]; % 1 x D [x, y, z]
walls_abs.west    = 1;
walls_abs.east    = 1;
walls_abs.south   = 1;       % west, east, south, north, floor, ceiling ...
walls_abs.north   = 1;       % as in pyroomacoustics
walls_abs.floor   = 0.1;
walls_abs.ceiling = 1;
room.walls_abs = walls_abs;
room.k_refl = 20;
room.do_diffusion = false;
room.c = 343;
% source
source.pos = [2,2,2;]; % N x D [x, y, z]
N = 1;

% anthenna
mic_bar = [1,1,0.2];
mic1 = mic_bar + [0.5,0.5,0];
mic2 = mic_bar - [0.3,0.3,0];
M = 2;
array.pos = [mic1; mic2]; % M x D [x, y, z]

% check dimension
assert(all(size(room.size) == [1,3]));
assert(all(size(source.pos)== [N,3]));
assert(all(size(array.pos) == [M,3]));

%% teoretical output
dist1 = norm(mic1-source.pos);
dist2 = norm(mic2-source.pos);
tau11 = round(room.Fs*dist1/room.c);
tau21 = round(room.Fs*dist2/room.c);
alpha1 = 1/(4*pi*dist1^2);
alpha2 = 1/(4*pi*dist2^2);

%% output
rirs = zeros(n_simulators,M,room.max_sim_smpl);

for i = 1:n_simulators
    current_sim = [PATH_TO_SIMs, MATLAB_RIRS_SIMs{i}];
    addpath(current_sim);

    rirs(i,:,:) = generate_rirs_with(MATLAB_RIRS_SIMs{i}, room, source, array);

    rmpath(current_sim)
end

%% plot result
n_max = floor(4*max(tau11,tau21));
for i = 1:n_simulators
    figure(i)
    plot(squeeze(rirs(i,1,1:n_max)))
    hold on
    plot(squeeze(rirs(i,2,1:n_max)))
    text(tau11,alpha1,'\leftarrow \tau_1')
    text(tau21,alpha2,'\leftarrow \tau_2')
    hold off
end

n_max = floor(4*max(tau11,tau21));
figure(n_simulators+1)
for i = 1:n_simulators
    subplot(211)
    plot(squeeze(rirs(i,1,1:n_max)))
    if i == 1; text(tau11,alpha1,'\leftarrow \tau_1'); end;
    hold on
    subplot(212)
    plot(squeeze(rirs(i,2,1:n_max)))
    if i == 1; text(tau21,alpha2,'\leftarrow \tau_2'); end;
    hold on
end
subplot(211)
legend(MATLAB_RIRS_SIMs)
subplot(212)
legend(MATLAB_RIRS_SIMs)
hold off
