function A=RT60toA(room_size,RT60)

% RT60toA Converts a given reverberation time into a single absorption
% coefficient for all surfaces
%
% A=RT60toA(room_size,RT60)
%
% Input:
% room_size: 3 x 1 vector containing the dimension of the room in m
% RT60: reverberation time in seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2003 Douglas R. Campbell
% Copyright 2008 Emmanuel Vincent
% This software is based on the Roomsim toolbox version 3.3 by Douglas R.
% Campbell , which was previously available at
% http://media.paisley.ac.uk/~campbell/Roomsim/
% It is distributed under the terms of the GNU Public License version 3
% (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Roomsimove, http://homepages.loria.fr/evincent/software/Roomsimove.zip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Norris-Eyring formula %%%
Lx=room_size(1); Ly=room_size(2); Lz=room_size(3);
V_room=Lx*Ly*Lz; % Volume of room m^3
Sxz=Lx*Lz; Syz=Ly*Lz; Sxy=Lx*Ly; S=2*(Sxz+Syz+Sxy); % Total area of shoebox room surfaces
A=1-exp(-0.1611*V_room/(S*RT60));

return