 
function [ERRm,ERRs,ERRa,ERRd] = estimate_error(GTm,GTs,GTw,Em,Es,Ew)
%This function evaluates the accuracy of reconstruction of microphones,
%sources and walls position.
%Procrustes alignment must be run in order to align estimated and ground
%truth data before error estimation.

%GTm and Em:  12 x 3 matrices of ground truth and estimated microphone positions; 
%GTs and Es:  17 x 3 matrices of ground truth and estimated source positions;
%GTw and Ew:   6 x 3 matrices of planar surfaces positions. Each row is a
%vector normal to the corresponding plane, whose length is equal to the
%distance between the plane and the center of coordinates.
%ERRm :  1 x 12 matrix of microphone position errors
%ERRs :  1 x 17 matrix of source position errors
%ERRa :  1 x 6 matrix of angle errors in planar surfaces positions
%ERRd :  1 x 6 matrix of errors in distance from the coordinate center to
%the planar surface position.

for i = 1 : 12 
    ERRm(i) = norm(GTm(i,:)-Em(i,:));
end
for i = 1 : 17 
    ERRs(i) = norm(GTs(i,:)-Es(i,:));
end
for i = 1 : 6
    GTvec = squeeze(GTw(i,:));
    Evec = squeeze(Ew(i,:));
    ERRa(i) = (180/pi)*acos((Evec/norm(Evec))*(GTvec/norm(GTvec))');
    ERRd(i) =  norm(GTvec-Evec);
end
