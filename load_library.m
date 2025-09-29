function load_library()
% load_library Add relevant library folders to the Matlab search path.
%
% SYNOPSIS:
%   load_library
%
% DESCRIPTION:
%   Adds relevant library folders to the Matlab search path such that they
%   may be called from the current folder.
%
% REQUIRED PARAMETERS:
% 
% OPTIONAL PARAMETERS:
%
% RETURNS:
%
% DEPENDENCIES:
%
% See also 
% 
% REFERENCES
% 
% CONTACT INFORMATION
%  tobk@dtu.dk
% 
% AUTHORS
% Tobias K. S. Ritschel

% Let the user know that the library is being loaded
fprintf('Loading distributed delay differential equation simulator library .. ');

% Add real thermodynamics functions
addpath(genpath(fullfile(pwd, './src')));

% Let the user know that the library is being loaded
fprintf('Done\n');