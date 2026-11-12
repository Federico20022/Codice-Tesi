% SETUPPATHS  Aggiunge al path le cartelle del progetto e definisce
% DATA_PATH e RES_PATH. Da lanciare dalla cartella principale.
root_path = fileparts(mfilename('fullpath'));
if isempty(root_path), root_path = pwd; end
addpath(fullfile(root_path, 'Functions'));
addpath(fullfile(root_path, 'Models'));
addpath(fullfile(root_path, 'Results'));
addpath(fullfile(root_path, 'Data'));
DATA_PATH = fullfile(root_path, 'Data');
RES_PATH  = fullfile(root_path, 'Results');
