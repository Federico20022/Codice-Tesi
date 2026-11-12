function rp = RunAll_path()
% Imposta i path del progetto e restituisce la root.
rp = fileparts(mfilename('fullpath'));
if isempty(rp), rp = pwd; end
addpath(fullfile(rp,'Functions'), fullfile(rp,'Models'), ...
        fullfile(rp,'Results'),   fullfile(rp,'Data'));
if ~exist(fullfile(rp,'Results'),'dir'), mkdir(fullfile(rp,'Results')); end
assignin('base','RES_PATH',  fullfile(rp,'Results'));
assignin('base','DATA_PATH', fullfile(rp,'Data'));
end
