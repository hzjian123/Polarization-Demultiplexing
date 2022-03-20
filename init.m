function init
% path management

% set MATLAB encoding character set
setenv('LANG','en_US.UTF-8');

% addpath
self_dir = fileparts(mfilename('fullpath'));
addpath(self_dir);
addpath(genpath(fullfile(self_dir, 'contrib')));
addpath(genpath(fullfile(self_dir, 'corefun')));
addpath(genpath(fullfile(self_dir, 'utils')));
addpath(genpath(fullfile(self_dir, 'tests')));
addpath(genpath(fullfile(self_dir, 'toolbox')));
addpath(genpath(fullfile(self_dir, 'examples')));
if exist(fullfile(self_dir, 'deprecated'), 'dir')
  addpath(genpath(fullfile(self_dir, 'deprecated')));
end