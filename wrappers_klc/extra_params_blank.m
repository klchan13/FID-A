% Makes a structure containing the extra parameters needed for the pulse
% sequence input into spectra_sim.m.
%
% Possible inputs:
% ----------------
% refocWaveform     = name of refocusing pulse waveform.
% editWaveform      = name of editing pulse waveform.
% editOnFreq        = frequency of edit on pulse[ppm]
% editOffFreq       = frequency of edit off pulse[ppm]
% refTp             = duration of refocusing pulses[ms]
% editTp            = duration of editing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% x                 = vector of X positions to simulate [cm]
% y                 = vector of y positions to simulate [cm]
% taus              = vector of pulse sequence timings  [ms]
% spinSys           = spin system to simulate 
% editPhCyc1        = vector of phase cycling steps for 1st editing pulse [degrees]
% editPhCyc2        = vector of phase cycling steps for 2nd editing pulse [degrees]
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]

clear more_params

% Include a list of all possible parameters
more_params.all_param_names = {'refocWaveform','disp_span', ...
    'refTp','Npts', 'sw', 'lw', 'Bfield', 'thkX', 'thkY', 'nx', 'ny', 'x', 'y', ...
    'TE', 'taus', 'centreFreq', 'editPhCyc1', 'editPhCyc2', 'refPhCyc1', 'refPhCyc2'...
    'editWaveform', 'editOnFreq', 'editOffFreq', 'editTp', 'editPhCyc1', 'editPhCyc2', 'pulse_type'};

% General Parameters
more_params.disp_span = []; % default: 1 ppm
more_params.refocWaveform = []; % default: sampleRefocPulse.pta
more_params.refTp = [];% default: 6.91
more_params.Npts = []; % default: 8192
more_params.sw = []; % default: 2000
more_params.lw = []; % default: 2.5
more_params.Bfield = []; % default: 3
more_params.thkX = []; % default 3
more_params.thkY= []; % default: 3
more_params.nx = []; % default: 19
more_params.ny = []; % default: 19
more_params.x=[]; % default:
more_params.y=[]; % default:
more_params.TE=[]; % default: 144
more_params.taus=[]; % default: 
more_params.centreFreq=[]; % default: 3.0 ppm
more_params.refPhCyc1=[]; % default: [0,90]
more_params.refPhCyc2=[]; % default: [0,90]

% MEGA-PRESS specific parameters
more_params.editWaveform = []; % default: sampleEditPulse.pta
more_params.editOnFreq = []; % default: 4.1 for lactate
more_params.editOffFreq = []; % default: 
more_params.editTp = []; % default: 14
more_params.editPhCyc1=[]; % default: [0 90]
more_params.editPhCyc2=[]; % default: [0 90 180 270]
more_params.pulse_type = []; % default: shaped

% Remove empty fields
all_field_names = fieldnames(more_params);
non_empty_fields = find(~structfun(@isempty,more_params));
cell_more_params = struct2cell(more_params);
more_params = cell2struct({cell_more_params{non_empty_fields}}, {all_field_names{non_empty_fields}}, 2);

clearvars -except more_params
save more_params.mat

% Reset the parameter generator to a blank
code = fileread('extra_params_blank.m');
fID = fopen('extra_params_gen.m', 'w');
fwrite(fID, code);

clear extra_params_gen.m