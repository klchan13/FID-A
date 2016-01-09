function [varargout] = spectra_sim(pulse_seq, spin_sys, obs_freq, varargin)
% Runs a simulation with the desired pulse sequence and spin system and
% parameters given by the structure given as the last input.  Can run
% either a single simulation at the given values or a series of simulations
% for multiple values for one parameter.
%
% 
% Inputs
% ------
% pulse_seq - function handle.  e.g. Using @sim_mega for pulse_seq will run
%             the MEGA-PRESS sequence.
% spin_sys - Spin system of interest.
% obs_freq - Observed frequency.  If not an empty array, the spectra will
%            be plotted in only a small range around the observed
%            frequency.
% varargin - Optional parameter that takes in a structure containing all
%            additional parameters.  If generated out of
%            extra_params_gen.m, this is a structure called "more_params'
% Outputs
% -------
% varargout - Output of a varying number.  If running a MEGA-PRESS
%             simulation, this produces 2 outputs - an on spectra and an off spectra.
%             If running a PRESS sequence, this produces only one output (spectra).

disp_span = 1; % Default display span around observed frequency if observed
               % obs_freq is given.

if nargin > 3
    more_params = varargin{1};
else
    more_params = [];
end

run check_param_fields

% Determine whether or not you want to do a parameter series
if isempty(find(param_len)) % No parameter series.
    
    % Run the simulation with the input pulse sequence.
    [spect] = pulse_seq(spin_sys, more_params);
    
    % Unpack the spectra out if the spect cell array, and prepare the
    % output.
    if length(spect) > 1
        ref_spect = spect{1};
        on_spect = spect{2};
        varargout{1} = ref_spect;
        varargout{2} = on_spect;
    else
        ref_spect = spect;
        spect = {spect};
        varargout{1} = ref_spect;
    end
    
    % Find the observed frequency and set the spectra range
    ppm_scale = ref_spect{1}{1}.ppm;
    if ~isempty(obs_freq)
        mid_ppm = find(round(ppm_scale*100)/100 == obs_freq);
        
        % Span of the metabolite's observed spectra [pts]
        disp_span_pts = floor((disp_span*ref_spect{1}{1}.Bo*42.577/ref_spect{1}{1}.spectralwidth)*ref_spect{1}{1}.n);
        pt_span = floor(disp_span_pts/2);
        obs_inds = (mid_ppm(1) - pt_span):(mid_ppm(1) + pt_span);
        varargout{1}{1}{1}.obs_inds = obs_inds;
    else
        obs_inds = 1:ref_spect.n;
    end

    figure 
    hold
    
    for s_idx=1:length(spect)
        if s_idx == 1
            color = 'k';
        elseif s_idx == 2
            color = 'r';
        end
        
        for n = 1:length(ref_spect{1}{1}.x)
            plot(ppm_scale(obs_inds),spect{s_idx}{n}{1}.specs(obs_inds)+5*n, color);
        end
    end

elseif length(find(param_len)) == 1  % Parameter series
    [spect_out] = parameter_series(pulse_seq, spin_sys, obs_freq, more_params);
    for out_idx = 1:length(spect_out)
        varargout{out_idx} = spect_out{out_idx};
    end
                                                 
elseif length(find(param_len)) > 1  % Sorry, only one parameter series can be run at a time.
    error('Only one parameter series can be run at a time.')
end
end

                                                              