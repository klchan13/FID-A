function [spect_out] = parameter_series(pulse_seq, spin_sys, obs_freq, more_params)
% Runs a simulation for each value in the parameter series and displays the
% spectra at the observed frequency for each value in a string.
%
% Inputs
% ------
% pulse_seq - function handle.  e.g. Using @sim_mega for pulse_seq will run
%             the MEGA-PRESS sequence.
% spin_sys - Spin system of interest.
% obs_freq - Observed frequency.  If not an empty array, the spectra will
%            be plotted in only a small range around the observed
%            frequency.
% more_params - A structure containing all additional parameters.  If
%               generated out of extra_params_gen.m, this is a structure 
%               called "more_params'.
% Outputs
% -------
% spect_out - Cell array output with each entry containing one spectra.
%             If running a MEGA-PRESS simulation, this produces 2 outputs -
%             an ON spectra and an OFF spectra.  If running a PRESS
%             sequence, this produces only one spectra.

tic

Npts = 8192; % Using default number of points until changed by more_params
disp_span = 1; % Using default number of points until changed by more_params

run check_param_fields

% Find which of the fields contains the parameter series.
all_field_names = fieldnames(more_params);
cell_more_params = struct2cell(more_params);
idx_series = find(param_len);
series = cell_more_params{idx_series};
series_name = all_field_names{idx_series};


progress_bar = waitbar(0, sprintf('Initializing...'));
for idx = 1:length(series)

    t2 = toc;
    
    progress_bar = waitbar((idx-0.5)/length(series), progress_bar,...
                           sprintf('Simulating spectra %d of %d %s.\nElapsed time is %3.3f minutes.',...
                           idx, length(series), series_name, t2/60));
    
    % Make a new parameter structure with only value from the parameter
    % series.
    
    if isa(series(idx), 'cell')
        cell_more_params{idx_series} = series{idx};    
    else
        cell_more_params{idx_series} = series(idx);
    end
    
    more_params = cell2struct(cell_more_params, all_field_names, 1);               
    
    % Run the pulse sequence simulation
    [spect] = pulse_seq(spin_sys, more_params);
    
    % Initialize the output arrays
    if idx == 1
        ref_spect_arr = zeros(length(series), Npts);
        if length(spect) > 1
            on_spect_arr = ref_spect_arr;
        end
    end
    
    % Unpack the spectra out of the spect cell array.
    if length(spect) > 1
        on_spect = spect{2};
        ref_spect = spect{1};
    else
        ref_spect = spect;
        spect = {spect};
    end
                       
    ref_spect_arr(idx, :) = ref_spect{1}{1}.specs;
    if exist('on_spect')
        on_spect_arr(idx, :) = on_spect{1}{1}.specs;
    end
end

close(progress_bar)
% Now plot a string of spectra around the observed frequency
ppm_scale = ref_spect{1}{1}.ppm;
if ~isempty(obs_freq)
    mid_ppm = find(round(ppm_scale*100)/100 == obs_freq);
    
    % Span of the metabolite's observed spectra [pts]
    disp_span_pts = floor((disp_span*ref_spect{1}{1}.Bo*42.577/ref_spect{1}{1}.spectralwidth)*ref_spect{1}{1}.n);
    pt_span = floor(disp_span_pts/2);
    obs_inds = (mid_ppm(1) - pt_span):(mid_ppm(1) + pt_span);
else
    obs_inds = 1:ref_spect.n;
end

% Prepare the output as a cell array.
spect_out{1} = ref_spect_arr;
if exist('on_spect')
    spect_out{2} = on_spect_arr;
end

figure
hold on

for s_idx=1:length(spect)
    if s_idx == 1
        color = 'k'; % OFF
    elseif s_idx == 2
        color = 'r'; % ON
        diff_spect = spect_out{2} - spect_out{1};
        spectra_string(real(diff_spect), obs_inds, 'b')
    end

    this_arr = spect_out{s_idx};
    spectra_string(real(this_arr), obs_inds, color)

end

end