% Takes the parameters given in the struct more_params, extracts the values
% from each field and assigns them to different variables with the same
% names as their field names.  Also produces a binary array indicating
% whether or not there are more than one value for each parameter.

if ~isempty(more_params)
    disp('Extracting additional parameters.  Default values will be used for the rest.')
    add_params = fieldnames(more_params);
    param_len = zeros(1,length(add_params));
    for p_idx = 1:length(add_params)
        switch add_params{p_idx}
            case 'refocWaveform'
                refocWaveform = getfield(more_params, 'refocWaveform');
                param_len(p_idx) = (isa(refocWaveform, 'cell') && length(refocWaveform) > 1);
            case 'refocWaveform1'
                refocWaveform1 = getfield(more_params, 'refocWaveform1');
                param_len(p_idx) = (isa(refocWaveform1, 'cell') && length(refocWaveform1) > 1);
            case 'refocWaveform2'
                refocWaveform2 = getfield(more_params, 'refocWaveform2');
                param_len(p_idx) = (isa(refocWaveform2, 'cell') && length(refocWaveform2) > 1);
            case 'editWaveform'
                editWaveform = getfield(more_params, 'editWaveform');
                param_len(p_idx) = (isa(editWaveform, 'cell') && length(editWaveform) > 1); 
            case 'excWaveform'
                excWaveform = getfield(more_params, 'excWaveform');
                param_len(p_idx) = (isa(excWaveform, 'cell') && length(excWaveform) > 1);
            case 'editOnFreq'
                editOnFreq = getfield(more_params, 'editOnFreq');
                param_len(p_idx) = length(editOnFreq) > 1;
            case 'editOffFreq'
                editOffFreq = getfield(more_params, 'editOffFreq');
                param_len(p_idx) = length(editOffFreq) > 1;
            case 'editTp'
                editTp = getfield(more_params, 'editTp');
                param_len(p_idx) = length(editTp) > 1;
            case 'excTp'
                excTp = getfield(more_params, 'excTp');
                param_len(p_idx) = length(excTp) > 1;
            case 'refTp'
                refTp = getfield(more_params, 'refTp');
                param_len(p_idx) = length(refTp) > 1;
            case 'refTp1'
                refTp1 = getfield(more_params, 'refTp1');
                param_len(p_idx) = length(refTp1) > 1;
            case 'refTp2'
                refTp2 = getfield(more_params, 'refTp2');
                param_len(p_idx) = length(refTp2) > 1;
            case 'Npts'
                Npts = getfield(more_params, 'Npts');
                param_len(p_idx) = length(Npts) > 1;
            case 'sw'
                sw = getfield(more_params, 'sw');
                param_len(p_idx) = length(sw) > 1;
            case 'lw'
                lw = getfield(more_params, 'lw');
                param_len(p_idx) = length(lw) > 1;
            case 'Bfield'
                Bfield = getfield(more_params, 'Bfield');
                param_len(p_idx) = length(Bfield) > 1;
            case 'thkX'
                thkX = getfield(more_params, 'thkX');
                param_len(p_idx) = length(thkX) > 1;
            case 'thkY'
                thkY = getfield(more_params, 'thkY');
                param_len(p_idx) = length(thkY) > 1;
            case 'nx'
                nx = getfield(more_params, 'nx');
                param_len(p_idx) = length(nx) > 1;
            case 'ny'
                ny = getfield(more_params, 'ny');
                param_len(p_idx) = length(ny) > 1;
            case 'x'
                x = getfield(more_params, 'x');
                param_len(p_idx) = (isa(x, 'cell') && length(x) > 1);
            case 'y'
                y = getfield(more_params, 'y');
                param_len(p_idx) = (isa(y, 'cell') && length(y) > 1);
            case 'Gz'
                Gz = getfield(more_params, 'Gz');
                param_len(p_idx) = length(Gz) > 1;
            case 'z'
                z = getfield(more_params, 'z');
                param_len(p_idx) = (isa(z, 'cell') && length(z) > 1);
            case 'TE1'
                TE1 = getfield(more_params, 'TE1');
                param_len(p_idx) = length(TE1) > 1;
            case 'TE'
                TE = getfield(more_params, 'TE');
                param_len(p_idx) = length(TE) > 1;
            case 'taus'
                taus = getfield(more_params, 'taus');
                param_len(p_idx) = (isa(taus, 'cell') && length(taus) > 1);
            case 'centreFreq'
                centreFreq = getfield(more_params, 'centreFreq');
                param_len(p_idx) = length(centreFreq) > 1;
            case 'editPhCyc1'
                editPhCyc1 = getfield(more_params, 'editPhCyc1');
                param_len(p_idx) = (isa(editPhCyc1, 'cell') && length(editPhCyc1) > 1);
            case 'editPhCyc2'
                editPhCyc2 = getfield(more_params, 'editPhCyc2');
                param_len(p_idx) = (isa(editPhCyc2, 'cell') && length(editPhCyc2) > 1);
            case 'refPhCyc1'
                refPhCyc1 = getfield(more_params, 'refPhCyc1');
                param_len(p_idx) = (isa(refPhCyc1, 'cell') && length(refPhCyc1) > 1);
            case 'refPhCyc2'
                refPhCyc2 = getfield(more_params, 'refPhCyc2');
                param_len(p_idx) = (isa(refPhCyc1, 'cell') && length(refPhCyc1) > 1);
            case 'disp_span'
                disp_span = getfield(more_params, 'disp_span');
                param_len(p_idx) = length(disp_span) > 1;
            case 'pulse_type'
                pulse_type = getfield(more_params, 'pulse_type');
                param_len(p_idx) = 0;
            case 'edit_angle'
                edit_angle = getfield(more_params, 'edit_angle');
                param_len(p_idx) = length(edit_angle) > 1;
            case 'k'
                edit_angle = getfield(more_params, 'k');
                param_len(p_idx) = length(k) > 1;
        end

    end
else
    disp('No additional parameter inputs.  Using default values.')
    param_len = 0;
end
    