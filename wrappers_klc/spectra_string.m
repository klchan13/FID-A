function [spect_string, spect_integral] = spectra_string(phased_spectra, disp_range, color)
% Take spectra from multiple parameter values and string the parts we're
% interested in together.

spect_string = zeros(size(phased_spectra,1)*length(disp_range), 1);
spect_integral = zeros(size(phased_spectra,1), 1);

for te_idx = 1:size(phased_spectra,1)
    spectra_range = (1 + (te_idx-1)*length(disp_range)):(te_idx*length(disp_range));
    spect_string(spectra_range) = phased_spectra(te_idx, disp_range);
    spect_integral(te_idx) = sum(phased_spectra(te_idx, disp_range));
end

plot(spect_string', color)

end