%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
% Function to calculate continuous N2 Raman signal of different           %
% Raman shift by convolution with voigt function                          %
%                                                                         %
% Input:                                                                  %
%      temperature : flame (N2) temperature                               %
%      offset : same as baseline(ax+b,but a = 0)                          %
%      shift : make the specturm shift to left or right                   %
%      wG : width of the Gaussian component (Half-width at half maximum)  %
%      wL : width of the Lorentzian component (Half-width at half maximum)%
%      stretch : expand or shrink of spectrum                             %
% Output:                                                                 %
%      Calculated continuous N2 Raman spectrum                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = Ramanspectra(temperature,offset,shift,wG,wL,stretch)

Wavelength_raman= 2200:1:2400;
% get discreted simulated N2 Raman spectrum
[nu_spectrum,intensity_spectrum] = newraman_spectrum(temperature);

% get wavelength after stretch and shift
Si =zeros(size(Wavelength_raman));
newWavelength_raman = Wavelength_raman + (2330-Wavelength_raman)*stretch  + shift;

% for each wavelength, convolution with vogit function, then add them
% together.
for i=1:1:length(nu_spectrum)
    Si = Si + voigt(newWavelength_raman,intensity_spectrum(i),nu_spectrum(i),wG,wL,offset);
end

% normalization
S = Si./max(Si);

end 


%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%