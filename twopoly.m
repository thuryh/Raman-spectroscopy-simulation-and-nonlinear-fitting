
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% 2  Degree polunomial for fitting baseline  %
% y = ax + b                                 %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = twopoly(a,b,c)
% get the Wavelength of the baseline
global Wavelength_offset


wavenumber = Wavelength_offset;
Yi = zeros(size(wavenumber));
for i = 1:length(Yi)
    Yi(i)= a * wavenumber(i)* wavenumber(i) + b* wavenumber(i)+c;
end

Y = Yi;
end

%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%