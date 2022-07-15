
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        %
% linear function for fitting baseline   %
% y = ax + b                             %
%                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = Linearfun(a,b)
global Wavelength_offset


wavenumber = Wavelength_offset;
Yi = zeros(size(wavenumber));
for i = 1:length(Yi)
    Yi(i)= a * wavenumber(i) + b;
end

Y = Yi;
end

%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%

