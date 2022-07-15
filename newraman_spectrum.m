%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
% Function to calculate discrete N2 Raman signal of different Raman shift %                                  
% with different temperature                                              %
%                                                                         %
%                                                                         %
% only for Raman spontaneous stokes scattering                            %
% selection rules for diatomic molecules: di = 1, -1. dj = 0, 2, -2       %
% here we only consider the stokes line, which means di = 1.              %
% Considering the influence of nuclear spin, Pauli exclusion principle,   %
% which cause a intensity alternation, even j : odd j = 1 : 2             %
% All the equations are from Raman Spectroscopy (D.A. Long) and Laser     %
% diagnostic for combustion temperature and speciec (Alan C. Eckbreth)    %
%                                                                         %
% Input:                                                                  %
%      temperature : flame (N2) temperature                               %
% Output:                                                                 %
%      nu_spectrum: wavelength of calculated discreted N2 Raman spectrum  %
%      intensity_spectrum: intensity of calculated discreted N2           %
%                          Raman spectrum                                 %
%      nu_intensity: 2D matrix of nu_spectrum and intensity_spectrum      %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [nu_spectrum,intensity_spectrum,nu_intensity] = newraman_spectrum(temperature)
global wavenumber rowavenumber


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% molecule constants %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

laserWaveLength = 355; % unit nm

% nm to cm-1
laserWaveNumber = 1e7 / laserWaveLength;

% molecular constants
h = 6.62607004e-34;   % Planck constant                  unit: m^2*kg*s^-1
c = 299792458;        % Speed of light                   unit: m/s
e0 = 8.85418782e-12;  % Permittivity of free space       unit: m^-3*kg^-1*s^4*A^2
kb = 1.38064852e-23;  % Boltzmann constant               Unit: m^2*kg*s^-2*K^-1

% Bv = be - alpha_e * (v+0.5)
alpha_e = 0.0173;  % cite(Jorgen Bendtsen, 1974) unit: cm-1
be = 1.99823;  % cite(D.A. Long, Raman Spectroscopy, page 156) unit: cm-1
ve = 2358.518;   % cite(CARSpy, Zhiyao Yin)  unit: cm-1
xe = 0.00606;   % cite(CARSpy, Zhiyao Yin), calculated by vexe/ve  unit: cm-1

% derived polarizability tensor α', λ'.
% cite(Vibration roration Raman spectroscopy of gas media, M.A. Buldakov, 1995)
alpha_tensor = 1.86; % derived mean polarizability  unit: 10^-30 m^3
gamma_tensor = 2.23; % derived anisotropy unit: 10^-30 m^3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocating, O Branch
wavenumber_O = zeros(1,2020);
intensity_O = zeros(1,2020);

% 20 vibrational energy states
% 100 rotational energy states
for vib = 1:20
    for ro = 2:102
        % Raman shift
        wavenumber_O_value = get_O_Branch_Shift(vib,ro);
        wavenumber_O(1,101*(vib-1)+ro-1) = wavenumber_O_value;
        
        % Raman signal intensity
        intensity_O_value=get_Intensity_wo_tensor(vib,wavenumber_O_value)*get_tensor_O(ro)*get_population_factor_new(vib,ro);
        intensity_O(1,101*(vib-1)+ro-1) = intensity_O_value;
    end
end

% preallocating, S Branch
wavenumber_S = zeros(1,2020);
intensity_S = zeros(1,2020);
for vib = 1:20
    for ro = 0:100
        wavenumber_S_value = get_S_Branch_Shift(vib,ro);
        wavenumber_S(1,101*(vib-1)+ro+1) = wavenumber_S_value;
        intensity_S_value=get_Intensity_wo_tensor(vib,wavenumber_S_value)*get_tensor_S(ro)*get_population_factor_new(vib,ro);
        intensity_S(1,101*(vib-1)+ro+1) = intensity_S_value;
    end
end

% preallocating, Q Branch
wavenumber_Q = zeros(1,2020);
intensity_Q = zeros(1,2020);
for vib = 1:20
    for ro = 0:100
        wavenumber_Q_value = get_Q_Branch_Shift(vib,ro);
        wavenumber_Q(1,101*(vib-1)+ro+1) = wavenumber_Q_value;
        intensity_Q_value=get_Intensity_wo_tensor(vib,wavenumber_Q_value)*get_tensor_Q(ro)*get_population_factor_new(vib,ro);
        intensity_Q(1,101*(vib-1)+ro+1) = intensity_Q_value;
    end
end


% merge Q,S,O branches
nu_spectrum=[wavenumber_Q,wavenumber_O,wavenumber_S];
intensity_spectrum=[intensity_Q,intensity_O,intensity_S];

% row transfer to column
nu_trans = nu_spectrum.';
intensity_trans = intensity_spectrum.';

% sorting
nu_intensity=[nu_trans,intensity_trans];
[sort_nu_intensity,~]=sortrows(nu_intensity,1);

% column transfer ro row
nu_intensity = sort_nu_intensity.';

% pick up the data we need
nu_spectrum = nu_intensity(1,:);
% pick up 2200 cm-1 < raman shift < 2400 cm -1
pos = find(nu_spectrum > 2200 & nu_spectrum < 2400);
start = pos(1);
finish = pos(end);
nu_spectrum = nu_intensity(1,start:finish);
intensity_spectrum=nu_intensity(2,start:finish);

%normalisiren
%intensity_spectrum=intensity_spectrum./max(intensity_spectrum);


%%%%%%%%%%%%%%%%%%%%%%
%%%   functions    %%%
%%%%%%%%%%%%%%%%%%%%%%

    % get wavnumber of rotational energy levels from global virables, saving time
    function [roNv] = get_F(vib,ro)
        % param vib: vibrational branches
        % param ro: rotational banches
        
        % return: wavenumber of 'ro' branch in the 'vib' energy state, unit: cm-1
        %         only rorational enerfgy state wavenumber for calculation
        %         of population factor
        
        roNv=rowavenumber(vib+1,ro+1);
    end
    

    % get wavnumber of vib+ro energy levels from global virables
    function [Nv] = get_G(vib,ro)
        % param vib: vibrational branches
        % param ro: rotational banches
        
        % return: wavenumber of 'ro' branch in the 'vib' energy state, unit: cm-1
        
        Nv = wavenumber(vib+1,ro+1);
    end
    

    % get q branch raman shift dv = 1, dj = 0
    function [qBranchShift] = get_Q_Branch_Shift(vib,ro)
        % param vib: vibrational branches
        % param ro: rotational banches
        
        % return: Raman shift of Q branch, unit: cm-1
            
        
        qBranchShift = get_G(vib,ro)-get_G(vib-1,ro);
    end
    
    % get o branch raman shift dv = 1, dj = -2
    function [oBranchShift] = get_O_Branch_Shift(vib,ro)
        % param vib: vibrational branches
        % param ro: rotational banches
        
        % return: Raman shift of O branch, unit: cm-1
        
        
        oBranchShift = get_G(vib,ro-2)-get_G(vib-1,ro);
    end

    % get s branch raman shift dv = 1, dj = 2
    function [sBranchShift] = get_S_Branch_Shift(vib,ro)
        
        % param vib: vibrational branches
        % param ro: rotational banches
        
        % return: Raman shift of S branch, unit: cm-1
        
        
        sBranchShift = get_G(vib,ro+2)-get_G(vib-1,ro);
    end
    

    % calculate Raman intensity for different raman shifts
    function [intensity_wo_tensor] = get_Intensity_wo_tensor(vib,ramanShift)
        
        % param vib: Vibrational branches
        % param ramanShift: calculated raman shift
        
        % return: Raman signal intensity of every signal energy level
        
        % detailed eqation see <Raman spectroscopy> (D.A. Long) and
        % <Laser diagnostics for combustion temperature and speciees> (Alan C. Eckbreth)
        intensity_scattering = (h*(laserWaveNumber*100-ramanShift*100)^4 /8/c/e0/e0/ramanShift/100);
        
        % boltzmann distribution law
        % Here is the simplified eqation
        % instneisy_number_density_vib1 = 1 - exp(-1*h*c*ve*100/kb/temperature); 
        % intensity_number_density_vib =(vib)*exp(-1*(vib-1)*h*c*ve*100/kb/temperature);
      
        vibsum = 0;
        for isum = 1:1:20

            vibsum = vibsum+exp(-1*((isum-0.5)*ve-(isum-0.5)^2*ve*xe)*100*h*c/kb/temperature);
        end

        
        population_factor_vib=vib*exp(-1*((vib-0.5)*ve-(vib-0.5)^2*ve*xe)*100*h*c/kb/temperature)/vibsum;

        
        intensity_wo_tensor = intensity_scattering*population_factor_vib;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Space-aceraged squares of transition polarizability tensor components
    % for diatomoic molecules for totation and vibration-rotation
    % transitions. cite:<Raman spectroscopy> (D.A. Long) CRS 13,14.
    
    % param ro: rotational branches

    % return: polarizability tensor
    function [tensor_O] = get_tensor_O(ro)
        b_jminus2_j=3*ro*(ro-1)/2/(2*ro+1)/(2*ro-1);
        tensor_O=(4/45)*(gamma_tensor)^2*b_jminus2_j;
    end

    function [tensor_S] = get_tensor_S(ro)
        b_jplus2_j =3*(ro+1)*(ro+2)/2/(2*ro+1)/(2*ro+3);
        tensor_S=(4/45)*gamma_tensor^2*b_jplus2_j;
    end

    function [tensor_Q] = get_tensor_Q(ro)
        b_j_j=ro*(ro+1)/(2*ro-1)/(2*ro+3);
        tensor_Q=(4/45)*gamma_tensor^2*b_j_j+alpha_tensor^2;
    end

    % population factor caused by Pauli exclusion principle times 

    % Here is the simplified equation
    % Rotational partition function
    %     function [population_factor] = get_population_factor(vib,ro)
    %         % param vib: Vibrational branches
    %         % param ro: rotational branches
    %         
    %         % return: population factor x Rotational partition function
    %         
    %         qr = kb*temperature/h/c/(be-alpha_e*(vib+0.5))/100;
    %         if mod(ro,2)==0
    %             gj = 6;
    %         else 
    %             gj = 3;
    %         end
    %         population_factor=gj*(2*ro+1)*exp(-1*get_F(vib-1,ro)*100*h*c/kb/temperature)/qr;
    %     end  
    function [population_factor] = get_population_factor_new(vib,ro)
        % param vib: Vibrational branches
        % param ro: rotational branches
        
        % return: population factor x Rotational partition function
        sum = 0;
        for i = 1:1:100
            if mod(i,2)==0
                gj = 6;
            else 
                gj = 3;
            end
            sum = sum+gj*(2*i+1)*exp(-1*get_F(vib-1,i)*100*h*c/kb/temperature);
        end
       if mod(ro,2)==0
                gj = 6;
            else 
                gj = 3;
       end
        
        population_factor=gj*(2*ro+1)*exp(-1*get_F(vib-1,ro)*100*h*c/kb/temperature)/sum;
    end  
end


%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%

