%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%                                                                  %
% Function to get all the N2 wavenumber of different v and J       %
%                                                                  %
% Input:                                                           %
%      number_vib_states : Vibrational energy state                %
%      number_ro_states  : Rotational energy state                 %
%                                                                  %
% Output:                                                          %
%      wavenumber : Vibrational and rotational wavenumber          %
%                   (G(v)+F(j))                                    %
%      rowavenumber : Pure Rotational wavenumber (F(j))            %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [wavenumber,rowavenumber] = get_all_nu(number_vib_states,number_ro_states)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Molecule constants    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bv = be - alpha_e * (v+0.5)
alpha_e = 0.0173;     % cite(Jorgen Bendtsen, 1974) unit: cm-1
be = 1.99823;         % cite(D.A. Long, Raman Spectroscopy, page 156) unit: cm-1

% G(v) = ve * (v+0.5) - ve *  xe * (v+0.5)^2
ve = 2358.518;        % cite(CARSpy, Zhiyao Yin)  unit: cm-1
xe = 0.00606;         % cite(CARSpy, Zhiyao Yin), calculated by vexe/ve  unit: cm-1

% Dv = de + beta_e * (v+0.5)
beta_e = 1.55e-8;     % cite(CARSpy, Zhiyao Yin)  unit: cm-1
de = 5.774e-6;        % cite(CARSpy, Zhiyao Yin)  unit: cm-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocating for saving time
wavenumberbox = zeros(number_vib_states+1,number_ro_states+3);
rowavenumberbox = zeros(number_vib_states+1,number_ro_states+3);

% number_ro_state+2, because when J = number_ro_state, according the
% selection rules J has two (number_ro_state +- 2) energy state

% get vib+ro wavenumber
for vib = 0:number_vib_states
    for ro = 0:number_ro_states+2
        wavenumber_row_value = get_Raman_Shift(vib,ro);
        wavenumberbox(vib+1,ro+1)=wavenumber_row_value;
    end
end

% get pure rotational wavenumber
for vib = 0:number_vib_states
    for ro = 0:number_ro_states+2
        rowavenumber_row_value = get_F(vib,ro);
        rowavenumberbox(vib+1,ro+1) = rowavenumber_row_value;
    end   
end

wavenumber = wavenumberbox;
rowavenumber = rowavenumberbox;

%%%%%%%%%%%%%%%%%%%%%%
%%%   functions    %%%
%%%%%%%%%%%%%%%%%%%%%%
    
    % function to get G        
    function [vibRaman] = get_G(vib)
        % param vib: vibrational branches

        % return: wavenumber of 'vib' branch, unit: cm-1

        % G(vib) = (vib + 1/2) * ve - (vib + 1/2)^2 * ve * xe   
        % higher order are neglected

        vibRaman=(vib+0.5)*ve-(vib+0.5)^2.*ve*xe;
    end
    
    % function to get F       
    function [roRaman] = get_F(vib,ro)
        % param vib: vibrational branches
        % param ro: rotational banches

        % return: wavenumber of 'ro' branch in the 'vib' energy state, unit: cm-1

        % F(v,ro) = (be - alpha_e * (v+0.5)) * ro * (ro + 1) - 
        % (de + beta_e * (v+0.5))* ro^2 * (ro+1)^2
        % higher order are neglected

        roRaman=(be-alpha_e*(vib+0.5))*ro*(ro+1)-(de+beta_e*(vib+0.5))*ro^2*(ro+1)^2;
    end
    
    function [wavenumber] = get_Raman_Shift(vib,ro)
        % param vib: vibrational branches
        % param ro: rotational banches
        
        % return : vib+ro wavenumber
        wavenumber = get_G(vib) + get_F(vib,ro);
    end
    
    
end

%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%
