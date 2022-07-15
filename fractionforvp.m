%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for mole fraction of major species (N2,O2,C2H4) by measured    %
% temperature and ideal gas law                                       %
% Hu Meng & Dr. Yihua Ren, 2022-07-14                                 %
% This programm is only for vertilcal polarized signals               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
global wavenumber rowavenumber
global Wavelength_offset 

h = 6.62607004e-34;   % Planck constant                  unit: m^2*kg*s^-1
cs = 299792458;        % Speed of light                   unit: m/s
kb = 1.38064852e-23;  % Boltzmann constant               Unit: m^2*kg*s^-2*K^-1
ve_N2 = 2358.518; % cite(CARSpy, Zhiyao Yin)  unit: cm-1
ve_O2 = 1580.193; % cite (Buldakov et al. Journal of Molecular Spectroscopy, 2003)
ve_C2H4_V2 = 1656.4; % cite (Martin et al. The Journal of Chemical physics, 1995)
ve_C2H4_V3 = 1372.3; % cite (Martin et al. The Journal of Chemical physics, 1995)

path = 'C:\Users\mengh\Desktop\matlab\daten\2021-12-14\ref only vp\';
Files = dir(strcat(path,'*.csv'));
LengthFiles = length(Files);
plotpath = 'C:\Users\mengh\Desktop\article\new boundary\ref_0\';

% Preallocating
% 16 HAB(files), each HAB has 5 frames. totally 80 frames
Datafororigin = zeros(80,1024);
HAB = strings(16,1);
Simulationforoffset = zeros(80,3);
HABnamecol = zeros(80,1);
simulationdataforoffset = zeros(80,427);
simulationdataforraman = zeros(80,201);
[wavenumber,rowavenumber] = get_all_nu(20,100);

% Get HAB name
for i = 1:LengthFiles
    [Wavelength,Data]=csvreaderforfraction(Files(i).name,path);
    filename = strsplit(Files(i).name,'-');
    plotname = filename{1,1};
    plotname1 = convertCharsToStrings(plotname);
    HABname = strsplit(plotname1,'mm');
    HABname0 = HABname{1,1};
    HABname1 = strsplit(HABname0,'HAB=');
    HABname2 = HABname1{1,2};
    
    
    for j = 1:length(Data(:,1))
        Datafororigin(5*(i-1)+j,:)=Data(j,:);
        HABnamecol(5*(i-1)+j,:)=str2double(HABname2);
    end

        HAB(i,:)=HABname2;
end

% % get averaged spectrum of each HAB for plot
% for coll = 1:1024
% 
%     spec = reshape(Datafororigin(:,coll),[5,16]);
%     Meanspec = mean (spec);
%     meanspectra(:,coll)=Meanspec';
% end

HAB_voll = [NaN;HABnamecol];
Datasammelung = [Wavelength;Datafororigin];
Datafinalsammelung = [HAB_voll,Datasammelung];

% Baseline selection
xdaten1 = Datafinalsammelung(:,2:76); % - 1300cm-1
% xdaten2 = Datafinalsammelung(:,137:259); % 1372 -1537cm-1
xdaten3 = Datafinalsammelung(:,404:534); % 1680 -1900cm-1
xdaten4 = Datafinalsammelung(:,620:678); % 2000 -2100cm-1
xdaten5 = Datafinalsammelung(:,751:816); % 2160 -2250cm-1
xdaten6 = Datafinalsammelung(:,930:1025); % 2360 - 

baselinedaten = [xdaten1,xdaten3,xdaten4,xdaten5,xdaten6];

Wavelength_offset = baselinedaten(1,:);

% baseline fitting with 2-degree polynomial
for n=2:length(baselinedaten(:,1))
    Iexp = baselinedaten(n,:);
    errorlsqfun = @(x)twopoly(x(1),x(2),x(3))-Iexp;
    initial0 = [0,0,0];
    ub = [1,5,30000];
    lb = [-1,-10,-1000];

    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxFunctionEvaluations',10^6);
    [Xfit] = lsqnonlin(errorlsqfun,initial0,lb,ub,options);
    nullarray = zeros(size(Wavelength_offset));

    % save fitting results
    Simulationforoffset(n-1,:)=Xfit;

    simulationdataforoffset(n-1,:) = twopoly(Xfit(1),Xfit(2),Xfit(3));
end

% preallocating
Dataforfraction = zeros(80,1024);
baseline = zeros(80,1024);

% baseline correction
for b = 1:length(Datafororigin(:,1))
    for c = 1:length(Wavelength(1,:))
        baseline(b,c) = Simulationforoffset(b,1) * Wavelength(c) *Wavelength(c)+ Simulationforoffset(b,2)*Wavelength(c)+Simulationforoffset(b,3);
       Dataforfraction(b,c) = Datafororigin(b,c) - baseline(b,c);  
    end
end

% get averaged baseline corrected spectrum of each HAB for plot
% for col = 1:1024
% 
%     Aftersub = reshape(Dataforfraction(:,col),[5,16]);
%     Meanafter = mean (Aftersub);
%     Aftersubstracting(:,col)=Meanafter';
% end

C2H4V2 = zeros(80,2);
C2H4V3 = zeros(80,2);
O2 = zeros(80,2);
N2 = zeros(80,2);
C2H4V3(:,1) = HABnamecol;
C2H4V2(:,1) = HABnamecol;
O2(:,1) = HABnamecol;
N2(:,1) = HABnamecol;
% preallocating
area=zeros(80,201);
newarea=zeros(80,201);
%load temperature, intensity, fitted result of Raman spectrum
load('C:\Users\mengh\Desktop\article\new boundary\0\0_expt_DMC——nostretch.mat','finalresult','meanvalue','Simulationforraman','Dataforramanoffset');
for i = 1:1:80
    
    % get fitted Raman spectrum
    simulationdataforraman(i,:) = Intensity_sim_wonormalization(Simulationforraman(i,1),Simulationforraman(i,2),Simulationforraman(i,3),Simulationforraman(i,4),Simulationforraman(i,5),Simulationforraman(i,6));

end
for q = 1:length(Dataforfraction(:,1))
    % intergration
    C2H4V3(q,2) = trapz(Dataforfraction(q,93:133));
    % intergration
    C2H4V2(q,2) = trapz(Dataforfraction(q,300:330));
    % intergration
    O2(q,2) = trapz(Dataforfraction(q,251:287)); 
    % reconstruction by fitted Raman spectrum
    area(q,:)=max(Dataforramanoffset(q,:)) *simulationdataforraman(q,:);
    % baseline correction of fitted Raman spectrum
    % newarea(q,:)=area(q,:)-mean(area(q,160:190));
    % intergration
    N2(q,2) =trapz(area(q,:));
end


% N2 FOR VP CALIBRATION HAB = 15 mm
N2_I0_vp = mean(N2(41:45,2));
N2_n0_vp = 0.79;
N2_T0_vp = mean(finalresult(41:45,2));
% O2 CALIBRATION HAB = 15 mm
O2_I0 = mean(O2(41:45,2));
O2_n0 = 0.21;
O2_T0 = mean(finalresult(41:45,2));
% C2H4V2 CALIBRATION HAB= 7 mm
C2H4V2_I0 = mean(C2H4V2(46:50,2));
C2H4V2_n0 = 0.4997;
C2H4V2_T0 = mean(finalresult(46:50,2));
% C2H4V3 CALIBRATION HAB= 7 mm
C2H4V3_I0 = mean(C2H4V3(46:50,2));
C2H4V3_n0 = 0.4997;
C2H4V3_T0 = mean(finalresult(46:50,2));



% START CALCULATE n_hab FOR C2H4V2
% preallocating
temp = zeros(80,2);
temp(:,1)=C2H4V2(:,1);
n_c2h4v2 = zeros(80,2);
n_c2h4v3 = zeros(80,2);
n_n2 = zeros(80,2);
n_o2 = zeros(80,2);
% save HAB to first col
n_c2h4v2(:,1) = C2H4V2(:,1);
n_c2h4v3(:,1) = C2H4V2(:,1);
n_n2(:,1) = C2H4V2(:,1);
n_o2(:,1) = C2H4V2(:,1);

% calculation
for hab_c2h4 = 1:length(C2H4V2(:,1))
    temp(hab_c2h4,2) = finalresult(hab_c2h4,2);
    n_c2h4v2(hab_c2h4,2) = (C2H4V2(hab_c2h4,2) * C2H4V2_n0 / C2H4V2_I0)* ((finalresult(hab_c2h4,2))/C2H4V2_T0)*((1 - exp(-1*h*cs*ve_C2H4_V2*100/kb/finalresult(hab_c2h4,2)))/(1 - exp(-1*h*cs*ve_C2H4_V2*100/kb/C2H4V2_T0)));
    n_c2h4v3(hab_c2h4,2) = (C2H4V3(hab_c2h4,2) * C2H4V3_n0 / C2H4V3_I0)* ((finalresult(hab_c2h4,2))/C2H4V3_T0)*((1 - exp(-1*h*cs*ve_C2H4_V3*100/kb/finalresult(hab_c2h4,2)))/(1 - exp(-1*h*cs*ve_C2H4_V3*100/kb/C2H4V3_T0)));
    n_n2(hab_c2h4,2) = (N2(hab_c2h4,2) * N2_n0_vp / N2_I0_vp)* ((finalresult(hab_c2h4,2))/N2_T0_vp)*((1 - exp(-1*h*cs*ve_N2*100/kb/finalresult(hab_c2h4,2)))/(1 - exp(-1*h*cs*ve_N2*100/kb/N2_T0_vp)));
    n_o2(hab_c2h4,2) = (O2(hab_c2h4,2) * O2_n0 / O2_I0)* (finalresult(hab_c2h4,2)/O2_T0)*((1 - exp(-1*h*cs*ve_O2*100/kb/finalresult(hab_c2h4,2)))/(1 - exp(-1*h*cs*ve_O2*100/kb/O2_T0)));
end

% preallocating for mean and standard value
mean_n_c2h4v2 = zeros(16,1);
std_n_c2h4v2 = zeros(16,1);

mean_n_c2h4v3 = zeros(16,1);
std_n_c2h4v3 = zeros(16,1);

mean_n_n2 = zeros(16,1);
std_n_n2 = zeros(16,1);

mean_n_o2 = zeros(16,1);
std_n_o2 = zeros(16,1);

% calculate mean and std
for index = 1:(length(n_c2h4v2(:,1))/5)
    matrixformean_c2h4v2 = n_c2h4v2(5*(index-1)+1:5*(index-1)+5,2);
    mean_n_c2h4v2(index,1)=mean(matrixformean_c2h4v2);
    std_n_c2h4v2(index,1)=std(matrixformean_c2h4v2);
    
    matrixformean_c2h4v3 = n_c2h4v3(5*(index-1)+1:5*(index-1)+5,2);
    mean_n_c2h4v3(index,1)=mean(matrixformean_c2h4v3);
    std_n_c2h4v3(index,1)=std(matrixformean_c2h4v3);
    
    
    matrixformean_n2 = n_n2(5*(index-1)+1:5*(index-1)+5,2);
    mean_n_n2(index,1)=mean(matrixformean_n2);
    std_n_n2(index,1)=std(matrixformean_n2);
    
    matrixformean_o2 = n_o2(5*(index-1)+1:5*(index-1)+5,2);
    mean_n_o2(index,1)=mean(matrixformean_o2);
    std_n_o2(index,1)=std(matrixformean_o2);
    
end

datenforexcel = [HAB,mean_n_c2h4v2,std_n_c2h4v2,mean_n_c2h4v3,std_n_c2h4v3,mean_n_n2,std_n_n2,mean_n_o2,std_n_o2];

datafilename = strcat(plotpath,'ref_fraction_0.79_final_test0.xlsx');
writematrix(n_c2h4v2,datafilename,'Sheet',1,'Range','A1');
writematrix(n_c2h4v3,datafilename,'Sheet',1,'Range','C1');
writematrix(n_n2,datafilename,'Sheet',1,'Range','E1');
writematrix(n_o2,datafilename,'Sheet',1,'Range','I1');
writematrix(datenforexcel,datafilename,'Sheet',2,'Range','A1');


matfilename = strcat(plotpath,'ref_fraction_0.79_final_test0.mat');
save(matfilename);
    


%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%


