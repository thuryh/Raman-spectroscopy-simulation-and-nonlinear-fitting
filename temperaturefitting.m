%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for postprecess of Raman spectra to get temperature by using   %
% lsq-fitting method.                                                 %
% Hu Meng & Dr. Yihua Ren, 2022-01-25                                 %
% This programm is only for vertilcal polarized signals               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc; close all;

global Wavelength_offset Wavelength_raman
global wavenumber rowavenumber

%%%%%%%%%%%%%%%%%%%%%%
%      read csv      %
%%%%%%%%%%%%%%%%%%%%%%

% get ready for reading data from csv files
path = 'C:\Users\mengh\Desktop\matlab\daten\2021-12-14\ref only vp\';
Files = dir(strcat(path,'*.csv'));
LengthFiles = length(Files);
plotpath = 'C:\Users\mengh\Desktop\article\no stretch\0_stretch0.001\';

% preallocating for saving time

% preallocating for the data for plots, raw data, 16 files (HABS), 
% each HAB 5 Frames, each frame 211 data.
frame = 5;
matrix_x = frame * LengthFiles;
Datafororigin = zeros(matrix_x,211);

% preallocating for HABs, 16 files
HAB = strings(LengthFiles,1);

% preallocating for baseline data
Simulationforoffset = zeros(matrix_x,2);

% preallocating for fitting result of each HABs
Simulationforraman = zeros(LengthFiles,6);

% preallocating for HABs, 80 frames
HABnamecol = zeros(matrix_x,1);

% preallocating for the data for baseline fitting
simulationdataforoffset = zeros(matrix_x,130);

% preallocating for the data for N2 spectra fitting
simulationdataforraman = zeros(matrix_x,53);
% preallocating for the 95% confidential interval for N2 spectra fitting
simulationdataforramanerror = zeros(matrix_x,53);
ciArray = zeros(matrix_x*2,7);
errArray = zeros(matrix_x,7);
relerrArray = zeros(matrix_x,7);


% get the calculated Raman shift, consist of Vibrational-rotational Raman 
% shift and pure rotational Raman shift(pure rotational Raman shift for 
% calculating population factor)
[wavenumber,rowavenumber] = get_all_nu(20,100);

%
% read csv files
%

for i = 1:LengthFiles
    % Wavelength 2230 - 2490 cm-1 
    [Wavelength,Data]=csvreader(Files(i).name,path);
    
    % get HAB name, char2string
    filename = strsplit(Files(i).name,'-');
    plotname = filename{1,1};
    plotname1 = convertCharsToStrings(plotname);
    HABname = strsplit(plotname1,'mm');
    HABname0 = HABname{1,1};
    HABname1 = strsplit(HABname0,'HAB=');
    HABname2 = HABname1{1,2};
    
    
    for j = 1:length(Data(:,1))
        %collecting data from each csv file to Matrix
        %each file has 5 frames
        Datafororigin(5*(i-1)+j,:)=Data(j,:);
        %HAB frames name string2double 
        HABnamecol(5*(i-1)+j,:)=str2double(HABname2);
    end
        %HAB files name
        HAB(i,:)=HABname2;
end

% put the wavelength,HAB and raw data together
HAB_voll = [NaN;HABnamecol];
Datasammelung = [Wavelength;Datafororigin];
Datafinalsammelung = [HAB_voll,Datasammelung];

% pick up the data that used to calculate the baseline
% 2250 cm-1 > Wavelength or Wavelength > 2350 cm-1 was chosen
% get the position of those chosen data
xdaten_pos_offset=find(Datasammelung(1,:)<2250|2350 < Datasammelung(1,:));

% preallocating for save wavelength and the data for baseline
Dataforoffset = zeros(matrix_x+1,length(xdaten_pos_offset));

% get the wavelength and signal intensity of the picked baseline by its
% position
for a = 1:length(xdaten_pos_offset)
    pos = xdaten_pos_offset(1,a);
    Dataforoffset(:,a)=Datasammelung(:,pos);
end
Wavelength_offset = Dataforoffset(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    baseline fitting      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=2:length(Dataforoffset(:,1))
    Iexp = Dataforoffset(n,:);
    % fiiting with linear function
    errorlsqfun = @(x)Linearfun(x(1),x(2))-Iexp;
    
    % initial values
    initial0 = [0,0];
    
    % upper boundary
    ub = [0,20000];
    
    % lower boundary
    lb = [-5,0];

    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxFunctionEvaluations',10^6);
    
    % fitting result
    [Xfit] = lsqnonlin(errorlsqfun,initial0,lb,ub,options);
    
    nullarray = zeros(size(Wavelength_offset));
    
%     % check fitting result by plots
%     figure(1)
%     subplot(2,1,1);
%     plot(Wavelength_offset,Iexp,Wavelength_offset,Linearfun(Xfit(1),Xfit(2)));
%     set(gcf,'visible','on');
%     xlim([2200 2550]);
%     ylim([0 20000]);
%     subplot(2,1,2);
%     plot(Wavelength_offset,Linearfun(Xfit(1),Xfit(2))-Iexp,Wavelength_offset,nullarray);hold on;
%     plotnameforsave = num2str(HABnamecol(n-1));
%     set(gcf,'visible','on');
%     xlim([2200 2550]);
%     ylim([-1 1]);
%     %saveas(gcf,strcat(plotpath,strcat(plotnameforsave,'_Frame_',num2str(mod(n-1,5)),'_',num2str(Xfit(1)),'.fig')));
%     %saveas(gcf,strcat(plotpath,strcat(plotnameforsave,'_Frame_',num2str(mod(n-1,5)),'_',num2str(Xfit(1)),'.jpg')));
%     delete(gcf);

    
    % save the fitting results
    Simulationforoffset(n-1,:)=Xfit;

    
    % save the fitted baseline
    simulationdataforoffset(n-1,:) = Linearfun(Xfit(1),Xfit(2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Raman fitting        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pick up the data that used to fitting calculated N2 Raman spectra
% 2280 cm-1 < Wavelength or Wavelength < 2345 cm-1 was chosen
% get the position of those chosen data
xdaten_pos_raman = find(2345>Datasammelung(1,:) & 2280<Datasammelung(1,:));
start = xdaten_pos_raman(1);
finish = xdaten_pos_raman(end);

% preallocating for save wavelength and the data for Raman fiiting
Dataforramanoffset = zeros(80,finish-start+1);
Dataforramanraw = Datasammelung(:,start:finish);
baseline = zeros(matrix_x,finish-start+1);

% subtract the baseline from raw data
for b = 2:length(Dataforramanraw(:,1))
    for c = 1:length(Dataforramanraw(1,:))
        
        %calaulate baseline
        baseline(b-1,c) = Simulationforoffset(b-1,1) * Dataforramanraw(1,c) + Simulationforoffset(b-1,2);
        
        % substract baseline
        Dataforramanoffset(b-1,c) = Dataforramanraw(b,c) - baseline(b-1,c);  
    end
end
% % calculated SN ration
%  signal = max(Dataforramanoffset(62,:));
% % noisearray = Dataforramanoffset(62,end-40:end);
% % noise = std(noisearray);
% % SR= signal/noise;

% fitting with calculated N2 Raman spectra
Wavelength_raman = Dataforramanraw(1,:);
for n=1:length(Dataforramanoffset(:,1))
    
    % normlization
    Iexp= Dataforramanoffset(n,:)./max( Dataforramanoffset(n,:));

    % x(1):temperature
    % x(2):offset % set as 0
    % x(3):shift
    % x(4):wG
    % x(5):lG
    % wG,lG for Voigt profile
    % x(6):stretch
    errorlsqfunraman = @(x)Intensity_sim(x(1),x(2),x(3),x(4),x(5),x(6))-Iexp;
    
    % initial values
    initial0 = [1000, 0, 0, 5, 1, 0.001];
    
    % upper boundary
    ub = [3000, 0, 10, 10, 2, 0.001];
    
    % lower boundary
    lb = [290 , 0, -20, 0, 0, 0.001];

    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxFunctionEvaluations',10^6);
    [Tfit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(errorlsqfunraman,initial0,lb,ub,options);
    ci = nlparci(Tfit,residual,'Jacobian',jacobian,'alpha',0.95);
    err = ((ci(:,2)-ci(:,1))./3.92);
    relerr = err./Tfit';
    nullarray = zeros(size(Wavelength_raman));


    % figure 
    figure(1)
    subplot(2,1,1);
    plot(Wavelength_raman,Iexp,Wavelength_raman,Intensity_sim(Tfit(1),Tfit(2),Tfit(3),Tfit(4),Tfit(5),Tfit(6)));
    xlim([2280 2345]);
    ylim([0 1]);
    subplot(2,1,2);
    plot(Wavelength_raman,Intensity_sim(Tfit(1),Tfit(2),Tfit(3),Tfit(4),Tfit(5),Tfit(6))-Iexp,Wavelength_raman,nullarray);hold on;
    plotnameforsave = num2str(HABnamecol(n));
    xlim([2280 2345]);
    ylim([-0.25 0.25]);
    saveas(gcf,strcat(plotpath,strcat(plotnameforsave,'_Frame_',num2str(mod(n,5)),'_',num2str(Tfit(1)),'.jpg')));
    delete(gcf);
    
    % save fitting results and confidential interval
    Simulationforraman(n,:)=Tfit;
    ciArray(2*n-1,1) = Tfit(1);
    ciArray(2*n-1:2*n,2:7) = ci';
    errArray(n,1) = Tfit(1);
    errArray(n,2:7) = err';
    relerrArray(n,1) = Tfit(1);
    relerrArray(n,2:7) = relerr';
    % save simulated N2 spectrum
    simulationdataforraman(n,:) = Intensity_sim(Tfit(1),Tfit(2),Tfit(3),Tfit(4),Tfit(5),Tfit(6));
    
    simulationdataforramanerror(n,:) = Intensity_sim(Tfit(1),Tfit(2),Tfit(3),Tfit(4),Tfit(5),Tfit(6))-Iexp;
end

% Merge HAB and simulated Raman signal intensity 
finalresult = [HABnamecol,Simulationforraman];

% calculation of meanvalue and standard value of temperature at each HAB
meanvalue = zeros(LengthFiles,1);
standardvalue = zeros(LengthFiles,1);
for index = 1:(length(Simulationforraman(:,1))/5)
    matrixformean = Simulationforraman(5*(index-1)+1:5*(index-1)+5);
    meanvalue(index,1)=mean(matrixformean);
    standardvalue(index,1)=std(matrixformean);
end

%%%%%%%%%%%%%%%%%%%
%   save data     %
%%%%%%%%%%%%%%%%%%%

simulationdatafororigin = [Wavelength_raman',simulationdataforraman'];
simulationdatafororiginerror = [Wavelength_raman',simulationdataforramanerror'];
datafilename = strcat(plotpath,'0_expt.xlsx');
writematrix(Datafinalsammelung,datafilename,'Sheet',1,'Range','A1');
writematrix(HAB_voll,datafilename,'Sheet',2,'Range','A1');
writematrix(Dataforramanoffset,datafilename,'Sheet',2,'Range','A2');
writematrix(ciArray,datafilename,'Sheet',3,'Range','A1');
writematrix(errArray,datafilename,'Sheet',4,'Range','A1');
writematrix(relerrArray,datafilename,'Sheet',5,'Range','A1');


datafilename = strcat(plotpath,'0_simu.xlsx');
writematrix(finalresult,datafilename,'Sheet',1,'Range','A1');
writematrix(HAB,datafilename,'Sheet',1,'Range','H2');
writematrix(meanvalue,datafilename,'Sheet',1,'Range','I2');
writematrix(standardvalue,datafilename,'Sheet',1,'Range','J2');
writematrix(HABnamecol',datafilename,'Sheet',2,'Range','B1');
writematrix(simulationdatafororigin,datafilename,'Sheet',2,'Range','A2');
writematrix(HABnamecol',datafilename,'Sheet',3,'Range','B1');
writematrix(simulationdatafororiginerror,datafilename,'Sheet',3,'Range','A2');
matfilename = strcat(plotpath,'0_expt_.mat');
save(matfilename);





%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%


