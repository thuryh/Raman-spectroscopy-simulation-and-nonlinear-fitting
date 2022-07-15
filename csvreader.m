%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
% read data from csv files                         %
% output nu_expt: Wavenumber                       %
% output intensity_expt: Raman signal intensity    %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nu_expt,intensity_expt] = csvreader(csv_filename,csv_folderpath)

% cd to the directory, where the given csv-file is located
old_folder = cd(csv_folderpath);

% read nu 
xdaten = csvread(csv_filename,2,1,[2,1,2,1024]);

% nm to cm-1
for i =1:length(xdaten)
    xdaten(i) = 1e7*(1/354.747 - 1/xdaten(i));
end

% pick up data, 2230 - 2490 cm-1 
xdaten_need = xdaten(xdaten>2230 & xdaten<2490);
xdaten_pos = find(xdaten>2230 & xdaten<2490);
start = xdaten_pos(1);
finish = xdaten_pos(end);

% preallcaoting
intensity = zeros(5,211);

% 40 rows for frame = 5, 85 for frame =10
% read intensity
 for row = 4:9:40
     ydaten = csvread(csv_filename,row,start,[row,start,row,finish]);
     intensity((row-4)/9 + 1,:)=ydaten;
 end
 
intensity_expt = intensity;

% average Frames
%intensity_expt = mean(intensity);

% get rid of blemish or error point
threshhold = 2000;
for j=1:length(intensity_expt)
    if j==1
        if abs(intensity_expt(j) - intensity_expt(j+1))>=threshhold || intensity_expt(j)<=0
            intensity_expt(j) = intensity_expt(j+1);
        end
    elseif j == length(intensity_expt)
        if abs(intensity_expt(j) - intensity_expt(j-1))>=threshhold || intensity_expt(j)<=0
            intensity_expt(j) = intensity_expt(j-1);
        end
    elseif intensity_expt(j)<=0
        intensity_expt(j) = 0.5*(intensity_expt(j-1)+intensity_expt(j+1));
    elseif abs(intensity_expt(j) - intensity_expt(j+1))>=threshhold && abs(intensity_expt(j) - intensity_expt(j-1))>=threshhold
        intensity_expt(j) = 0.5*(intensity_expt(j-1)+intensity_expt(j+1));
    end
end

% output
nu_expt = xdaten_need;

% normalisiren
% intensity_expt = intensity_expt./max(intensity_expt);

% Change back to the function's directory
cd(old_folder);

end


%%%%%%%%%%%%%%%
%%%%% END %%%%%
%%%%%%%%%%%%%%%

