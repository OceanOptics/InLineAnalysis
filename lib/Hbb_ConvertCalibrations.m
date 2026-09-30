% Hbb_ConvertCalibrations will convert .mat plaque and temperature
% calibration files (obsolete) to a .hbb_cal and .hbb_tcal files, respectivly.
% The new files (text based) are used with Hyper-bb software version 2.0
% and greater.
%
% INPUTS: 'PlaqueCal' and 'TemperatureCal' are a string path the plaque and
% temperature .mat calibration files.
%
% OUTPUTS: The new files will be saved in the same location as the 
% provided .mat calibration files.
%
% Sequoia Scientific, Inc. - 10/15/2024

function Hbb_ConvertCalibrations(PlaqueCal, TemperatureCal)

plaque_cal = load(PlaqueCal);
plaque_cal = plaque_cal.cal;

temp_cal = load(TemperatureCal);
temp_cal = temp_cal.cal_temp;

if isfield(temp_cal,'serialNumber')
    SN = temp_cal.serialNumber;
else 
    prompt = {'Enter Instrument Serial Number'};
    dlgtitle = 'Serial Number';
    fieldsize = [1 40];
    answer = inputdlg(prompt,dlgtitle,fieldsize);
    
    if isempty(answer)
        return
    end
    
    SN = str2double(answer{1});
end

if ~isfield(plaque_cal,'processingVersion')
    plaque_cal.processingVersion = 1.0;
end

%% Plaque Cal

[filepath,name,~] = fileparts(PlaqueCal); 
outputFilename = fullfile(filepath,[name '.hbb_cal']);

if isfile(outputFilename)
  % Ask user if they want to overwrite the file.
  promptMessage = sprintf('This file already exists:\n%s\nDo you want to overwrite it?', PlaqueCal);
  titleBarCaption = 'Overwrite?';
  buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
  if strcmpi(buttonText, 'No')
    % User does not want to overwrite. 
    return
  end
end

tempCoeff_mu = GetTemperatureCoefficients(temp_cal,plaque_cal.muWavelengths,plaque_cal.muLedTemp);

fileID = fopen(outputFilename,'w');
startingPosition = ftell(fileID);
fwrite(fileID, 'Sequoia Hyper-bb Cal','char');
fwrite(fileID, zeros(1, 4),'uint8');
fwrite(fileID, 0, 'uint16'); % skip cal byte length field, fill in at the end
fwrite(fileID, plaque_cal.processingVersion * 100, 'uint16');
fwrite(fileID, SN, 'uint16');
fwrite(fileID, year(plaque_cal.timestamp) - 1900, 'uint8');
fwrite(fileID, month(plaque_cal.timestamp), 'uint8');
fwrite(fileID, day(plaque_cal.timestamp), 'uint8');
fwrite(fileID, hour(plaque_cal.timestamp), 'uint8');
fwrite(fileID, minute(plaque_cal.timestamp), 'uint8');
fwrite(fileID, second(plaque_cal.timestamp), 'uint8');
fwrite(fileID, year(temp_cal.timestamp) - 1900, 'uint8');
fwrite(fileID, month(temp_cal.timestamp), 'uint8');
fwrite(fileID, day(temp_cal.timestamp), 'uint8');
fwrite(fileID, hour(temp_cal.timestamp), 'uint8');
fwrite(fileID, minute(temp_cal.timestamp), 'uint8');
fwrite(fileID, second(temp_cal.timestamp), 'uint8');
fwrite(fileID, plaque_cal.pmtRefGain, 'uint16');
fwrite(fileID, plaque_cal.pmtGamma, 'float');
fwrite(fileID, 0, 'float'); % PMT gamma RMSE not computed
fwrite(fileID, round(plaque_cal.gain12 * 1000), 'uint16');
fwrite(fileID, 0, 'float'); % gain STD not computed;
fwrite(fileID, round(plaque_cal.gain23 * 1000), 'uint16');
fwrite(fileID, 0, 'float'); % gain STD not computed;
fwrite(fileID, 0, 'uint16'); % mu Factor PMT gain not recorded
fwrite(fileID, 31.4 * 100, 'uint16'); % Transmit/Receive beam distance (mm)
fwrite(fileID, 1.1 * 100, 'uint8'); % Plaque reflectivity (rho)
fwrite(fileID, numel(plaque_cal.muWavelengths), 'uint8');
fwrite(fileID, numel(plaque_cal.darkCalWavelength), 'uint8');
fwrite(fileID, numel(plaque_cal.darkCalPmtGain), 'uint8');
fwrite(fileID, round(plaque_cal.muWavelengths * 10), 'uint16');
fwrite(fileID, plaque_cal.muFactors, 'float');
fwrite(fileID, tempCoeff_mu, 'float');
fwrite(fileID, round(plaque_cal.muLedTemp * 100), 'uint16');  
fwrite(fileID, round(plaque_cal.darkCalWavelength * 10), 'uint16');
fwrite(fileID, plaque_cal.darkCalPmtGain, 'uint16');
fwrite(fileID, plaque_cal.darkCalScat1, 'float');
fwrite(fileID, plaque_cal.darkCalScat2, 'float');
fwrite(fileID, plaque_cal.darkCalScat3, 'float');

endingPosition = ftell(fileID);

% fill in cal byte length
fseek(fileID,24,'bof');
fwrite(fileID, endingPosition - startingPosition, 'uint16');

fclose(fileID);

msgbox([sprintf('Converted plaque calibration saved!\n\n') outputFilename]);

%% Temperature Cal
 
[filepath,name,~] = fileparts(TemperatureCal); 
outputFilename = fullfile(filepath,[name '.hbb_tcal']);

if isfile(outputFilename)
  % Ask user if they want to overwrite the file.
  promptMessage = sprintf('This file already exists:\n%s\nDo you want to overwrite it?', TemperatureCal);
  titleBarCaption = 'Overwrite?';
  buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
  if strcmpi(buttonText, 'No')
    % User does not want to overwrite. 
    return
  end
end

fileID = fopen(outputFilename,'w');
fwrite(fileID, 'Sequoia Hyper-bb T Cal','char')
fwrite(fileID, zeros(1, 2),'uint8');
fwrite(fileID, plaque_cal.processingVersion * 100, 'uint16');
fwrite(fileID, SN, 'uint16');
fwrite(fileID, year(temp_cal.timestamp) - 1900, 'uint8');
fwrite(fileID, month(temp_cal.timestamp), 'uint8');
fwrite(fileID, day(temp_cal.timestamp), 'uint8');
fwrite(fileID, hour(temp_cal.timestamp), 'uint8');
fwrite(fileID, minute(temp_cal.timestamp), 'uint8');
fwrite(fileID, second(temp_cal.timestamp), 'uint8');
fwrite(fileID, round(temp_cal.normalizedTemp * 100), 'uint16');
fwrite(fileID, 2, 'uint8') % polynomial order
fwrite(fileID, numel(temp_cal.wl), 'uint8')
fwrite(fileID, round(temp_cal.wl * 10), 'uint16');
fwrite(fileID, temp_cal.coeff, 'float');
fclose(fileID);

msgbox([sprintf('Converted temperature calibration saved!\n\n') outputFilename]);


end


function tempCoeff = GetTemperatureCoefficients(cal_temp,wavelength,temperature) 

% Generate temperature correction grid
LEDTempRange = [min(temperature):0.1:max(temperature)+0.1]; % need to make sure the max value is included
TempCorrGrid = NaN(length(cal_temp.wl),length(LEDTempRange));
for n = 1:length(cal_temp.wl)
   TempCorrGrid(n,:) = polyval(cal_temp.coeff(n,:),LEDTempRange);
end

% 2D interpolate (wavelength x temperature) to find temperature
% correction factor
tempCoeff = interp2(LEDTempRange, cal_temp.wl, TempCorrGrid, temperature, wavelength, 'linear');

end