function [ data ] = importLISST200x_raw( filename, verbose )
% Import LISST200X data logged internally
% Author: Guillaume Bourdin
% Date: 2022-10-06
%%
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end
[param.zsc, param.fzs, param.dcal, raw_rata, ~, ~, param.config, param.housek] = parse200XBinary(filename);
% Convert into table
data = array2table(raw_rata);
% merge ring values variable
data = mergevars(data,1:36);

% assign variable names
data.Properties.VariableNames = {
'RingValues', ... % 1:36
'LaserTransmission', ...
'SupplyVoltage', ...
'ExternalAnalogInput1', ...
'LaserReference', ...
'Depth', ...
'Temperature', ...
'Year', ...
'Month', ...
'Day', ...
'Hour', ...
'Minute', ...
'Second', ...
'ExternalAnalogInput2', ...
'SauterMeanDiameter', ...
'TotalVolumeConcentration', ...
'RelativeHumidity', ...
'AccelerometerX', ...
'AccelerometerY', ...
'AccelerometerZ', ...
'RawPressureMSB', ...
'RawPressureLSB', ...
'AmbientLight', ...
'ExternalAnalogInput3' ...
};

% add datetime variable
data = addvars(data, datenum(data.Year, data.Month, data.Day, data.Hour, data.Minute, ...
  data.Second), 'NewVariableNames', 'dt', 'Before', 'RingValues');
f = fieldnames(param);
c = struct2cell(param);
data = addprop(data, f, repmat({'table'}, size(f)));
for k = 1:numel(f)
    data.Properties.CustomProperties.(f{k}) = c{k};
end

if verbose; fprintf('Done\n'); end
end

function [zsc,fzs,dcal,data,Tv,Ta,config,housek] = parse200XBinary(datafile)
  % Sequoia Scientific LISST200X import function
  % Date: 2022-03-22
  % Author: Thomas Leeuw

  % record IDs
  DCAL_RECORD_ID = 44778;
  TV1_RECORD_ID = 49391;
  TV2_RECORD_ID = 49394;
  TA1_RECORD_ID = 44047;
  TA2_RECORD_ID = 44274;
  ZSCAT_RECORD_ID = 47820;
  FZSCAT_RECORD_ID = 64428;
  CONFIG_RECORD_ID = 19529;
  HOUSEK_RECORD_ID = 52928;
  DATA_RECORD_ID = 56026;
  fid = fopen(datafile,'r','b'); % open file for reading using big endian format
  fseek(fid, 0, 'eof');
  fileSize = ftell(fid); % get the file size
  fseek(fid, 0, 'bof');
  recordID = fread(fid,1,'uint16'); % read the first record ID
  % calculate number of records in the file
  RecordSize = 120; % 120 bytes per record
  numRecords = fileSize/RecordSize;

  % Number of data records should be an integer
  if floor(numRecords) ~= numRecords
     warning('File contains incomplete data records');
     numRecords = floor(numRecords);
  end
  % calculate the number of variables in each record
  num16 = (RecordSize/2)-1;
  num32 = floor((RecordSize-2)/4);
  [zsc,fzs,dcal,Tv,Ta,config,housek] = deal([]);
  data = NaN(numRecords,num16); % preallocate data array for big speed increase
  % loop through the file and read in the data according to the record ID
  for recordNumber = 1:numRecords
      switch recordID
          case ZSCAT_RECORD_ID
              zsc = fread(fid,num16,'uint16')';
          case DATA_RECORD_ID
              data(recordNumber,:) = fread(fid,num16,'uint16')';
          case DCAL_RECORD_ID
              dcal = fread(fid,num16,'uint16')';
              dcal = dcal(2:37)./dcal(1);
          case FZSCAT_RECORD_ID
              fzs = fread(fid,num16,'uint16')';
          case CONFIG_RECORD_ID
              fseek(fid,-2,'cof');
              config.name = fread(fid,20,'*char')';
              config.serialNumber = fread(fid,1,'uint16');
              config.firmwareVer = fread(fid,1,'uint16').*0.001;
              config.VCC = fread(fid,1,'uint32');
              config.fullPath = fread(fid,1,'uint16').*0.01;
              config.effPath = fread(fid,1,'uint16').*0.01;
              config.bioBlock = fread(fid,1,'uint8');
              config.sTube = fread(fid,1,'uint8');
              config.analogConcScale = fread(fid,1,'uint16');
              config.endcap = fread(fid,1,'uint16');
              config.startCond = fread(fid,1,'uint16');
              config.startCondData = fread(fid,20,'*char')';
              config.stopCond = fread(fid,1,'uint16');
              config.stopCondData = fread(fid,20,'*char')';
              config.measurementAve = fread(fid,1,'uint16');
              config.sampleInterval = fread(fid,1,'uint16');
              config.sampleMode = fread(fid,1,'uint16');
              config.burstSamples = fread(fid,1,'uint16');
              config.burstInterval = fread(fid,1,'uint16');
              config.transmitRaw = fread(fid,1,'uint16');
              config.lifetimeSamples = fread(fid,1,'uint32');
              config.lifetimeLaserOn = fread(fid,1,'uint32');
              config.supportBoard = fread(fid,1,'uint16');
              config.ambientLight = fread(fid,1,'uint16');
              config.uvOn = fread(fid,1,'uint16');
              config.uvOff = fread(fid,1,'uint16');
              config.autoStart = fread(fid,1,'uint8');
              config.auxPower = fread(fid,1,'uint8');
              config.sensorWarmup = fread(fid,1,'uint8');
          case HOUSEK_RECORD_ID
              housek=fread(fid,num32,'float')';
          case TV1_RECORD_ID
              Tv(1:num32)=fread(fid,num32,'float')';
          case TV2_RECORD_ID
              Tv(30:36)=fread(fid,7,'float')';
          case TA1_RECORD_ID
              Ta(1:num32)=fread(fid,num32,'float')';
          case TA2_RECORD_ID
              Ta(30:36)=fread(fid,7,'float')';
          otherwise
              warning('Unrecognized data record ID found in .RBN file')
      end
      fseek(fid,recordNumber * RecordSize,'bof'); % go the location of the next record ID
      recordID = fread(fid,1,'uint16');        % read the next record ID
  end
  fclose(fid);
  % remove NaN rows from data matrix that correspond to header data rows
  data(all(isnan(data),2),:) = [];
  % negative ring values are possible, data must be corrected
  data(data(:,1:36)>40950) = data(data(:,1:36)>40950) - 65536;
  fzs(fzs(:,1:36)>40950) = fzs(fzs(:,1:36)>40950) - 65536;
  zsc(zsc(:,1:36)>40950) = zsc(zsc(:,1:36)>40950) - 65536;
end

