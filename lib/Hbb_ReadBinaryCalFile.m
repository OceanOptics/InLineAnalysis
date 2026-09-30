% Reads Hyper-bb binary calibration file (.hbb_cal)

function cal = Hbb_ReadBinaryCalFile(filename)

fid = fopen(filename);

fseek(fid, 0, 'eof');
EOF = ftell(fid);
fseek(fid,0,'bof');

idx = 0;

while ftell(fid) ~= EOF
    
idx = idx + 1;

cal(idx).ID =                 char(fread(fid,24,'char')');
cal(idx).calLengthBytes =          fread(fid,1,'uint16');
cal(idx).processingVersion =       fread(fid,1,'uint16') / 100;
cal(idx).serialNumber =            fread(fid,1,'uint16');
cal(idx).date =           datetime(fread(fid,6,'uint8')') + calyears(1900);
cal(idx).tempCalDate =    datetime(fread(fid,6,'uint8')') + calyears(1900);
cal(idx).PMTReferenceGain =        fread(fid,1,'uint16');
cal(idx).PMTGamma =                fread(fid,1,'float');
cal(idx).PMTGammaRMSE =            fread(fid,1,'float');
cal(idx).gain1_2 =                 fread(fid,1,'uint16') / 1000;
cal(idx).gain1_2_std =             fread(fid,1,'float');
cal(idx).gain2_3 =                 fread(fid,1,'uint16') / 1000;
cal(idx).gain2_3_std =             fread(fid,1,'float');
cal(idx).muFactorPMTGain =         fread(fid,1,'uint16');
cal(idx).transmitReceiveDistance = fread(fid,1,'uint16') / 100;
cal(idx).plaqueReflectivity =      fread(fid,1,'uint8')  / 100;
numMuFactorWl =                    fread(fid,1,'uint8');
numDarkOffsetWl =                  fread(fid,1,'uint8'); 
numDarkOffsetPMTGain =             fread(fid,1,'uint8');
cal(idx).muFactorWl =              fread(fid,numMuFactorWl,'uint16') / 10;
cal(idx).muFactors =               fread(fid,numMuFactorWl,'float');
cal(idx).muFactorTempCorr =        fread(fid,numMuFactorWl,'float');
cal(idx).muFactorLEDTemp =         fread(fid,numMuFactorWl,'uint16') / 100;
cal(idx).darkOffsetWl =            fread(fid,numDarkOffsetWl,'uint16') / 10;
cal(idx).darkOffsetPMTGain =       fread(fid,numDarkOffsetPMTGain,'uint16');
darkScat1 =                        fread(fid,numDarkOffsetWl * numDarkOffsetPMTGain, 'float');
darkScat2 =                        fread(fid,numDarkOffsetWl * numDarkOffsetPMTGain, 'float');
darkScat3 =                        fread(fid,numDarkOffsetWl * numDarkOffsetPMTGain, 'float');

cal(idx).darkOffsetScat1 = reshape(darkScat1, numDarkOffsetWl, numDarkOffsetPMTGain);
cal(idx).darkOffsetScat2 = reshape(darkScat2, numDarkOffsetWl, numDarkOffsetPMTGain);
cal(idx).darkOffsetScat3 = reshape(darkScat3, numDarkOffsetWl, numDarkOffsetPMTGain);

end

fclose(fid);

