% Reads Hyper-bb binary temperature calibration file (.hbb_tcal)

function cal_temp = Hbb_ReadBinaryTempCalFile(filename)

fid = fopen(filename);

fseek(fid, 0, 'eof');
EOF = ftell(fid);
fseek(fid, 0, 'bof');

idx = 0;

while ftell(fid) ~= EOF

idx = idx + 1;

cal_temp(idx).ID =                 char(fread(fid,24,'char')');
cal_temp(idx).processingVersion =       fread(fid,1,'uint16') / 100;
cal_temp(idx).serialNumber =            fread(fid,1,'uint16');
cal_temp(idx).date =           datetime(fread(fid,6,'uint8')') + calyears(1900);
cal_temp(idx).normalizedTemp =          fread(fid,1,'uint16') / 100;
cal_temp(idx).polynomialOrder =         fread(fid,1,'uint8');
numWLs =                           fread(fid,1,'uint8');
cal_temp(idx).wl =                      fread(fid,numWLs,'uint16') / 10;
coeffs =                           fread(fid,numWLs * (cal_temp(idx).polynomialOrder + 1), 'float');

cal_temp(idx).coeff = reshape(coeffs, numWLs, (cal_temp(idx).polynomialOrder + 1));

end

fclose(fid);
