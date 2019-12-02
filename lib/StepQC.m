function [data, Nbad_a, Nbad_c] = StepQC (data, lambda_a, lambda_c, fudge_factor_a, fudge_factor_c)
% author: Guillaume Bourdin
% created: Nov 13, 2019

%Detect filter events on ACS and BB3 inline
%
% INPUT:
%   - data_in: ACS data <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectrum <NxM double> in column
%         - c spectrum <NxM double> in column
%   - lambda_a: wavelenghts for absorption <1xM double> 
%   - lambda_c: wavelenghts for attenuation <1xM double> 
%   - fudge_factor_a: threshold for automatic QC of absorption spectrum (varies between ACS, must be >= 3)
%       (default = 3: Step = 3x > mean(diff))
%   - fudge_factor_c: threshold for automatic QC of attenuation spectrum (varies between ACS, must be >= 3)
%       (default = 3: Step = 3x > mean(diff))
%
% OUTPUT:
%   - data_in: quality filtered ACS data <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectrum <NxM double> in column
%         - c spectrum <NxM double> in column
%   - Nbad_a: <double> percentage of absorption spectrum deleted
%   - Nbad_c: <double> percentage of attenuation spectrum deleted

if fudge_factor_a < 3 | fudge_factor_c < 3
    error('QC threshold too low (minimum 3)');
end
if nargin < 3
    error('input attenuation wavelenghts missing');
end
if nargin < 2
    error('input absorption wavelenghts missing');
end
if nargin < 1
    error('input data in missing');
end

% lambda_a = ila.instrument.ACS091.lambda_a; lambda_c = ila.instrument.ACS091.lambda_c;
% lambda_ref = ila.instrument.ACS091.lambda_ref;
% data = ila.instrument.ACS091.data; 
% dataini = ila.instrument.ACS091.data;

wl_a = lambda_a(lambda_a > 500 & lambda_a < 600);
wl_c = lambda_c(lambda_c > 500 & lambda_c < 600);

diff_a = [diff(data.a(:,lambda_a > 500 & lambda_a < 600),[],2) NaN(size(data,1),1)];
diff_c = [diff(data.c(:,lambda_c > 500 & lambda_c < 600),[],2) NaN(size(data,1),1)];

bad_a = max(abs(diff_a(:,wl_a > 560 & wl_a < 600)),[],2)...
          > fudge_factor_a*mean(abs(diff_a(:,wl_a > 500 & wl_a < 550)),2);
bad_c = max(abs(diff_c(:,wl_c > 560 & wl_c < 600)),[],2)...
          > fudge_factor_c*mean(abs(diff_c(:,wl_c > 500 & wl_c < 550)),2);

% Nbad_tot = sum(bad_a & bad_c)/size(data,1)*100;
Nbad_a = sum(bad_a)/size(data,1)*100;%-Nbad_tot;
Nbad_c = sum(bad_c)/size(data,1)*100;%-Nbad_tot;
data.a(bad_a,:) = NaN;
data.c(bad_c,:) = NaN;
data(bad_a & bad_c,:)=[];
end

% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.a(100000:101000,:), false, 'Wavelength', false, 72); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.c(100000:101000,:), false, 'Wavelength', false, 73); zlabel('c_p (m^{-1})');
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.a(100000:121000,:), false, 'Wavelength', false, 74); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.c(100000:121000,:), false, 'Wavelength', false, 75); zlabel('c_p (m^{-1})');
