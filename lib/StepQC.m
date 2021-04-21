function [data, Nbad] = StepQC (data, lambda, fudge_factor, bb_dark, bb_threshold)
% author: Guillaume Bourdin
% created: Nov 13, 2019
%
% StepQC Auto QC ACs/AC9/BB spectrum
%
% INPUT:
%   - data_in: ACS data <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectrum <NxM double> in column
%         - c spectrum <NxM double> in column
%     or for BB data
%         - bb spectrum <NxM double> in column
%   - lambda: structure with three fields:
%       - a: wavelenghts for absorption <1xM double> 
%       - c: wavelenghts for attenuation <1xM double>
%       - bb: wavelenghts of BB sensor <1xM double>
%   - fudge_factor: structure with three fields:
%       - a: threshold for automatic QC of absorption spectrum (varies between ACS, must be >= 3)
%           (default = 3: Step = 3x > mean(diff))
%       - c: threshold for automatic QC of attenuation spectrum (varies between ACS, must be >= 3)
%           (default = 3: Step = 3x > mean(diff))
%       - bb: threshold for automatic QC of beta spectrum (must be >= 3)
%           (default = 3: Step = beta(lambda(i)) > 3 * mean(beta(lambda(~i))
%   - bb_dark: beta dark measurements <1xM double>
%   - bb_threshold: threshold for automatic QC of beta spectrum saturated (default = 4100)
%           (bb sensor saturation = 4130 counts)
%
% OUTPUT:
%   - data_in: quality filtered ACS data <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectrum <NxM double> in column
%         - c spectrum <NxM double> in column
%     or for BB data
%         - bb spectrum <NxM double> in column
%   - Nbad: structure with three fields:
%       - a: <double> percentage of absorption spectrum deleted
%       - c: <double> percentage of attenuation spectrum deleted
%       - bb: <1xM double> percentage of beta vales deleted

if nargin < 2
    error('input wavelenghts missing');
end
if nargin < 1
    error('input data in missing');
end
if nargin < 3
    fudge_factor.a = 3;
    fudge_factor.c = 3;
    warning('fudge_factor missing, set to default (3)');
end

if any(strcmp(data.Properties.VariableNames, 'a') | ...
        strcmp(data.Properties.VariableNames, 'c'))
    if any(fudge_factor.a < 3 | fudge_factor.c < 3)
        warning('QC threshold might be too low, data might be lost');
    end

    % lambda.a = ila.instrument.ACS091.lambda_a; lambda.c = ila.instrument.ACS091.lambda_c;
    % lambda.ref = ila.instrument.ACS091.lambda_ref;
    % data = ila.instrument.ACS091.data; 
    % dataini = ila.instrument.ACS091.data;

    wl_a = lambda.a(lambda.a > 500 & lambda.a < 600);
    wl_c = lambda.c(lambda.c > 500 & lambda.c < 600);

    diff_a = [diff(data.a(:,lambda.a > 500 & lambda.a < 600),[],2) NaN(size(data,1),1)];
    diff_c = [diff(data.c(:,lambda.c > 500 & lambda.c < 600),[],2) NaN(size(data,1),1)];

    bad_a = max(abs(diff_a(:,wl_a > 560 & wl_a < 600)),[],2)...
              > fudge_factor.a*mean(abs(diff_a(:,wl_a > 500 & wl_a < 550)),2);
    bad_c = max(abs(diff_c(:,wl_c > 560 & wl_c < 600)),[],2)...
              > fudge_factor.c*mean(abs(diff_c(:,wl_c > 500 & wl_c < 550)),2);

    % Nbad_tot = sum(bad_a & bad_c)/size(data,1)*100;
    Nbad.a = sum(bad_a)/size(data,1)*100;%-Nbad_tot;
    Nbad.c = sum(bad_c)/size(data,1)*100;%-Nbad_tot;
    data.a(bad_a,:) = NaN;
    data.c(bad_c,:) = NaN;
    data(bad_a & bad_c,:)=[];
    
elseif any(strcmp(data.Properties.VariableNames, 'beta'))
    if nargin < 5
        bb_threshold = 4100;
        warning('beta threshold missing, set to default (4100 counts)');
    end
    
    if nargin < 4
        error('beta dark missing (4th argument)');
    end
    Nbad.bb = NaN(1, size(lambda.bb,2));
    for ii = 1:size(lambda.bb,2)
        Nbad.bb(ii) = sum(data.beta(:,ii) >= bb_threshold);
        data.beta(data.beta(:,ii) >= bb_threshold, ii) = NaN;
    end
    
    if any(fudge_factor.bb < 3)
        error('QC threshold BB3 too low (minimum 3)');
    end
    bad_bb = NaN(size(data,1), size(lambda.bb, 2));
    for ii = 1:size(lambda.bb, 2)
        other = 1:size(lambda.bb, 2); other(ii)=[];
        bad_bb (:,ii) = data.beta(:,ii) - bb_dark(ii) > fudge_factor.bb * (nanmean(data.beta(:,other),2) - nanmean(bb_dark(other),2));
        Nbad.bb(ii) = Nbad.bb(ii) + sum(bad_bb(:,ii));
    end
    data.beta(logical(bad_bb)) = NaN;
    Nbad.bb = Nbad.bb / size(data.beta(:,ii),1) * 100;
end
end


% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.a(100000:101000,:), false, 'Wavelength', false, 72); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.c(100000:101000,:), false, 'Wavelength', false, 73); zlabel('c_p (m^{-1})');
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.a(100000:121000,:), false, 'Wavelength', false, 74); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.c(100000:121000,:), false, 'Wavelength', false, 75); zlabel('c_p (m^{-1})');
