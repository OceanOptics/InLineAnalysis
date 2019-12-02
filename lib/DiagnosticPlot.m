function DiagnosticPlot(data, i, level, dt)
if contains(i,'BB')
    wl = data.lambda;
elseif contains(i,'AC')
    wla = data.lambda_a;
    wlc = data.lambda_c;
else
    error('Intrument not supported')
end

for j = 1:length(level)
    if strcmp(level{j},'prod')
        if isempty(data.(level{j}).p)
            error('%s %s data missing', i, (level{j}))
        end  
        p = data.(level{j}).p;
        if cell2mat(dt(2)) < min(p.dt) | cell2mat(dt(1)) > max(p.dt)
            error('Check "dt" argument (out of boundary)')
        end  
        idp = p.dt >= cell2mat(dt(1)) & p.dt <= cell2mat(dt(2));
        if size(p(idp,:),1) > 400000
            error('Data too large, please reduce date & time interval (dt)')
        end
        if contains(i,'BB')
            visProd3D(wl, p.dt(idp,:), p.bb(idp,:), false, 'Wavelength', false, j*2+1); zlabel([(level{j}) ' bb p (m^{-1})']); %, 'Wavelength', true
        elseif contains(i,'AC')
            visProd3D(wl, p.dt(idp,:), p.c(idp,:), false, 'Wavelength', false, j*2+1); zlabel([(level{j}) ' c p (m^{-1})']); %, 'Wavelength', true
            visProd3D(wl, p.dt(idp,:), p.a(idp,:), false, 'Wavelength', false, j*3+10); zlabel([(level{j}) ' a p (m^{-1})']); %, 'Wavelength', true
        end
    else
        if isempty(data.(level{j}).tsw) | isempty(data.(level{j}).fsw)
            error('%s %s data missing', i, (level{j}))
        end  
        tot = data.(level{j}).tsw;
        filt = data.(level{j}).fsw;
        if cell2mat(dt(2)) < min(tot.dt) | cell2mat(dt(1)) > max(tot.dt)
            error('Check "dt" argument (out of boundary)')
        end  
        idtot = tot.dt >= cell2mat(dt(1)) & tot.dt <= cell2mat(dt(2));
        idfilt = filt.dt >= cell2mat(dt(1)) & filt.dt <= cell2mat(dt(2));
        if size(tot(idtot,:),1) > 400000 | size(filt(idfilt,:),1) > 400000
            error('Data too large, please reduce date & time interval (dt)')
        end
        if contains(i,'BB')
            visProd3D(wl, tot.dt(idtot,:), tot.bb(idtot,:), false, 'Wavelength', false, j*2+1); zlabel([(level{j}) ' bb tot (m^{-1})']); %, 'Wavelength', true
            visProd3D(wl, filt.dt(idfilt,:), filt.bb(idfilt,:), false, 'Wavelength', false, j*2+2); zlabel([(level{j}) ' bb filt(m^{-1})']); %, 'Wavelength', true
        elseif contains(i,'AC')
            visProd3D(wlc, tot.dt(idtot,:), tot.c(idtot,:), false, 'Wavelength', false, j*2+1); zlabel([(level{j}) ' c tot (m^{-1})']); %, 'Wavelength', true
            visProd3D(wlc, filt.dt(idfilt,:), filt.c(idfilt,:), false, 'Wavelength', false, j*2+2); zlabel([(level{j}) ' c filt (m^{-1})']); %, 'Wavelength', true
            visProd3D(wla, tot.dt(idtot,:), tot.a(idtot,:), false, 'Wavelength', false, j*3+10); zlabel([(level{j}) ' a tot (m^{-1})']); %, 'Wavelength', true
            visProd3D(wla, filt.dt(idfilt,:), filt.a(idfilt,:), false, 'Wavelength', false, j*3+11); zlabel([(level{j}) ' a filt (m^{-1})']); %, 'Wavelength', true
        end
    end
end
end
