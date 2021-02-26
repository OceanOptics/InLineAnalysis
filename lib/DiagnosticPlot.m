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
        if all(cell2mat(dt(2)) < min(p.dt) | cell2mat(dt(1)) > max(p.dt))
            error('Check "dt" argument (out of boundary)')
        end  
        idp = p.dt >= cell2mat(dt(1)) & p.dt <= cell2mat(dt(2));
        if size(p(idp,:),1) > 400000
            error('Data too large, please reduce date & time interval (dt)')
        end
        if contains(i,'BB')
            visProd3D(wl, p.dt(idp,:), p.bbp(idp,:), false, 'Wavelength', false, j*2+1);
            zlabel([(level{j}) ' bb p (m^{-1})']);
            xlabel('lambda (nm)');
            ylabel('time');
        elseif contains(i,'AC')
            visProd3D(wl, p.dt(idp,:), p.c(idp,:), false, 'Wavelength', false, j*2+1);
            zlabel([(level{j}) ' c p (m^{-1})']);
            xlabel('lambda (nm)');
            ylabel('time');
            visProd3D(wl, p.dt(idp,:), p.a(idp,:), false, 'Wavelength', false, j*3+10);
            zlabel([(level{j}) ' a p (m^{-1})']);
            xlabel('lambda (nm)');
            ylabel('time');
        end
    else
        if all(isempty(data.(level{j}).tsw) & isempty(data.(level{j}).fsw))
            error('%s %s data missing', i, (level{j}))
        end
        if any(cell2mat(dt(2)) < min(data.(level{j}).tsw.dt) | cell2mat(dt(1)) > max(data.(level{j}).tsw.dt))
            error('Check "dt" variable (out of boundary)')
        end
        if ~isempty(data.(level{j}).tsw)
            if size(data.(level{j}).tsw(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:),1) > 400000
                error('Data too large, please reduce date & time interval (dt)')
            end
        end
        if ~isempty(data.(level{j}).fsw)
            if size(data.(level{j}).fsw(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:),1) > 400000
                error('Data too large, please reduce date & time interval (dt)')
            end
        end
        if contains(i,'BB')
            if isempty(data.(level{j}).tsw)
                visProd3D(wl, data.(level{j}).fsw.dt(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).fsw.beta(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+2);
                zlabel([(level{j}) ' bb filt(m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
            elseif isempty(data.(level{j}).fsw)               
                visProd3D(wl, data.(level{j}).tsw.dt(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).tsw.beta(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+1);
                zlabel([(level{j}) ' bb tot (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
            else
                visProd3D(wl, data.(level{j}).tsw.dt(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).tsw.beta(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+1);
                zlabel([(level{j}) ' bb tot (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
                visProd3D(wl, data.(level{j}).fsw.dt(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).fsw.beta(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+2);
                zlabel([(level{j}) ' bb filt(m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
            end
        elseif contains(i,'AC')
            if isempty(data.(level{j}).tsw)
                visProd3D(wlc, data.(level{j}).fsw.dt(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).fsw.c(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+2);
                zlabel([(level{j}) ' c filt (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
                visProd3D(wla, data.(level{j}).tsw.dt(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).tsw.a(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*3+10);
                zlabel([(level{j}) ' a tot (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
            elseif isempty(data.(level{j}).fsw)
                visProd3D(wlc, data.(level{j}).tsw.dt(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).tsw.c(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+1);
                zlabel([(level{j}) ' c tot (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
                visProd3D(wla, data.(level{j}).fsw.dt(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).fsw.a(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*3+11);
                zlabel([(level{j}) ' a filt (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
            else
                visProd3D(wlc, data.(level{j}).tsw.dt(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).tsw.c(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+1);
                zlabel([(level{j}) ' c tot (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
                visProd3D(wlc, data.(level{j}).fsw.dt(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).fsw.c(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*2+2);
                zlabel([(level{j}) ' c filt (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
                visProd3D(wla, data.(level{j}).tsw.dt(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).tsw.a(data.(level{j}).tsw.dt >= cell2mat(dt(1)) & data.(level{j}).tsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*3+10);
                zlabel([(level{j}) ' a tot (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
                visProd3D(wla, data.(level{j}).fsw.dt(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), ...
                    data.(level{j}).fsw.a(data.(level{j}).fsw.dt >= cell2mat(dt(1)) & data.(level{j}).fsw.dt <= cell2mat(dt(2)),:), false, 'Wavelength', false, j*3+11);
                zlabel([(level{j}) ' a filt (m^{-1})']);
                xlabel('lambda (nm)');
                ylabel('time');
            end
        end
    end
end
end
