%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
cd ccf_fits/
file_list   = dir('*.fits');
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);
MJD         = zeros(N_FILE, 1);

for n = 1:N_FILE
    FILE    = file_name{n};
    spec    = fitsread(char(FILE));                                         % spectra of all orders
    H       = fitsinfo(char(FILE));                                         % fits information including header
    header  = H.PrimaryData.Keywords;                                       % get header only
    MJD(n)  = header{30, 2};
end

cd ../

ORDER           = 5;                                                        % Highest Hermite order 
coeff           = zeros((ORDER+1), N_FILE);
coeff_rvc       = zeros((ORDER+1), N_FILE);
RV_gauss        = zeros(N_FILE,1);


grid_size       = 0.1;
v               = (-4.7 : grid_size : 4.7)';                                % km/s
% VEL             = importdata('GL479_vels.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
h               = waitbar(0,'Please wait...');

for n = 1:N_FILE
    
    i           = n - 1;
    filename    = ['./ccf_dat/ccf', num2str(i), '.dat'];
    A           = importdata(filename);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0, 0, 2, 0] );
    b           = f.b;  % shift
    RV_gauss(n) = b;    

    for order = 0:ORDER
        temp                    = hermite_nor(order, v) * grid_size;
        temp_rvc                = hermite_nor(order, v - b) * grid_size;
        coeff(order+1, n)         = sum(A .* temp);  
        coeff_rvc(order+1, n)     = sum(A .* temp_rvc); 
    end
    
    waitbar( n / N_FILE )
end
close(h)  

for order = 0:ORDER
    
    [pxx,f] = plomb(coeff(order+1, :), MJD - min(MJD), 0.5);
    [pmax,lmax] = max(pxx);
    f0 = f(lmax);
    disp(['T_planet: ', num2str(1/f0)]);

    [pxx_rvc,f_rvc] = plomb(coeff_rvc(order+1, :), MJD - min(MJD), 0.5);
    [pmax_rvc,lmax_rvc] = max(pxx_rvc);
    f0_rvc = f_rvc(lmax_rvc);
    disp(['T_activity: ', num2str(1/f0_rvc)]);
    
    if 0
        [pxx_v,f_v] = plomb(VEL(:,2), MJD - min(MJD), 0.5);
        [pmax_v,lmax_v] = max(pxx_v);
        f0_v = f(lmax_v);
    end
    
    h = figure;
        hold on
        
        % plot(f_v, pxx_v / max(pxx_v), 'k--')
        plot(f, pxx / max(pxx), 'r')
        plot(f_rvc, -pxx_rvc / max(pxx_rvc), 'b')
        %plot(f0_v, pmax_v / max(pxx_v), 'ko')
        %text(f0_v, pmax_v / max(pxx_v), ['\leftarrow', num2str(1/f0_v)]);        
        plot(f0, pmax / max(pxx), 'ro')
        text(f0, pmax / max(pxx), ['\leftarrow', num2str(1/f0)]);
        plot(f0_rvc, -pmax_rvc / max(pxx_rvc), 'bo')
        text(f0_rvc, -pmax_rvc / max(pxx_rvc), ['\leftarrow', num2str(1/f0_rvc)])
        legend('RV', 'Rest frame', 'Observed frame', 'Location', 'Best')
        xlabel('Frequency')
        ylabel('Normalized Power')
        title_name = ['Order', num2str(order)];
        title(title_name);
    hold off

    out_eps = [title_name, '.eps'];
    print(out_eps, '-depsc')
    close(h);
end