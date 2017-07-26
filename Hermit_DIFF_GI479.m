% Branch from Hermit_coeff_nor.m
% DIFF;
% Suppose introducing a shift would not affect the diff; so I introduced
% the shift b, to move all line profiles to the centre. 

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
cd ccf_fits/
file_list   = dir('*.fits');
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);
MJD         = zeros(N_FILE, 1);

% read MJD %
for n = 1:N_FILE
    FILE = file_name{n};
    spec    = fitsread(char(FILE));                                         % spectra of all orders
    H       = fitsinfo(char(FILE));                                         % fits information including header
    header  = H.PrimaryData.Keywords;                                       % get header only
    MJD(n)  = header{30, 2};
end

cd ../

ORDER           = 5;                                                       % Highest Hermite order 
array_order     = 0:ORDER;
idx_even        = mod(0:ORDER, 2) == 0;
order_even      = array_order(idx_even);
order_odd       = array_order(~idx_even);

SHIFT           = 11;
Her_Spe         = zeros(N_FILE, ORDER);
Her_Spe_rvc     = zeros(N_FILE, ORDER);
coeff           = zeros((ORDER+1), SHIFT);
coeff_rvc       = zeros((ORDER+1), SHIFT);
RV_gauss        = zeros(N_FILE,1);

V_BEGIN         = -40;                                                      % m/s                                                     
V_END           = -V_BEGIN;
V_GRID          = (V_BEGIN : (V_END - V_BEGIN)/(SHIFT-1) : V_END) / 1000;   % km/s

grid_size       = 0.1;
v               = (-4.7 : grid_size : 4.7)';                                % km/s


%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
h               = waitbar(0,'Please wait...');
for n = 1:N_FILE
    
    i           = n - 1;
    filename    = ['./ccf_dat/ccf', num2str(i), '.dat'];
    A           = importdata(filename);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    b           = f.b;  % shift
    RV_gauss(n) = b;    

    for order = 1:ORDER+1
        for shift = 1:SHIFT
            temp                    = hermite_nor(order-1, v - V_GRID(shift)) * grid_size;
            temp_rvc                = hermite_nor(order-1, v - b - V_GRID(shift)) * grid_size;
            coeff(order, shift)     = sum(A .* temp);  
            coeff_rvc(order, shift) = sum(A .* temp_rvc); 
        end
    end

    % fitting %
    for order = order_even
        p = polyfit(V_GRID * 1000, coeff(order+1, :), 2);
        Her_Spe(n, order+1) = p(1);
        
        p_rvc = polyfit(V_GRID * 1000, coeff_rvc(order+1, :), 2);
        Her_Spe_rvc(n, order+1) = p_rvc(1);       
    end

    for order = order_odd
        p = polyfit(V_GRID * 1000, coeff(order+1, :), 1);
        Her_Spe(n, order+1) = p(1);
        
        p_rvc = polyfit(V_GRID * 1000, coeff_rvc(order+1, :), 2);
        Her_Spe_rvc(n, order+1) = p_rvc(1);           
    end
    
    waitbar( n / N_FILE )
end
close(h)  

for n_hermite = 0:ORDER
    
    [pxx,f] = plomb(Her_Spe(:, n_hermite+1), MJD - min(MJD), 0.5);
    [pmax,lmax] = max(pxx);
    f0 = f(lmax);
    disp(['T_planet: ', num2str(1/f0)]);

    [pxx_rvc,f_rvc] = plomb(Her_Spe_rvc(:, n_hermite+1), MJD - min(MJD), 0.5);
    [pmax_rvc,lmax_rvc] = max(pxx_rvc);
    f0_rvc = f_rvc(lmax_rvc);
    disp(['T_activity: ', num2str(1/f0_rvc)]);

    h = figure;
        hold on
        plot(f, pxx / max(pxx), 'r')
        plot(f_rvc, pxx_rvc / max(pxx_rvc), 'b')
        plot(f0, pmax / max(pxx), 'ro')
        
        text(f0, pmax / max(pxx), ['\leftarrow', num2str(1/f0)]);
        plot(f0_rvc, pmax_rvc / max(pxx_rvc), 'bo')
        text(f0_rvc, pmax_rvc / max(pxx_rvc), ['\leftarrow', num2str(1/f0_rvc)])
        legend('Rest frame', 'Observed frame', 'Location', 'Best')
        xlabel('Frequency')
        ylabel('Normalized Power')
        title_name = ['Order', num2str(n_hermite)];
        title(title_name);
    hold off

    out_eps = [title_name, '.eps'];
    print(out_eps, '-depsc')
    close(h);    
end

if 0
    %%%%%%%%%%%%
    % Plotting %
    %%%%%%%%%%%%
    for n = 1:N_FILE    

        h1 = figure;
        RELATIVE = Her_Spe(:, order_odd+1) - repmat(Her_Spe(1, order_odd+1), N_FILE, 1);
        bar(order_odd, Her_Spe(n, order_odd+1) - Her_Spe(1, order_odd+1), 'r')    
        ylim([min(min(RELATIVE)) max(max(RELATIVE))])
        title_name = num2str(MJD(n));
        title(title_name)
        % out_eps = ['ODD_', title_name, '.eps'];
        % print(out_eps, '-depsc')
        out_jpg = ['ODD_', title_name, '.jpg'];
        print(out_jpg, '-djpeg')    
        close(h1)

        h2 = figure;
        RELATIVE = Her_Spe(:, order_even+1) - repmat(Her_Spe(1, order_even+1), N_FILE, 1);
        bar(order_even, Her_Spe(n, order_even+1) - Her_Spe(1, order_even+1), 'b')    
        ylim([min(min(RELATIVE)) max(max(RELATIVE))])
        title_name = num2str(MJD(n));
        title(title_name)
        out_jpg = ['EVEN_', title_name, '.jpg'];
        print(out_jpg, '-djpeg')
        % out_eps = ['EVEN_', title_name, '.eps'];
        % print(out_eps, '-depsc')
        close(h2)

    end
end




if 0

    N_CCF           = 74;                               % number of CCF files
    N_hermite       = 5;                                % Highest Hermite order 
    coeff           = zeros((N_hermite+1), (N_CCF+1));
    coeff_noise     = zeros((N_hermite+1), (N_CCF+1));
    coeff_rvc       = zeros((N_hermite+1), (N_CCF+1));
    coeff_noise_rvc = zeros((N_hermite+1), (N_CCF+1));
    grid_size       = 0.1;
    v               = (-20 : grid_size : 20)';          % km/s
    RV              = importdata('RV.dat') / 1000;      % activity induced RV [km/s]
    RV_gauss        = zeros(size(0:N_CCF))';

    idx             = (v > -10) & (v < 10);
    v               = v(idx);

    % template %
    A_tpl           = 1 - importdata('CCF_tpl.dat');
    A_tpl           = A_tpl(idx);
    f_tpl           = fit( v, A_tpl, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    b_tpl           = f_tpl.b;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Coefficient %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    for n_CCF = 0:N_CCF

        v_planet    = 5 * sin(n_CCF/25/0.618*2*pi + 1) * 0.001; % km/s

        filename    = ['CCF_dat/CCF', num2str(n_CCF), '.dat'];
        A           = 1 - importdata(filename);
        A           = A(idx);
        f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
        b           = f.b;  % shift
        RV_gauss(n_CCF+1) = b;

        disp([n_CCF, b*1000, (b-b_tpl)*1000])

        for n_hermite = 0:N_hermite
            temp                                = A .* hermite_nor(n_hermite, v - b_tpl + v_planet) * grid_size;
            coeff(n_hermite+1, n_CCF+1)         = sum(temp);                                % Without noise
            coeff_noise(n_hermite+1, n_CCF+1)   = sum(temp .* (1+normrnd(0, A.^0.5/SN)));   % With noise

            temp_rvc                                = A .* hermite_nor(n_hermite, v - b) * grid_size;
            coeff_rvc(n_hermite+1, n_CCF+1)         = sum(temp_rvc);                                % Without noise
            coeff_noise_rvc(n_hermite+1, n_CCF+1)   = sum(temp_rvc .* (1+normrnd(0, A.^0.5/SN)));   % With noise
        end

    end     
    RV_gauss = RV_gauss - b_tpl;


    for n_hermite = 0:N_hermite

        [pxx_noise,f_noise] = plomb(coeff_noise(n_hermite+1, :), 0:N_CCF);
        [pmax_noise,lmax_noise] = max(pxx_noise);
        f0_noise = f_noise(lmax_noise);
        disp(['T_planet: ', num2str(1/f0_noise)]);

        [pxx_noise_rvc,f_noise_rvc] = plomb(coeff_noise_rvc(n_hermite+1, :), 0:N_CCF);
        [pmax_noise_rvc,lmax_noise_rvc] = max(pxx_noise_rvc);
        f0_noise_rvc = f_noise_rvc(lmax_noise_rvc);
        disp(['T_activity: ', num2str(1/f0_noise_rvc)]);

        h = figure;
            hold on
            plot(f_noise, pxx_noise / max(pxx_noise), 'r')
            plot(f_noise_rvc, pxx_noise_rvc / max(pxx_noise_rvc), 'b')
            plot(f0_noise, pmax_noise / max(pxx_noise), 'ro')

            text(f0_noise, pmax_noise / max(pxx_noise), ['\leftarrow', num2str(1/f0_noise)]);
            plot(f0_noise_rvc, pmax_noise_rvc / max(pxx_noise_rvc), 'bo')
            text(f0_noise_rvc, pmax_noise_rvc / max(pxx_noise_rvc), ['\leftarrow', num2str(1/f0_noise_rvc)])
            legend('Rest frame', 'Observed frame', 'Location', 'Best')
            xlabel('Frequency')
            ylabel('Normalized Power')
            title_name = ['Order', num2str(n_hermite)];
            title(title_name);
        hold off

        out_eps = [title_name, '.eps'];
        print(out_eps, '-depsc')
        close(h);
    end
end

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficient Plotting %
    %%%%%%%%%%%%%%%%%%%%%%%%
    x = 0:N_CCF;

    for n_hermite = 0:N_hermite

        h = figure;
        hold on
            plot(x, coeff(n_hermite+1, :), 'r')
            plot(x, coeff_noise(n_hermite+1, :), 'rs--', 'MarkerSize',6, 'MarkerFaceColor', 'r')
            plot(x, coeff_rvc(n_hermite+1, :), 'b')
            plot(x, coeff_noise_rvc(n_hermite+1, :), 'bo--', 'MarkerSize',6, 'MarkerFaceColor', 'b')  
            legend('Rest frame', 'Rest frame w n', 'Shifted frame', 'Shifted frame w n', 'Location', 'Best')
            xlabel('Observation Number','Interpreter','latex')
            Y_label = ['a', num2str(n_hermite)];
            ylabel(Y_label,'Interpreter','latex')
            TITLE1 = ['SN', num2str(SN)];
            title_name = ['Order', num2str(n_hermite), '--', TITLE1];
            title(title_name)
        hold off

        out_eps = [title_name, '_nor.eps'];
        print(out_eps, '-depsc')
        close(h);
    end


    %%%%%%%%%%%%%%%%%%%%%
    % Correlation Plots %
    %%%%%%%%%%%%%%%%%%%%%

    for n_hermite = 0:N_hermite

        h = figure;
        hold on 
            plot(RV_gauss * 1000, coeff(n_hermite+1, :), 'r')
            plot(RV_gauss * 1000, coeff_rvc(n_hermite+1, :), 'b')
            plot(RV_gauss * 1000, coeff_noise(n_hermite+1, :), 'rs', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
            plot(RV_gauss * 1000, coeff_noise_rvc(n_hermite+1, :), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b')
            xlabel('RV_a (m/s)')
            ylabel(['a', num2str(n_hermite)])
            legend('Ideal', 'Observed Frame, S/N = 10000', 'Location', 'northeast')
            legend('Rest frame', 'Observed Frame', 'Location', 'southeast')
            title_name = ['Order', num2str(n_hermite), ' -- ', TITLE1];
            title(['a', num2str(n_hermite)])
        hold off

        out_eps = [title_name, '-nor.png'];
        print(out_eps, '-dpng')
        close(h);
    end   
end



