%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
cd ccf_fits/
file_list   = dir('*.fits');
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);
MJD         = zeros(N_FILE, 1);

for n = 1:N_FILE
    FILE = file_name{n};
    spec    = fitsread(char(FILE));                                         % spectra of all orders
    H       = fitsinfo(char(FILE));                                         % fits information including header
    header  = H.PrimaryData.Keywords;                                       % get header only
    MJD(n)  = header{30, 2};
end

cd ../

ORDER           = 9;                                                       % Highest Hermite order 
array_order     = 0:ORDER;
idx_even        = mod(0:ORDER, 2) == 0;
order_even      = array_order(idx_even);
order_odd       = array_order(~idx_even);

SHIFT           = 11;
Her_Spe         = zeros(N_FILE, ORDER);
coeff           = zeros((ORDER+1), SHIFT);

V_BEGIN         = -40;                                                      % m/s                                                     
V_END           = -V_BEGIN;
V_GRID          = (V_BEGIN : (V_END - V_BEGIN)/(SHIFT-1) : V_END) / 1000;   

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

    parfor order = 1:ORDER+1
        for shift = 1:SHIFT
            temp                    = hermite_nor(order-1, v - V_GRID(shift)) * grid_size;
            coeff(order, shift)     = sum(A .* temp);         
        end
    end

    % fitting %
    for order = order_even
        p = polyfit(V_GRID * 1000, coeff(order+1, :), 2);
        Her_Spe(n, order+1) = p(1);
    end

    for order = order_odd
        p = polyfit(V_GRID * 1000, coeff(order+1, :), 1);
        Her_Spe(n, order+1) = p(1);
    end
    
    waitbar( n / N_FILE )
end
close(h)  

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

array = [MJD'; Her_Spe(:, 1+1)'];
fileID = fopen('Her_spe_1.txt','w');
fprintf(fileID,'%f %1.20f\n',array);
fclose(fileID);

% array = [MJD', Her_Spe(:, 1+1)'];

% VEL = importdata('GL479_vels.txt');
% fprintf(fileID, '%f %f \n', MJD', Her_Spe(:, 1+1)');
% fclose(fileID);

% dlmwrite('Her_spe_1.txt', [MJD, Her_Spe(:, 1+1)]);