clear all
clc

% Master file

% Compare CVX + Mosek to Rmosek
% Read the contains genrated data: 100 replications
% Focus on innerloop only with beta = [.9; .9]

%% DATA GENENRATION
% We included data simulated as csv files.
% One can simulate new set of data by uncomment the follow code block

% [y, X, Z] = dgpLinearIV(ones(2,1), 120*100, 80);
% Data_120_80 = [y, X, Z];

% [y, X, Z] = dgpLinearIV(ones(2,1), 120*100, 160);
% Data_120_160 = [y, X, Z];

% [y, X, Z] = dgpLinearIV(ones(2,1), 240*100, 80);
% Data_240_80 = [y, X, Z];

% [y, X, Z] = dgpLinearIV(ones(2,1), 240*100, 160);
% Data_240_160 = [y, X, Z];

% csvwrite('Data_120_80.csv', Data_120_80)
% csvwrite('Data_120_160.csv', Data_120_160)
% csvwrite('Data_240_80.csv', Data_240_80)
% csvwrite('Data_240_160.csv', Data_240_160)

%% Start

cvx_quiet(true)
b = 0.9 * ones(2,1);

L_Result = zeros(100, 4);
T_Result = zeros(100, 4);
c = 0;
for n = [120, 240]
    for m = [80, 160]
        c = c + 1; 
        
        file_name = strcat('Data_', num2str(n), '_', num2str(m), '.csv');
        D = csvread(file_name);
        
        tau = 0.5 * sqrt( log(m) / n );
        
        for r = 1:100
            disp([n,m,r]);
            
            id = ((r-1)*n+1):(r*n);
            y = D(id, 1);
            X = D(id, 2:3);
            Z = D(id, 4:end);
            tic
            [L, p] = REL_inner(b, y, X, Z, tau);
            t = toc;
            T_Result(r, c) = t;
            L_Result(r, c) = L;
        end
    end
end

Result = [L_Result, T_Result];
Time_Sum = sum(T_Result)';

save('REL_CVX_Result.mat', 'Result')
csvwrite('REL_time_compare_matlab.csv', Time_Sum)
