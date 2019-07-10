% The main file using PLS for China's provincial total GDP estimation
global p K_max
cvx_solver mosek

tol = 0.0001;
R = 80;

%% Import data: Norminal Case
data=csvread('Data_for_CLasso_0007.csv');
year=data(:,1);
code=data(:,2);
ngdp=data(:,3);   % log term of norminal provincial GDP
light=data(:,4);  % log term of night light
ntax=data(:,5);   % log term of national tax
expor=data(:,6);  % log term of export
impor=data(:,7);  % log term of import
ecos=data(:,8);   % log term of electricity consumption
rail=data(:,9);   % log term of railway cargo volumn

%% Choose variables

% 1)This is the all variable case
% 2)Drop light from X and re-run the code for case 2: no light
% 3)Drop light and ntax from X and re-run the code for case 3: no light&tax
X = [light,ntax,expor,impor,ecos];
y = ngdp;

p = size(X, 2);
T = max(year);
N = max(code);

K_max = 4;
lamb.grid = 10;

lamb.min  = 0.001;
lamb.max  = 0.01;
lamb_const = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) );
numlam = length(lamb_const);

index = dataset(code, year, y, X );
index.Properties.VarNames = {'N'  'T'  'y'  'X'};

y_raw = y;
X_raw = X;

%%
for i = 1:N
    yi = y(index.N == i);
    mean_yi = mean(yi);
    yi = bsxfun(@minus, yi, mean(yi) );
    y(index.N == i) = yi/std(yi, 1);
    y_raw(index.N==i) = y(index.N == i) + mean_yi;
    
    Xi = X(index.N == i, : );
    mean_Xi = mean(Xi);
    Xi = bsxfun(@minus, Xi, mean(Xi) );
    X(index.N == i, :) = Xi./repmat( std(Xi, 1), [T 1] ) ;
    X_raw(index.N == i, :) = X(index.N == i, :) + repmat( mean(Xi), [T 1]);
end

ds = dataset( code, year, y, X, y_raw, X_raw );
ds.Properties.VarNames = {'N'  'T'  'y'  'X' 'y_raw' 'X_raw'};
%% initial values
beta_hat0 = zeros(N, p);
for i = 1:N
    yi = ds.y(ds.N == i );
    Xi = ds.X(ds.N == i, : );
    beta_hat0(i,:) = regress( yi , Xi );
end

%% estimation
TT = T;
IC_total = ones(K_max, numlam );
Time = zeros(K_max, numlam);

for ll = 1:numlam
    

    a = ds.X \ ds.y; 
    bias = SPJ_PLS(T,ds.y_raw, ds.X_raw);
    a_corr = 2 * a - bias;
    IC_total(1, :) = mean( ( y - X*a_corr ).^2 );


    for K = 2:K_max
        
        disp([ll, K])
        
        Q = 999*zeros(K,1);

        lam = lamb_const(ll)*var(y) * T^(-1/3);
        tic
        [b_K, hat.a] = PLS_est(N, TT, y, X, beta_hat0, K, lam, R, tol); % estimation
        time = toc;
        Time(K, ll) = time;
        [~, H.b, ~, group] = report_b( b_K, hat.a, K );
        sum(group)            

        post_b = zeros(N, p);
        post_a = zeros(K, p);
        if K >=2
            for i = 1:K
                NN = 1:N;
                H.group = logical(group);
                this_group = group(:,i);
                if sum(this_group) > 0
                    g_index = NN(this_group);
                    g_data = ds( ismember(ds.N, g_index), : );

                    post = post_est_PLS_dynamic(T, g_data);

                    e = g_data.y - g_data.X * post.post_a_corr ;
                    Q(i) = sum( e.^2 );
                    post_b(this_group,:) = repmat(post.post_a_corr', [sum(this_group), 1] );
                end
            end
        end


        IC_total(K , ll) = sum(Q) / (N*T) ;

    end
end
%calculate the IC
pen = 2/3 * (N*T)^(-.5) * p .* repmat( (1:K_max)', [1 numlam]);
IC_final = log(IC_total) + pen;
disp(IC_final)


%% PLS estimation
[K,lamb_index]=find(IC_final==min(min(IC_final)));

if length(K)>1
    K=K(end);
    lamb_index=lamb_index(end);
end

lam = lamb_const(lamb_index) *var(y) * T^(-1/3);

[b_K, a] = PLS_est(N, T, y, X, beta_hat0, K, lam, R, tol);
[~, b, ~ , group] = report_b( b_K, a, K );