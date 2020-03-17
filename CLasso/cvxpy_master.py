#%%
import numpy as np
from scipy.linalg import block_diag
import cvxpy as cp
import pandas as pd
import datetime as dt
import itertools

#%%
# pen_generate
def pen_generate (b,a,N,p,K,kk):
    a_out_exp = np.transpose(np.array([a]*N),(1,0,2))  
    p_temp = (b - a_out_exp)**2
    p_norm = np.sqrt(np.apply_along_axis(np.sum,2,p_temp).T)
    
    ind = [int (i)-1 for i in set(list(np.linspace(1,K,K))).symmetric_difference(set([kk]))]
    
    pen = np.apply_along_axis(np.prod,1,p_norm[:,ind]).reshape(N,1)
    
    return pen

# criterion 

def criterion (a_old, a_new,b_old,b_new,tol):
    
    d = False
        
    a_cri = np.sum(np.absolute(a_old - a_new)) / (np.sum(np.absolute(a_old)) + 0.0001)
    b_cri = np.mean(np.absolute(b_old - b_new)) / (np.mean(np.absolute(b_old)) + 0.0001)
    
    if (a_cri < tol) and (b_cri < tol):
        d = True
        
    return d


# post_lasso

def extend_ind(group_ind, TT):
    c = []
    for i in list(group_ind):
        c.append( list( range((i-1)*TT, i*TT) ) )
    c = list( itertools.chain.from_iterable(c) )
    return c


def post_lasso(group_est,y,X,K,p,N,TT):
    
    a_out = np.zeros(K*p).reshape(K,p)
    
    for k in range(K):
        group_ind = np.where(group_est == k)[0] + 1  
        data_ind = extend_ind(group_ind, TT)
        yy = y[np.r_[data_ind]]
        XX = X[np.r_[data_ind]]
        
        a_out[k,] = (np.linalg.inv(XX.T @ XX) @ (XX.T @ yy)).T     
        
    return a_out

# summarize simualtion results

def group_coerce(group_est, a_out, group0, a0, N, N_frac, K, p):
    a_est = np.array([[0.0]*2]*3)
    group = np.array([None]* N)
    
    for k in range(K):
        a_out_k = np.array( list(a_out[k,:]) * 3 ).reshape((3,2))
        gg = np.argmin( np.sum( (a_out_k - a0)**2 , 1) )
        a_est[gg,:] = np.array( a_out[k,:] )
        group_id = list(np.where(group_est == k)[0])
        group[group_id] = gg
         
    if sum(group == None) > 0:
        print("Error: Some individuals are not classfied.")
        return -1
    
    # Correct ratio
    ratio = sum( group == group0 ) / N
    # weighted se
    weighted_se = sum(( a_est[:, 1] - a0[:, 1] )**2 * N_frac )
    
    return ratio, weighted_se, a_est, group

#%% Estimation function
def PLS_cvxpy(N,TT,y,X,K,lambda_,R,tol = 0.0001):
    
    p = X.shape[1] 

    beta0 = np.zeros(N*p).reshape(N,p)

    for i in range (1,N+1):      
        ind = []
        ind.append((i-1)*TT)     
        ind.append(i*TT)            
        yy = y[ind[0]:ind[1]]    
        XX = X[ind[0]:ind[1]]
        beta0[i-1,] = (np.linalg.inv(XX.T @ XX) @ (XX.T @ yy)).T  

    b_out = np.array([beta0]*K)
    a_out = np.zeros(K*p).reshape(K,p)

    b_old = np.ones(N*p).reshape(N,p)
    a_old = np.ones(p).reshape(1,p)   
    

    for r in range (1,R+1): 
    
        for k in range (1,K+1):
        
            gamma = pen_generate(b_out,a_out,N,p,K,k)
        
            X_list = []
            for i in range (1,N+1):
                ind = []
                ind.append((i-1)*TT)      
                ind.append(i*TT)
                id_ = []
                id_.append((i-1)*p)
                id_.append(i*p)
                X_list.append(X[ind[0]:ind[1]])
            
            b = cp.Variable((p, N))
            a = cp.Variable((p,1))
            A = np.ones(N).reshape(1,N)
            obj1 = cp.norm(b - (a @ A),2,axis = 0) @ gamma     
        
            XX = block_diag(*X_list)
        
            obj = cp.Minimize(cp.sum_squares(y-cp.reshape((XX @ cp.vec(b)),(N*TT,1)))/(N * TT) + obj1 * (lambda_/N))
            prob = cp.Problem(obj)
            try:
                prob.solve(solver=cp.MOSEK)
            except:
                prob.solve(solver=cp.ECOS)
            
            a_out[k-1] = a.value.reshape(1,p)
            b_out[k-1] = b.value.T        
   
        a_new = np.copy(a_out[K-1])
        b_new = np.copy(b_out[K-1])
            
        if (criterion (a_old, a_new,b_old,b_new,tol) == True):
            break
    
        a_old = np.copy(a_out[K-1])
        b_old = np.copy(b_out[K-1])

    a_out_exp = np.transpose(np.array([a_out]*N),(1,0,2)) 
    d_temp = (b_out - a_out_exp)**2
    dist = np.sqrt(np.apply_along_axis(np.sum,2,d_temp).T)
    group_est = np.apply_along_axis(np.argmin,1,dist)   

    count = np.unique(group_est,return_counts = True)   
    if np.count_nonzero(count[1] > p) == K:
        a_out = post_lasso(group_est,y,X,K,p,N,TT)

    b_est = np.empty(N*p).reshape(N,p)
    for i in range (1,N+1):
        group = group_est[i-1]
        b_est[i-1,] = a_out[group ,]

    return b_est, a_out, group_est

#%%
# The main part

a0 = np.array([[0.4, 1.6], [1,1],[1.6, 0.4]])
p = 2
N_cut = np.array([0.3, 0.6, 1])
N_frac = np.array([0.3, 0.3, 0.4])


Rep = 30
MaxIter = 500
tol = 1e-4
K = 3

p = 2

correct_ratio = np.array([[None] * Rep] * 6)
se_record = np.array([[None] * Rep] * 6)
time_record = np.array([[None] * Rep] * 6)

case = -1
for N in [100, 200]:
    
    N_group_num = N*N_frac
    group0 = np.repeat( np.array(range(K)), N_group_num.astype(int))
    
    for T in [15, 25, 50]:
        
        case = case + 1
        
        for r in range(Rep):
            
            print([N,T,r])
            
            file_name = "./simu_data/pls_data_N_" + str(N) + "_T_" + str(T) + "_r_" + str(r+1) + ".csv"
            data = pd.read_csv(file_name)
            X_data = np.array(data[["V2","V3"]])
            y_data = np.array(data[["V1"]])
            
            # np.var returns the variance without ajustment of the degree of freedom
            # For comparison to R
            lambda_NT = 0.5 * sum( y_data**2 )/(y_data.size - 1) / (T**(1/3))
            
            t_start = dt.datetime.now()
            b_est, a_out, group_est = PLS_cvxpy(N = N, TT = T, y=y_data, X = X_data, K=K, lambda_ = lambda_NT, R = MaxIter, tol = tol)
            t_end = dt.datetime.now() - t_start
            
            time_record[case, r] = t_end.seconds + t_end.microseconds / 1e6
            
            ratio, se, a_est, group = group_coerce(group_est, a_out, group0, a0, N, N_frac, K, p)
            
            correct_ratio[case, r] = ratio
            se_record[case, r] = se
            
# %%
rmse = np.sqrt( np.mean(se_record, 1).astype(float) )
cr = np.mean(correct_ratio, 1)
time_sum = np.sum(time_record, 1)

result = pd.DataFrame(data = {'RMSE': rmse, 'Correct_Ratio': cr, 'Time': time_sum})
result.to_csv('Python_result.csv')
