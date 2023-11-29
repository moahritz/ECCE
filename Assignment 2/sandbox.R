T_end = 1 #samo
dt = 1/252 
NSim = 10000

r = 0.01
k = 2.80
theta = 0.06
eta = 0.50
rho = -0.70

S0 = 1
V0 = theta

W_X = sqrt(T_end) * rnorm(NSim)
W_V = rho * W_1 + sqrt(1 - rho^2) * sqrt(T_end) * rnorm(NSim)

S_h = numeric(T_end/dt)
S_h[1] = 1

for (j in 1:NSim){
  for (i in 2:nsteps){
    
    S_h[i] = S_h[i-1] + r * S_h[i-1] *dt + sigma * S_h[i-1] * W_X[i-1]
    
    V[i] = V[i-1] + k * (theta - V[i-1]) * dt + sigma * sqrt(V[i-1]) * W_V[i-1] 
  }
  
}