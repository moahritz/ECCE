### Exercise 3 ###
#1.
# params:
S0 = 100
mu = 0.1
sigma = 0.2
T_end = 1 #aka a year
dt = 1/252 # daily frequency, saw this number in class (is it always 252 trading days?)
r = 0.05

#European Option Pricing (closed form solution)
K = 110

d2 = (log(S0/K) + (r- 0.5 * sigma^2) * T_end) / (sigma * sqrt(T_end))
d1 = d2 + sigma * sqrt(T_end)

C = S0 * pnorm(d1) - K * exp(-r * T_end) * pnorm(d2)


#Monte Carlo Simulation
set.seed(1234)
NSim = c(1000, 10000, 100000, 50000)
WT = list()
ST = list()
MC_pay = list()
C_MC = list()

#check again if the way sqrt(T_end) enters into the formula is actually correct, or does it have to be smt else. Concerning T_end jsut being 1, but daily frequency means in a year it's actually 252. This is confusing..

#Definitely go back in the script to actually check for the proper definition of T, although everything makes sense, and it seems to be the same as in Luca's Matlab code....=> T_End = 1 is actually FINE, because 1 = 252 / 252 !!!!

for (i in 1:(length(NSim)-1)){
  WT[[i]] = rnorm(NSim[i])
  ST[[i]] = S0 * exp((r - 0.5 * sigma^2) * T_end + sigma * sqrt(T_end) * WT[[i]])
  
  MC_pay[[i]] = exp(-r * T_end) * pmax(ST[[i]]- K, 0) 
  
  
  a = mean(MC_pay[[i]])
  b = sd(MC_pay[[i]])
  c = a - qnorm(0.975) * (b/sqrt(NSim[i]))
  d = a + qnorm(0.975) * (b/sqrt(NSim[i]))
  
  cat(NSim[i], "simulations \nThe mean: ", a,", sd: ", b," and 95% CI:  [", c, ",", d, "]
      \n")
}
cat("The closed form price is", C,"\n")


#2.
#For the Asian Option we need the entire path and the take the mean of it!
#So we have a sequence of 252 days -> take the mean, put it into the max -> then we have one of the 50k simulations 

MC_A = c()
for (i in 1:50000){
  WT = rnorm(T_end/dt)
  ST = S0 * exp((r-0.5 * sigma^2) * T_end) + sigma * cumsum(sqrt(T_end) * WT)
  
  MC_A[i] = exp(-r * T_end) * max(mean(ST) - K, 0)
  
}


cat("For a Strike Price K =", K,": out of 50.000 Simulations, \nthe mean of the stock's price path was larger than K:\n", length(MC_A[MC_A > 0]),"time(s), i.e,",(length(MC_A[MC_A > 0])/50000),"% of the time, \nby amount(s) of (discounted):", MC_A[MC_A > 0],"\n
    The MC Estimation produces an average discounted payoff of:", mean(MC_A),"\n")

### Exercise 4 ###

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

# dont forget about Absorbtion: V(t) = max(V(t), 0)



#simulate N = 10 000 times XT using both the Heston and the Black-Scholes model (with σ = η)


W_X = rnorm(NSim)
W_V = rho * W_1 + sqrt(1 - rho^2) * rnorm(NSim)







