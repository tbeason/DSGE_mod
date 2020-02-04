var V $V$
    y $y$
    c $c$
    k $k$
    invest $i$
    z $z$
    s $s$
    E_t_SDF_plus_1 ${E_t(SDF_{t+1}}$
    E_t_R_k ${E_t(R^k_{t+1})}$
    R_f ${R^f}$
    q $q$
    d $d$
    w $w$
    p_e $p_e$
    cf_cons $dC$
    ;


varexo e $\varepsilon$
    ;

parameters 
    beta $\beta$
    gamma $\gamma$
    delta $\delta$
    alpha $\alpha$
    psi $\psi$
    a $a$
    b $b$
    xi $\xi$
    mu $\mu$
    sigma ${\sigma}$
    ;

beta = 0.998;
gamma = 5;
psi = 1.5;
delta = 0.025;
alpha = 0.36;
mu = 0.004;
sigma = 0.04;
xi = 3;
a = (exp(mu)-1+delta)^(1/xi);
b = (exp(mu)-1+delta)/(1-xi);

model;
#theta = (1 - gamma)/(1 - (1/psi));
// Define Value function
V = ((1-beta)*c^((1-gamma)/theta) + beta*s^(1/theta))^(theta/(1-gamma));

// Define an auxiliary variable s that captures E_t[V(+1)^sigma]
s = (exp(z(+1))*V(+1))^(1-gamma);


// Euler equation
1 = beta*(exp(z(+1))*c(+1)/c)^(-1/psi) * ((exp(z(+1))*V(+1))/(s^(1/(1-gamma))))^(1/psi-gamma) * (a*(invest/k(-1))^(-1/xi))*(((alpha-1)*y(+1)+c(+1))/k + ((a/(1-1/xi)*(invest(+1)/k)^(1-1/xi)+b) + 1 - delta)/(a*(invest(+1)/k)^(-1/xi)));

//define net return to capital
E_t_R_k = (a*(invest/k(-1))^(-1/xi))*(((alpha-1)*y(+1)+c(+1))/k + ((a/(1-1/xi)*(invest(+1)/k)^(1-1/xi)+b) + 1 - delta)/(a*(invest(+1)/k)^(-1/xi)));

//define expected value of stochastic discount factor
E_t_SDF_plus_1=beta*(exp(z(+1))*c(+1)/c)^(-1/psi) * ((exp(z(+1))*V(+1))/(s^(1/(1-gamma))))^(1/psi-gamma);
//E_t_SDF_plus_1=beta*exp(z(+1))^(-gamma)*(c(+1)/c)^(-1/psi) * (V(+1)/(s^(1/(1-gamma))))^(1/psi-gamma);

//define net risk-free rate
R_f=(1/E_t_SDF_plus_1-1);
    

//Budget constraint
//c +invest = y;

// Law of motion of capital
k = exp(-z(+1))*((1-delta)*k(-1) + k(-1)*(a/(1-1/xi)*(invest/k(-1))^(1-1/xi)+b));

// Technology shock
z = mu + sigma*e;

// Output definition
y = k(-1)^(alpha);

q=1/(a*(invest/k(-1))^(-1/xi));
d = alpha*y-invest;
c + p_e = w + d + p_e;
p_e = beta*exp(z(+1))^(-1/psi)*(c(+1)/c)^(-1/psi) * (exp(z(+1))*V(+1)/(s^(1/(1-gamma))))^(1/psi-gamma) * exp(z(+1))* (d(+1) + p_e(+1));
w = (1-alpha)*y;
cf_cons = log(c) - log(c(-1)) + z;
end;





steady_state_model;
    k = ((1/beta - 1 + delta)/alpha)^(1/(alpha-1));
    invest = (exp(mu)-1+delta)*k;
    y=k^alpha;
    c=y-invest;
    z=mu;
    V=((1-beta)/(1-beta*exp(z*(1-1/psi))))^(1/(1-1/psi))*c;
    s=(exp(z)*V)^(1-gamma);
    E_t_SDF_plus_1=beta*exp(-gamma*z);
    R_f=(1/E_t_SDF_plus_1-1);
    E_t_R_k = ((alpha-1)*y+c)/k + (exp(mu)-1+delta) + 1 - delta;
    q=1;
    d=alpha*y-invest;
    w = (1-alpha)*y;
    p_e = d/(1-beta*exp((1-1/psi)*z));
    cf_cons = z;
end;

shocks;
var e; stderr 1;
end;

steady;//(nocheck);

stoch_simul(order=3,periods=100000,drop=1000,irf=0) c k y E_t_R_k R_f;

