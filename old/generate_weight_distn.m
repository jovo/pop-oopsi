rate=5;
x=1/3;  % (peak epsp size)/(threshold - rest) 
dt=0.01;
b=log(-log(1-rate*dt)/dt); % baseline firing rate
w=log(-log(exp(-exp(b)*dt)-x)/dt)-b;
