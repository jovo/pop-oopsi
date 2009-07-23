function v=lograte2voltage(z,rate,dt)
%converts log-rate units to voltage
v=sign(z).*(exp(-rate*dt)-exp(-exp(abs(z))*rate*dt));
