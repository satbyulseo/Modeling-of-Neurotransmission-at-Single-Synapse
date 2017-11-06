%ODE's for M-mode Kinectics
function dpdt = M_mode(t,p,C)
dpdt=zeros(7,1); %a column vector for the system
a = 1e-3; %if the unit of inc is ms, then a should be 1e-3, but if unit of inc is s, a should be 1;
dpdt(1,:) = a*(58*p(2) - 39*C*p(1));
dpdt(2,:) = a*((116*p(3) - 19*C*p(2)) - (58*p(2) - 39*C*p(1)));
dpdt(3,:) = a*(-(116*p(3) - 19*C*p(2)) + (173*p(4) - 150*p(3)));
dpdt(4,:) = a*((2412*p(5) - 902*p(4)) - (173*p(4) - 150*p(3)));
dpdt(5,:) = a*(-(2412*p(5) - 902*p(4)) + (1283*p(6) - 4467*p(5)));
dpdt(6,:) = a*((526*p(7) - 4630*p(6)) - (1283*p(6) - 4467*p(5)));
dpdt(7,:) = a*(-(526*p(7) - 4630*p(6)));
