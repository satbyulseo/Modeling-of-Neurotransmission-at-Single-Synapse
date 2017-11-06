%Low group
function dpdt = L_mode(t,p,C)
% C is the concentration of glutamate at the given time step.
% First set of seven equations is for NMDA Receptors
dpdt=zeros(7,1); %a column vector for the system
a = 1e-3;  %if the unit of inc is ms, then a should be 1e-3, but if unit of inc is s, a should be 1;
dpdt(1,:) = a*(60*p(2) - 34*C*p(1));
dpdt(2,:) = a*((120*p(3) - 17*C*p(2)) - (60*p(2) - 34*C*p(1)));
dpdt(3,:) = a*(-(120*p(3) - 17*C*p(2)) + (161*p(4) - 127*p(3)));
dpdt(4,:) = a*((2610*p(5) - 580*p(4)) - (161*p(4) - 127*p(3)));
dpdt(5,:) = a*(-(2610*p(5) - 580*p(4)) + (2167*p(6) - 2508*p(5)));
dpdt(6,:) = a*((662*p(7) - 3449*p(6)) - (2167*p(6) - 2508*p(5)));
dpdt(7,:) = a*(-662*p(7) + 3449*p(6));


 