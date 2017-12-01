function [Time,probM,totM]=NMDA_kinetics(C,Ndt,av,inc);

%,varargout
%We take 50,000 values from the diffusion simulation and average the values
%of [G] over 20 microsecond intervals. That is, every 10 entries, since
%dt=2 microseconds. Thus we are left with a vector with 5,000 entries
%spanning an interval of 10 ms. Since all the concenctrations are near
%zeros at this point, we use zeros for the remaining 50,000 entries of the
%concentration. Thus our total time interval is 100 ms.

SS = size(C,1); %Initial size of concentration vector 
%inc = 2e-5*av;
%inc=2e-3;
p0L=[1 0 0 0 0 0 0];p0M=[1 0 0 0 0 0 0];
%outM=p0M; outL=p0L;
tfinal = Ndt;
probL = zeros(tfinal,1); probM = zeros(tfinal,1);
totL = zeros(tfinal,1); totM = zeros(tfinal,1);
scale = (6.022e-4)^(-1);

%G is the concentration of Glutamate averaged over a period of 10 
%time steps. 
G=zeros(tfinal,1);
for i= 1:SS/av;        %only need to do for av(1), after that conc = 0
    G(i)=(sum(C(av(1)*i-(av(1)-1):av(1)*i))/av(1))*scale;
    %G(i)=sum(C(((av*i-av+1):av*i), 16))/av*scale;
end
%Multiply by 1.66e3 since [G] is in mM, and concentration units in ODE file
%is microMolar (uM)
time = [0 inc]; %initial interval for solving ODE, 
%Time=(1:1:Ndt);

Time=(1:1:100);

for i=1:tfinal
    %i=1:Ndt
    [T,Y] = ode23s(@M_mode,time,p0M,[],G(i));
    %outM=cat(1, outM, Y(size(Y,1),:));
    p0M=Y(size(Y,1),:);
    time(1) = time(2);
    time(2) = inc + time(2);
    probM(i)=p0M(6)+p0M(7);
    totM(i) = sum(p0M);
end

%time = [0 inc]; %initial interval for solving ODE, 

%for i=1:tfinal
    %i=1:10:tfinal
 %   [T,Y] = ode23s(@L_mode,time,p0L,[],G(i));
    %outL=cat(1, outL, Y(size(Y,1),:));
 %   p0L=Y(size(Y,1),:);
 %   time(1) = time(2);
 %   time(2) = inc + time(2);
 %   probL(i)=p0L(6)+p0L(7);
 %   totL(i) = sum(p0L);
%end


