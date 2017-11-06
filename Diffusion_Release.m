%% Parameters for diffusion of glutamate

%first.m, for initializing;
%synapse_size = 600; % 200, 400, 600
%cleft_geo='c'; % narrow near edge
%pos=9; %5 for center release and 9 for edge release
%R='i';
%Dglut2=0.1;% High affinity near edge \n - zone2
%Dglut3=0.4;% high affinity at the center? \n - zone1
%Dglut2=0.4;% Dglut value for outer region of cleft? \n - zone2
%tin=0;
%tf=10000;




%% Size of Synapse
synapse_size = input('Size of synapse, 200 or 600 nm? \n');

if synapse_size == 600
    i = 60; j = 101; k = 101;
    cleft.i = 29:31; cleft.j = 21:81; cleft.k = 21:81;
    zone.i = 29:31; zone.j = 34:67; zone.k = 34:67;
    %Receptor Locations
    R_Loc = zeros(16,2);
    R_Loc(1,:)=[26,76];R_Loc(2,:)=[42,76];R_Loc(3,:)=[60,76];R_Loc(4,:)=[76,76];
    R_Loc(5,:)=[26,60];R_Loc(6,:)=[42,60];R_Loc(7,:)=[60,60];R_Loc(8,:)=[76,60];
    R_Loc(9,:)=[26,42];R_Loc(10,:)=[42,42];R_Loc(11,:)=[60,42];R_Loc(12,:)=[76,42];
    R_Loc(13,:)=[26,26];R_Loc(14,:)=[42,26];R_Loc(15,:)=[60,26];R_Loc(16,:)=[76,26];
    
elseif synapse_size == 200
    i = 60; j = 61; k = 61;
    cleft.i = 29:31; cleft.j = 21:41; cleft.k = 21:41;
    zone.i = 29:31; zone.j = 26:36; zone.k = 26:36;
    %Receptor Locations
    R_Loc = zeros(16,2);
    R_Loc(1,:)=[23,39];R_Loc(2,:)=[28,39];R_Loc(3,:)=[34,39];R_Loc(4,:)=[39,39];
    R_Loc(5,:)=[23,34];R_Loc(6,:)=[28,34];R_Loc(7,:)=[34,34];R_Loc(8,:)=[39,34];
    R_Loc(9,:)=[23,28];R_Loc(10,:)=[28,28];R_Loc(11,:)=[34,28];R_Loc(12,:)=[39,28];
    R_Loc(13,:)=[23,23];R_Loc(14,:)=[28,23];R_Loc(15,:)=[34,23];R_Loc(16,:)=[39,23];
    
else
    disp('Size must either be 200 or 600 nm. \n');
end
%%

%% Geometry of cleft, reduced edge or middle?
cleft_geo = input('Reduce the cleft on outer edge (e), center (c), or flat (f)?\n','s');
%%

%% Release location
pos = input('Release over which receptor, numeric value only. \n');
if (pos == 6 && cleft_geo == 'c') || (pos == 16 && cleft_geo == 'e')
    Release_loc(1) = cleft.i(1)+1; 
    Release_loc(2) = R_Loc(pos,1); 
    Release_loc(3) = R_Loc(pos,2);
else
    Release_loc(1) = cleft.i(1); 
    Release_loc(2) = R_Loc(pos,1); 
    Release_loc(3) = R_Loc(pos,2);
end
            
%%

%Glutamate diffusion with multiple values within the synaptic cleft. 
%Cleft Glutamate Concentration Simulations

OK = 1;
global S C

R = input('Instantaneous, I, or Vesicle release, V?\n','s');
switch lower(R)
    case 'i'
        initial = 4000;
    case 'v'
        
        Release = vect_sum(Exit_2nm,20);
        %Release = vect_sum(Exit_10nm,20); 
        %Release = vect_sum(Exit_2nm_3,20); 
        Release = Release'; %Make Release a column vector
        initial =Release(1);
    otherwise
        OK = 0;
        disp('Enter V or I');
end

%% Diffusion coefficients
%Set Diffusion constants, Dglut1 is outside the cleft and Dglut2 is inside

%the cleft. dt = 0.2 microseconds. 
Dglut1=0.75;    %Outside the cleft
Dglut2 = input('Dglut value for outer region of cleft? \n');
%Dglut2=0.40;    %Outer region of cleft
Dglut3 = input('Dglut value for inner region of cleft? \n');
%Dglut3=0.2;    %Innermsot region of cleft

dx=0.01;        %0.01 micrometer or 10 nm  
dt=2e-5;        %2e-5 ms or 0.02 microseconds

D1=(Dglut1*dt)/(dx)^2;
D2=(Dglut2*dt)/(dx)^2;
D3=(Dglut3*dt)/(dx)^2; % for stability D1, D2, & D3 must be less than 0.125 (1/2^3) for CFL condition

Ddiff12=Dglut1/Dglut2;
Ddiff23=Dglut2/Dglut3;

Ndt=50000;
av=10;
inc=2e-3; %ms

%%
tin=input('Initial Timestep \n');
tf=input('Final Timestep \n');
t = tin;

if tin >= tf
    disp('Initial timestep must be less than final');
    OK = 0;
elseif tin > 1 %add additional space to existing variables
    Ctemp = C; Stemp = S;
    %global S C
    C = Ctemp; S = Stemp;
    Cleft_Conc=cat(1,Cleft_Conc,zeros(tf-tin,1));
    Cleft_Zone=cat(1,Cleft_Zone,zeros(tf-tin,1));
    Cconc=cat(1,Cconc,zeros(tf-tin,1));
    R_dual=cat(1,R_dual,zeros(tf-tin,16));
else %start from scratch
    %global S C
    %pos = input('Release over which receptor, numeric value only \n');
    S = Smat(i,j,k,cleft_geo,cleft,zone);
    C = Cmat(i,j,k,initial,Release_loc);
    R_dual=zeros(tf,16);
    [i_n j_n k_n]=size(C);
    Ctemp=zeros(i_n,j_n,k_n);
    Cconc=zeros(tf+2,1);
    Cconc(2)=sum(sum(sum(C)));
    Cleft_Conc=zeros(tf+2,1);
    Cleft_Conc(2)=sum(sum(sum(C(29:31,21:81,21:81))));
    Cleft_Zone=zeros(tf+2,1);
    Cleft_Zone(2)=sum(sum(sum(C(29:31,34:67,34:67))));
end

tic;

while OK == 1
    
%Start with the exterior
i=[2:i_n-1]; %same value for all for whole exterior region

%Region 1 and 2
k=[2:k_n-1]; j=[2:20,82:j_n-1]; 
    Ctemp(i,j,k) = diffALL(i,j,k,D1);
%Region 3 and 4
k=[2:20,82:k_n-1]; j=[21:81];
    Ctemp(i,j,k) = diffALL(i,j,k,D1);
    
%Interior, we first start with the 'outer' region of the cleft
%Boundary between exteriors and interior, the corners are dealt with
%separately.
i=29:31;    %i is constant for this region
j=21;  k=[22:80];
    Ctemp(i,j,k) = diffW(i,j,k,D2,Ddiff12); %West Boundary
j=81; 
    Ctemp(i,j,k) = diffE(i,j,k,D2,Ddiff12); %East Boundary
j=[22:80]; k=21;  
    Ctemp(i,j,k) = diffS(i,j,k,D2,Ddiff12); %South Boundary
k=81;
    Ctemp(i,j,k) = diffN(i,j,k,D2,Ddiff12); %North Boundary

%Corners with free space and cleft
j=21; k=21; 
    Ctemp(i,j,k) = diffSW(i,j,k,D2,Ddiff12); %SW Corner
j=81; k=21; 
    Ctemp(i,j,k) = diffSE(i,j,k,D2,Ddiff12); %SE Corner
j=21; k=81; 
    Ctemp(i,j,k) = diffNW(i,j,k,D2,Ddiff12); %NW Corner
j=81; k=81; 
    Ctemp(i,j,k) = diffNE(i,j,k,D2,Ddiff12); %NE Corner

%Outer area in cleft
%1st and 2nd outer region
j=[22:33,68:80]; k=[22:80];
    Ctemp(i,j,k) = diffALL(i,j,k,D2);
%3rd and 4th outer region    
j=[34:67];       k=[22:33,68:80];
    Ctemp(i,j,k) = diffALL(i,j,k,D2);

%Inner Area of cleft
j=[35:66]; k=[35:66];
    Ctemp(i,j,k) = diffALL(i,j,k,D3);
    
%Boundary between inner and outer area of synaptic cleft
j=34; k=[35:66];
    Ctemp(i,j,k) = diffW(i,j,k,D3,Ddiff23); %West boundary
j=67;  
    Ctemp(i,j,k) = diffE(i,j,k,D3,Ddiff23); %East boundary
j=35:66; k=34;  
    Ctemp(i,j,k) = diffS(i,j,k,D3,Ddiff23); %South boundary
k=67; 
    Ctemp(i,j,k) = diffN(i,j,k,D3,Ddiff23); %North boundary
    
%Now consider the four corners of the inner/outer cleft
j=34; k=34; 
    Ctemp(i,j,k) = diffSW(i,j,k,D3,Ddiff23); %SW Corner
j=67; k=34; 
    Ctemp(i,j,k) = diffSE(i,j,k,D3,Ddiff23); %SE Corner
j=34; k=67; 
    Ctemp(i,j,k) = diffNW(i,j,k,D3,Ddiff23); %NW Corner
j=67; k=67; 
    Ctemp(i,j,k) = diffNE(i,j,k,D3,Ddiff23); %NE Corner

%Updating C with Ctemp
C=Ctemp;

R_dual(t+1,1)=C(31,26,75);R_dual(t+1,2)=C(31,42,75);R_dual(t+1,3)=C(31,59,75);
R_dual(t+1,4)=C(31,75,75);R_dual(t+1,5)=C(31,26,59);R_dual(t+1,6)=C(31,42,59);
R_dual(t+1,7)=C(31,59,59);R_dual(t+1,8)=C(31,75,59);R_dual(t+1,9)=C(31,26,42);
R_dual(t+1,10)=C(31,42,42);R_dual(t+1,11)=C(31,59,42);R_dual(t+1,12)=C(31,75,42);
R_dual(t+1,13)=C(31,26,26);R_dual(t+1,14)=C(31,42,26);R_dual(t+1,15)=C(31,59,26);
R_dual(t+1,16)=C(31,75,26);


%Cconc will be a function of t, that gives the total concentration of the
%3-dimensional array to check for conservation of glutamat.
Cconc(t+2)=sum(sum(sum(C)));
Cleft_Conc(t+2)=sum(sum(sum(C(29:31,21:81,21:81))));
Cleft_Zone(t+2)=sum(sum(sum(C(29:31,34:67,34:67))));

if mod(t,10000) == 0
    disp('Total Concentration is')
    disp(Cconc(t+2))
    disp('Current Cleft Concentration is')
    disp(Cleft_Conc(t+2))
    disp('At timestep (dt)')
    disp(t)
    disp('R6 Concentration is')
    disp(R_dual(t+1,6));
    etime=toc/60;
    disp('Elapsed time (minutes) is:')
    disp(etime)
    disp('.'),disp('.')
end

%Add release to cleft
t=t+1;
if initial ~= 4000
    if pos == 6 && t <= size(Release,1)  %%#ok<AND2>
        C(29,42,59) = C(29,42,59) + Release(t);
    elseif pos == 16 && t <= size(Release,1)
        C(30,75,26) = C(30,75,26) + Release(t);
    end
end


if t == tf
    OK = 0;
end

end
conctime=toc/60;

%Add release to cleft
t=t+1;
if initial ~= 4000;
    if pos == 6 && t <= size(Release,1)  %%#ok<AND2>
        C(29,42,59) = C(29,42,59) + Release(t);
    elseif pos == 16 && t <= size(Release,1)
        C(29,75,26) = C(29,75,26) + Release(t);
    end
end


if t == tf
    OK = 0;
end
conctime=toc/60;
toc

% Total concentration plot
plot(Cconc)
set(gca,'XTickLabel',[0 0.02 0.04 0.06 0.08 0.1 0.12]);
% Concetrations on each receptor
figure
subplot(4,4,1);
plot(R_dual(:,1));
title('R1');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,2);
plot(R_dual(:,2));
title('R2');
set(gca,'XTickLabel',[0 0.05 0.1])

subplot(4,4,3);
plot(R_dual(:,3));
title('R3');
set(gca,'XTickLabel',[0 0.05 0.1])

subplot(4,4,4);
plot(R_dual(:,4));
title('R4');
set(gca,'XTickLabel',[0 0.05 0.1])

subplot(4,4,5);
plot(R_dual(:,5));
title('R5');
set(gca,'XTickLabel',[0 0.05 0.1])

subplot(4,4,6);
plot(R_dual(:,6));
title('R6');
set(gca,'XTickLabel',[0 0.05 0.1])

subplot(4,4,7);
plot(R_dual(:,7));
title('R7');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,8);
plot(R_dual(:,8));
title('R8');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,9);
plot(R_dual(:,9));
title('R9');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,10);
plot(R_dual(:,10));
title('R10');
set(gca,'XTickLabel',[0 0.1]);

subplot(4,4,11);
plot(R_dual(:,11));
title('R10');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,12);
plot(R_dual(:,12));
title('R12');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,13);
plot(R_dual(:,13));
title('R13');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,14);
plot(R_dual(:,14));
title('R14');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,15);
plot(R_dual(:,15));
title('R15');
set(gca,'XTickLabel',[0 0.05 0.1]);

subplot(4,4,16);
plot(R_dual(:,16));
title('R16');
set(gca,'XTickLabel',[0 0.05 0.1]);


