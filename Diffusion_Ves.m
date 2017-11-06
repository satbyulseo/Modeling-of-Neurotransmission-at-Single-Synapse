%Glutamate diffusion with multiple values within the synaptic cleft. 
%Cleft Glutamate Concentration Simulations with the addition of a 40nm by
%40 nm vesicle at the release site. 

%2 nm fusion pore!!!!
tic
%create S and C matrices 
OK = 1;
s = input('Fusion neck, 10 nm or 2 nm? \n');
if s == 2
    p = [21 23];
elseif s == 10
    p = [17 27];
else
    OK = 0;
    disp('Either input 2 or 10.');
end

%Set Diffusion constants, Dglut1 is outside the cleft and Dglut2 is inside
%the cleft. dt = 0.2 microseconds. 
% raise the diffusion constant
%Dglut1=0.3; %Within vesicle
%Dglut2=0.1; %Fusion pore


%Dglut1=0.3; %Within vesicle
%Dglut2=0.3; %Fusion pore
Dglut1=0.15; %Within vesicle
%Dglut2=0.15; %Fusion pore
Dglut2=0.0375; %Fusion pore
%Dglut1=0.4; %Within vesicle
%Dglut2=0.4; %Fusion pore
dx=0.001; %0.001 micrometer or 1 nm  (since we have a smaller mesh)l
dt=1e-6; %1e-6 ms or 0.01 microseconds
D1=(Dglut1*dt)/(dx)^2;
D2=(Dglut2*dt)/(dx)^2;
Ddiff21=Dglut2/Dglut1;
tic;

tin=input('Initial Timestep \n');
tf=input('Final Timestep \n');
t = tin;

if tin >= tf
    disp('Initial timestep must be less than final');
    OK = 0;
elseif tin > 1
    Ctemp = C; Stemp = S;
    global S C
    C = Ctemp; S = Stemp;
    Exit=cat(1,Exit,zeros(tf-tin,1));
    Ves_Conc=cat(1,Ves_Conc,zeros(tf-tin,1));
    Ves_Conc(tin)=sum(sum(sum(C)));
else
    clear global S
    clear global C
    global S C
    S = Svesmat(p);
    C = Cvesmat;
    [i_n j_n k_n]=size(C);
    Ctemp = zeros(i_n,j_n,k_n);
    Exit=zeros(tf,1);
    Ves_Conc=zeros(tf+1,1);
end


while OK == 1
    
        %Diffusion within vesicle, excluding South boundary layer, since
        %it has exit with fusion pore that we deal with separately.
        i=2:41; j=2:42; k=2:42;
            Ctemp(i,j,k) = diffALL(i,j,k,D1);
            
        %Diffusion along Southern boundary of vesicle (exluding entry into
        %fusion pore.
        i=42; j=2:42; k=[2:p(1)-1,p(2)+1:42];
                Ctemp(i,j,k) = diffALL(i,j,k,D1);
         
        k=p(1):p(2); j=[2:p(1)-1,p(2)+1:42];
                Ctemp(i,j,k) = diffALL(i,j,k,D1);
            
        %Diffusion at entry of fusion pore to vesicle, i.e. North plate.
        i = 42; j=p(1):p(2); k=p(1):p(2);
                Ctemp(i,j,k) = diffBot(i,j,k,D1,Ddiff21);
        
        %Diffusion along fusion pore, except at endplates. 
        i=43:50;
            Ctemp(i,j,k) = diffALL(i,j,k,D2);
            
        
        %Here we don't allow for particles to flow back up since they
        %exit just below.
        i=51;
            Ctemp(i,j,k) = diff_noN(i,j,k,D2);
            
        %Diffusion at exit of fusion pore to presynaptic membrane, i.e.
        %Southern plate.
        i = 52;
            Ctemp(i,j,k) = diffSExit(i,j,k,D2);
        
        %Ves_Cconc(t+1)=sum(sum(sum(C)));
        %Exit(t)=sum(sum(sum(C(41,20:22,20:22))));


            
        C = Ctemp;
        %Exit(t+1)=sum(sum(sum(C(41,20:22,20:22))));
        Exit(t) = sum(sum(C(52,p(1):p(2),p(1):p(2))));
        Ccheck=C(50:53,p(1):p(2),22);
        ExitTot = sum(Exit);
        C(52,p(1):p(2),p(1):p(2)) = 0; %account for these molecules leaving the system.
        Ves_Conc(t+1) = sum(sum(sum(C)));
        
        
        if (Ves_Conc(t+1) <= 0.5 || t == tf)
            OK = 0;
        end
        if mod(t,10000) == 0 
            disp('Exit Total is')
            disp(ExitTot)
            disp('Current Vesicle Concentration is')
            disp(Ves_Conc(t+1))
            disp('At timestep (dt)')
            disp(t)
            etime=toc/60;
            disp('Elapsed time (minutes) is:')
            disp(etime)
            Ccheck;
        end
    t = t+1;
end
vestime=toc/60;        
