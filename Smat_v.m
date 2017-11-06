%Creating S Matrix
%These are the correct dimensions for the 
%paper and scaled down so that dx=10 nm as in the Nielsen paper. 
function S_out = Smat_v(i,j,k,cleft_geo,cleft,zone)

%cleft is defined outside;
%cleft = evalin('base', 'cleft');

%First make entire domain ones, then we make geometry of zeros for synapse.
S_out = ones(i,j,k);
        
%Top half of synapse
x=1:cleft.i(1)-1; y=cleft.j(1)+10:cleft.j(end)-10; z=cleft.k(1)+10:cleft.k(end)-10;
x_n=size(x,2); y_n=size(y,2); z_n=size(z,2);
S_out(x,y,z)=zeros(x_n,y_n,z_n);
        
%Bottom half of synapse
x=cleft.i(end)+1:i; x_n=size(x,2);
S_out(x,y,z)=zeros(x_n,y_n,z_n);
        
if cleft_geo == 'c'
    %Add more zeros here for the decrease in the center of the synapse
    x=cleft.i(1); y=zone.j; z=zone.k;
    y_n=size(y,2); z_n=size(z,2);
    S_out(x,y,z) = zeros(y_n,z_n); 
elseif cleft_geo == 'e'
    %This case reduces the edge of the synapse
    %First add another row of zeros for entire top half of cleft
    x=cleft.i(1); y=cleft.j; z=cleft.k;
    y_n=size(y,2); z_n=size(z,2);
    S_out(x,y,z)=zeros(y_n,z_n);
    
    %Now remove interior of synapse by changing from zeros to ones.
    y=zone.j; z=zone.k;
    y_n=size(y,2); z_n=size(z,2);
    S_out(x,y,z)=ones(y_n,z_n);
end %end of if

end %end of function