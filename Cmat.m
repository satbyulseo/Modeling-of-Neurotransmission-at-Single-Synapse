%Creating  C Matrix
%With initial values of Glutamate at 4000 at release sites
%Build individual blocks natrices 'Cb'
%Matrix is constructed similar to S matrix, with no concentrations
%These are the correct dimensions for the paper and scaled down so that
%dx=10 nm as in the Nielsen paper. 
function C_out = Cmat(i,j,k,initial,Release_loc)
C_out = zeros(i,j,k);
x = Release_loc(1); 
y = Release_loc(2);
z = Release_loc(3);

C_out(x,y,z) = initial;


