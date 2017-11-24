function Ctemp = diffALL(i,j,k,D)
global C S
Ctemp=C(i,j,k)+D.*(S(i+1,j,k).*C(i+1,j,k)+S(i-1,j,k).*C(i-1,j,k)... 
                     +S(i,j+1,k).*C(i,j+1,k)+S(i,j-1,k).*C(i,j-1,k) ...
                     +S(i,j,k+1).*C(i,j,k+1)+S(i,j,k-1).*C(i,j,k-1) ...
                     -(S(i+1,j,k).*C(i,j,k)+S(i-1,j,k).*C(i,j,k)...
                      +S(i,j+1,k).*C(i,j,k)+S(i,j-1,k).*C(i,j,k) ...
                      +S(i,j,k-1).*C(i,j,k)+S(i,j,k+1).*C(i,j,k)));
end