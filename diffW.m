function Ctemp = diffW(i,j,k,D,Ddiff)
global S C
Ctemp = C(i,j,k)+D.*(S(i+1,j,k).*C(i+1,j,k)+S(i-1,j,k).*C(i-1,j,k)... 
                     +S(i,j+1,k).*C(i,j+1,k)+Ddiff.*S(i,j-1,k).*C(i,j-1,k) ...
                     +S(i,j,k+1).*C(i,j,k+1)+S(i,j,k-1).*C(i,j,k-1) ...
                     -(S(i+1,j,k).*C(i,j,k)+S(i-1,j,k).*C(i,j,k)...
                     +S(i,j+1,k).*C(i,j,k)+Ddiff.*S(i,j-1,k).*C(i,j,k) ...
                     +S(i,j,k+1).*C(i,j,k)+S(i,j,k-1).*C(i,j,k)));
end
