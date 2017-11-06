function Ctemp = diffBot(i,j,k,D,Ddiff)
global S C
Ctemp = C(i,j,k)+D.*(Ddiff.*S(i+1,j,k).*C(i+1,j,k)+S(i-1,j,k).*C(i-1,j,k)... 
                     +S(i,j+1,k).*C(i,j+1,k)+S(i,j-1,k).*C(i,j-1,k) ...
                     +S(i,j,k+1).*C(i,j,k+1)+S(i,j,k-1).*C(i,j,k-1) ...
                     -(Ddiff.*S(i+1,j,k).*C(i,j,k)+S(i-1,j,k).*C(i,j,k)...
                     +S(i,j+1,k).*C(i,j,k)+S(i,j-1,k).*C(i,j,k) ...
                     +S(i,j,k+1).*C(i,j,k)+S(i,j,k-1).*C(i,j,k)));
end
%This gives diffusion between different diffusion constants, where D is the
%original constant inputed and Ddiff is the scaling factor from D1 to D2 for
%instance. Here we have D as diffusion constant everywhere except in the
%south direction (i+1), where we have D*Ddiff = other constant as defined
%in the calling .m file.