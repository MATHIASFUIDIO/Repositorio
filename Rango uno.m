##
% Actualizacion de Rango 1
% Entrda: L - Triangular Superior , V- Vector de actualizacion
% Salida: B - Triangular superior de actualizacion
%
% Funciones auxiliaries para calcular los valores de  la matriz B
##
function [B] = CholeskyModifyB(L, V)
n = length(V);
B = zeros(n,n);
Vnueva = zeros(1,n);
 for i = 1 : n
      Laux = L(i,i);
      Vaux = V(i);
    [Caux , Saux ,Laux] = Compute( Laux, Vaux);
      B(i,i) =Laux;
    for j = (i+1) : n
     Laux = L(i,j);
     Vaux = V(j);
     [Laux,Vaux] = Apply(Caux , Saux, Laux ,Vaux );
     B(i,j)= Laux;
     V(j)  = Vaux;
    endfor
  endfor  
endfunction


      


#// helper function: 
function[Caux , Saux , Laux] = Compute( Laux , Vaux )
  w     = sqrt( Laux^2 + Vaux^2); 
  Caux  = w/Laux;
  Saux  = Vaux/Laux;
  Laux  = w;
endfunction

#// helper function:

function [Laux, Vaux] =  Apply(Caux , Saux, Laux , Vaux )
  Laux  = ( Laux + Saux*Vaux ) / Caux;
  Vaux    = Caux*Vaux - Saux * Laux;
endfunction





