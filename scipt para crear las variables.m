## Script auxiliar  del experimento
% Se crear la funcion F , Q, b, y a de dimension n.
%
%% Como el experimento hay que repetirno para distintos n, creo el vector con 
% los distintos taman;os y luego  hago las pruebas en el FOR
%
%% En resultados voy a guardarn: Largo, ErrAbs, ||X||, iteraciones, Tiempo , ||F(x)|| 

Ej        = [ 100 ; 200; 300; 400; 500; 1000; 2000 ; 3000; 4000; 5000]
largo     = size(Ej);
Resultado = zeros(largo ,6 );

for j = 1:largo
    Resultado(j,1) =  Ej(j);   
    
    n = Ej(j);      ## Defino dimension
   
    Q = zeros(n,n);      
    b = zeros(n,1);
    uno = ones(1, n);
    x =  linspace(1/n, 1/n , n);
    a =  linspace(sqrt( exp( sum(x) ) ), sqrt( exp( sum(x) ) ) , n);
    
    for i = 1:n
    
       if i < n
        Q(i,i+1) = ( (-1)^i )/ (3*i);
        Q(i+1,i) = ( (-1)^i )/ (3*i);
        Q(i,i)   = 2*i - 1;
        b(i)     = ((-1)^i )* i;
       else
        Q(i,i)   = 2*i - 1;
        b(i)     = ( (-1)^i )* i;
       endif
      
     endfor 
     
  F = @(x) (Q*x' - b) + uno'*e.^sum(x);
  
    
    
   [ XSol, Fxsol ,IterTOT ,Tiempo ] = NR2 (F,x,Q,a , 20 , 0.001 );
    
    SolReal = fsolve(F,x);
    Resultado(j,2) = norm(XSol - SolReal);
    Resultado(j,3) = norm(XSol);
    Resultado(j,4) = IterTOT;
    Resultado(j,5) = Tiempo;
    Resultado(j,6) = norm(Fxsol);
endfor


% Salida de resultados obtenidos luego del experimento
Resultado



%% Para armar grafica de resultados obtenidos
variables = Resultado(:,1);
tiempo    = Resultado(:,5);
NormaF    = Resultado(:,6);
logTiempo = log(tiempo);



[p,ErrorEstimado] = polyfit(variables,tiempo,2);  % Modelo la curva
Test              = polyval(p,variables,ErrorEstimado); % Estimo los datos


plot(variables , Test , '-' , variables , tiempo , '+' )
  ylabel( 'Tiempo' );
  xlabel( 'n');
   title( 'Tiempo en funcion de n');
    legend('Modelo polinomico','Datos Observados','Location','NorthWest');


   
   
% Escala Logaritmica   

[p2,ErrorEstimado2] = polyfit(variables,logTiempo,2);  % Modelo la curva
Test2              = polyval(p2,variables,ErrorEstimado); % Estimo los datos



plot(variables , Test2 , '-' , variables , logTiempo , '+' )
  ylabel( 'Log Tiempo ');
  xlabel( 'n');
   title( 'Tiempo en funcion de n en escala logaritmica');
    legend('Modelo polinomico','Datos Observados','Location','NorthWest');


   