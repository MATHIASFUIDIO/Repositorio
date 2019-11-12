##Ejercicio 2.a PMCL

#función g1
x=[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];
g1=[3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.7,0.32,0.4,0.26,0.32,0.25];
X=log(x);
Y=log(g1);
Xcuad=X.^2;
A=[length(x),sum(X);sum(X),sum(Xcuad)];
B=[sum(Y),sum(X.*Y)];
Sol=linsolve(A,B');
c=e^Sol(1);
p=-Sol(2);

#función g2
g2=[15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.2];
Y2=log(g2);
B2=[sum(Y2),sum(X.*Y2)];
Sol2=linsolve(A,B2');
d=e^Sol2(1);
q=-Sol2(2);

##Ejercicio 2.a QR

[Q, R] = gsC(A);

#g1
SolQR =inv(R)*inv(Q)*B';
cQR  = exp(SolQR(1));
pQR = -SolQR(2);

#g2
Sol2QR =inv(R)*inv(Q)*B2';
dQR  = exp(Sol2QR(1));
qQR = -Sol2QR(2);


##Ejercicio 2.b PMCNL

x       =[0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.1; 1.2; 1.3; 1.4; 1.5; 1.6; 1.7; 1.8; 1.9; 2.0];
g1_real =[3.89 ;2.75 ;2.01 ;1.61 ;1.21 ;0.89 ;0.69 ;0.63 ;0.44 ;0.42 ;0.7  ;0.32 ;0.4  ;0.26 ;0.32 ;0.25];
g2_real =[15.96;9.45 ;5.75 ;3.82 ;2.89 ;2.17 ;1.22 ;1.05 ;0.86 ;0.63 ;0.69 ;0.40 ;0.44 ;0.29 ;0.43 ;0.2];


% Defino g1 y g2, las variables que queremos hayar son (c,p) y (d,q)
g1 = @(y , x) y(1)./( x.^(y(2)) )
g2 = @(y , x) y(1)./( x.^(y(2)) )


% Defino la Jacobiana
Jg1 = @(y,x) [1./( x.^(y(2)) ) , (y(1).*log(x) ) ./ ( -x.^(y(2)) )]
Jg2 = @(y,x) [1./( x.^(y(2)) ) , (-y(1).*log(x) ) ./ ( x.^(y(2)) )]

%%% g1
y1 = [2;2];           % initial guess
n = 0;
while 1
  b = g1(y1,x)- g1_real ;          % evaluate f
  A = Jg1(y1,x);         % evaluate Jacobian
  m = -A\b;          % solve linear least squares problem norm(A*d+b)=min
  y1 = y1 + m;
  n= n+1;
  if norm(m)<=1e-15  % stop iteration if norm(d)<=StepTolerance
    break
  end
end

%%% g1
y2 = [2;2];           % initial guess
n = 0;
while 1
  b = g2(y2,x)- g2_real ;          % evaluate f
  A = Jg2(y2,x);         % evaluate Jacobian
  m = -A\b;          % solve linear least squares problem norm(A*d+b)=min
  y2 = y2 + m;
  n= n+1;
  if norm(m)<=1e-15  % stop iteration if norm(d)<=StepTolerance
    break
  end
end

##Probamos las funciones obtenidas por PMCL y por PMCNL y comparamos cuál aproxima
##mejor los valores dados.

x=[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];
g1=[3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.7,0.32,0.4,0.26,0.32,0.25];
g2=[15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.2];

#PMCL

#g1
g1PMCL=[];
for i=1:length(x)
    g1PMCL=[g1PMCL,(0.9614*(x(i)^(-2)))];
endfor

#g2
g2PMCL=[];
for i=1:length(x)
    g2PMCL=[g2PMCL,(1.9598*(x(i)^(-3)))];
endfor

#PMCNL

#g1
g1PMCNL=[];
for i=1:length(x)
    g1PMCNL=[g1PMCNL,(0.9675*(x(i)^(-2)))];
endfor

#g2
g2PMCNL=[];
for i=1:length(x)
    g2PMCNL=[g2PMCNL,(1.9949*(x(i)^(-3)))];
endfor

#Errores para g1
errores_g1_PMCL=[];
errores_g1_PMCNL=[];
for i=1:length(x)
  errores_g1_PMCL=[errores_g1_PMCL,abs(g1(i)-g1PMCL(i))];
  errores_g1_PMCNL=[errores_g1_PMCNL,abs(g1(i)-g1PMCNL(i))];
endfor
comparacion_g1=[errores_g1_PMCL',errores_g1_PMCNL']

#Errores para g2
errores_g2_PMCL=[];
errores_g2_PMCNL=[];
for i=1:length(x)
  errores_g2_PMCL=[errores_g2_PMCL,abs(g2(i)-g2PMCL(i))];
  errores_g2_PMCNL=[errores_g2_PMCNL,abs(g2(i)-g2PMCNL(i))];
endfor
comparacion_g2=[errores_g2_PMCL',errores_g2_PMCNL']

#Ejercicio 3

x=[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];

##Forma analítica
Fy = @(x) (-0.90602*(exp(0.9614/x))) + (2.03849/x) + 2.12033
F_analit=[];
for i=1:length(x)
  F_analit=[F_analit,Fy(x(i))];
endfor

g1 = @(x) 0.9614*x^(-2)
g2 = @(x) 1.9598*x^(-3)

##Euler hacia adelante
function [X]=euler_pa_del(X0,n,a,b,g1,g2)
X=zeros(n,2);
h = abs(b-a)/(n-1);
X(1,:)=X0;

for i = 1:n-1;
  k =(g2(X(i,1)) - (g1(X(i,1))*X(i,2)));
  X(i+1,2)= X(i,2) + (h*k); 
  X(i+1,1)= X(i,1) + h;
endfor
return
endfunction

X0 = [0.5, 0];
XX = euler_pa_del(X0,length(x),0.5,2,g1,g2);


##Euler hacia atrás
function [X]=euler_pa_atr(X0,n,a,b,g1,g2)
X=zeros(n,2);
X(1,:)=X0;
h = abs(b-a)/(n-1);

for i=1:n-1;
  X(i+1,1)= X(i,1) + h;
endfor

for i = 1:n-1;
  X(i+1,2)=(X(i,2)+(h*g2(X(i+1,1))))/(1+(h*g1(X(i+1,1))));
endfor
return
endfunction

X0 = [0.5, 0];
XX2 = euler_pa_atr(X0,length(x),0.5,2,g1,g2);




##Runge-Kutta
Fprima =  @(t, y)  g2(t ) - ( g1( t )*y )
Rec = [0.5,2]
X0 = [ 0];

[t, y] = ode45(Fprima,Rec, X0);


hold on;
plot(x,F_analit);
plot(XX(:,1),XX(:,2));
plot(XX2(:,1),XX2(:,2));
plot(t,y);
title('Aproximación de "y" con 15 pasos');
xlabel('x');
ylabel('y');
legend({'y = Resultado analítio','y = Euler hacia adelante','y = Euler hacia atrás','y = Runge-Kutta'},'Location','southeast')
hold off;

#Los Euler con 200 pasos
#X0 = [0.5, 0];
#XX = euler_pa_del(X0,200,0.5,2,g1,g2);
#XX2 = euler_pa_atr(X0,200,0.5,2,g1,g2);
#hold on;
#plot(x,F_analit);
#plot(XX(:,1),XX(:,2));
#plot(XX2(:,1),XX2(:,2));
#plot(t,y);
#title('Aproximación de "y" con 200 pasos');
#xlabel('x');
#ylabel('y');
#legend({'y = Resultado analítio','y = Euler hacia adelante','y = Euler hacia atrás','y = Runge-Kutta'},'Location','southeast')
#hold off;


#Estudio de errores
#analítica evaluada en x, E adelante, E atrás, error adelante, error atrás
[F_analit',XX(:,2),XX2(:,2),abs(F_analit'-XX(:,2)),abs(F_analit'-XX2(:,2))]
Fy = @(x) (-0.90602*(exp(0.9614/x))) + (2.03849/x) + 2.12033
F_analit_rk=[];
for i=1:length(t)
  F_analit_rk=[F_analit_rk,Fy(t(i))];
endfor
#analítica evaluada en t, RK, error de RK
[F_analit_rk',y,abs(F_analit_rk'-y)]


#Ejercicio 4

%% INTERPOLACION LINEAL

xint = 0.5:0.01:2

vq1 = interp1(t , y , xint);


%% SPLINE CUBICO


vq2 = interp1(t,y,xint,'spline');



hold on;
plot(xint,vq2,':.');
plot(xint,vq1,':.' )
plot(t,y,'o')
plot(x, F_analit , 'x' );
  title('Interpolacion');
  xlabel('x');
  ylabel('y');
  legend({'Spline Cubico', 'Lineal a trozo', 'Runge Kutta', 'Estimacion analitica' },'Location','southeast')
hold off;


##Errores de interpolación
vq1 = interp1(t , y , x);
vq2 = interp1(t,y,x,'spline');
#analítica evaluada en x,interp trozos en x, interp cubico en x, error trozos,error cúbicos
[F_analit',vq1',vq2',abs(F_analit'-vq1'),abs(F_analit'-vq2')]

