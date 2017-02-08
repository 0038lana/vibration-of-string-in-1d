l = 1;
n = 200;

h = l/n; %korak

x = 0:h:l;
A = zeros(n+1, n+1);
M = zeros(n+1, n+1);

for i=2:n
    A(i,i) = 1/(2*h) * (a(x(i-1)) + 2*a(x(i)) + a(x(i+1))) + h * b(x(i));
    M(i,i) = h * rho(x(i));
end

for i=1:n
    A(i,i+1) = -1/(2*h) * (a(x(i+1)) + a(x(i)));
    A(i+1,i) = -1/(2*h) * (a(x(i+1)) + a(x(i)));
end

%Dirichelt u lijevom rubu
A(1,1)=1000;
A(1,2)=0;
A(2,1)=0;
M(1,1)=10;

%Dirichelt u desnom rubu
% A(n+1,n+1)=1000;
% A(n,n+1)=0;
% A(n+1,n)=0;
% M(n+1,n+1)=10;

%Neumann u lijevom rubu
% A(1,1) = 1/(2*h) * (a(x(2))+a(x(1))) + h/2 * b(x(1));
% M(1,1) = h/2 * rho(x(1));

%Neumann u desnom rubu
A(n+1,n+1) = 1/(2*h) * (a(x(n)) + a(x(n+1))) + h/2 * b(x(n+1));
M(n+1,n+1) = h/2 * rho(x(n+1));

%solve generalized eigenvalues equation -> A*X = lambda*M*X
[X,D] = eig(A,M);

%rjesavamo T'' + lambdaT = 0
A = zeros(n+1,1);
B = zeros(n+1,1);

for k=1:n+1
    for i=1:n
        A(k) = A(k) + h/2 * (rho(x(i+1))*u_0(x(i+1))*X(i+1,k) + rho(x(i))*u_0(x(i))*X(i,k));
        B(k) = B(k) + 1/sqrt(D(k,k)) * h/2 * (rho(x(i+1))*u_1(x(i+1))*X(i+1,k) + rho(x(i))*u_1(x(i))*X(i,k));
    end
end

%diskretizacija t na [0,T]
TIME = 1;
t = 0:0.01:TIME;

T = zeros(n+1,length(t));

for k=1:n+1
    for i=1:length(t)
        T(k,i) = A(k) * cos(sqrt(D(k,k))*t(i)) + B(k) * sin(sqrt(D(k,k))*t(i));
    end
end

u = X*T;

F1 = zeros(n+1, length(t));

for k=1:n+1
    for i=1:length(t)
        for j=1:n
            F1(k,i) = F1(k,i) + h/2 * (f1(x(j+1),t(i))*X(j+1,k) + f1(x(j),t(i))*X(j,k));
        end
    end
end

T1 = zeros(n+1,length(t));

for k=1:n+1
    for i=1:length(t)-1
        temp=0;
        TIME=t(i);
        for j=1:i-1
            temp = temp + 1/sqrt(D(k,k)) * h/2 * (F1(k,j+1)*sin(sqrt(D(k,k))*(TIME - t(i+1))) + ...
                F1(k,j)*sin(sqrt(D(k,k))*(TIME - t(i))));
        end
         T1(k,i) = temp + A(k) * cos(sqrt(D(k,k))*t(i)) + B(k) * sin(sqrt(D(k,k))*t(i));
    end
end

u1 = X*T1;

plot_times = [1 45 75 95];

% for i = 1:4
%     
%   figure(i)
%   k = plot_times(i);
%   plot(x,u(:,k), 'b--', x, u1(:,k), 'r:','linewidth',2);
%   grid on;
%   axis([min(x) max(x) -3 3]);
%   xlabel('X axis','fontSize',14);
%   ylabel('Wave Amplitude','fontSize',14);              
%   titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(t(k)), 'second'];
%   title(titlestring ,'fontsize',14);   
%   legend('f=0','f\neq0');
%   h=gca; 
%   get(h,'FontSize'); 
%   set(h,'FontSize',14);
%   fh = figure(i);
%   set(fh, 'color', 'white'); 
%  
% end

figure(5)

for j = 1:length(t)-1  
  plot(x,u(:,j), 'b--', x, u1(:,j), 'r:','linewidth',2);
  grid on;
  axis([min(x) max(x) -3 3]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
  title(titlestring ,'fontsize',14);
  legend('f=0','f\neq0');
  h=gca; 
  get(h,'FontSize'); 
  set(h,'FontSize',14);
  F=getframe;
            
end


