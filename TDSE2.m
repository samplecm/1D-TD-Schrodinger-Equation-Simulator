delta_t = 0.05; % what interval of time for each iteration?
tEnd = 1800;%How many time iterations? %9000 for low energy
E = 5.1;

sigma = 20;
k=sqrt(2*(E-(1/(8*sigma^2))));
x0 = 150;
x_int = 3/25;
xPoints = 5000;

%above all the constants are defined

psi = zeros(tEnd,xPoints);
x = zeros(xPoints);
V = zeros(1,xPoints);
%Get the initial gaussian packet ready, and the potential function.
for i = 1:xPoints
  x(i) = -300 + i*x_int; %-60 for down potential..
      psi(1,i) = Gaussian(x(i),sigma,k,x0);
      V(1,i) = Potential(x(i));
end
%Get the initial normalzation factor
P_tot0 = 0;
for i = 1:(xPoints-1)
absPsi = psi.*conj(psi);
P_tot0 = P_tot0 + x_int*((absPsi(1,i)+absPsi(1,i+1))/2);
    
    
end


    


%Get the Ds:
D = zeros(xPoints);
phi = zeros(xPoints);
%Get the phis
mu = -sqrt(-1)*delta_t/(4*x_int^2);
for i = 1:xPoints
   D(i) = 1+(sqrt(-1)*delta_t/2)*(V(1,i)+1/x_int^2);    
end


%
%
%
%
%TTime Loop:
%


t = 0;
for z = 1:tEnd
    t = t + delta_t;
%
%
%
%
%

for i = 2:(xPoints-1)
       phi(i) = psi(z,i)-sqrt(-1)*(delta_t/2)*(-((psi(z,i-1)-2*psi(z,i)+psi(z,i+1))/(2*x_int^2))+(V(1,i)*psi(z,i)));
end




%Now the final algorithm:
d = zeros(xPoints);
f = zeros(xPoints);
%1.)
d(1) = D(1);
f(1) = 0; %because phi(1) must be zero?
 for j = 2:xPoints 
    
     f(j) = phi(j) - mu*f(j-1)/d(j-1);
     d(j) = D(j) - mu*mu/d(j-1); 
 end
 
 psi(z+1,xPoints) = f(xPoints)/d(xPoints);
 for j = (xPoints-1):-1:1

     psi(z+1,j) = (f(j)-mu*psi(z+1,j+1))/d(j);
 end
%  for j = 1:xPoints  %This loop adds on the time dependence of the SE (e^-iwt)
%      psi(z+1,j) = psi(z+1,j)*exp(-sqrt(-1)*E*t);    %exp(-0.5*sqrt(-1)*k^2*t); %w = hk^2/2m
%  end
end
absPsi = psi.*conj(psi);
%Finding probabilities: 
%Originally, 

P_R = 0;
temp = 0;
temp2 = 0;
P_T = 0;
for i = 1:((xPoints/2)-1)
temp = (temp + x_int*((absPsi(tEnd,i)+absPsi(tEnd,i+1))/2));
P_R = temp/P_tot0;    
end
%for i = (xPoints/2):((xPoints-1))
%temp2 = temp2 + x_int*((absPsi(tEnd,i)+absPsi(tEnd,i+1))/2);
%P_T = temp2/P_tot0;    
    
%end


load chirp
sound(y,1/2*Fs)



%Plotting

for i = 1:200:tEnd
hold off;

plot(x,V(1,:))
axis([-240 240 0 1]);
hold on;%set the potential 
plot(x,6*absPsi(i,:))
 legend('V(x)','|\Psi|^{2}')
pause(0.001)


end