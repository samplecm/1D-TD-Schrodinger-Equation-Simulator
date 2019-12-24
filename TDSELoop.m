delta_t = 0.2; % what interval of time for each iteration?
tEnd = 1600;%How many time iterations? %9000 for low energy
P_R1 = zeros(1,10);

for w = 0:4
E = 0.2+(w/5);

sigma = 30;%spread of the wavepacket.
k=sqrt(2*(E-(1/(8*sigma^2))));
x0 = -120;
x_int = 2/6.25;
xPoints = 3000;

%above all the constants are defined

psi = zeros(tEnd,xPoints);
x = zeros(xPoints);
V = zeros(1,xPoints);
%Get the initial gaussian packet ready, and the potential function.
for i = 1:xPoints
  x(i) = -480 + i*x_int; %%here -480 is my initial starting x value. adjust this 
                         %accoridng to your x_int so that at the end of the
                         %loop you are at a point symmetrically far from
                         %the y axis.
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
 for j = 1:xPoints  %This loop adds on the time dependence of the SE (e^-iwt)
     psi(z+1,j) = psi(z+1,j)*exp(-sqrt(-1)*E*t);    %exp(-0.5*sqrt(-1)*k^2*t); %w = hk^2/2m
 end
end
absPsi = psi.*conj(psi);
%Finding probabilities: 
%Originally, 
% 
% P_R = 0;
temp = 0;
temp2 = 0;
% P_T = 0;
for i = 1:((xPoints/2)-1)
temp = (temp + x_int*((absPsi(tEnd,i)+absPsi(tEnd,i+1))/2));
P_R1(w+3) = temp/P_tot0;    
end
%for i = (xPoints/2):((xPoints-1))
%temp2 = temp2 + x_int*((absPsi(tEnd,i)+absPsi(tEnd,i+1))/2);
%P_T = temp2/P_tot0;      
%end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





delta_t = 0.2; % what interval of time for each iteration?
tEnd = 1600;%How many time iterations? %9000 for low energy
P_R2 = zeros(1,10);

for w = 0:4
E = 0.2+(w/5);

sigma = 1;
k=sqrt(2*(E-(1/(8*sigma^2))));
x0 = -120;
x_int = 1/6.25;
xPoints = 6000;

%above all the constants are defined

psi = zeros(tEnd,xPoints);
x = zeros(xPoints);
V = zeros(1,xPoints);
%Get the initial gaussian packet ready, and the potential function.
for i = 1:xPoints
  x(i) = -480 + i*x_int; %-60 for down potential..
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
 for j = 1:xPoints  %This loop adds on the time dependence of the SE (e^-iwt)
     psi(z+1,j) = psi(z+1,j)*exp(-sqrt(-1)*E*t);    %exp(-0.5*sqrt(-1)*k^2*t); %w = hk^2/2m
 end
end
absPsi = psi.*conj(psi);
%Finding probabilities: 
%Originally, 
% 
% P_R = 0;
temp = 0;
temp2 = 0;
% P_T = 0;
for i = 1:((xPoints/2)-1)
temp = (temp + x_int*((absPsi(tEnd,i)+absPsi(tEnd,i+1))/2));
P_R2(w+1) = temp/P_tot0;    
end
%for i = (xPoints/2):((xPoints-1))
%temp2 = temp2 + x_int*((absPsi(tEnd,i)+absPsi(tEnd,i+1))/2);
%P_T = temp2/P_tot0;    
    
end

% fileID = fopen('PvsE.txt', 'w');
% for i = 1:numel(E)
%    fprintf(fileID, '%f %f\n', E(i), PRi(i));
% end
% fclose(fileID);