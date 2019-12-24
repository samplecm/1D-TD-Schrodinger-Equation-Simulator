function psiG = GaussianR(x, sigma, k,x0)
  i = sqrt(-1);
  psiG = ((2*pi*sigma^2)^-0.25)*exp(-i*k*x)*exp(-(x-x0)^2/(4*sigma^2));
  

  
end
