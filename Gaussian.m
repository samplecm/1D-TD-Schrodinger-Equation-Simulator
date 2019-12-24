function psiG = Gaussian(x, sigma, k,x0)

  psiG = ((2*pi*sigma^2)^-0.25)*exp(-sqrt(-1)*k*x)*exp(-(x-x0)^2/(4*sigma^2));
  

  
end
