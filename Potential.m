function pot = Potential(x)

pot = 0;
if (x>-5)&&(x<-2.5)
    pot = 2;
elseif (x>-2.5)&&(x<0)
    pot = 5;
elseif (x>0)&&(x<2.5)
    pot = -3;
elseif (x>2.5) 
    pot = 3;
end
% pot = -4/(1+exp(-x/5))+2;
% pot = 0;
% if (x>0) && (x < 2)
%     pot = 5;
% elseif (x>=2)&& (x<3)
%     pot = -2;
% elseif (x>=3)
%     pot = 3;
% end

% if (x>0)&& (x<1)
%     pot = -5;
% elseif (x>=1)&& (x<3)
%     pot = 2;
% elseif (x>=3)
%     pot = -3;
% end

% pot = 0;
% if (x>0) && (x<2)
%     pot = 1;
% end
% if x > 2
%     pot = 0.5;
% end



 
% pot = 0;
% if (x>0)
%     pot = 1;
% end
 
 
 
 
%Weird Step left (*2)
% pot = 0;
% if (x<-2) && (x>-4) 
%     pot = 2;
% end
% 
% if (x>-2) && (x < 0)
% pot = 1;
% end
% if (x>0) && (x<2)
%     pot = -2;
% end
% if (x>2) && (x<6)
%    pot = 0;
% end
% 
% if (x>6)
%     pot = 0.5;
% end
%First Run - E = 1.5, sigma = 10: PR = .7369





%Weird Step Right
% function pot = Potential(x)
% pot = 0;
% if (x<-1) && (x>-3) 
%     pot = -0.5;
% end
% 
% if (x>-1) && (x < 0)
% pot = -2.5;
% end
% if (x>0) && (x<1)
%     pot = 0.5;
% end
% if (x>1) && (x<2)
%    pot = 1.5;
% end
% 
% if (x>2)
%     pot = -0.5;
% end


% if x <= 0
%   pot = exp(4*x);
%   else
%   pot = 2-exp(-4*x);
%  end
% 

%pot = -1/(1+exp(-2*x/3));