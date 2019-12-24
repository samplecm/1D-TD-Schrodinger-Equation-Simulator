 x = -20:0.05:20;
 %V = zeros(1,401);
% for i = 1:401
%V(1,i) = Potential(x(i));
pot = zeros(1,numel(x));
for i=1:numel(x)
    pot(1,i) = Potential(x(i));
end
% end

%plot(x,V(1,:))

fileID = fopen('expPot.txt', 'w');
for i = 1:numel(x)
   fprintf(fileID, '%f %f \n', x(i),pot(1,i));
end
fclose(fileID);