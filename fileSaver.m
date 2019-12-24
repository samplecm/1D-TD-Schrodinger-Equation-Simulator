


% Store data in 'xy.txt'
fileID = fopen('rightProp.txt', 'w');
for i = 1:5000
   fprintf(fileID, '%f %f\n', x(i,1), absPsi(1750,i));
end
fclose(fileID);