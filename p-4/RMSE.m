error = 0;
sum = 0;
wref = 16;
for i=1:3001
    sum = sum + (speed(i)/wref - model(i)/wref)^2;
end
mean = sum / 3001;
square = mean^0.5