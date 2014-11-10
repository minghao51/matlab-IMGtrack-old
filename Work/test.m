% Create some example input data.
x = 1:10
y = cumsum( randn(1,10) );
lower = y - ( rand(1,10) );
upper = y + ( rand(1,10) );

% Convert absolute lower and upper bounds into the relative values 
% values that are expected by the errorbar function.
L = y - lower;
U = upper - y;

figure(1);
clf;
hold('on');
plot( x, y, 'b-' );
errorbar( x, y, L, U, 'r', 'Marker', 'none', 'LineStyle', 'none' );