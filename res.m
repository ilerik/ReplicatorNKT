A = load('u.txt');

G = [ 0.8 1.1 ; 1.2 0.9 ];

k = 1;

ind1 = (mod(1:1:size(A,1), 2*k) == 1);
ind0 = (mod(1:1:size(A,1), 2*k) == 0);

t = load('t.txt');
t = t(:, mod(1:1:length(t), k) == 0);
x = load('x.txt');
    
surf(x, t, A(ind1, :));