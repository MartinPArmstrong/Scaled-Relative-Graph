close all
clear
tic

% Define a 2x2 complex matrix A
A = [1 0 -1 0; 
     0 2  0 1;
     1 1  0 0;
     0 0  1 1];

% Compute the Beltrami-Klein transform of A
Z = BeltramiKlein(A);

% Numerical Range
figure;
Nr = fv(Z, 1, 48);
grid on;
axis equal;
title('Numerical Range');


% Compute G using the inverse mapping
G = BeltramiKleinInverse(Nr);

% Plot the complex numbers in vector G in the complex plane and fill the area
figure;
% Use 'fill' for points in G as well
fill(real(G), imag(G), [0.8, 0.8, 1], 'EdgeColor', 'none');
hold on;
fill(real(G), imag(-G), [0.8, 0.8, 1], 'EdgeColor', 'none');
grid on;
axis equal;
title('Scaled Relative Graph');
toc

%% Beltrami Klein Functions

function result = BeltramiKlein(A)
    % Define identity matrix of appropriate size
    I = eye(size(A));
    
    % Compute the matrix (I + A'*A)
    M = I + A'*A;
    
    % Use matrix square root and then invert to get (I + A'*A)^(-0.5)
    M_inv_sqrt = inv(sqrtm(M));
    
    % Compute the result using the provided expression
    result = M_inv_sqrt * (A' - 1i*I) * (A - 1i*I) * M_inv_sqrt;
end

function result = BeltramiKleinInverse(z)
    % Compute the inverse mapping using element-wise operations.
    % Here, z is assumed to be a complex number or an array of complex numbers.
    numerator   = imag(z) - 1i*sqrt(1 - z.*conj(z));
    denominator = real(z) - 1;
    
    result = numerator ./ denominator;
end

