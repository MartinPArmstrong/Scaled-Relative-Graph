close all

% Define matrix A
A = [1 0 -1 0; 
     0 2  0 1;
     1 1  0 0;
     0 0  1 1];

% Define a range of c values
c_values = linspace(-5, 5, 100);
num_c = length(c_values);

% Preallocate arrays for singular values
min_singular = zeros(size(c_values));
max_singular = zeros(size(c_values));

% Precompute theta and the unit circle using complex exponentials
theta_vals = linspace(0, 2*pi, 200);  % Increase resolution if desired
unit_circle = exp(1i * theta_vals);  % Complex unit circle for plotting circles

% Prepare figure
figure;
hold on;
xlim([-5 5])  
ylim([-5 5])  
grid on;
xlabel('Real axis');
ylabel('Imaginary axis');
title('Scaled Relative Graph');

% Storage for circles data
circles_center = zeros(1, num_c);
circles_min_radius = zeros(1, num_c);
circles_max_radius = zeros(1, num_c);

% Plot circles iteratively and store parameters
for k = 1:num_c
    c = c_values(k);
    M = A - c*eye(size(A));
    
    % Compute singular values
    s = svd(M);
    min_singular(k) = min(s);
    max_singular(k) = max(s);
    
    % Store center and radi for envelope calculation later
    circles_center(k) = c;
    circles_min_radius(k) = min(s);
    circles_max_radius(k) = max(s);
    
    % Plot circles using complex transformation
    z_min_circle = c + min_singular(k) * unit_circle;
    z_max_circle = c + max_singular(k) * unit_circle;
    plot(real(z_min_circle), imag(z_min_circle), 'b');
    plot(real(z_max_circle), imag(z_max_circle), 'r');
    
    drawnow;
    pause(0.01);
end