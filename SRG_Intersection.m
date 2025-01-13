close all
clear
tic

% Define the matrix L
L = [1 0 -1 0; 
     0 2  0 1;
     1 1  0 0;
     0 0  1 1];

% Define the range of c values
c_values = linspace(-10, 10, 1000); % Replace -10 and 10 with your desired range, 100 points

% Initialize arrays to store singular values
min_singular_values = zeros(size(c_values));
max_singular_values = zeros(size(c_values));

% Loop over c values to calculate singular values
for i = 1:length(c_values)
    c = c_values(i);
    M = L - c * eye(size(L)); % Compute L - cI
    singular_values = svd(M); % Singular value decomposition
    min_singular_values(i) = min(singular_values);
    max_singular_values(i) = max(singular_values);
end

% Initialize arrays for (x, y) points
x_min = zeros(1, length(c_values)-1);
y_min = zeros(1, length(c_values)-1);
x_max = zeros(1, length(c_values)-1);
y_max = zeros(1, length(c_values)-1);

% Compute (x_min, y_min) and (x_max, y_max) values
for i = 1:length(c_values)-1
    % For minimum singular values
    x_min(i) = ( ...
        min_singular_values(i)^2 - min_singular_values(i+1)^2 + ...
        c_values(i+1)^2 - c_values(i)^2 ...
    ) / ( 2*( c_values(i+1) - c_values(i) ) );

    y_min(i) = sqrt( ...
        max(0, min_singular_values(i)^2 - ( x_min(i) - c_values(i) )^2)); % Ensure no negative values under square root

    % For maximum singular values
    x_max(i) = ( ...
        max_singular_values(i)^2 - max_singular_values(i+1)^2 + ...
        c_values(i+1)^2 - c_values(i)^2) / ( 2*( c_values(i+1) - c_values(i) ) );

    y_max(i) = sqrt( ...
        max(0, max_singular_values(i)^2 - ( x_max(i) - c_values(i) )^2)); % Ensure no negative values under square root
end

% Plot the results
figure;
hold on;

% Define the x and y coordinates for the filled area
fill_x = [x_min, x_max]; % Combine x_min and reversed x_max
fill_y = [y_min, y_max]; % Combine y_min and reversed y_max

% Fill the area inside the points
fill(fill_x, fill_y, [0.8, 0.8, 1], 'EdgeColor', 'none'); % Light blue fill
fill(fill_x, -fill_y, [0.8, 0.8, 1], 'EdgeColor', 'none'); % Light blue fill

% Plot singular value curves (optional for reference)
%plot(c_values, min_singular_values, 'b-', 'LineWidth', 2); % Min singular values
%plot(c_values, max_singular_values, 'r-', 'LineWidth', 2); % Max singular values

% Customize the plot
title('Scaled Relative Graph');
axis equal;
grid on;
hold off;
toc