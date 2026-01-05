close all;
clear all;

mu = 500;  % airfoil density (kg/m3)
c = 1;  % chord length (m)
x_f = 0.45 * c; % Location of elastical axis measured from leading edge of airfoil
a = -0.05 * c;  % Location of elastic axis measured from center of airfoil
t = 0.01;   % Airfoil thickness (m)
b = c / 2;  % Airfoil mid chord (m)



omega_h = 5 * 2 * pi;    % Plunge natural frequency (rad/s)
omega_alpha = 15 * 2 * pi;   % Pitch natural frequency (rad/s)

m = mu * c * t; % Mass per unit span of airfoil (kg/m)

S = -m * a;  

K_h = m * (omega_h ^ 2);    % Plunge stiffness

I_c = (1 / 12) * m * (c ^ 2);   % Moment of inertia of airfoil about center of airfoil
I_alpha = I_c + (m * (a ^ 2));  % Parallel axis theorem to obtain moment of inertia about elastic axis

K_alpha = I_alpha * (omega_alpha ^ 2);  % Pitch stiffness

rho = 1.225; % Density of air (kg/m3)

V1 = [];    %plunge eigenvalue
d1 = [];    % plunge damping
V2 = [];    % Pitch eigenvalue
d2 = [];    % Pitch damping
d_all = [];

U_min = 0;  % (m/s)
U_max = 80; % (m/s)
dU = 0.1; % (m/s)

tolerance = 1e-3;


U_range = U_min:dU:U_max;   % Range of velocity values from 0 to 80 m/s in steps of 0.1 m/s


for U = U_range

    ii = 3; % Eigenvalue index to extract

    for omega = [sqrt(K_h / m) sqrt(K_alpha / I_alpha)] % Initial estimate of flutter frequency

        error = 10; % Arbitrary number to determine convergence to stop p-k method iteration

        while error > tolerance     % While loop for eigenvalue iteration, stops when eigenvalue converges

            k = (omega * b) / U;
            C_k = 1 - (0.165 ./ (1 - 1i * (0.0455 ./ k))) - (0.335 ./ (1 - 1i * (0.30 ./ k)));  % Reduced frequency k using Theodorsen Function approximation

            M_matrix = [m, S; S, I_alpha];
            K_matrix = [K_h, 0; 0, K_alpha];
            D_matrix = rho * pi * (b ^ 2) * [-1, a; a, -((a ^ 2) + ((b ^ 2) / 8))];
            E_k_matrix = rho * pi * b * U * [-2 * C_k, -b + (2 * C_k * (a - (b / 2))); b - (b - (2 * a + b) * C_k), (b - (2 * a + b) * C_k) * (a - (b / 2))];
            F_k_matrix = rho * pi * b * (U ^ 2) * [0, -2 * C_k; 0, b - (b - (2 * a + b) * C_k)];
       
            A_k = [zeros(2), eye(2); (M_matrix - D_matrix) \ (F_k_matrix - K_matrix), (M_matrix - D_matrix) \ E_k_matrix];

            Eigenvalues = eig(A_k); % Eigenvalue analysis

            [Eigenvalues_Imag, Index] = sort(imag(Eigenvalues), 'ascend');  % Arranging the eigenvalues in ascending order by their imaginary parts

            Eigenvalue_Imag_1 = Eigenvalues_Imag(ii);   % Extracting imaginary part of largest/2nd largest eigenvalue for checking against initial assumed frequency (omega)
            error = abs(Eigenvalue_Imag_1 - omega);     % Calculate error between imaginary part of eigenvalue against initial assumed frequency (omega)
            omega = Eigenvalue_Imag_1;  % Initial assumed frequency is updated to the new eigenvalue imaginary part

        end

        Eigenvalues_sorted = Eigenvalues(Index);

        % Section to obtain and store frequencies and damping

        % plunge values
        if ii == 3

            V1 = [V1 abs(Eigenvalues_sorted(ii))];
            d1 = [d1 -1 * real(Eigenvalues_sorted(ii)) / abs(Eigenvalues_sorted(ii))];

        end

        % pitch values
        if ii == 4

            V2 = [V2 abs(Eigenvalues_sorted(ii))];

            d2 = [d2 -1 * real(Eigenvalues_sorted(ii)) / abs(Eigenvalues_sorted(ii))];


        end

        d_all = [d_all Eigenvalues_sorted];

        ii = ii + 1;

    end

end

% Loop to find the flutter speed and frequency at the point where damping
% goes negative
for i = 1:1:length(U_range)
    
    if d1(i)<0 
        Flutter_speed = U_range(i);
        Flutter_frequency = V1(i)/2/pi;
        Flutter_damping = d1(i);
        break;
    end
    
    if d2(i)<0 
        Flutter_speed = U_range(i);
        Flutter_frequency = V2(i)/2/pi;
        Flutter_damping = d2(i);
        break;
    end
    
end

fprintf("Flutter Speed: %.2f m/s \n", Flutter_speed);
fprintf("Flutter Frequency: %.3f Hz \n", Flutter_frequency);

figure(1)
title("Natural Frequency vs Airspeed")
xlabel("U (m/s)")
ylabel("Natural Frequency (Hz)")
hold on
plot(U_range,V1/2/pi,"DisplayName","P-K Plunge");
plot(U_range,V2/2/pi,"DisplayName","P-K Pitch");
plot(Flutter_speed, Flutter_frequency,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
legend show
hold off
grid on

figure(2)
title("Damping Ratio vs Airspeed")
xlabel("U (m/s)")
ylabel("Damping Ratio")
hold on
plot(U_range,d1,"DisplayName","P-K Plunge");
plot(U_range,d2,"DisplayName","P-K Pitch");
plot(Flutter_speed, Flutter_damping,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
legend show
hold off
grid on