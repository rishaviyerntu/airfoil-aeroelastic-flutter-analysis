close all;
clear all;

syms tt 

mu = 500;  % airfoil density (kg/m3)
c = 1;  % chord length (m)
x_f = 0.45 * c; % Location of elastical axis measured from leading edge of airfoil
a = -0.05 * c;  % Location of elastic axis measured from center of airfoil
t = 0.01;   % Airfoil thickness (m)
b = c / 2;  % Airfoil mid chord (m)
e = x_f - 0.25 * c;   % Elastic axis distance fraction from aerodynamic centre

omega_h = 5 * 2 * pi;    % Plunge natural frequency (rad/s)
omega_alpha = 15 * 2 * pi;   % Pitch natural frequency (rad/s)

m = mu * c * t; % Mass per unit span of airfoil (kg/m)

S = -m * a;  

K_h = m * (omega_h ^ 2);    % Plunge stiffness

I_c = (1 / 12) * m * (c ^ 2);   % Moment of inertia of airfoil about center of airfoil
I_alpha = I_c + (m * (a ^ 2));  % Parallel axis theorem to obtain moment of inertia about elastic axis

K_alpha = I_alpha * (omega_alpha ^ 2);  % Pitch stiffness

rho = 1.225; % Density of air (kg/m3)

% Wagner function parameters
psi_1 = 0.165;
psi_2 = 0.335;
epsilon_1 = 0.0455;
epsilon_2 = 0.3;

U_min = 0;  % (m/s)
U_max = 80; % (m/s)
dU = 0.1; % (m/s)

U_range = U_min:dU:U_max;   % Range of velocity values from 0 to 80 m/s in steps of 0.1 m/s

EV = [];
Frequency = [];
Damping_Ratio = [];

Q = zeros(8);

for U = U_range

    phi = 1 - (psi_1 * exp((-epsilon_1 * U * tt) / b)) - (psi_2 * exp((-epsilon_2 * U * tt) / b));    % Wagner Function
    phi_dot = diff(phi, tt);

    phi_0 = subs(phi, tt, 0);
    phi_dot_0 = subs(phi_dot, tt, 0);

    M = [m + (rho * pi * (b ^ 2)), S - (rho * pi * (b ^ 2) * (x_f - (c / 2))); S - (rho * pi * (b ^ 2) * (x_f - (c / 2))), I_alpha + (rho * pi * (b ^ 2)) * (((x_f - (c/2)) ^ 2) + ((b ^ 2) / 8))];
    C = rho * pi * U * c * [phi_0, (c / 4) + phi_0 * (0.75 * c - x_f); -e * c * phi_0, (0.75 * c - x_f) * (0.25 * c - (e * c * phi_0))];
    K = [K_h + (rho * pi * U * c * phi_dot_0), rho * pi * U * c * ((U * phi_0) + (0.75 * c - x_f) * phi_dot_0); -rho * pi * U * e * (c ^ 2) * phi_dot_0, K_alpha - (rho * pi * U * e * (c ^ 2)) * ((U * phi_0) + (0.75 * c - x_f) * phi_dot_0)];
    W = 2 * rho * pi * (U ^ 3) *[-psi_1 * (epsilon_1) ^ 2 / b, -psi_2 * (epsilon_2) ^ 2 / b, psi_1 * epsilon_1 * (1 - epsilon_1 * (1 - 2 * e)), psi_2 * epsilon_2 * (1 - epsilon_2 * (1 - 2 * e)); e * c * psi_1 * (epsilon_1) ^ 2 / b, e * c * psi_2 * (epsilon_2) ^ 2 / b, -e * c * psi_1 * epsilon_1 * (1 - epsilon_1 * (1 - 2 * e)), -e * c * psi_2 * epsilon_2 * (1 - epsilon_2 * (1 - 2 * e))];
    
    B = [1 0; 1 0; 0 1; 0 1];
    G = [-epsilon_1 * U / b, 0, 0, 0; 0, -epsilon_2 * U / b, 0, 0; 0, 0, -epsilon_1 * U / b, 0; 0, 0, 0, -epsilon_2 * U / b];

    Q(1:2,1:2) = -M \ C;
    Q(1:2,3:4) = -M \ K;
    Q(1:2,5:8) = -M \ W;
    Q(3:4,1:2) = eye(2);
    Q(5:8,3:4) = B;
    Q(5:8,5:8) = G;
    
    Eigenvalues = eig(Q);   % Eigenvalue analysis
    [Eigenvalues_imag, Index] = sort(imag(Eigenvalues), 'ascend');
    Eigenvalues_sorted = Eigenvalues(Index, :);

    EV = [EV; Eigenvalues_sorted(7) Eigenvalues_sorted(8)];

end

found = false;

tolerance = 1e-6;

for i = 1:1:size(U_range,2)

    
    Frequency(i,1) = abs(EV(i,1));
    Frequency(i,2) = abs(EV(i,2));
    Damping_Ratio(i,1) = -real(EV(i,1))/abs(EV(i,1));
    Damping_Ratio(i,2) = -real(EV(i,2))/abs(EV(i,2));
    
    if Damping_Ratio(i,1)<0 && found == false && abs(Damping_Ratio(i,1)) > tolerance
        Flutter_speed = U_range(i);
        Flutter_frequency = Frequency(i,1)/2/pi;
        Flutter_damping = Damping_Ratio(i,1);
        found = true;
    end
    
    if Damping_Ratio(i,2)<0 && found == false && abs(Damping_Ratio(i,2)) > tolerance
        Flutter_speed = U_range(i);
        Flutter_frequency = Frequency(i,2)/2/pi;
        Flutter_damping = Damping_Ratio(i,2);
        found = true;
    end
    
end

fprintf("Flutter Speed: %.1f m/s \n", Flutter_speed);
fprintf("Flutter Frequency: %.3f Hz \n", Flutter_frequency);

figure(1)
title("Natural Frequency vs Airspeed")
xlabel("U (m/s)")
ylabel("Natural Frequency (Hz)")
hold on
plot(U_range, Frequency(:,1)/2/pi,"DisplayName","Wagner Plunge");
plot(U_range,Frequency(:,2)/2/pi,"DisplayName","Wagner Pitch");
plot(Flutter_speed, Flutter_frequency,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
legend show
hold off
grid on

figure(2)
title("Damping Ratio vs Airspeed")
xlabel("U (m/s)")
ylabel("Damping Ratio")
hold on
plot(U_range,Damping_Ratio(:,1),"DisplayName","Wagner Plunge");
plot(U_range,Damping_Ratio(:,2),"DisplayName","Wagner Pitch");
plot(Flutter_speed, Flutter_damping,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
legend show
hold off
grid on

figure(3)
title("Natural Frequency and Damping Ratio vs Airspeed")
xlabel("U (m/s)")
hold on
yyaxis left
ylabel("Natural Frequency (Hz)")
plot(U_range,Frequency(:,1)/2/pi,"DisplayName","Plunge");
plot(U_range,Frequency(:,2)/2/pi,"DisplayName","Pitch");
plot(Flutter_speed, Flutter_frequency,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
yyaxis right
ylabel("Damping Ratio")
plot(U_range,Damping_Ratio(:,1),"DisplayName","Plunge");
plot(U_range,Damping_Ratio(:,2),"DisplayName","Pitch");
plot(Flutter_speed, Flutter_damping,"x","DisplayName", "Flutter Point", "Color", "green", "MarkerSize", 15)
legend show
hold off
grid on
