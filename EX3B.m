% coded in matlab r2025a
% exercise 1.2
close all
clearvars
format short e

% beam constants
% let x denote the distance along the span of the beam

L = 1; % length of the beam - m
A = 10^-4; % cross section - m^2
rho = 7850; % mass per unit volume - kg/m^3
E = 2.1e11; % Young's modulus - Pascals
I = (10^-8)/12;  % Moment of inertia - m^4 

% convert to mass per unit length
m = rho*A;

nel = 100; % number of elements
nno = nel + 1; % number of nodes
nbc = 2;  % number of boundary conditions (used for error detection)

lel = L/nel; % element length

%% construct mass and stiffness matrices
Mel = (rho*A*lel/420).*[    156   22*lel      54  -13*lel;
                         22*lel  4*lel^2  13*lel -3*lel^2;
                             54   13*lel     156  -22*lel;
                        -13*lel -3*lel^2 -22*lel  4*lel^2]; % element mass matrix

Kel = (E*I/lel^3).* [   12   6*lel    -12   6*lel;
                     6*lel 4*lel^2 -6*lel 2*lel^2;
                       -12  -6*lel     12  -6*lel;
                     6*lel 2*lel^2 -6*lel 4*lel^2]; % element stiffness matrix

% Initialize global mass matrix
M = zeros(nno*2); % 2 degrees of freedom per node
K = zeros(nno*2);
% Assemble the global mass matrix
for e = 1:nel
    idx = [2*e-1, 2*e, 2*e+1, 2*e+2]; % global index for the element
    M(idx, idx) = M(idx, idx) + Mel; % add element mass matrix to global
    K(idx, idx) = K(idx, idx) + Kel;
end

% Apply boundary conditions (pin at first node, vertical slider at second
% node)
rowColIdxs = 2:2*nno-1; % everything but the first and last row/col
Mbc = M(rowColIdxs,rowColIdxs);
Kbc = K(rowColIdxs,rowColIdxs);
rowColIdxs = 2:2*nno-1;  % you used this earlier

lastNode = nno;
full_dof_last_disp = 2*lastNode - 1;    % should be 201
% find reduced index
reduced_index = find(rowColIdxs == full_dof_last_disp);
fprintf('full DOF for last node disp = %d, reduced index = %d\n', full_dof_last_disp, reduced_index);
fprintf('respDof (size(Mbc,1)-1) = %d\n', size(Mbc,1)-1);
% check if size matches with expected size based on the number of boundary
% conditions
assert(size(Mbc,1) == 2*nno-nbc, "The size of the matrix after applying boundary" + ...
    " conditions does not match with the number of boundary conditions specified: nbc = %d",nbc)
%%
% create  C and D as in slide 30 "SD2 Numerical modal analysis.pdf"
zeroM = zeros(size(Mbc));
Cbc = [zeroM Mbc; Mbc zeroM];
Dbc = [Kbc zeroM; zeroM -Mbc];

%% question a
% calculate the six lowest eigenfrequencies (in Hz) of Finite Element models
% of the beam using the MATLAB command 'eig'
[eigenVectorsa, eigenValuesa] = eig(Cbc, Dbc);
eigenValuesa = imag(diag(eigenValuesa));
validEvs = eigenValuesa > 0;

% filter negatives
eigenValuesa = eigenValuesa(validEvs); % Filter out negative eigenvalues
eigenVectorsa = eigenVectorsa(1:size(Mbc,1), validEvs);

% sort eigenvalues
eigenValuesa = sort(eigenValuesa);

% Extract eigenfrequencies from the eigenvalues matrix
eigenfrequenciesa = eigenValuesa / (2 * pi);
eigenfrequenciesa = eigenfrequenciesa(1:6); % Select the six lowest frequencies
% Display the calculated eigenfrequencies
disp('The six lowest eigenfrequencies calculated with ''eig'' (in Hz) are:');
disp(eigenfrequenciesa);


% Calculate the mode shapes corresponding to the eigenfrequencies
modeShapesa = eigenVectorsa(:, 1:6);
% add the bc columns/rows back to the eigenvectors (all zeros)
modeShapesaFull = zeros(2*nno,6);
modeShapesaFull(2:end-1,:) = modeShapesa;


% Display the calculated mode shapes
tiledlayout(2,3);
title("Mode shapes calculated with ''eig''")
subtitle(num2str(nel) + " number of elements.")
xno = 1:1:nno;
dispDOFs = 1:2:nno*2;
for p = 1:6
    nexttile;
    plot(xno, imag(modeShapesaFull(dispDOFs,p)))
end

%% Exercise 3a 
% Method 1

f = 0.2:0.2:500;  % Frequency range in Hz

% Dynamic stiffness matrix 
for i = 1:length(f)
    K_dynamic(:,:,i) = inv(Kbc - (2*pi*f(i))^2 * Mbc); % Computing the inverse dynamic stiffness matrix
    lastDiag(i) = K_dynamic(end,end,i); % Getting the last diagonal components
end

% FRF
figure;

% Magnitude
subplot(2,1,1);              
loglog(f, abs(lastDiag), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FRF Magnitude for 100 elements');
grid on;

% Phase
subplot(2,1,2);              
semilogx(f, angle(lastDiag)*180/pi, 'LineWidth', 1.5); % convert phase to degrees
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('FRF Phase for 100 elements');
grid on;

%% Method 2
% Inputs:
% omega  : frequency (scalar)
% u0     : matrix of mode shapes, each column u0(:,k) corresponds to u0k
% m      : vector of modal masses, m(k)
% omega0 : vector of natural frequencies, omega0(k)

omega = f*2*pi
omega0 = eigenfrequenciesa*2*pi


% Number of modes
n = length(omega0);

% Initialize H
H = zeros(size(u0,1));

% Sum over modes
for k = 1:n
    H = H + (u0(:,k) * u0(:,k)') / (m(k) * (omega0(k)^2 - omega^2));
end
 
%% B


L  = 1.0;              % length [m]
E  = 2.1e11;           % Young's modulus [N/m^2]
A  = 1e-4;             % cross-sectional area [m^2]
rho = 7850;            % density [kg/m^3]
I  = 1e-8/12;          % second moment of area [m^4]
m  = rho*A;            % mass per unit length

% convert to mass per unit length
m = rho*A;
nev = 6; % number of modes/eigenvalues to analyze

nel = 100; % number of elements
nno = nel + 1; % number of nodes
nbc = 2;  % number of boundary conditions (used for error detection)

lel = L/nel; % element length

solver = "eig";

 Mel = (rho*A*lel/420).*[    156   22*lel      54  -13*lel;
                         22*lel  4*lel^2  13*lel -3*lel^2;
                             54   13*lel     156  -22*lel;
                        -13*lel -3*lel^2 -22*lel  4*lel^2]; % element mass matrix

Kel = (E*I/lel^3).* [   12   6*lel    -12   6*lel;
                     6*lel 4*lel^2 -6*lel 2*lel^2;
                       -12  -6*lel     12  -6*lel;
                     6*lel 2*lel^2 -6*lel 4*lel^2]; % element stiffness matrix

% Initialize global mass matrix
M = zeros(nno*2); % 2 degrees of freedom per node
K = zeros(nno*2);
% Assemble the global mass matrix
for e = 1:nel
    idx = [2*e-1, 2*e, 2*e+1, 2*e+2]; % global index for the element
    M(idx, idx) = M(idx, idx) + Mel; % add element mass matrix to global
    K(idx, idx) = K(idx, idx) + Kel;
end

% Apply boundary conditions (pin at first node, vertical slider at second
% node)
rowColIdxs = 2:2*nno-1; % everything but the first and last row/col
Mbc = M(rowColIdxs,rowColIdxs);
Kbc = K(rowColIdxs,rowColIdxs);

lastNode = nno;
full_dof_last_disp = 2*lastNode - 1;    % should be 201
% find reduced index
reduced_index = find(rowColIdxs == full_dof_last_disp);
fprintf('full DOF for last node disp = %d, reduced index = %d\n', full_dof_last_disp, reduced_index);
fprintf('respDof (size(Mbc,1)-1) = %d\n', size(Mbc,1)-1);
% check if size matches with expected size based on the number of boundary
% conditions
assert(size(Mbc,1) == 2*nno-nbc, "The size of the matrix after applying boundary" + ...
    " conditions does not match with the number of boundary conditions specified: nbc = %d",nbc)

% create  C and D as in slide 30 "SD2 Numerical modal analysis.pdf"
zeroM = zeros(size(Mbc));
Cbc = [zeroM Mbc; Mbc zeroM];
Dbc = [Kbc zeroM; zeroM -Mbc];


% calculate the six lowest eigenfrequencies (in Hz) of Finite Element models
% of the beam using the MATLAB command specified in solver
if solver == "eig"
    t0 = cputime; % start timing
    [eigenVectors, eigenValues] = eig(Dbc, -Cbc);
elseif solver == "eigs"
    % make matrices sparse
    spaDbc = sparse(Dbc);
    spaCbc = sparse(Cbc);
    t0 = cputime; % start timing
    [eigenVectors, eigenValues] = eigs(spaDbc, -spaCbc, nev*2, "smallestabs");
else
    error('Specified solver not recognized.')
end

t1 = cputime; % stop timing
reqTime = t1 - t0;


eigenValues = imag(diag(eigenValues));
validEvs = eigenValues > 0;

% filter negatives
eigenValues = eigenValues(validEvs); % Filter out negative eigenvalues
eigenVectors = eigenVectors(1:size(Mbc,1), validEvs);

% sort eigenvalues
sortedEigenValues = sort(eigenValues, 'ascend');
sortedEigenVectors = zeros(size(eigenVectors,1),nev);
for i = 1:nev
    smallEigenValue = sortedEigenValues(i);
    idx = find(eigenValues == smallEigenValue);
    sortedEigenVectors(:,i) = eigenVectors(:,idx);
end

% Extract eigenfrequencies from the eigenvalues matrix
eigenfrequencies = sortedEigenValues ./ (2 * pi);
eigenfrequencies = eigenfrequencies(1:6); % Select the six lowest frequencies
% Display the calculated eigenfrequencies
disp('The six lowest eigenfrequencies calculated with ''eig'' (in Hz) are:');
disp(eigenfrequencies);  % Hz

omega_ok = (eigenfrequencies.^2 / L^2) * sqrt(E*I/m);  % rad/s
freq_n  = omega_ok / (2*pi);                      % Hz

disp('Natural frequencies (Hz):');
disp(freq_n.');

xi = 0.02;

mu   = -xi * omega_ok;
omega_k = omega_ok .* sqrt(1 - xi^2);
lambda_k = mu + 1i*omega_k;

disp('Damped eigenvalues (rad/s):');
disp(lambda_k.');
% respDof as you used it
respDof = size(Mbc,1)-1;   % last node displacement in reduced system

% show the mode shape entries at the response DOF
phi_at_resp = sortedEigenVectors(respDof, :);    % 1 x 6
disp('Mode shape values at response DOF (respDof):');
disp(phi_at_resp);

disp('Absolute values:');
disp(abs(phi_at_resp));

% Also print indices where abs>eps
threshold = 1e-8;
nonzero_idx = find(abs(phi_at_resp) > threshold);
fprintf('Modes with |phi| > %.1e at response DOF: %s\n', threshold, mat2str(nonzero_idx));

excDof  = respDof;         % same node for excitation
f = 0.2:0.2:500; % Hz 
w = 2*pi*f; % rad/s
mk = zeros(nev,1);
for k = 1:nev
    phi_k = sortedEigenVectors(:,k);
    mk(k) = phi_k' * Mbc * phi_k;
end
H_total = zeros(length(w),1);           % total FRF
H_modes = zeros(length(w),6);           % contributions from each mode

for k = 1:6
    phi_r = sortedEigenVectors(:,k)';
    phi_q = sortedEigenVectors(:,k);
    Hk = (phi_r * phi_q) ./ ...
        (mk(k) * (omega_ok(k)^2 - w.^2 + 2j*xi*omega_k(k)*w));  % FRF of mode k
    
    H_modes(:,k) = Hk;    % store individual mode contribution
    H_total = H_total + Hk;
end


figure;

% Magnitude plot
subplot(2,1,1);
loglog(f, abs(H_modes), 'LineWidth', 2);  % total FRF in black
hold on;
for k = 1:6
    loglog(f, abs(H_modes(:,k)), 'LineWidth', 1);  % each mode
end
hold off;
xlabel('Frequency [Hz]'); ylabel('|H| [m/N]');
title('FRF (modal superposition, 6 modes, \xi=0.02)');
grid on;

% Phase plot
subplot(2,1,2);
semilogx(f, angle(H_modes)*180/pi, 'LineWidth', 2);  % total FRF phase
hold on;
for k = 1:6
    semilogx(f, angle(H_modes(:,k))*180/pi, 'LineWidth', 1);  % each mode
end
hold off;
xlabel('Frequency [Hz]'); ylabel('Phase [deg]');
grid on;
