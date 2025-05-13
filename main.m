close all
clear
clc

rng('default');

% system 9 in Table 1 and its definition is in table 5
A = [0.3, 0.4, 0.2, 0.2;
     0.2, 0.3, 0.3, 0.2;
     0.2, 0.2, 0.4, 0.4;
     0.4, 0.2, 0.2, 0.4];
B = eye(4); 
Q = [3, -1, 0, -1;
     -1, 3, -1, 0;
     0, -1, 3, -1;
     -1, 0, -1, 3]; 
R = eye(4); 
% subsystem total count
N = 4;

% neighbors
neighbors = {[1, 2], [2, 3], [3, 4], [4, 1]};

% data collection
l = 4; % sample
T = 25; % total time step in one round 
x0 = [1; 2; 3; 4]; % initial state
X = zeros(4, T); % state record in one round
U = zeros(4, T); % input record in one round
K_rec = zeros(4, 2, T); % K matrix record

X(:, 1) = x0;

total_round = 10; % total round
X_rec = zeros(total_round, 4, T);

for m = 1: total_round
    for k = 1:T-1
        % use random K to get input u
        if k < l
            K_rec(:, :, k) = randn(4, 2);
        end
        for i = 1:N
            neighbor_idx = neighbors{i};
            X_Ni_temp = X(neighbor_idx, k);
            U(i, k) = K_rec(i, :, k)*X_Ni_temp; % random input
        end

        X(:, k+1) = A * X(:, k) + B * U(:, k);

        if k >= l
            % each subsystem
            for i = 1:N
                % get neighbors input and state
                neighbor_idx = neighbors{i};
                X_Ni = X(neighbor_idx, k-l+1:k+1);
                U_Ni = U(neighbor_idx, k-l+1:k+1);

                % construct matrix D_Ni and X_Ni 
                D_Ni = [X_Ni(:, 1:l); U_Ni(:, 1:l)];
                X_Ni_next = X_Ni(:, 2:l+1);

                % define matrix Lambda_Ni
                Q_Ni = Q(neighbor_idx, neighbor_idx);
                R_Ni = R(neighbor_idx, neighbor_idx);
                Lambda_Ni = blkdiag(Q_Ni, R_Ni);

                % use CVX to solve the problem
                cvx_begin sdp
                variable H_Ni_11(size(Q_Ni, 1), size(Q_Ni, 1)) symmetric
                variable H_Ni_12(size(Q_Ni, 1), size(R_Ni, 1))
                variable H_Ni_22(size(R_Ni, 1), size(R_Ni, 1)) symmetric
                variable W_Ni(size(Q_Ni, 1), size(Q_Ni, 1)) symmetric

                minimize -(trace(W_Ni))
                subject to
                % First LMI
                [H_Ni_11 - W_Ni, H_Ni_12;
                    H_Ni_12', H_Ni_22] >= 0;

                % Second LMI
                [X_Ni_next' * H_Ni_11 * X_Ni_next - D_Ni' * ([H_Ni_11, H_Ni_12; H_Ni_12', H_Ni_22] - Lambda_Ni) * D_Ni, X_Ni_next' * H_Ni_12;
                    H_Ni_12' * X_Ni_next, H_Ni_22] >= 0;
                cvx_end

                % get K
                K_Ni = -inv(H_Ni_22) * H_Ni_12';
                disp(['Subsystem ', num2str(i), ' controller gain:']);
                disp(K_Ni);
                K_rec(i, :, k+1) = K_Ni(1, :);
            end
        end
    end
    X_rec(m, :, :) = X; % store data
    X = zeros(4, T); % refresh record arrays
    U = zeros(4, T); 

    % refresh
    X(:, 1) = x0;
end