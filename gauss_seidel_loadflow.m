% =========================================================================
% Gauss-Seidel Load Flow Solver
% Builds Y-bus from line data, accepts Slack/PV/PQ bus data, and
% iteratively solves node voltages with optional over-relaxation and
% PV Q-limit handling.
% =========================================================================

clc; clear; close all;

%% --- Input section: system size and topology ---
nbus  = input('Enter total number of buses: ');
nline = input('Enter total number of transmission lines: ');

fprintf('\nEnter line data [From To R X]:\n');
linedata = zeros(nline,4);
for k = 1:nline
    fprintf('Line %d: ', k);
    linedata(k,:) = input('');
end

%% --- Bus data input ---
% Bus | Type (1=Slack, 2=PV, 3=PQ) | V | delta | Pgen | Qgen | Pload | Qload | Qmin | Qmax
fprintf('\nEnter bus data [Bus Type V delta Pgen Qgen Pload Qload Qmin Qmax]:\n');
busdata = zeros(nbus,9);
for k = 1:nbus
    fprintf('Bus %d: ', k);
    busdata(k,:) = input('');
end

%% --- Forming the Y-bus (short-line model) ---
Y = zeros(nbus,nbus);
for k = 1:nline
    i = linedata(k,1); j = linedata(k,2);
    R = linedata(k,3); X = linedata(k,4);
    Z = R + 1i*X; y = 1/Z;
    Y(i,j) = -y; Y(j,i) = -y;
end
for i = 1:nbus
    for k = 1:nline
        if linedata(k,1)==i
            Z = linedata(k,3)+1i*linedata(k,4);
            Y(i,i) = Y(i,i) + 1/Z;
        elseif linedata(k,2)==i
            Z = linedata(k,3)+1i*linedata(k,4);
            Y(i,i) = Y(i,i) + 1/Z;
        end
    end
end

%% --- Initialization of voltages and powers ---
V = busdata(:,3).*exp(1i*busdata(:,4)*pi/180);
Psp = (busdata(:,5)-busdata(:,7))/100; % P specified
Qsp = (busdata(:,6)-busdata(:,8))/100; % Q specified

%% --- Convergence and relaxation parameters ---
tol = 1e-6;
maxIter = 100;
alpha = 1.6; % acceleration factor

%% --- Bus type identification ---
slack = find(busdata(:,2)==1);
pv    = find(busdata(:,2)==2);
pq    = find(busdata(:,2)==3);

%% --- Iteration loop (Gauss-Seidel core) ---
for iter = 1:maxIter
    V_prev = V;
    for i = 1:nbus
        if i == slack
            continue; % skip slack bus
        end

        % Calculate power injection
        sumYV = 0;
        for j = 1:nbus
            if j ~= i
                sumYV = sumYV + Y(i,j)*V(j);
            end
        end

        % Update voltage based on bus type
        if ismember(i,pq)  % PQ bus
            S = Psp(i) + 1i*Qsp(i);
            V_new = (1/Y(i,i))*((S/conj(V(i)))-sumYV);
            V(i) = V_new;

        elseif ismember(i,pv) % PV bus
            % Estimate reactive power first
            Qcalc = -imag(conj(V(i))*(Y(i,:)*V));
            Qsp(i) = Qcalc;

            % Q-limit check
            if Qsp(i) < busdata(i,8)/100
                Qsp(i) = busdata(i,8)/100;
                busdata(i,2) = 3; % convert to PQ
            elseif Qsp(i) > busdata(i,9)/100
                Qsp(i) = busdata(i,9)/100;
                busdata(i,2) = 3; % convert to PQ
            end

            % Update voltage (magnitude fixed)
            S = Psp(i) + 1i*Qsp(i);
            V_new = (1/Y(i,i))*((S/conj(V(i)))-sumYV);
            V(i) = busdata(i,3)*exp(1i*angle(V_new)); % keep |V| fixed
        end
    end

    % Acceleration (optional)
    V = V_prev + alpha*(V - V_prev);

    % Check convergence
    if max(abs(V - V_prev)) < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
end

%% --- Output reporting ---
disp(' ');
disp('Bus Voltages in Polar Form:');
for i = 1:nbus
    fprintf('Bus %d: |V| = %.4f pu, angle %.3f deg\n', i, abs(V(i)), angle(V(i))*180/pi);
end

fprintf('\nReactive powers at PV buses:\n');
for i = pv'
    fprintf('Bus %d: Q = %.4f pu\n', i, Qsp(i));
end
