% =========================================================================
% Gauss-Seidel Load Flow Solver -- DEMO / NON-INTERACTIVE VERSION
%
% Same solver logic as gauss_seidel_loadflow.m, but with the interactive
% input() prompts replaced by a hardcoded sample 3-bus test system, so it
% can run unattended (e.g. from the command line / CI) and produce
% repeatable output files for the repo.
%
% Sample system (100 MVA base):
%   Bus 1: Slack, V = 1.05 pu, delta = 0 deg
%   Bus 2: PV,    V = 1.00 pu, Pgen = 200 MW, Q limits [-100, 300] MVAr
%   Bus 3: PQ,    Pload = 200 MW, Qload = 100 MVAr
%
%   Lines [From To R X] (pu):
%     1-2: R=0.02   X=0.04
%     1-3: R=0.01   X=0.03
%     2-3: R=0.0125 X=0.025
% =========================================================================

clc; clear; close all;

diary('gauss_seidel_output.txt');
diary on;

%% --- System size and topology (hardcoded sample data) ---
nbus  = 3;
nline = 3;

fprintf('\nLine data [From To R X]:\n');
linedata = [ ...
    1 2 0.0200 0.0400; ...
    1 3 0.0100 0.0300; ...
    2 3 0.0125 0.0250];
disp(linedata);

%% --- Bus data ---
% Bus | Type (1=Slack, 2=PV, 3=PQ) | V | delta | Pgen | Qgen | Pload | Qload | Qmin | Qmax
fprintf('\nBus data [Bus Type V delta Pgen Qgen Pload Qload Qmin Qmax]:\n');
busdata = [ ...
    1 1 1.05 0   0    0  0    0    0     0;   ...
    2 2 1.00 0   200  0  0    0    -100  300; ...
    3 3 1.00 0   0    0  200  100  0     0];
disp(busdata);

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

fprintf('\nY-bus matrix:\n');
disp(Y);

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
iterHistory = zeros(maxIter,1);
converged = false;
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

    maxDelta = max(abs(V - V_prev));
    iterHistory(iter) = maxDelta;

    % Check convergence
    if maxDelta < tol
        fprintf('Converged in %d iterations.\n', iter);
        converged = true;
        break;
    end
end
iterHistory = iterHistory(1:iter);

if ~converged
    fprintf('Did NOT converge within %d iterations (last max change = %.3e).\n', maxIter, iterHistory(end));
end

%% --- Output reporting ---
disp(' ');
disp('Bus Voltages in Polar Form:');
Vmag = zeros(nbus,1);
Vang = zeros(nbus,1);
for i = 1:nbus
    Vmag(i) = abs(V(i));
    Vang(i) = angle(V(i))*180/pi;
    fprintf('Bus %d: |V| = %.4f pu, angle %.3f deg\n', i, Vmag(i), Vang(i));
end

fprintf('\nReactive powers at PV buses:\n');
for i = pv'
    fprintf('Bus %d: Q = %.4f pu\n', i, Qsp(i));
end

diary off;

%% --- Save results to files for the repo ---
% 1) Bus results as CSV
fid = fopen('bus_results.csv','w');
fprintf(fid, 'Bus,Type,Vmag_pu,Vangle_deg\n');
typeNames = {'Slack','PV','PQ'};
for i = 1:nbus
    fprintf(fid, '%d,%s,%.6f,%.6f\n', i, typeNames{busdata(i,2)}, Vmag(i), Vang(i));
end
fclose(fid);

% 2) Iteration convergence history as CSV
fid = fopen('convergence_history.csv','w');
fprintf(fid, 'Iteration,MaxVoltageChange\n');
for k = 1:length(iterHistory)
    fprintf(fid, '%d,%.10e\n', k, iterHistory(k));
end
fclose(fid);

% 3) Full workspace as a .mat file
save('gauss_seidel_results.mat', 'V', 'Vmag', 'Vang', 'Y', 'busdata', ...
     'linedata', 'iterHistory', 'converged', 'iter');

fprintf('\nSaved: gauss_seidel_output.txt, bus_results.csv, convergence_history.csv, gauss_seidel_results.mat\n');
