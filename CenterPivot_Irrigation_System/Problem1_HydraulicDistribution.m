%% Center-pivot Irrigation - Problem 1 - Hydraulic Distribution
% Author: Amir Foruzanfard
% Version: 1.0
% Copyright: © 2025, Amir Foruzanfard
% Licensed under the BSD 3-Clause License. 
% This code is licensed under the BSD 3-Clause License. See the LICENSE file for details.
% For details, visit: https://opensource.org/licenses

% Description: Hydraulic distribution by constructing and solving systems of nonlinear algebraic equations.

% clear workspace
clc; clear; close all


% pre-allocate memory for matrices
section_in = zeros(1, 5);
pipe_flows = zeros(14, 5);
headloss = zeros(14, 5);
k_factor = zeros(4,7);

% define constants
Qin = 263047/(3600*1000);                % convert from l/h to m^3/s 
inlet_pressure = 65 * 6894.757;          % convert from psig to Pa
section_inlet_pressure = inlet_pressure; % used later to keep track of the pressure between sections
pipe_length = [6 4 4 6 4.472 2.828 2.828 2.828 2.828 2.828 2.828 4.472 20 0]';          % (m)
pipe_diameter = [0.05 0.05 0.05 0.05 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.1 0.1]'; % (m)
rho = 998;              % (kg/m^3)
mu = 0.0013076;         % (Pa*s)
g = 9.81;               % (m/s^2)

S(22) = 0.02 * Qin;     % Final sprinkler uses 2% of flow

% Assume first sprinkler covers first 10m radius section. Flowrate % is 
% determined by the ratio between central section area and the total area
S(1) = 0.98 * Qin * (pi * 10^2)/(pi * 110^2); 

% create vector that contains the distance of each sprinkler from the center
A = 10:4:110;
no_sprink = 10:20:110;
sprink_pos = A(~ismember(A, no_sprink)); 

sprink_circ = sprink_pos * 2 * pi;      % length of path traveled by each sprinkler

% calculate flow of each sprinkler based on the ratio of the path length they cover
S(2:21) = (Qin - S(1) - S(22)) * (sprink_circ / sum(sprink_circ)); 


% determine the the inlet flowrate of each section
for i = 1:5
    b = (i*4)-3;                     % index of final sprinkler in i'th section
    section_in(i) = Qin-sum(S(1:b)); % Qin - sprinkler flows of prior sections 
end


for i = 1:5     % loops once for each of the 5 sections.

    count = 0;                              % iteration counter
    precision = 0.000001;                   % the maximum allowable difference between q and Q
    difference = 1;                         % any value above precision will work
    sprink_section = S((i*4)-2:(i*4)+1);    % stores the flow rates of the 4 sprinklers in the current section
    q = (ones(14, 1)/150);                  % initial guess (q = 1/150 chosen so convergence graphs look nice)
    
    figure;                         % initialise figure     
    hold on                         % hold figure so points can be added
    grid on                         
    xlabel('Iteration Count');
    ylabel('Value of Q (m^3/s)');
    title(['Iterative Evolution of Section ' num2str(i)]); % title will increment for each loop
    colours = lines(14);            % creates 14x3 matrix representing 14 distinct RGB colours
    qData = zeros(50, 14);          % create matrix to store a maximum of 50 iterations
    qData(1, :) = q;                % store inital guess in qData

    while difference > precision    % loops until (difference between q and Q) < precision
    
        Re = (4 * q * rho) ./ (pi * mu * pipe_diameter); % Reynolds number
        fric_coef = colebrook(abs(Re), 0);               % use colebrook function to get fi

        % colebrook() only accepts positive Re, so abs(Re) is used
        % sign(Re) is used when calculating K to re-introduce the sign of Re
        K = (8 * sign(Re) .* fric_coef .* pipe_length) ./ (pi^2 * g * (pipe_diameter.^5));
        D = 2 * K .* q;
        Dprime = -K .* q .^2;
        
        % refer to portfolio to see how the following system of equations was set up 
        system_matrix = [...
            1    0    0    0    1    0    0    0    0    0    0    0    1    0;...
            1   -1    0    0    0    1   -1    0    0    0    0    0    0    0;...
            0    1   -1    0    0    0    0    1   -1    0    0    0    0    0;...
            0    0    1   -1    0    0    0    0    0    1   -1    0    0    0;...
            0    0    0    1    0    0    0    0    0    0    0    1    1   -1;...
            0    0    0    0    1   -1    0    0    0    0    0    0    0    0;...
            0    0    0    0    0    0    1   -1    0    0    0    0    0    0;...
            0    0    0    0    0    0    0    0    1   -1    0    0    0    0;...
            0    0    0    0    0    0    0    0    0    0    1   -1    0    0;...
            D(1) 0    0    0 -D(5) -D(6)  0    0    0    0    0    0    0    0;...
            0   D(2)  0    0    0    0 -D(7) -D(8)  0    0    0    0    0    0;...
            0    0   D(3)  0    0    0    0    0 -D(9) -D(10) 0    0    0    0;...
            0    0    0   D(4)  0    0    0    0    0    0 -D(11) -D(12) 0   0;...
            D(1) D(2) D(3) D(4) 0    0    0    0    0    0    0    0  -D(13) 0];
        
        output_matrix = [...
            section_in(i);...
            0;...
            0;...
            0;...
            0;...
            sprink_section(1);...
            sprink_section(2);...
            sprink_section(3);...
            sprink_section(4);...
            Dprime(5)+Dprime(6)-Dprime(1);...
            Dprime(7)+Dprime(8)-Dprime(2);...
            Dprime(9)+Dprime(10)-Dprime(3);...
            Dprime(11)+Dprime(12)-Dprime(4);...
            Dprime(13)-Dprime(1)-Dprime(2)-Dprime(3)-Dprime(4)];
        
        Q = gaussPiv(system_matrix, output_matrix); % system of equations solved using Gaussian with pivoting

        difference = max(abs(Q-q));        % calculate the maximum difference between individual values of Q and q
        q = Q;
        count = count+1;                   % increment the iteration counter.
        qData(count+1, :) = q;             % stores q in qData. count+1 used because count(1) is the initial guess

        for h = 1:14
            plot(0:count, qData(1:count+1, h), '-', 'Color', colours(h, :), 'LineWidth', 1.5);
                    % plot the values of Q on the y axis and number of iterations on the x
                    % the colour of each line is set to the RGB values in colours()
        end
        drawnow;    % updates the graph after every iteration so it is drawn in real time
    end

    hold off;       % ends editing of graph so a new graph can be made for the next section

    headloss(:,i) = (rho * g * K .* q.^2); % pressure loss in each pipe is calculated (Pa)
    pipe_flows(:, i) = 1000*q;             % final values of q are stored in i'th column of final_results
    
    % pressures at the 4 sprinkler nodes calculated by subtracting headloss from section inlet pressure
    node_pressure(1:4) = [section_inlet_pressure-headloss(5,i), section_inlet_pressure-sum(headloss(5:7,i)),...
                          section_inlet_pressure-sum(headloss(5:9,i)), section_inlet_pressure-sum(headloss(5:12,i))]; 
    
    node_pressure = node_pressure./6894.757;                                    % node pressures converted to psig
    k_factor(:, i+1) = (15850.3 * S(i*4-2:i*4+1) ./ (node_pressure.^(1/2)))';   % k factors are calculated
    section_inlet_pressure = section_inlet_pressure - sum(headloss(1:4, i));    % inlet pressure of next section
end


% calculate k factor of the first and last sprinkler
k_factor(1,1) = 15850.3 * S(1) ./ (inlet_pressure/6894.757)^(1/2);        
k_factor(1,7) = 15850.3 * S(22) / (section_inlet_pressure/6894.757)^(1/2);


%% create tables to neatly display the k factors and pipe flows with headings

pipe_names = compose("P%d", 1:14);              % create array containing pipe labels P1-P14
sprink_names = compose("S%d", 1:4);             % create array containing sprinkler labels S1-S4
section_names = compose("Section %d", 1:5);     % create array containing section labels

% use array2table() to convert pipe flows and k factor matrices into tables with the headings created above
pipe_flow_table = array2table(pipe_flows, 'VariableNames',section_names,'RowNames',pipe_names);
k_factor_table = array2table(k_factor, 'VariableNames', ['First Sprinkler',section_names,'Final Sprinkler'],...
'RowNames',sprink_names);

format short g
display(pipe_flow_table)       % output the table of pipe flowrates
display(k_factor_table)        % output the table of nozzle k factors
