%% Center-pivot Irrigation - Problem 2 - SprinklerCamData_CurveFitting
% Author: Amir Foruzanfard
% Version: 1.0
% Copyright: © 2025, Amir Foruzanfard
% Licensed under the BSD 3-Clause License. 
% This code is licensed under the BSD 3-Clause License. See the LICENSE file for details.
% For details, visit: https://opensource.org/licenses

% Description: Curve fitting to the sprinkler modifiers and camera theta angles data.

% clear workspace
clc; clear; close all

load('SprinklerModifier_Theta_datafile_problem2.mat');  % load data file
headloss = zeros(1000, 14);                    % preallocate variable size
pipe_flows = zeros(1000, 14);                  % preallocate variable size
smooth_mods = zeros(size(Sprinkler_modifier)); % preallocate variable size

%%
tic             % start timer
for i = 1:4     % repeat once for each sprinkler
    
    % fit() function uses the curve fitting toolbox to fit a smoothing spline to the modifier data
    % Different parameters were attempted manually in the curve fitter first to find the most suitable ones
    fit_result = fit(Theta', Sprinkler_modifier(:, i), 'smoothingspline', 'SmoothingParam', 0.001);

    % smoothed modifiers at each point in theta are stored in the ith column
    smooth_mods(:, i) = fit_result(Theta'); 
    end

for i = 1:1000  
    results = Problem1FlowDistribFunc(smooth_mods(i, :));  % i'th row of smoothed modifiers input into Problem1FlowDistribFunc()
    pipe_flows(i,:) = results(1,:);                % the two Problem1FlowDistribFunc outputs stored in seperate variables
    headloss(i,:) = results(2,:);
end
toc             % end timer


%% initialise figure and maximise window to fill the screen
figure('WindowState', 'maximized');

% create first panel, give it a title, and define its position from edges of figure window
panel1 = uipanel('Title', 'Volume Flow Rate', 'FontSize', 12, 'Position', [0.01 0.05 0.48 0.95]);

% create a subplot in the first panel and define its position from edges
subplot('Position', [0.1, 0.7, 0.87, 0.25], 'Parent', panel1)
hold on;                                             % hold the plot to overlay data
plot(Theta, pipe_flows(:,13));                       % plot pipe flows of P13 against theta
plot(Theta, pipe_flows(:,14));                       % plot pipe flows of P14 against theta
hold off;
ylabel('Flow Rate (l/s)'); 
title('Bus Pipes');
legend('P13', 'P14', 'Location','northeastoutside'); % create a legend and define its location
xlim([0, 360]);                                      % set x-axis limit of the plot
grid on

% create another subplot in the first panel and define its position from edges
subplot('Position', [0.1, 0.375, 0.87, 0.25], 'Parent', panel1)
hold on;
for i = 1:4                 
    plot(Theta, pipe_flows(:,i));                    % plot flows of P1-P4 against theta
end
hold off;
ylabel('Flow Rate (l/s)');
title('Horizontal Pipes');
% create a legend, create data labels using compose() and define its location
legend(compose('P%d', 1:4), 'Location','northeastoutside');
xlim([0, 360]);                             % set x-axis limit of the plot
grid on

% create third subplot in the first panel and define its position from edges
subplot('Position', [0.1, 0.05, 0.87, 0.25], 'Parent', panel1)
hold on;

for i = 5:12                                % iterate through pipes 5 to 12
    if i == 12
        % use dashed line for pipe 12 as there are only 6 colours in matlab
        plot(Theta, pipe_flows(:, i), '--');  
    else
        plot(Theta, pipe_flows(:, i));
    end
end
hold off;
ylabel('Flow Rate (l/s)');
title('Sprinkler Pipes');

% create a legend, create data labels using compose() and define its location
legend(compose('P%d', 5:12), 'Location','northeastoutside');
xlim([0, 360]);
grid on                                     % set x-axis limit of the plot

% create second panel, give it a title, and define its position from edges of figure window
panel2 = uipanel('Title', 'Head Loss', 'FontSize', 12, 'Position', [0.5 0.05 0.48 0.95]);

% plotting headloss() instead of pipe_flows()
subplot('Position', [0.1, 0.7, 0.87, 0.25], 'Parent', panel2)
hold on;
plot(Theta, headloss(:,13));
plot(Theta, headloss(:,14));
hold off;
ylabel('Head-Loss (m)');
title('Bus Pipes');
legend('P13', 'P14', 'Location','northeastoutside');
xlim([0, 360]);
grid on
