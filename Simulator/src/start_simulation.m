% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

close all
clear all 
clc

%% LOAD DATA

run('config.m');

%% START THE CHOSEN SIMULATION
% T = vector of time used by ODE, [s] also for Tf Ta
% Y = State = ( x y z | u v w | p q r | q0 q1 q2 q3 ) also for Ya,Yf corresponding to T

tic
    
% Checking if stochastic or standard simulation needed
if settings.ballistic
    if settings.stoch.N > 1
        if settings.wind.model 
            fprintf('Stochastic Ballistic Simulation With Wind Model Started...\n\n');
        else
            fprintf('Stochastic Ballistic Simulation With Random Wind Started...\n\n');
        end
            
        [LP,X,ApoTime] = stoch_run_bal(settings);
    else
        fprintf('Standard Ballistic Simulation Started...\n\n');
        [T,Y,Ta,Ya,bound_value] = std_run_ballistic(settings);
    end
else
    if settings.stoch.N > 1
        if settings.wind.model 
            fprintf('Stochastic Simulation With Wind Model Started...\n\n');
        else
            fprintf('Stochastic Simulation With Random Wind Started...\n\n');
        end
        [LP,X,ApoTime] = stoch_run(settings);
    else
        fprintf('Standard Simulation Started...\n\n');
        [T,Y,Ta,Ya,bound_value] = std_run(settings);
    end
end

toc

%% NOT-STOCHASTIC SIMULATIONS (N=1)

if settings.stoch.N == 1 
    
    N = length(Y(:,1));
    n = length(Ya(:,1));
    
    % POSITIONS
    
    x = Y(:,1);
    y = Y(:,2);
    z = -Y(:,3);
    X = [x, y, z];
    
    z_a  = z(1:n);
    
    [apogee, i_apogee] = max(-Y(:,3)); % position, index of position at apogee
    
    
    % VELOCITIES
    
    
    u = Y(:,4);
    v = Y(:,5);
    w = -Y(:,6);
    V = [u, v, w];
    
    % ACCELERATIONS
    
    % main derivatives
    ax = (u(3:N)-u(1:N-2))./(T(3:N)-T(1:N-2));
    ay = (v(3:N)-v(1:N-2))./(T(3:N)-T(1:N-2));
    az = (w(3:N)-w(1:N-2))./(T(3:N)-T(1:N-2));
    
    % add derivative at the boundaries
    ax = [u(2)/T(2); ax; (u(end)-u(end-1))/(T(end)-T(end-1))];
    ay = [v(2)/T(2); ay; (v(end)-v(end-1))/(T(end)-T(end-1))];
    az = [w(2)/T(2); az; (w(end)-w(end-1))/(T(end)-T(end-1))];
    
    A = [ax, ay, az];
    
    clear('ax', 'ay', 'az', 'vx', 'vy')
    
    % MAXIMUM POSITIONS, VELOCITIES AND ACCELERATIONS
    
    % pre-allocation
    abs_X = zeros(length(X), 1);
    abs_V = abs_X;
    abs_A = abs_X;
    
    % determine the module of every row element
    for k = 1:length(X)
        abs_X(k) = norm(X(k,:));
        abs_V(k) = norm(V(k,:));
        abs_A(k) = norm(A(k,:));
    end
    
    abs_Va = abs_V(1:n);
    abs_Aa = abs_A(1:n);
    
    
    [max_dist, imax_dist] = max(abs_X);
    [max_v, imax_v] = max(abs_V);
    [max_a, imax_a] = max(abs_A);
    iexit = find(abs_X <= settings.lrampa);  % checking where the missile is undocked from the hook of the launch pad
    iexit = iexit(end);
    
    % TEMPERATURE AND MACH NUMBER
    
    % pre-allocation
    M = zeros(length(X), 1);
    Tamb = M;
    Ttot = Tamb;
    for n = 1:length(M)
        h = z(n);
        if (h < 0)
            h = 0;
        end
        [Tamb(n), a, ~, ~] = atmosisa(h);
        M(n) = abs_V(n)/a;
        Ttot(n) = Tamb(n)*(1+0.2*M(n)^2);
    end
    % determine the maximum Mach number
    [max_M, imax_M] = max(M);
    
    % DATA RECORD (display)
    
    disp(' ')
    disp('DATA RECORD:')
    fprintf('apogee reached: %g [m] \n', apogee);
    fprintf('apogee velocity: %g [m/s] \n', abs_V(i_apogee));
    fprintf('time: %g [sec] \n\n', Ta(end))
    
    fprintf('max speed reached: %g [m/s] \n', max_v)
    fprintf('altitude: %g [m] \n', z(imax_v))
    fprintf('Mach: %g [-] \n', M(imax_v))
    fprintf('time: %g [sec] \n\n', T(imax_v))
    
    fprintf('max Mach reached: %g [-] \n', max_M)
    fprintf('altitude: %g [m] \n', z(imax_M))
    fprintf('veocity: %g [m/s] \n', abs_V(imax_M))
    fprintf('time: %g [sec] \n\n', T(imax_M))
    
    fprintf('max acceleration reached: %g [m/s2] = %g [g] \n', max_a, max_a/9.80665)
    fprintf('veocity: %g [m/s] \n', abs_V(imax_a))
    fprintf('time: %g [sec] \n\n', T(imax_a))
    
    fprintf('run on launch pad: %g [m] \n', abs_X(iexit))
    fprintf('speed at launch pad exit: %g [m/s] \n', abs_V(iexit))
    fprintf('time: %g [sec] \n\n', T(iexit))
    
    
%% STOCHASTIC SIMULATIONS (N>1)

else
    
    % CHECKING BAD SIMULATION
    
    if numel(LP(LP(:,3) < -10,:))
        fprintf(['Some landing points might be incorrect ' ...
            'please check parameters!\n']);
    end
    
    % PRINTING VALUES
   
    % Mean Apogee Time
    ApoTimem = mean(ApoTime);
    
    % Std. Deviation Apogee Time
    ApoTimestd = std(ApoTime);
    
    % Std. Deviation Altitude
    zstd = std(X(:,3));
    
    % Mean Apogee Points
    xapom = mean(X(:,1));
    yapom = mean(X(:,2));
    zapom = mean(X(:,3));
    
    if settings.ao
        text =['Mean Apogee Point:X:%3.3f m, Y:%3.3f m\n',...
            'Mean Altitude: %3.3f m || STD: %3.3f m\n',...
            'Mean Apogee Time: %3.3f s || STD: %3.3f s\n'];
        fprintf(text,xapom,yapom,zapom,zstd,ApoTimem,ApoTimestd);
        
    else
        xm = mean(LP(:,1));
        ym = mean(LP(:,2));
        text =['Mean Landing Point:X:%3.3f m, Y:%3.3f m\n',...
            'Mean Altitude: %3.3f m || STD: %3.3f m\n',...
            'Mean Apogee Time: %3.3f s || STD: %3.3f s\n'];
        fprintf(text,xm,ym,zapom,zstd,ApoTimem,ApoTimestd);
    end
    
    delete(gcp('nocreate'))
    delete('parfor_progress.txt')
    
end


%% PLOTS

if settings.plots
    run('plots.m')
end

if settings.stoch.N == 1
    delete('ascent_plot.mat')    
end

clearvars -except ascent descent_bal descent_para T Ta Y Ya
