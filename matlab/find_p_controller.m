function [Kp] = find_p_controller()
%% FIND_P_CONTROLLER Find Kp value of controller
%   Detailed explanation goes here

s = tf('s');
G = (s + 1)/(s*(s^2 + 4*s + 5));

% Perfomance spec
spec = struct(...
    'RiseTime',     {1,0.5}, ... seconds
    'SettlingTime', {4,nan}, ... seconds
    'SettlingMin',  {nan,nan}, ...
    'SettlingMax',  {nan,nan}, ...
    'Overshoot',    {5,20}, ... percent 
    'Undershoot',   {nan,nan}, ... percent
    'Peak',         {nan,nan}, ...
    'PeakTime',     {nan,nan}, ... seconds
    'SSE',          {nan,5} ... percent
);

Kp = nan(3, length(spec));
Kp_values = cell(size(spec));
for i = 1:length(spec)    
    % Construct targets for this controller
    results = struct();
    targets = struct();
    fields = fieldnames(spec);
    for j = 1:length(fields)
        if ~isnan(spec(i).(fields{j}))
            targets.(fields{j}) = spec(i).(fields{j});
        end
    end
    
    for Gc = 5:1e-1:20
        % Analyse step response of compensated system
        sys = feedback(series(Gc, G), 1);
        info = stepinfo(sys);
        info.SSE= 100*abs(dcgain(sys) - 1);
        
        % Check results against specification
        fields = fieldnames(targets);
        pass = true;
        for j = 1:length(fields)
            results.(fields{j}) = info.(fields{j});
            pass = pass & results.(fields{j}) <= targets.(fields{j});
        end
        if pass
            % Save this gain
            Kp_values{i}(end+1) = Gc;
        end
    end
    if ~isempty(Kp_values{i})
        Kp(:,i) = [min(Kp_values{i})
                   round(mean(Kp_values{i}))
                   max(Kp_values{i})];
        
        % Plot pole locations.
        f = figure();
        ax = axes(f);
        rlocusplot(ax, G, 'k');
        ax.XLim = [-2.5 0];
        ax.YLim = [-5 5];
        hold('on');
        rlocusplot(ax, G, Kp([1 end],i), 'r-+');
        rlocusplot(ax, G, Kp(2,i), 'bs');
        hold('off');
        
        % Plot step response.
        f = figure();
        ax = axes(f);
        stepplot(ax, ...
            feedback(series(Kp(1,i), G), 1), 'r--',...
            feedback(series(Kp(3,i), G), 1), 'r--',...
            feedback(series(Kp(2,i), G), 1), 'k');
        ax.YLim = [0 1.4];
        hold('on');
        
        % Upper envelope
        ylim = 1 + spec(i).Overshoot/100;
        y = [ylim
             ylim];
        t = [ax.XLim(1)-1
             ax.XLim(end)+1];
        if ~isnan(spec(i).SettlingTime)
            if spec(i).Overshoot > 2 
                y = [y
                     1.02
                     1.02];
                t = [ax.XLim(1)-1
                     spec(i).SettlingTime
                     spec(i).SettlingTime
                     ax.XLim(end)+1];
            end
        end
        area([ax.XLim(1)-1 ax.XLim(end)+1],...
            [ax.YLim(end)+1 ax.YLim(end)+1],...
            'FaceColor', [0.9 0.9 0.5],...
            'BaseValue', ax.YLim(1)-1,...
            'LineStyle', '-.');
        area(t, y,...
            'FaceColor', 'w',...
            'BaseValue', ax.YLim(1)-1,...
            'LineStyle', '-.');
        
        % Lower envelope
        y = [0.98
             0.98];
        t = [spec(i).SettlingTime
             ax.XLim(end)+1];
        area(t, y,...
            'FaceColor', [0.9 0.9 0.5],...
            'BaseValue', ax.YLim(1)-1,...
            'LineStyle', '-.');
        
        % Rise time condition
        [Y,T] = step(feedback(series(Kp(1,i), G), 1));
        t01 = interp1(Y,T,0.1);
        t02 = t01 + spec(i).RiseTime;
        y = [0.1 0.9
             0.0 0.0];
        t = [t01 t02
             t01 t02];
        plot(t, y, 'k-.');
        scatter(t(1,:), y(1,:), 15, 'k', 'filled');
       
        hold('off');
    end
end
end

