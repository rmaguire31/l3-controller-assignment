function [Kp] = find_p_controller(Gc, dirout, suffix)
%% FIND_P_CONTROLLER Find Kp value of controller
%   Detailed explanation goes here

s = tf('s');
G = (s + 1)/(s*(s^2 + 4*s + 5));
if nargin >= 1
    sys = series(Gc, G);
else
    sys = G;
end
if nargin < 3
    suffix = [];
end

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
    
    for Kc = 5:1e-1:400
        % Analyse step response of compensated system
        sys0 = feedback(series(Kc, sys), 1);
        info = stepinfo(sys0);
        info.SSE= 100*abs(dcgain(sys0) - 1);
        
        % Check results against specification
        fields = fieldnames(targets);
        pass = true;
        for j = 1:length(fields)
            results.(fields{j}) = info.(fields{j});
            pass = pass & results.(fields{j}) <= targets.(fields{j});
        end
        if pass
            % Save this gain
            Kp_values{i}(end+1) = Kc;
        end
    end
    if ~isempty(Kp_values{i})
        Kp(:,i) = [min(Kp_values{i})
                   round(mean(Kp_values{i}))
                   max(Kp_values{i})];
        
        % Plot pole locations.
        f = figure();
        ax = axes(f);
        rlocusplot(ax, sys, 'k');
        ax.XLim = [-2.5 0];
        ax.YLim = [-5 5];
        hold('on');
        rlocusplot(ax, sys, Kp([1 end],i), 'r-+');
        rlocusplot(ax, sys, Kp(2,i), 'bs');
        hold('off');
        
        % Plot step response.
        f = figure();
        ax = axes(f);
        stepplot(ax, ...
            feedback(series(Kp(1,i), sys), 1), 'r--',...
            feedback(series(Kp(3,i), sys), 1), 'r--',...
            feedback(series(Kp(2,i), sys), 1), 'k');
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
        [Y,T] = step(feedback(series(Kp(1,i), sys), 1));
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

if nargin >= 2
    mkdir(dirout);
    
    % Write root locus to file for use in pgf plots
    [~,k] = rlocus(sys);
    [r,k] = rlocus(sys, sort([k Kp(:)']));
    rloc = ones(size(k,2), size(k,1) + 2*size(r,1));
    rloc(:,1:size(k,1)) = k';
    header = 'K';
    for j = 1:size(r,1)
        header = sprintf('%s\tre%d', header, j);
        rloc(:,size(k,1)+2*j-1) = real(r(j,:))';
        header = sprintf('%s\tim%d', header, j);
        rloc(:,size(k,1)+2*j) = imag(r(j,:))';
    end
    fname = [dirout '/rlocus' suffix '.dat'];
    fid = fopen(fname, 'wt');
    fprintf(fid, '%s\n', header);
    fclose(fid);
    dlmwrite(fname, rloc, '-append', 'delimiter', sprintf('\t'));

    % Write step responses to file for use in pgfplots
    t = 0:1e-2:20;
    Kc = sort([1; Kp(:)]);
    stepres = ones(length(t),1 + length(Kc));
    stepres(:,1) = t';
    header = 't';
    for j = 1:length(Kc)
        header = sprintf('%s\ty%.1f', header, Kc(j));
        y = step(feedback(series(Kc(j), sys), 1), t);
        stepres(:,1+j) = y';
    end
    fname = [dirout '/stepres' suffix '.dat'];
    fid = fopen(fname, 'wt');
    fprintf(fid, '%s\n', header);
    fclose(fid);
    dlmwrite(fname, stepres, '-append', 'delimiter', sprintf('\t'));
end
end

