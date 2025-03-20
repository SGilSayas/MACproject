%% Figure settings
fi = gcf;
fi.InvertHardcopy = 'off';
fi.Color = 'white';
fi.PaperPositionMode = 'auto';



%% Axes settings
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 11;
ax.Color    = [0.9,0.9,0.9];

%% Grid
ax.GridAlpha  = 1;
ax.GridColor  = 'white';
ax.XGrid      = 'on';
ax.YGrid      = 'on';
drawnow
ax.XGridHandle.LineWidth = 1; % Undocumented!
ax.YGridHandle.LineWidth = 1; % Undocumented!

%% Minor
if strcmpi(ax.XMinorGrid,'on');
    ax.MinorGridAlpha = 1;
    ax.MinorGridColor = 'white';
    ax.XGridHandle.MinorGridLineStyle = '-';
end
ax.GridAlpha  = 1;
ax.GridColor  = 'white';
ax.XGrid      = 'on';
ax.YGrid      = 'on';
drawnow
%ax.XGridHandle.LineWidth = 1; % Undocumented!
%ax.YGridHandle.LineWidth = 1; % Undocumented!

%% Rulers
ax.XRuler.Color = [0.4,0.4,0.4]; % Undocumented!
ax.YRuler.Color = [0.4,0.4,0.4]; % Undocumented!
ax.YRuler.Axle.Visible = 'off';  % Undocumented!
ax.XRuler.Axle.Visible = 'off';  % Undocumented!
ax.XLabel.Color = 'black';
ax.YLabel.Color = 'black';

%% Remove major ticks
try % Old MATLAB
    ax.XRuler.MajorTicks.Visible = 'off'; % Undocumented!
    ax.YRuler.MajorTicks.Visible = 'off'; % Undocumented!
catch 
    try % Newer MATLABs
        ax.XAxis.MajorTickChild.Visible = 'off'; % Undocumented!
        ax.YAxis.MajorTickChild.Visible = 'off'; % Undocumented!
    catch
        ax.TickLength = [0 0]; % Works but leaves a small mark
    end
end

%% Remove minor ticks
try % Old MATLAB
    ax.XRuler.MinorTicks.Visible = 'off'; % Undocumented!
catch 
    try % Newer MATLABs
        ax.XAxis.MinorTickChild.Visible = 'off'; % Undocumented!
    catch
        ax.TickLength = [0 0]; % Works but leaves a small mark
    end
end
try % Old MATLAB
    ax.YRuler.MinorTicks.Visible = 'off'; % Undocumented!
catch 
    try % Newer MATLABs
        ax.YAxis.MinorTickChild.Visible = 'off'; % Undocumented!
    catch
        ax.TickLength = [0 0]; % Works but leaves a small mark
    end
end