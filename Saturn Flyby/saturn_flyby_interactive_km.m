function saturn_flyby_interactive_km
    % ===========================
    % COSTANTI (in km, km^3/s^2)
    % ===========================
    Mus = 3.7931187e7;      % Saturn GM [km^3/s^2]
    MuS = 1.32712440018e11; % Sun GM [km^3/s^2]
    rs  = 58232;            % Saturn radius [km]
    Rs  = 9.537e8;          % Saturn orbit radius [km]
    vp  = sqrt(MuS/Rs);     % Saturn orbital velocity [km/s]
    r_soi = Rs*(Mus/MuS)^(2/5); % Saturn SOI radius [km]

    % Default input values
    rp0 = rs + 3000;  % periapsis radius [km]
    v_minus0 = 12;    % incoming heliocentric speed [km/s]
    delta_i0 = 30;    % [deg]

    % Create figure
    fig = figure('Name','Saturn Flyby Interactive','NumberTitle','off',...
                 'Position',[100 100 800 600]);

    % Axes
    ax = axes('Parent',fig,'Position',[0.08 0.25 0.85 0.7]);
    axis equal; hold on; grid on;
    xlabel('x [km]'); ylabel('y [km]');
    title('Saturn Flyby Trajectory');
    
    % Sliders
    slider_rp = uicontrol('Style','slider','Min',rs-1000000,'Max',rs+1000000,...
        'Value',rp0,'Units','normalized','Position',[0.15 0.15 0.7 0.03]);
    text_rp = uicontrol('Style','text','Units','normalized',...
        'Position',[0.15 0.19 0.7 0.03],...
        'String',sprintf('r_p = %.0f km',rp0));
    
    slider_v = uicontrol('Style','slider','Min',5,'Max',20,...
        'Value',v_minus0,'Units','normalized','Position',[0.15 0.10 0.7 0.03]);
    text_v = uicontrol('Style','text','Units','normalized',...
        'Position',[0.15 0.14 0.7 0.03],...
        'String',sprintf('v_{-} = %.1f km/s',v_minus0));

    slider_d = uicontrol('Style','slider','Min',-180,'Max',180,...
        'Value',delta_i0,'Units','normalized','Position',[0.15 0.05 0.7 0.03]);
    text_d = uicontrol('Style','text','Units','normalized',...
        'Position',[0.15 0.09 0.7 0.03],...
        'String',sprintf('\\delta_i = %.1f°',delta_i0));

    % Plot update function
    function update_plot(~,~)
        cla(ax); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
        xlabel(ax,'x [km]'); ylabel(ax,'y [km]');
        
        % Read sliders
        rp = get(slider_rp,'Value');
        v_minus = get(slider_v,'Value');
        delta_i = get(slider_d,'Value');

        % Update labels
        set(text_rp,'String',sprintf('r_p = %.0f km',rp));
        set(text_v,'String',sprintf('v_{-} = %.1f km/s',v_minus));
        set(text_d,'String',sprintf('\\delta_i = %.1f°',delta_i));

        % Flyby calculations (tutto in km e km/s)
        v_inf = sqrt(v_minus^2 + vp^2 - 2*v_minus*vp*cosd(delta_i));
        a_hyp = -Mus/v_inf^2;
        e_hyp = 1 + (rp*v_inf^2)/Mus;
        theta_star = acos(-1/e_hyp);
        xi = 2*theta_star - pi;

        % Hyperbolic trajectory
        theta = linspace(-theta_star, theta_star, 2000);
        r = (a_hyp*(e_hyp^2 - 1)) ./ (1 + e_hyp*cos(theta));

        % FILTRO: elimina valori non validi
        r(~isreal(r)  | r == 0) = NaN;

        x = r .* cos(theta);
        y = r .* sin(theta);

        plot(ax, x, y, 'b', 'LineWidth', 1.5); % tutto già in km


        % Periapsis
        plot(ax,-rp,0,'ro','MarkerFaceColor','r');

        % Saturn (filled circle)
        fill(rs*cos(linspace(0,2*pi,200)),...
             rs*sin(linspace(0,2*pi,200)),[0.9 0.8 0.6]);

        % SOI for context
        plot((r_soi)*cos(linspace(0,2*pi,300)),...
             (r_soi)*sin(linspace(0,2*pi,300)),'--k');

        title(ax,sprintf('Deflection angle ξ = %.1f° | v_∞ = %.2f km/s',...
            rad2deg(xi),v_inf));
    end

    % Link sliders to update function
    addlistener(slider_rp,'Value','PostSet',@update_plot);
    addlistener(slider_v,'Value','PostSet',@update_plot);
    addlistener(slider_d,'Value','PostSet',@update_plot);

    % Initial plot
    update_plot();
end
