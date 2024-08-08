function ellipse_stats = process_ellipse(tvec, N , x_all, y_all, x0, y0, left_wall_list, A)

    [~,isort] = sort(x0);
    iskip = 1;

    ellipse_stats = zeros(length(isort(1:iskip:end)), 4);
    nn_list = isort(1:iskip:N);
    
    % random_plot_sample = randsample(nn_list, 10);
    % debug_counter = 1;
    for nn_counter = 1:length(nn_list)
        % waitbar(nn_counter/length(nn_list)); % COMMENTED OUT BECAUSE: On hamming waitbar not supported under the -nojvm startup option.
        nn = nn_list(nn_counter);
        if(~left_wall_list(nn))
            x_temp = x_all(nn,:);
            y_temp = y_all(nn,:);
            if length(unique(x_temp))>100
    
                % find index when particles rises to .5 its max amplitude
                i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');
    
                % Calculate the amplitude as the Euclidean distance from the initial position (x0(nn), y0(nn)).
                amp = sqrt((x_temp-x0(nn)).^2+(y_temp-y0(nn)).^2);
    
                % Calculate the angle from the initial position to each point in Cartesian coordinates.
                ang = atan2((y_temp-y0(nn)),(x_temp-x0(nn)));
                t_temp = tvec(i0:end);
    
                medamp = median(amp);
                t = 1:length(amp);
    
                % Normalize the x and y coordinates by subtracting the initial position and scaling by A.
                X = (x_temp(100:end)-x0(nn))/A; % why does this start at 100? Should it start at t_temp then plot t_temp, X? Aka time we mark as the rise above .5 max amplitude? Otherise you get a long period of no oscillation
                Y = (y_temp(100:end)-y0(nn))/A;
    
                % Why are these set to their same values? 
                X = X; 
                Y = Y;
    
                % Compute the radius squared for each point relative to the origin.
                R = X.^2+Y.^2;
    
                % Find the maximum radius squared from the array R. This should become the semi-major axis
                R_max = max(R);
    
                % Find the index of the first occurrence where R equals R_max.
                max_ind = find(R == R_max,1,"first");
    
                % Extract the X and Y coordinates corresponding to the maximum radius.
                X_max = X(max_ind);
                Y_max = Y(max_ind);
                af = [];
    
                % *** NEEDS fit_ellipse FUNCTION *** Fit an ellipse to the latter half of the data points, returning a structure with ellipse parameters.
                ellipse_t = fit_ellipse(X(round(end/2):end)',Y(round(end/2):end)');
                % test_x = (x_temp(i0:end)-x0(nn))/A;
                % test_y = (y_temp(i0:end)-y0(nn))/A;
                % ellipse_t = fit_ellipse(test_x',test_y');
                % nn_counter
    
                % Check if the ellipse fitting structure 'ellipse_t' is not empty and contains  parameter 'a' that comes from the fit_ellipse function
                if ~isempty(ellipse_t) && ~isempty(ellipse_t.a)
                    % if nn_counter==336
                    
                    % end
    
                    % Assign the semi-major axis, semi-minor axis, and negative rotation angle from 'ellipse_t' to 'af'
                    af(1) = ellipse_t.a; % sub axis (radius) of the X axis of the non-tilt ellipse
                    af(2) = ellipse_t.b; % sub axis (radius) of the Y axis of the non-tilt ellipse
                    af(3) = -ellipse_t.phi; % orientation in radians of the ellipse (tilt)
                else
                    % Define initial parameter estimates for fitting based on maximum radius and the angle at maximum point
                    a0 = [R_max R_max atan(Y_max/X_max)];
                    c = [0 0]; % Center of transformation set to origin
    
                    % Define the function 'f' to minimize an ellipse fit error function parameterized by 'a'
                    f = @(a) ((((X-c(1))*cos(a(3))+(Y-c(2))*sin(a(3))).^2)/a(1).^2 + (((X-c(1))*sin(a(3))-(Y-c(2))*cos(a(3))).^2)/a(2).^2 -1);
                    %[af fval] = fminsearch(f,a0);
    
                    % Use nonlinear least squares to optimize parameters 'a' starting from 'a0'
                    af = lsqnonlin(f, a0);%, [], [], options);
                end
    
                if ~isempty(af)

                    % Make sure to consider only the abs of the ellipse. for the next step
                    af(1:2) = abs(af(1:2));
    
                    % If the semi-major axis (af(1)) is less than the semi-minor axis (af(2))
                    if af(1)<af(2)
    
                        % Rotate the ellipse by -pi/2 to correct the orientation
                        af(3) = af(3)-pi/2;
    
                        % Swap the semi-major and semi-minor axes to maintain convention
                        temp = af(2);
                        af(2) = af(1);
                        af(1) = temp;
                    end
    
                    % Normalize the rotation angle af(3) within the range -pi/2 to pi/2
                    af(3) = (mod(af(3)+pi/2,pi)-pi/2);
    
    
                    % Store the fitted ellipse parameters (semi-major axis, semi-minor axis, rotation angle, initial position from oscillation) in the ellipse_stats matrix for the current particle.
                    ellipse_stats(nn_counter,:) = [af, x0(nn)];

                    % debug_aspect_ratio = af(2)./af(1)

                    % if ismember(nn_counter, random_plot_sample)
                    % if debug_aspect_ratio > 1
                    %     debug_counter = debug_counter +1
                    %     % if debug_counter>5
                    %     %     return
                    %     % end
                    %     figure(nn_counter), clf, plot((x_temp(100:20:end)-x0(nn))/A,(y_temp(100:20:end)-y0(nn))/A,'ro'), axis equal, hold on
                    %     % plot(X_max,Y_max,'ro','markerfacecolor','r','markersize',20)
                    %     debug_counter = debug_counter +1
                    %     % plot(c(1), c(2), 'r*')
                    %     axis equal
                    %     theta_rot = af(3);
                    %     v = linspace(0, 1);
                    %     rx = af(1);
                    %     ry = af(2);
                    %     x = rx*cos(2*pi*v);
                    %     y = ry*sin(2*pi*v);
                    %     xx = x*cos(theta_rot) - y*sin(theta_rot);
                    %     yy = x*sin(theta_rot) + y*cos(theta_rot);
                    %     % figure
                    %     plot(xx, yy,'k-','linewidth',2)
                    %     plot([0 af(1)*cos(af(3))],[0 af(1)*sin(af(3))],'k-')
                    %     plot([0 -af(2)*sin(af(3))],[0 af(2)*cos(af(3))],'k-')
                    %     plot([0 af(1)*cos(af(3))],[0 0],'k--')
                    %     text(0.4, 0.1, '$a$','Interpreter','latex','fontsize',14)
                    %     text(-0.1, 0.02, '$b$','Interpreter','latex','fontsize',14)
                    %     text(0.8,0, '$\theta$','Interpreter','latex','fontsize',14)
                        
                    %     xlim([1.05*min((x_temp(100:20:end)-x0(nn))/A), 1.05*max((x_temp(100:20:end)-x0(nn))/A)])
                    %     ylim([1.05*min((y_temp(100:20:end)-y0(nn))/A), 1.05*max((y_temp(100:20:end)-y0(nn))/A)])
                    %     xlabel('$\Delta x / A$','Interpreter','latex','fontsize',18)
                    %     ylabel('$\Delta y / A$','Interpreter','latex','fontsize',18)
                    %     title(['$\theta=$' num2str(af(3))], "Interpreter", "latex")
                    %     % grid
                    %     % axis([-1  1    -1  1])
                    %     % axis('equal')
                    %     % t=0:pi/10:2*pi;
                    %     % plot(c(1) + af(1)*cos(t), c(2) + af(2)*sin(t), 'r')
                    %     % axis([-2 2 -2 2])

                    % end
                end
            end
        end
    end
end
