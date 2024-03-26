function movie(Y, L1, L2, titleText)
    x1=L1*sin(Y(:,1));
    y1=-L1*cos(Y(:,1));
    x2=x1+L2*sin(Y(:,3));
    y2=y1-L2*cos(Y(:,3));
    
    frames = [];
    
    % Create a figure for the animation
    animationFig = figure;
    
    % Set manual axes limits for animation
    animationLimits = [-3 3 -3 1];  % Adjust as needed
    
    % Set initial axes limits
    axis(animationLimits);
    
    for i = 1:length(Y)
        % Update the position of the pendulum bobs and rods in the animation plot
        plot([0, x1(i)], [0, y1(i)], 'b', 'LineWidth', 2);
        hold on;
        plot(x1(i), y1(i), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');  % Mass for Pendulum 1
        plot([x1(i), x2(i)], [y1(i), y2(i)], 'r', 'LineWidth', 2);  % Rod for Pendulum 2
        plot(x2(i), y2(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Mass for Pendulum 2
        
        xlabel('X Position (m)');
        ylabel('Y Position (m)');
        title(titleText)
        
        % Set initial y-axis limits
        ylim(animationLimits(3:4));
        xlim(animationLimits(1:2));
    
        % Capture the current frame for the animation
        frame = getframe(animationFig);
        frames = [frames, frame];
    
        % Clear the previous frame in the animation plot
        if i < length(Y)
            cla(animationFig);
        end
    end
    close(animationFig);
end