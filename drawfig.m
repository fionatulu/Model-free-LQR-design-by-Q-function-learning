x_axis = 1:25;
%mean_X_rec = squeeze(mean(X_rec, 1));
for j = 1:N
    %y_axis = mean_X_rec(j, :);
    y_axis = squeeze(X_rec(1, j, :));
    plot(x_axis, y_axis);
    hold on
end

hold off
legend('x1', 'x2', 'x3', 'x4');