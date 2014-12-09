file = fopen('ground_truth_score_i.txt', 'r');
A  = fscanf(file, '%d\t%f', [2 Inf]);
y0 = A(2,:);
x  = 1:(size(y0, 2));


file = fopen('naive_score_i.txt', 'r');
A  = fscanf(file, '%d\t%f', [2 Inf]);
y1 = A(2,:);

file = fopen('pl_score_i.txt', 'r');
A  = fscanf(file, '%d\t%f', [2 Inf]);
y2 = A(2,:);

figure();
grid on
hold on
plot(x, y0, '-r', 'LineWidth', 1);
plot(x, y1, '-b', 'LineWidth', 1);
plot(x, y2, '-m', 'LineWidth', 1);
plot(x, 0.5 * ones(size(x,2)), '--k', 'LineWidth', 1);
%plot(x, y2, '-r', 'LineWidth', 1);
legend('Ground-truth', 'Naive','PLM');
xlabel('Ranking', 'FontSize', 20);
ylabel('Estimated Score', 'FontSize', 20);
set(gca, 'FontSize', 20);


figure();
grid on
hold on
% inds = 1:20; %find(y0>0.5)
% y0 = y0(inds)
% y1 = y1(inds)
% y2 = y2(inds)
% x  = 1:(size(y0, 2));
plot(x, log(y1) - log(y0) - (log(y1(1)) - log(y0(1))), '-b', 'LineWidth', 1);
plot(x, log(y2) - log(y0) - (log(y2(1)) - log(y0(1))), '-m', 'LineWidth', 1);
%plot(x, y2, '-r', 'LineWidth', 1);
%legend('Naive','PLM');
xlabel('Ranking', 'FontSize', 20);
ylabel('Estimated Score Ratio', 'FontSize', 20);
set(gca, 'FontSize', 20);