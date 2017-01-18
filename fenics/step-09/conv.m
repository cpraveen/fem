load conv1.dat;
load conv2.dat;

figure(1)
loglog(conv1(:,1), conv1(:,3), 'o-', conv2(:,1), conv2(:,3), '*--', ...
       'LineWidth', 1.5)
xlabel('ndof', 'FontSize', 16)
ylabel('L2 error', 'FontSize', 16)
leg=legend('Uniform','Adaptive');
set(leg, 'FontSize', 16)
set(gca, 'FontSize', 16)
