load conv1.dat;
load conv2.dat;

figure(1)
loglog(conv1(:,2), conv1(:,3), 'o-', conv2(:,2), conv2(:,3), '*-')
xlabel('h')
ylabel('L2 error')
legend('P1','P2')

figure(2)
loglog(conv1(:,1), conv1(:,3), 'o-', conv2(:,1), conv2(:,3), '*-')
xlabel('ndof')
ylabel('L2 error')
legend('P1','P2')
