% What happens here? An assymetric domain is added to the normal nearest neighbour interaction domain, which grows in size (e.g .
% more and more nodes interact with each other).
% The size of the sync basin should grow towards 1 (the values in c(1,i)), while all the others towards zero.
n = 20;
x = zeros([1 n]);
dt = 0.005; % size of time-step
T = 30000; % number of iterations; total time is dt*T
Test = 10;
c = zeros(5,Test+1);
NRuns = 1000; % number of runs
for lambda = 0:Test
    for Stat = 1:NRuns
        x = 2*pi*(rand(n,1)-0.5); % random configuration
        x = dyKuramoto(x,lambda,dt,T); % simulation
        q = 0;
        q = mod(pi+x(1)-x(n),2*pi)-pi;
        for i = 1:n-1
            q = q+mod(pi+x(i+1)-x(i),2*pi)-pi; % detection
        end
        q = q/(2*pi);
        for j = 0:4
            if abs(abs(q)-j)<0.4
                c(j+1,lambda+1) = c(j+1,lambda+1)+1;
            end
        end
        calt = c./NRuns;
        calt(:,lambda+1) = c(:,lambda+1)./Stat;
        disp(calt);
    end
end
disp(c);
hold 'off';
for i = 1:5
    plot(c(i,:));
    hold 'on';
end
hold 'off';

function y = dyKuramoto(x,lambda,dt,T)
    yalt = [x x];
    yneu = yalt;
    n = length(x);
    for t = 1:T
        for i = 1:n
            yneu(i) = yneu(i)+dt*(sin(yalt(i+1)-yalt(i))+sin(yalt(i+n-1)-yalt(i)));
            for j = 1:lambda
                yneu(i) = yneu(i)+dt*(sin(yalt(i+j)-yalt(i))-sin(yalt(i+n-j)-yalt(i)));
            end
        end
        yalt = [yneu(1:n) yneu(1:n)];
    end
    y = yneu(1:n);
end
