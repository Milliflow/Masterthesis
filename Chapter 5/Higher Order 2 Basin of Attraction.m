n = 20;
x = zeros([1 n]);
dt = 0.005; % size of time-step
T = 30000; % number of runs; total time is dt*T
Test = 10;
c = zeros(6,Test+1); % chimera states in c(6,:)
NRuns = 10000; % number of runs
for lambda = 0:Test
    for Stat = 1:NRuns
        x = 2*pi*(rand(n,1)-0.5); % Random configuration
        x = dyKuramoto(x,lambda,dt,T); % Simulation
        qvec = zeros([1,n]);
        q = 0;
        qvec(n) = mod(pi+x(1)-x(n),2*pi)-pi;
        for i = 1:n-1
            qvec(i) = q+mod(pi+x(i+1)-x(i),2*pi)-pi; % Detection
        end
        q = sum(qvec)/(2*pi);

        % chimera detection
        chimera = false;
        for i = 1:n
            if abs(q-n*qvec(i)/(2*pi))>0.1
                chimera = true;
            end
        end

        % twist detection
        for j = 0:4
            if abs(abs(q)-j)<0.4 && chimera == false
                c(j+1,lambda+1) = c(j+1,lambda+1)+1;
            end
        end
        if chimera == true
            c(6,lambda+1) = c(6,lambda+1)+1;
        end
        calt = c./NRuns;
        calt(:,lambda+1) = c(:,lambda+1)./Stat;
        disp(calt);
    end
end
disp(c);
hold 'off';
for i = 1:length(c)
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
            yneu(i) = yneu(i)+dt*(sin(yalt(i+1)-yalt(i))+sin(yalt(i+n-1)-yalt(i))+lambda*(2*sin(yalt(i+1)-yalt(i))+2*sin(yalt(i+n-1)-yalt(i))+sin(2*yalt(i+1)-2*yalt(i))+sin(2*yalt(i+n-1)-2*yalt(i))));
        end
        yalt = [yneu(1:n) yneu(1:n)];
    end
    y = yneu(1:n);
end
