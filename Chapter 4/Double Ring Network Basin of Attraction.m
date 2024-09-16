n = 20;
x = zeros(n,2);
dt = 0.005; % size of time-step
T = 30000; % number of runs; total time is dt*T
Test = 10;
c = zeros(5,Test);
NRuns = 10000; % number of runs
for lambda = 1:Test % for lambda = 1 all rings will be independent -> random twists -> programm not suited for this
    for Stat = 1:NRuns
        x = 2*pi*(rand(n,2)-0.5); % random configuration
        x = dyKuramoto(x,lambda,dt,T); % simulation
        q = zeros([0,2]);
        for j = 1:2
            q(j) = mod(pi+x(1,j)-x(n,j),2*pi)-pi;
            for i = 1:n-1
                q(j) = q(j)+mod(pi+x(i+1,j)-x(i,j),2*pi)-pi; % detection
            end
        end
        q = q./(2*pi);
        for j = 0:4
            qtrue = true;
            for i = 1:2
                if abs(abs(q(i))-j)>0.4
                    qtrue = false;
                end
            end
            if qtrue
                c(j+1,lambda) = c(j+1,lambda)+1;
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
    sz = size(x);
    n = sz(1);
    yalt = [x; x];
    yneu = yalt;
    for t = 1:T
        for i = 1:n
            for j = 1:2
                yneu(i,j) = yneu(i,j)+dt*(sin(yalt(i+1,j)-yalt(i,j))+sin(yalt(i+n-1,j)-yalt(i,j)));
            end
            for j = 1:lambda
                for l = 1:2
                    for m = 1:2
                        if m ~= l
                            yneu(i,l) = yneu(i,l)+dt*(sin(yalt(i+j,m)-yalt(i,l))+sin(yalt(i+n-j,m)-yalt(i,l)));
                        end
                    end
                end
            end                                                 
        end
        for i = 1:2
            yalt(:,i) = [yneu(1:n,i); yneu(1:n,i)];
        end
    end
    y = yneu(1:n,1);
    for i = 1:2
        y = [y yneu(1:n,i)];
    end
end
