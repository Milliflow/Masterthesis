n = 20;
x = zeros([1 n]);
dt = 0.005; % size of time-steps
T = 30000; % number of time-steps; total times is dt*T
Test = 10;
c = zeros(6,Test+1); % stores sizes of the basins of attraction
d = c; % stores return time
errorstore = 0;
warning = false;
for lambda = 0:Test
    retTime = 0;
    for Stat = 1:10000
        xstart = 2*pi*(rand(n,1)-0.5);
        x = xstart;
        x = dyKuramoto(x,lambda,dt,T);%Simulation
        q = 0;
        q = mod(pi+x(1)-x(n),2*pi)-pi;
        for i = 1:n-1
            q = q+mod(pi+x(i+1)-x(i),2*pi)-pi;%Detektion
        end
        q = q/(2*pi);
        for j = 0:4
            if abs(abs(q)-j)<0.4
                retTime = dyKuramoto2(xstart,q,lambda,dt,T);
                if retTime == T
                    warning = true;
                    errorstore = xstart;
                end
                c(j+1,t+1) = c(j+1,lambda+1)+1;
                d(j+1,t+1) = d(j+1,lambda+1)*(c(j+1,lambda+1)-1)/c(j+1,lambda+1)+dt*retTime/c(j+1,lambda+1);
            end
        end
        disp(d);
    end
end
disp(d);
hold 'off';
for i = 1:5
    plot(d(i,:));
    hold 'on';
end
hold 'off';

function y = dyKuramoto(x,lambda,dt,T)
    yalt = [x x];
    yneu = yalt;
    n = 20;
    for t = 1:T
        for i = 1:length(x)
            yneu(i) = yneu(i)+dt*(sin(yalt(i+1)-yalt(i))+sin(yalt(i+n-1)-yalt(i))+lambda*(2*sin(yalt(i+1)-yalt(i))+2*sin(yalt(i+n-1)-yalt(i))+0*sin(yalt(i+1)+yalt(i+n-1)-2*yalt(i))+sin(2*yalt(i+1)-2*yalt(i))+sin(2*yalt(i+n-1)-2*yalt(i))));
        end
        yalt = [yneu(1:length(x)) yneu(1:length(x))];
    end
    y = yneu(1:length(x));
end

function ret = dyKuramoto2(x,q,lambda,dt,T)
    yalt = [x x];
    b = false;
    n = length(x);
    for ret = 1:T
        yneu = yalt;
        for i = 1:n
            yneu(i) = yneu(i)+dt*(sin(yalt(i+1)-yalt(i))+sin(yalt(i+n-1)-yalt(i))+lambda*(2*sin(yalt(i+1)-yalt(i))+2*sin(yalt(i+n-1)-yalt(i))+0*sin(yalt(i+1)+yalt(i+n-1)-2*yalt(i))+sin(2*yalt(i+1)-2*yalt(i))+sin(2*yalt(i+n-1)-2*yalt(i))));
        end 
        yalt = [yneu(1:length(x)) yneu(1:length(x))];
        b = true;
        for tr = 2:n
            if abs(mod(pi+yalt(tr)-yalt(1)-2*pi*(tr-1)*q/n,2*pi)-pi)>0.4
                b = false;
            end
        end
        if b == true
            break;
        end
    end
end
