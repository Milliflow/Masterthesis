n = 20;
x = zeros(n,2);
dt = 0.005; % size of time-step
T = 30000; % number of time-steps; total time is dt*T
Test = 10;
c = zeros(5,Test); % stores the size of the basins of attraction
d = c; % stores the return times
errorstore = 0;
warning = false;
for lambda = 1:Test
    retTime = 0;
    for Stat = 1:10000
        xstart = 2*pi*(rand(n,2)-0.5);
        x = xstart;
        x = dyKuramoto(x,lambda,dt,T); % simulation
        q = zeros([0,2]);
        for j = 1:2
            q(j) = mod(pi+x(1,j)-x(n,j),2*pi)-pi;
            for i = 1:n-1
                q(j) = q(j)+mod(pi+x(i+1,j)-x(i,j),2*pi)-pi; % detection
            end
        end
        q = q./(2*pi);
        retTime = 0;
        retTime = dyKuramoto2(xstart,q,lambda,dt,T);
        if retTime == T
            warning = true;
            errorstore = xstart;
            Stat = 10001;
            lambda = Test+1;
            break;
        end
        for j = 0:4
            qtrue = true;
            for i = 1:2
                if abs(abs(q(i))-j)>0.4
                    qtrue = false;
                end
            end
            if qtrue
                c(j+1,lambda) = c(j+1,lambda)+1;
                d(j+1,lambda) = d(j+1,lambda)*(c(j+1,lambda)-1)/c(j+1,lambda)+dt*retTime/c(j+1,lambda);
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

function y = dyKuramoto(x,lambda,dt,T) % returns the attractor
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
    y = [y yneu(1:n,2)];
end

function ret = dyKuramoto2(x,q,lambda,dt,T) % returns the time it needs to get close enough to the attractor
    sz = size(x);
    n = sz(1);
    yalt = [x; x];
    yneu = yalt;
    b = false;
    for ret = 1:T
        yneu = yalt;
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
            for i = 1:2
            yalt(:,i) = [yneu(1:n,i); yneu(1:n,i)];
            end
        end
        y = yneu(1:n,1);
        y = [y yneu(1:n,2)];
        b = true;
        for i = 1:2
            for tr = 2:n
                if abs(mod(pi+yalt(tr,i)-yalt(1,i)-2*pi*(tr-1)*q/n,2*pi)-pi)>0.4
                    b = false;
                end
            end
        end
        if b == true
            break;
        end
    end
end
