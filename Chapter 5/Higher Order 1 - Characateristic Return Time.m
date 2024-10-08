q = 1; % twist
r = pi/10; % should translate to the discretized model with 20 nodes
for i = 1:400 % combined with the next code line this translates to the intervall [0,4]
    x = i/100;
    k = 1;
    a(i) = (1+2*pi*x*V(r,q))*C1(q,1,r); % stores the highest eigenvalue for the q-twist
    b(i) = (1+2*pi*x*V(r,0))*C1(0,1,r); % stores the highest eigenvalue for the 0-twist
    for k = 2:5 % iterates through several k's to find the highest eigenvalue
        achallenger = (1+2*pi*x*V(r,q))*C1(q,k,r);
        bchallenger = (1+2*pi*x*V(r,0))*C1(0,k,r);
        if a(i) < achallenger
            a(i) = achallenger;
        end
        if b(i) < bchallenger
            b(i) = bchallenger;
        end
    end
end
plot(a./b); % plots the quotient; for the other graphs, plot -1./a and -1./b
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt/100)
xlabel("\lambda");
ylabel("quotient");
title("Higher Order 1 Characteristic Return Time quotient");

function c1 = C1(q,k,r)
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function v = W(r,q,k)
% V Calculates the fourier coefficient for the given values
    v = V(r,q)*V(r,k)/2;
end

function v = V(r,k)
% V Calculates the fourier coefficient for the given values
    v = 2;
    for i=1:length(k)
        if k(i) == 0
            v_i = r(i)/pi;
        else
            v_i = sin(k(i)*r(i))/(pi*k(i));
        end
        v = v*v_i;
    end
end
