q = 1; % twist
m = 1000;
A = zeros(m,m);
n = 10;
B = A;
criticalValue = 0;
criticalk = 0;
flag = false;

% pre-processing: only r_1 determines the stable zone
for i = 1:m
    if flag
        break;
    end
    for k = 1:n
        if flag
            break;
        end
        r = pi*i/m;
        z = c1Calc(q,k,r);
        if z>0
            criticalValue = i;
            criticalk = k;
            flag = true;
            break;
        end
    end
end

for i = criticalValue:criticalValue+4
    for j = 1:m
        k = criticalk;
        r = [pi*i/m pi*j/m];
        % bifurcation formula
        B(i,j) = -real(-1/2*(Vexp(q+k,r)+Vexp(-q+k,r)-3*(Vexp(q,r)+Vexp(-q,r))+3*(Vexp(q-k,r)+Vexp(-q-k,r))-(Vexp(q-2*k,r)+Vexp(-q-2*k,r))));
        B(i,j) = B(i,j)+real(1/4*(Vexp(q-2*k,r)-Vexp(-q-2*k,r)-2*(Vexp(q-k,r)-Vexp(-q-k,r))+Vexp(q,r)-Vexp(-q,r))*(-Vexp(q-2*k,r)+Vexp(-q-2*k,r)+Vexp(q-k,r)-Vexp(-q-k,r)+Vexp(q,r)-Vexp(-q,r)-Vexp(q+k,r)+Vexp(-q+k,r))*1/(pi*(Vsin(q-k,r)+Vsin(q+k,r))-C1(q,2*k,r(1))-pi/2*(Vsin(q-2*k,r)+Vsin(q+2*k,r))));
        
    end
end

% post-processing
for j = 1:m
    for i = 1:criticalValue-1
        B(i,j) = 1;
    end
    for i = criticalValue:criticalValue+4
        if B(i,j) < 0
            B(i,j) = 2;
        else
            B(i,j) = 4;
        end
    end
    for i = criticalValue+5:m
        B(i,j) = 3;
    end
end
B = B';

cmap = zeros(10, 3);
cmap = [0, 1, 0; ...   % Blue for 1
  1, 1, 0; ...       % Black for 2
  1, 0, 0];       % Green for 3


figure('pos',[100 100 1000 800]);
ax = subplot(1,1,1);
imagesc([0,pi],[0,pi],B);
xlabel("r_1");
ylabel("r_2");
ax.YDir = "normal";
title('Bifurcation Portrait');
%C = colormap();
C = cmap();
L = size(C,1);
Gs = round(interp1(linspace(min(B(:)),max(B(:)),L),1:L,B));
H = reshape(C(Bs,:),[size(Bs) 3]);
subplot(1,2,2);
imshow(B, cmap);

function v = Vexp(k,r)
% V Calculates the exp-fourier coefficient for the given values
    if k ~ 0
        v = 2*sin(k*r(1))/k+(exp(1i*k*r(2))+exp(-1i*k*r(2))-2)/(1i*k);
    else
        v = 2*r(2);
    end
end
function v = Vsin(k,r)
% V Calculates the sin-fourier coefficient for the given values
    if k ~ 0
        v = (2-2*cos(k*r(2)))/(pi*k);
    else
        v = 0;
    end
end

function c1 = C1(q,k,r)
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function v = V(r,k)
% V Calculates the fourier coefficient for the given values
    v = 2;
    for i = 1:length(k)
        if k(i) == 0
            v_i = r(i)/pi;
        else
            v_i = sin(k(i)*r(i))/(pi*k(i));
        end
        v = v*v_i;
    end
end
