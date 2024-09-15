q = 1; % twist
m = 1000;
A = zeros(m,m);
n = 20;
B = A;
for i = 1:m
    for l = 1:m
        for k1 = 0:n
            k = [k1];
            r = [pi*i/m];
            lambda = l/100-3;
            z = (1+2*pi*lambda*V(r,q))*C1(q,k,r);
            if z > 0 && A(i,l) == 0
                A(i,l) = 1;
                c = C5(q,k,r) + lambda*(3*pi^2/4*(2*(W(r,q,q-2*k)-W(r,q,q-k)+W(r,q,q)-W(r,q,q+k)+W(r,q,q+2*k))-W(r,k-q,q-k)+W(r,k+q,q-k)+W(r,k-q,k+q)-W(r,q+k,q+k)-W(r,k-q,q+2*k)+W(r,k+q,q+2*k)+W(r,k-q,q-2*k)-W(r,k+q,q-2*k)));
                B(i,l) = c - C2(q,k,r)*C3(q,2*k,k,r)/((1+2*pi*lambda*V(r,q))*C1(q,2*k,r));
            end
        end
    end
end
for i = 6:m-5
    for lambda = 1:m-5
        if A(i,lambda) ~ 0
            if [A(i-1,lambda+1)==0 ||A(i-2,lambda+2)==0||A(i-3,lambda+3)==0||A(i-4,lambda+4)==0||A(i-5,lambda+5)==0||A(i+1,lambda+1)==0 ||A(i+2,lambda+2)==0||A(i+3,lambda+3)==0||A(i+4,lambda+4)==0||A(i+5,lambda+5)==0] && B(i,lambda)>0
                B(i,lambda) = 4;
            elseif [A(i-1,lambda+1)==0 ||A(i-2,lambda+2)==0||A(i-3,lambda+3)==0||A(i-4,lambda+4)==0||A(i-5,lambda+5)==0||A(i+1,lambda+1)==0 ||A(i+2,lambda+2)==0||A(i+3,lambda+3)==0||A(i+4,lambda+4)==0||A(i+5,lambda+5)==0] && B(i,lambda)<0
                B(i,lambda) = 2;
            else
                B(i,lambda) = 0;
            end
        end
    end
end
for i = 1:m
    for lambda = 1:m
        if B(i,lambda) == 2
            B(i,lambda) = 2;
        elseif B(i,lambda) == 4
            B(i,lambda) = 4;
        elseif A(i,lambda) == 0
            B(i,lambda) = 1;
        else
            B(i,lambda) = 3;
        end
    end
end
B=B';

cmap = zeros(10, 4);
cmap = [0, 1, 0; ...   % Blue for 1
  1, 1, 0; ...       % Yellow for 2
  1, 0, 0; ...       % Green for 3
  0, 1, 1];         %Lightblue

figure('pos',[100 100 1000 800]);
ax = subplot(1,1,1);
imagesc([0,pi],[-2,8],B);
xlabel("r");
ylabel("\lambda");
ylim([-1.5,7.5]);
yticklabels({'-2','-1','0','1','2','3','4','5','6'});
ax.YDir = "normal";
title('Bifurcation Portrait');
%C = colormap();
C = cmap();
L = size(C,1);
Gs = round(interp1(linspace(min(B(:)),max(B(:)),L),1:L,B));
H = reshape(C(Bs,:),[size(Bs) 3]);
subplot(1,2,2);
imshow(B, cmap);
colormap(cmap);

function c1 = C1(q,k,r)
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function c2 = C2(q,k,r)
    c2 = (2*pi)^length(r)*(-V(r,q-2*k)+2*V(r,q-k)-2*V(r,q+k)+V(r,q+2*k))/8;
end

function c3 = C3(q,k,j,r)
    c3 = (2*pi)^length(r)*(-V(r,q-k)+V(r,q-k+j)+V(r,q-j)+V(r,q+k)-V(r,q+k-j)-V(r,q+j))/8;
end

function c5 = C5(q,k,r)
    c5 = (2*pi)^length(r)*(V(r,q-2*k)-4*V(r,q-k)+6*V(r,q)-4*V(r,q+k)+V(r,q+2*k))/16;
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
