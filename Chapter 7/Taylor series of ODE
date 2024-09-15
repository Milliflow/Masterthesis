tay = [0 -1 0 1/6 0 -1/120 0 1/5040 0 -1/362880 0 1/39916800];
A = zeros(length(tay)); % stores the taylor coefficients of the derivatives
A(1,:) = tay;
for i = 2:length(tay)
    z = A(i-1,:);
    z = circshift(z,-1);
    for n = 1:length(z)-1
        z(n) = n*z(n);
    end
    z(length(z)) = 0;
    z = fold(z,tay);
    A(i,:) = z;
end
B = zeros(length(tay)); % stores the taylor coefficients of the products of the derivatives

B(1,:) = fold(fold(A(1,:),A(1,:)),A(1,:)); % taylor coefficients of x'^3
B(2,:) = fold(A(1,:),A(2,:)); % taylor coefficients of x'*x''
B(3,:) = A(3,:); % taylor coefficients of x'''
% those taylor represenations have partially to be shifted too the right, due too the multiplication with x^n
% the first coefficients of all representations should be in the same row; nullifying as much as possible
% gives a system of linear equations

function f = fold(a,b)
% Calculates the Cauchy product of two taylor series a and b
    f = zeros(length(a),1);
    for i = 1:length(a)
        for j   =1:i
            f(i) = f(i)+a(j)*b(i-j+1);
        end
    end
end
