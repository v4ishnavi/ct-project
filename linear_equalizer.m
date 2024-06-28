function Y = linear_equalizer(received_signal, desired_signal, mu, num_coeff)
    N = length(received_signal);
    % disp('Total Samples: ' + string(N));
    w = 0.005 * randn(num_coeff, 1); 
    received_signal = received_signal';
    desired_signal = desired_signal';
    Y_arr = zeros(N, 1);

    for n = num_coeff:N
        x_new = received_signal(n:-1:n-num_coeff+1);
        y = w' * x_new;
        Y_arr(n) = y;
        e = desired_signal(n) - y;
        e_conj = conj(e);

        w = w + mu * e_conj * x_new;

        if any(isnan(w))
            disp('NaN detected in filter coefficients at iteration ' + string(n));
            break; 
        end
    end
    Y = Y_arr';
end
