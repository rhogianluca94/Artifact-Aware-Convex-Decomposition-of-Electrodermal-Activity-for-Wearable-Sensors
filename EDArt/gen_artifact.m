function irf = gen_artifact(t0, n, fs, artType)

    % artType -> 1: mano aperta; 2: mano chiusa; 3: pressione elettrodo; 
    %            4: tensione filo elettrodo/mano aperta/mano chiusa
    %            (falangi prossimali);
    %            5: 1,2,3 lungo (fase in salita)
    %            6: 1,2,3 lungo (fase in discesa)

    switch artType
        case 1
            ampl = 2.3542 * .1919;
%             ampl = gamrnd(2.3542, .1919);
%             a = gamrnd(3.8202, 2.3378);
%             b = gamrnd(3.3680, 2.4481);
%             a = gamrnd(5.8202, 2.3378);
%             b = gamrnd(5.3680, 2.4481);
            b = (4.3680) * 2.4481;
            a = (4.8202) * 2.3378;
            T = 1/fs;
            ar = [1, ...
                  -2*exp(-b/fs)-2*exp(-a/fs), ...
                  4*exp(-(a+b)/fs) + exp(-2*b/fs)+exp(-2*a/fs),...
                  -2*exp(-(a+2*b)/fs)-2*exp(-(2*a+b)/fs), ...
                  exp(-2*(a+b)/fs)];
            A = exp(-a/fs)/fs;
            B = exp(-b/fs)/fs;
            ma = [1 ...
                  A+B ...
                  -2*A*exp(-b*T)-2*B*exp(-a*T) ...
                  exp(-2*b*T) + exp(-2*a*T)] ./(2000) * ampl;%./(3636.11) * ampl;
        case 2
            ampl = 2.3542 * .1919;
%             ampl = gamrnd(2.3542, .1919);
%             a = gamrnd(3.8202, 2.3378);
%             b = gamrnd(3.3680, 2.4481);
%             a = gamrnd(5.8202, 2.3378);
%             b = gamrnd(5.3680, 2.4481);
            b = (4.3680) * 2.4481;
            a = (4.8202) * 2.3378;
            T = 1/fs;
            ar = [1, ...
                  -2*exp(-b/fs)-2*exp(-a/fs), ...
                  4*exp(-(a+b)/fs) + exp(-2*b/fs)+exp(-2*a/fs),...
                  -2*exp(-(a+2*b)/fs)-2*exp(-(2*a+b)/fs), ...
                  exp(-2*(a+b)/fs)];
            A = exp(-a/fs)/fs;
            B = exp(-b/fs)/fs;
            ma = -[1 ...
                  A+B ...
                  -2*A*exp(-b*T)-2*B*exp(-a*T) ...
                  exp(-2*b*T) + exp(-2*a*T)] ./(2000) * ampl;        
        case 3
            ampl = gamrnd(4.711, 0.197);
            a1 = 20.319;
            a2 = 54.869;
            w = round(unifrnd(a1, a2));
            ma = [1 1];
            ar = [1 -1];
        case 4
            ampl = gamrnd(4.711, 0.197);
            b1 = 0.312;
            b2 = 2.927;
            A = 1;
            b = unifrnd(b1, b2);
            ma = A.*[0 -1+exp(-b/fs)];
            ar = [1 -1-exp(-b/fs) exp(-b/fs)];
        case 5
            ampl = gamrnd(4.711, 0.197);
            a1 = 20.319;
            a2 = 54.869;
            w = round(unifrnd(a1, a2));
            ma = [-1 -1];
            ar = [1 -1];
        case 6
            ampl = gamrnd(4.711, 0.197);
            b1 = 0.312;
            b2 = 2.927;
            A = 1;
            b = unifrnd(b1, b2);
            ma = A.*[0 1-exp(-b/fs)];
            ar = [1 -1-exp(-b/fs) exp(-b/fs)];
    end

    % generate artifact IRF
    imp = zeros(n,1);
    imp(t0,1) = fs; % unitary area pulse
    yd = filter(ma, ar, imp);
    if artType == 3 || artType == 5
        movAvg = ones(1,w)/w;
        fc = .05;
        [bButter, aButter] = butter(1, fc, 'low');
        yd = filter(movAvg, 1, yd);
        yd = filter(bButter, aButter, yd);
    end
%     irf_norm = (yd-min(abs(yd)))/(max(abs(yd))-min(abs(yd)));
%     irf = irf_norm * ampl;
    irf = yd;
end