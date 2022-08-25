function [s, t] = get_synthetic_signal(name, n)

t = linspace(0, 1, n)';

if strcmpi(name, 'blocks')
    h = [4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 5.1, -4.2]';
    tt = [.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81]';
    s = zeros(n, 1);
    for ii=1:length(h)
        s = s + h(ii) * (1 + sign(t - tt(ii))) / 2;
    end
    
    return
end

if strcmpi(name, 'bumps')
    tt = [.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81]';
    h = [4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2]';
    w = [.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005]';
    s = zeros(n, 1);
    for ii=1:length(h)
        s = s + h(ii) ./ (1 + abs((t - tt(ii))/w(ii)).^4);
    end
    return
end


if strcmpi(name, 'parabolas')
    p = @(x)(8*x.^2 - 2*x + 2);
    s = p(t);
    s(t>0.5) = s(t>0.5)+ 7;
end
    
    
if strcmpi(name, 'cubic')
    p = @(x)(50*(x - 0.1).*(x - 0.5).* (x - 0.7));
    s = p(t);
end

if strcmpi(name, 'sine')
    p = @(x)(sin(2./(x + 0.05)));
    s = p(t);
end
    
