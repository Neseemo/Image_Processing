function [x] = soft_thr(x, lambda)
    x(x>lambda) = x(x>lambda) - lambda;
    x(x<-lambda) = x(x<-lambda) + lambda;
    x( abs(x) < lambda) = 0;
end