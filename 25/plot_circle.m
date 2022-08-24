function  plot_circle(x_c,y_c,r, col)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x_c;
yunit = r * sin(th) + y_c;
plot(xunit, yunit, 'Color', col);
