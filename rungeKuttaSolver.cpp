/* Runge Kutta 4th order solver for 2nd order ODE

http://lampx.tugraz.at/~hadley/num/ch8/rk4ode2.php

Any second order differential equation can be written as two coupled first order equations,
dx1/dt = f1(x1, x2, t)
dx2/dt = f2(x1, x2, t).

These coupled equations can be solved numerically using a fourth order Runge - Kutta routine.

The relevant equations for an electron moving inside a magnetic bottle take the form:
ds/dt = v					--->	f1(x1, x2, t) = x2
dv/dt = - mue * gradB / me	--->	f2(x1, x2, t) = - mue * gradB(x1) / me

with intial conditions x1(t = 0) = 0, x2(t = 0) = 1.
*/
double RungeKutta4(double xo, double vo, double t, double dt) {

}

<script>
var x1 = 0;
var x2 = 1;
var t = 0;
var dt = 0.01;

function f1(x1, x2, t) {
	dx1dt = x2;
	return dx1dt;
}

function f2(x1, x2, t) {
	dx2dt = -x2 - Math.sin(x1) + Math.sin(t);
	return dx2dt;
}

for (j = 0; j < 20; j++) {
	document.write("t=" + t + "  x1=" + x1 + "  x2=" + x2 + "<\/br>");
	k11 = dt * f1(x1, x2, t);
	k21 = dt * f2(x1, x2, t);
	k12 = dt * f1(x1 + 0.5 * k11, x2 + 0.5 * k21, t + 0.5 * dt);
	k22 = dt * f2(x1 + 0.5 * k11, x2 + 0.5 * k21, t + 0.5 * dt);
	k13 = dt * f1(x1 + 0.5 * k12, x2 + 0.5 * k22, t + 0.5 * dt);
	k23 = dt * f2(x1 + 0.5 * k12, x2 + 0.5 * k22, t + 0.5 * dt);
	k14 = dt * f1(x1 + k13, x2 + k23, t + dt);
	k24 = dt * f2(x1 + k13, x2 + k23, t + dt);
	x1 += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
	x2 += (k21 + 2 * k22 + 2 * k23 + k24) / 6;
	t += dt;
}
< / script>