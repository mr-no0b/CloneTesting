# Equation-Solver
Done as a part of CSE 2208 project. Capable of solving linear, non linear and differential equations.

To operate the console application:
1. Choose the desired mode of operation by using the key value of each mode as input.
2. Then similarly choose the desired method to solve the problem (if any).

Description:

LINEAR SYSTEM:

   Jacobi Iterative:  Uses the iterative equations used in Jacobian iteration to approximately obtain the solutions to a system of linear equations.

   Gauss-Seidel Iterative: Similar to Jacobian, uses the iterative equations used in Gauss-Seidel iteration to approximately obtain the solutions to a system of linear equations.
 
   Gauss elimination: An algorithm starts by performing row operation on the rows of the matrix to convert the matrix into an upper triangular matrix. Then the matrix is passed on to the solver method that solves the equations using substitution method.

   Gauss Jordan elimination: The algorithm first performs gaussian elimination to turn the matrix into an upper triangular matrix and then performs more row operations to eliminate the elements above the diagonal to create a diagonal matrix. A solver method then solves the system by dividing the final values by the diagonal values.

   LU Decomposition:  The initial matrix A from the relation AX = B is decomposed into an upper and a lower triangular matrix namely L and U such that LUX = B. Then the we take UX = y so that LY = B; Solving this relation yields the value of Y which can then be used to obtain value of X using the relation UX = Y.

NON LINEAR SYSTEM:

   Bisection Method: This algorithm first take two value a and b from user so that f(a)*f(b)<0. Then we obtain a probable value  of x using x=(a+b)/2 formula. If f(x) and f(a) one of the same sign, we take the interval (x,b) for the iteration. Otherwise (a,x) is taken. After multiple iterations we reach a state when f(xi)=0.
   
   False Position Method: this method also firstly take two value a and b such that f(a)*f(b)<0. But here formula is x=(af(b)-bf(a))/(f(b)-f(a)). If f(a).f(x)>0, the interval (x,b) is chosen for the next iteration. Otherwise (a,x) is chosen as the next interval. After multiple iterations f(x) tends to 0.

   Secant Method: Here two initial value x0,x1 are take from user to find a root. In each iteration newer value can be found x(n+1)= xn- f(xn)*((xn-x(n-1))/(f(xn)-f(xn-1)). Then the next two values for the next iterations becomes xn and x(n+1), while x(n-1) is discarded. After multiple iterations we reach an approximate solution  very close to the actual root.

   Newton-Raphson Method: This method starts with an initial point x=x0 in each iteration. A better approximation is obtained using the following equation: x(n+1)=xn-f(xn)/f'(xn). In each iteration xn reaches mean the root and after multiple iterations a value very close to the root can be obtained.


 DIFFERENTIAL EQUATIONS:

    Runge-Kutta method: This algorithm effectively solves Ordinary differential equations of the form dy/dx=f(x,y) with initial condition y(x0)=y0. Fourth order Runge-Kutta method(RK4) provides good accuracy with less computation than higher ordered methods in most case. To implement RK4 method here we used 4 equations.

MATRIX INVERSION:
    
  At first we added an identity matrix to form augmented matrix AI and then performed the required row operations for gaussian elimination on the new augmented matrix. Then more row operations were performed so that the newly obtained upper triangular form so that only the diagonal values remain. Then the first half of the matrix was reduced to row echelon form thus resulting in an identity matrix. The augmented part now holds the inverse matrix of the given matrix.
