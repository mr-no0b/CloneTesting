#include <bits/stdc++.h>
using namespace std;
vector<double> coeffs,ans;

double fx(double x)
{
    double res = 0.0;
    int n=coeffs.size();
    for (int i = 0; i < n; ++i)
        res+=coeffs[i]*pow(x,n-i-1);
    return res;
}

double fdx(double x)
{
    double res = 0.0;
    int n=coeffs.size();
    for (int i = 0; i < n - 1; ++i)
        res+= coeffs[i]*(n-i-1)*pow(x,n-i-2);
    return res;
}

double nr_method(double x0)
{
    double tol = 1e-6;
    int maxIt = 1000;
    double x = x0;
    for (int i = 0; i < maxIt; i++)
    {
        double f = fx(x);
        double fPrime = fdx(x);
        if (fabs(fPrime) < 1e-10)
        {
            cout << "Method failed for too small derivative." << endl;
            return x;
        }
        double x1=x- f /fPrime;
        if (fabs(x1 - x) < tol)
            return x1;
        x = x1;
    }

    cout << "Method not converge within the maximum iterations." << endl;
    return x;
}

double fp_method(double a, double b)
{
    double c = a;
    int maxIt = 1000;
    double tol = 0.00001;
    for (int i = 0; i < maxIt; i++)
    {
        double fa = fx(a);
        double fb = fx(b);
        c = (a * fb - b * fa) / (fb - fa);
        double fc = fx(c);
        if (fabs(fc) < tol or fabs(b - a) < tol)
            return c;
        if (fa * fc < 0)
            b = c;
        else
            a = c;
    }
    return c;
}


vector<double> def(double root)
{
    int n = coeffs.size();
    vector<double> deft(n - 1);
    deft[0] = coeffs[0];
    for (int i = 1; i < n - 1; ++i)
    {
        deft[i] = coeffs[i];
        deft[i]+= deft[i - 1] * root;
    }
    return deft;
}
int false_position(int n)
{
    cout<<endl<<"**********FALSE POSITION METHOD********"<<endl;
    while(n--)
    {
        double an = coeffs[0];
        double an_1 = coeffs[1];
        double an_2 = coeffs[2];

        double xmax = sqrt(pow(an_1 / an, 2) - 2 * (an_2 / an));

        double a = -fabs(xmax), b = fabs(xmax);
        double root = fp_method(a, b);
        ans.push_back(root);
        coeffs = def(root);
    }
}

int newton_raphson(int n)
{
    cout<<endl<<"**********NEWTON-RAPHSHON METHOD********"<<endl;
    double x0=0;
    while(n--)
    {
        double root = nr_method(x0);
        ans.push_back(root);
        coeffs = def( root);
    }
}

vector <int> * gaus(vector <int> ara[], int n){
    int i = n-1, j = 0;
    //Start:
    while(1){
        if(ara[i][j] != 0) {
            int k = 1;
            /*
            if(ara[i][i] == 0){
                for(int k = 0; ara[i][k] != 0 && ara[k][i] != 0;k++){
                    swap(ara[i],ara[k]);
                    goto Start;
                }
            }
            */
            while(ara[i-k][j] == 0) {
                k++;
            }
            int temp1 = ara[i][j];
            for(int l = 0; l <= n; l++) {
                    ara[i][l] *= ara[i-k][j];
            }
            for(int l = 0; l <= n; l++) {
                ara[i-k][l] *= temp1;
            }
            for(int x = 0; x <= n;x++) {
                ara[i][x] -= ara[i-k][x];
            }
        }
        i--;
        if(i == j) {
            j++; i = n-1;
        }
        if(j == n-1) break;
    }
    return ara;
}
vector <int> * Jordan(vector <int> ara[], int n){
    int i = 0, j = n-1;
    while(1){
        if(ara[i][j] != 0){
            int k = 1;
            while(ara[i+k][j] == 0) {
                k++;
            }
            int temp1 = ara[i][j];
            for(int l = 0; l <= n; l++) {
                    ara[i][l] *= ara[i+k][j];
            }
            for(int l = 0; l <= n; l++) {
                ara[i+k][l] *= temp1;
            }
            for(int x = 0; x <= n;x++) {
                ara[i][x] -= ara[i+k][x];
            } 
        }
        i++;
        if(i == j) {
            j--; i = 0;
        }
        if(j == 0) break;
    }
    return ara;
}
void Jordan_solver(vector <int> ara[], int n){
    Jordan(ara, n);
    for(int i = 0; i < n; i++){
        ara[i][n] = ara[i][n]/ara[i][i];
        ara[i][i] = ara[i][i]/ara[i][i];
    }
    for(int i = 0; i < n; i++){
        cout << "x" << i+1 << " = " << ara[i][n] << endl;
    }
}
void gaus_solver(vector <int> ara[], int n){
    int ans[n] = {};
    ans[n-1] =ara[n-1][n] / ara[n-1][n-1];
    for(int i = n - 2; i >= 0; i--){
        ans[i] = ara[i][n];
        for(int j = i+1; j < n; j++){
            ans[i] -= ara[i][j]*ans[j];
        }
        ans[i] /= ara[i][i];
    }
    for(int i = 0 ; i < n; i++){
        cout << "x"<< i+1 << " = " << ans[i] << endl;
    }
}
double fx(double x, double y){
    return x-y;
}
double ffx(double x){
    return x = x - 1.0 + 2.0/exp(x);
}
void rungekutta(){
    vector <double> arax;
    vector <double> aray;
    vector <double> araye;
    int i = 0;
    double x = 0, y = 1, step = 0.1;
    araye.push_back(1);
    while(x <= 10){
        cout << i++ << endl;
        arax.push_back(x);
        aray.push_back(y);
        double k1 = step*fx(x,y);
        double k2 = step*fx(x+step/2.0,y+k1/2.0);
        double k3 = step*fx(x+step/2.0,y+k2/2.0);
        double k4 = step*fx(x+step,y+k3); 
        y = y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        x +=  step;
        araye.push_back(ffx(x));
    }
    for(auto it: arax){
        cout << it << endl;
    }
    for(auto it: aray){
        cout << it << endl;
    }
    
}

void jacobiMethod(vector<int> vec[], int n)
{
    vector<double> x(n,0.0),x1(n, 0.0);
    int maxIt=1000;
    double tol=0.00001;
    cout<<endl<<"**********JACOBI ITERATIVE METHOD********"<<endl;
        for (int it = 0; it < maxIt; it++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = vec[i][n];
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sum -= vec[i][j] * x[j];
            }
            x1[i] = sum / vec[i][i];
        }

        double err = 0.0;
        for (int i = 0; i < n; i++)
            err += fabs(x1[i] - x[i]);

        err /= n;

        if (err< tol)
        {
            cout<<endl<<"The result: "<<endl;
            for (int i = 0; i < n; i++)
                cout << "x" << i + 1 << " = " << x1[i] << "\n";
            cout<<endl<<endl;
            return;
        }
        x = x1;
    }
    cout << "Did not converge within the maximum iterations.\n";
}

void gauss_seidel(vector<int>vec[], int n)
{
    vector<double> x(n, 0.0);
    int maxIt=1000;
    double tol=0.00001;
    cout<<endl<<"**********GAUSS-SEIDEL METHOD********"<<endl;
    for (int it = 0; it < maxIt; it++)
    {
        double err = 0.0;
        for (int i = 0; i < n; i++)
        {
            double sum = vec[i][n];
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sum -= vec[i][j] * x[j];
            }
            double x1 = sum / vec[i][i];
            err += fabs(x1 - x[i]);
            x[i] = x1;
        }
        err /= n;
        if (err < tol)
        {
            cout<<endl<<"The result: "<<endl;
            for (int i = 0; i < n; i++)
                cout << "x" << i + 1 << " = " << x[i] << "\n";
            cout<<endl<<endl;
            return;
        }
    }
    cout << "Did not converge within the maximum iterations.\n";
}

double secantMethod(double x0, double x1)
{
    int maxIter=1000;
    double tol=0.00001;
    double x2;
    for (int i = 0; i < maxIter; ++i)
    {
        double f_x0 = fx(x0);
        double f_x1 = fx(x1);

        if (fabs(f_x1 - f_x0) < 1e-10)
        {
            cout << "Division by zero error in method." << endl;
            return x1;
        }
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);

        if (fabs(x2 - x1) < tol)
        {
            return x2;
        }

        x0 = x1;
        x1 = x2;
    }

    cout << "Method did not converge within the maximum iterations." << endl;
    return x2;
}


int main()
{
    int OPcode = 0, method = 0;
    printf("Please choose your desired operation mode:\n\n");
    printf("  0 - Solution of System of Linear equations:\n");
    printf("  1 - Solution of non linear equations:\n");
    printf("  2 - Solution of differential equations:\n");
    printf("  3 - Matrix inversion:\n");
    cin >> OPcode;
    if(OPcode == 0){
        printf("Choose the desired method: \n");
        printf("\t0 - Jacobi Iteration\n");
        printf("\t1 - Gauss-Seidel Iterative method\n");
        printf("\t2 - Gauss elimination\n");
        printf("\t3 - Gauss-Jordan elimination\n");
        printf("\t4 - LU Factorization\n\n");
        cin >> method;
        printf("Please enter the number of equations: ");
        int n, temp;
        cin >> n;
        vector <int> ara[n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j <= n; j++){
                cin >> temp;
                ara[i].push_back(temp);
            }
        }
        if(method == 0){
          jacobiMethod(ara,n);
        }
        else if(method == 1){
            gauss_seidel(ara,n);
        }
        else if(method == 2){
            gaus(ara, n);
            gaus_solver(ara, n);
            
        }
        else if(method == 3){
            gaus(ara, n);
            Jordan(ara, n);
            Jordan_solver(ara, n);
        }
        else if(method == 4){

        }
        for(int i = 0; i < n; i++){
            for(int j = 0; j <= n; j++){
                cout << ara[i][j] << " ";
            }
            cout << endl;
        }
    }
    if(OPcode  == 1){
          int n;
    cout << "Enter degree of equation: ";
    cin >> n;
    coeffs.resize(n + 1);
    cout << "Enter the coefficients: ";
    for (int i=0; i<=n; i++)
        cin>>coeffs[i];
        printf("Choose the desired method: \n");
        printf("\t5 - Bisection method\n");
        printf("\t6 - False position method\n");
        printf("\t7 - Secant method\n");
        printf("\t8 - Newton Raphson method\n\n");
        cin >> method;
        if(method == 5){

        }
        else if(method == 6){
         false_position(n);
        }
        else if(method == 7){
                for (int k = 0; k<n; k++)
    {
        double x0, x1;
        cout << "Enter two initial guesses to find root "<<k+1<<": ";
        cin >> x0 >> x1;

        double root = secantMethod(x0, x1);
        ans.push_back(root);

        coeffs = def( root);
    }

        }
        else if(method == 8){
          newton_raphson(n);
        }
     cout << endl<<"Roots found:" << endl;
    for (int i=0; i<ans.size(); i++)
    {
        cout <<"x"<<i+1<<" = "<<ans[i]<<endl;
    }
        
    }
    return 0;
}
