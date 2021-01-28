#include <iostream>
#include <cmath>

//euler declaration
const double EulerConstant = std::exp(1.0);

//function declarations
double function(double x);
void newton(double p0, double tol, double max_iter);
double dfunction(double x);

//newton's method implementation
void newton(double p0, double tol, double max_iter){
    int i = 1;
    static double p = 0.0;
    static double p_0 = p0;
    while (i < max_iter){
        p = p_0-(function(p_0)/dfunction(p_0));
        
        std::cout << 'P' << i << " = " << p << std::endl;
        
        double absolute = abs(p-p_0);
        if (absolute < tol){
            std::cout << "Finished in " << i << " iterations" << std::endl;
            return;
        }
        i++;
        p_0 = p;
    }
    std::cout << "Method Failed" << std::endl;
}

//function
double function(double x){
    
  //  double efun = pow(EulerConstant, x);
  //  double xpow = pow(2, -x);
  //  double cosin = 2*cos(x);
  //  double constant = -6;
  //  double square = pow(x, 2);
  //  double squareRoot = sqrt(x);
  //  double variable = 0;
     
  //  double result = efun + xpow + cosin + constant; //e^x + 2^-x + 2cosx - 6
  //  std::cout << result << std::endl;
    
    double cube = pow(x, 3);
    double cosin = cos(x);

    double result = -cube - cosin;
    return result;
    
}

//derivative
double dfunction(double x){
    
  //  double efun = pow(EulerConstant, x);
  //  double xpow = pow(2, -x);
  //  double cosin = 2*cos(x);
  //  double constant = -6;
  // double cube = pow(x, 3);
  // double cosin = cos(x);
  // double squareRoot = sqrt(x);
  // double variable = 0;
  // double result = efun + xpow + cosin + constant; //e^x + 2^-x + 2cosx - 6
  // std::cout << result << std::endl;
    
    double sinf = sin(x);
    double square = pow(x, 2);

    double result = -3*square + sinf;
    return result;
}

int main(int argc, const char * argv[]) {
    //startin p0
    double p0 = -1.0;
    
    //tolerance
    double tol = pow(10.0, -4.0);
    
    //iterations
    double max_iter = 100.0;
    
    //function call
    newton(p0, tol, max_iter);

}
