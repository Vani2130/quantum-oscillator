#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


// Constants
const double hbar = 1.0;  // Reduced Planck constant
const double m = 1.0;     // Mass of the particle
const double omega = 1.0; // Angular frequency of the oscillator

// Function to define the Schrödinger equation
vector<double> schrodingerEquation(double x, const vector<double>& y, double E) {
    // y[0]: Psi(x)
    // y[1]: Psi'(x)
    double psi = y[0];
    double psi_prime = y[1];

    // Calculate second derivative of Psi
    double psi_double_prime = -(2 * m / hbar) * (E - 0.5 * m * omega * omega * x * x) * psi;

    // Return the derivatives as a vector
    return {psi_prime, psi_double_prime};
}

// Runge-Kutta method to solve the ODE
vector<double> rungeKutta(double x, double step, vector<double> y, double E) {
    // Calculate the derivatives at the current point
    vector<double> k1 = schrodingerEquation(x, y, E);
    vector<double> y1 = {y[0] + 0.5 * step * k1[0], y[1] + 0.5 * step * k1[1]};

   vector<double> k2 = schrodingerEquation(x + 0.5 * step, y1, E);
   vector<double> y2 = {y[0] + 0.5 * step * k2[0], y[1] + 0.5 * step * k2[1]};

    vector<double> k3 = schrodingerEquation(x + 0.5 * step, y2, E);
    vector<double> y3 = {y[0] + step * k3[0], y[1] + step * k3[1]};

    vector<double> k4 = schrodingerEquation(x + step, y3, E);

    // Calculate the final solution
    y[0] += (step / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    y[1] += (step / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

    return y;
}

// Function to solve the Schrödinger equation using the Runge-Kutta method
void solveSchrodinger(double E, double x0, double x1, double step) {
    // Initial conditions
    vector<double> y = {0.0, 1.0}; // Psi(x0) = 0, Psi'(x0) = 1

    // Solve the equation from x0 to x1
    for (double x = x0; x <= x1; x += step) {
        y = rungeKutta(x, step, y, E);
    }

    // Print the final Psi at x1
    cout << "Psi(" << x1 << ") = " << y[0] << endl;
}

// Main function
int main() {
    // Define the energy eigenvalue (guess)
    double E = 0.5 * hbar * omega; // Guess the ground state energy

  // Prompt user to input the range and step size
    double x0, x1, step;
    cout << "Enter the starting value of x (x0): ";
    cin >> x0;
    cout << "Enter the ending value of x (x1): ";
    cin>>x1;
    cout << "Enter the step size (step): ";
    cin >> step;



    // Solve the Schrödinger equation
    solveSchrodinger(E, x0, x1, step);

    return 0;
}

