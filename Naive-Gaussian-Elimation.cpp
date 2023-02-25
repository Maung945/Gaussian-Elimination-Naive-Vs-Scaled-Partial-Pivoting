#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <cmath>

using namespace std;

// Forward Elimination with Scaled Partial Pivoting
void SPPFwdElimination(vector<vector<double>>& coeff, vector<double>& constant, vector<double>& ind) {
    int n = coeff.size();
    vector<double> scaling(n);

    // Initialize index and scaling vectors
    for (int i = 0; i < n; i++) {
        double smax = 0.0;
        for (int j = 0; j < n; j++) {
            smax = max(smax, abs(coeff[i][j]));
        }
        scaling[i] = smax;
    }

    for (int k = 0; k < n - 1; k++) {
        double rmax = 0.0;
        int maxInd = k;
        // Finding the max absolute value
        for (int i = k; i < n; i++) {
            double r = abs(coeff[ind[i]][k]) / scaling[ind[i]];
            if (r > rmax) {
                rmax = r;
                maxInd = i;
            }
        }
        swap(ind[maxInd], ind[k]);
        // Naive Elimination
        for (int i = k + 1; i < n; i++) {
            double mult = coeff[ind[i]][k] / coeff[ind[k]][k];
            for (int j = k + 1; j < n; j++) {
                coeff[ind[i]][j] -=  mult * coeff[ind[k]][j];
            }
            constant[ind[i]] -= mult * constant[ind[k]];
        }
    }
}
// Back Substitution with Scaled Partial Pivoting
void SPPBackSubst(vector<vector<double>>& coeff, vector<double>& constant, vector<double>& sol, vector<double>& ind) {
    int n = coeff.size();
    sol[n - 1] = constant[ind[n - 1]] / coeff[ind[n - 1]][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        double sum = constant[ind[i]];
        for (int j = i + 1; j < n; j++) {
            sum -= coeff[ind[i]][j] * sol[j];
        }
        sol[i] = sum / coeff[ind[i]][i];
    }
}

// SPP Gaussian Algorithm
void SPPGaussian(vector<vector<double>>& coeff, vector<double>& constant, vector<double>& sol) {
    int n = coeff.size();
    vector<double> ind(n);
    for (int i = 0; i < n; i++) {
        ind[i] = i;
    }
    SPPFwdElimination(coeff, constant, ind);
    SPPBackSubst(coeff, constant, sol, ind);
}

// Forward Elimination- NAIVE GAUSSIAN ELIMINATIN 
void FwdElimination(vector<vector<double>> &coeff, vector<double> &consts) {		// n^3 runtime complexity
    int n = coeff.size();
    for (int k = 0; k < n-1; k++) {
        for (int i = k+1; i < n; i++) {
            double mult = coeff[i][k] / coeff[k][k];
            for (int j = k+1; j < n; j++) {
                coeff[i][j] -= mult * coeff[k][j];
            }
            consts[i] -= mult * consts[k];
        }
    }
}

// Back Substitution- NAIVE GAUSSIAN ELIMINATIN 
void BackSubst(vector<vector<double>> &coeff, vector<double> &consts, vector<double> &sol) {
    int n = coeff.size();
    sol[n-1] = consts[n-1] / coeff[n-1][n-1];
    for (int i = n-2; i >= 0; i--) {
        double sum = consts[i];
        for (int j = i+1; j < n; j++) {
            sum -= coeff[i][j] * sol[j];
        }
        sol[i] = sum / coeff[i][i];
    }
}

// Naive Gaussian algorithm
void NaiveGaussian(vector<vector<double>> &coeff, vector<double> &consts, vector<double> &sol) {
    FwdElimination(coeff, consts);
    BackSubst(coeff, consts, sol);
}

int main() {
    ofstream outputFile;
    ifstream infile("sys1.lin");                                                // open input file
    double sizeofMatrix;
    infile >> sizeofMatrix;                                                     // read first line
    vector<vector<double>> coeff(sizeofMatrix, vector<double>(sizeofMatrix));   // create matrix
    vector<double> constant(sizeofMatrix);                                      // create vector to store last line
    vector<double> sol(sizeofMatrix);
    // read matrix elements and last line
    for (int i = 0; i < sizeofMatrix + 1; i++) {
        if (i == sizeofMatrix) {
            for (int j = 0; j < sizeofMatrix; j++) {
                infile >> constant[j];                                          // read last line
            }
        } else {
            for (int j = 0; j < sizeofMatrix; j++) {
                infile >> coeff[i][j];                                          // read matrix element
            }
        }
    }
    infile.close();                                                             // close input file
    cout << endl;
    cout << "Size: " << sizeofMatrix << endl;

    // print matrix and last line
    cout << "The Coefficient Matrix:" << endl;;
    for (int i = 0; i < sizeofMatrix; i++) {
        for (int j = 0; j < sizeofMatrix; j++) {
            cout << setw(9) << coeff[i][j];
        }
        cout << endl;
    }
    cout << "Constants: " << endl;
    for (int i = 0; i < sizeofMatrix; i++) {
        cout << setw(9) << constant[i];
    }
    cout << endl;
    int choice;
    while (true) {
        cout << "\n1-Scaled Partial Pivoting (OR) 2-Naive Gaussian algorithm ?: ";
        cin >> choice;
        
        if(choice == 1) {
            outputFile.open("extension.sol", ios::app);
            auto start = chrono::high_resolution_clock::now();
            NaiveGaussian(coeff, constant, sol);
            auto end = chrono::high_resolution_clock::now();
            auto naive_time_taken = chrono::duration_cast<chrono::microseconds>(end - start);
            cout << "Time Taken by Naive Gaussian Algorithm: " << naive_time_taken.count() << "-Microseconds" << endl << endl; 
            outputFile << "Solution Naive Gaussian Algorithm: " << endl;
            for (double x : sol) {
                outputFile << setw(12) << x;
            }
            outputFile << endl;
            outputFile << "Time Taken by Naive Gaussian Algorithm: " << naive_time_taken.count() << "-Microseconds" << endl << endl;
            outputFile.close();
            cout << endl;
            
            break;
        } else if (choice == 2) {
            outputFile.open("extension.sol", ios::app);
            auto start_SSP = chrono::high_resolution_clock::now();
            SPPGaussian(coeff, constant, sol);
            auto end_SSP = chrono::high_resolution_clock::now();
            auto SSP_time_taken = chrono::duration_cast<chrono::microseconds>(end_SSP - start_SSP);
            cout << "Time Taken by SPP Gaussian Algorithm: " << SSP_time_taken.count() << "-Microseconds" << endl;
            outputFile << "Solution Scaled Partial Pivoting: " << endl;
            for (double x : sol) {
                outputFile << setw(12) << x;
            }
            outputFile << endl;
            outputFile << "Time Taken by SPP Gaussian Algorithm: " << SSP_time_taken.count() << "-Microseconds.\n" << endl;
            outputFile.close();
            
            break;
        } else {
            cout << "Chosse 1 or 2!" << endl;
            false;
        }
    }
    return 0;
}
