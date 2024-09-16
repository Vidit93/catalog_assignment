#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>  
#include "json.hpp"  
#include <cmath>    

using namespace std;
using json = nlohmann::json;

// Function to decode y value from a given base
long long decodeBase(const string& value, int base) {
    long long decodedValue = 0;
    int len = value.size();
    for (int i = 0; i < len; ++i) {
        char digit = value[i];
        int num;
        if (digit >= '0' && digit <= '9') {
            num = digit - '0';
        } else if (digit >= 'A' && digit <= 'Z') {
            num = digit - 'A' + 10;  // For bases greater than 10
        } else {
            cerr << "Invalid character in base representation: " << digit << endl;
            exit(EXIT_FAILURE);
        }
        if (num >= base) {
            cerr << "Digit " << digit << " is out of range for base " << base << endl;
            exit(EXIT_FAILURE);
        }
        decodedValue = decodedValue * base + num;
    }
    return decodedValue;
}

// Gaussian elimination to solve system of linear equations
vector<double> gaussianElimination(vector<vector<double>>& matrix, int n) {
    for (int i = 0; i < n; i++) {
        // Find pivot element
        int pivotRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(matrix[k][i]) > fabs(matrix[pivotRow][i])) {
                pivotRow = k;
            }
        }
        // Swap rows if necessary
        if (pivotRow != i) {
            swap(matrix[i], matrix[pivotRow]);
        }
        
        // Make the diagonal element 1 and eliminate below
        double pivot = matrix[i][i];
        for (int j = 0; j <= n; j++) {
            matrix[i][j] /= pivot;
        }
        for (int k = i + 1; k < n; k++) {
            double factor = matrix[k][i];
            for (int j = 0; j <= n; j++) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }

    // Back substitution
    vector<double> result(n);
    for (int i = n - 1; i >= 0; i--) {
        result[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) {
            result[i] -= matrix[i][j] * result[j];
        }
    }
    return result;
}

int main() {
    // Open and read JSON file
    ifstream inputFile("testcase1.json");
    if (!inputFile.is_open()) {
        cerr << "Failed to open file." << endl;
        return EXIT_FAILURE;
    }

    // Parse the JSON data from file
    json jsonData;
    try {
        inputFile >> jsonData;
    } catch (const json::exception& e) {
        cerr << "Failed to parse JSON: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    // Close the file
    inputFile.close();

    // Extract n and k from the JSON
    int n = jsonData["keys"]["n"];  
    int k = jsonData["keys"]["k"];  

    // Map to store the (x, y) points
    map<int, long long> points;

    // Loop through the JSON data to extract x and y values
    for (auto& el : jsonData.items()) {
        if (el.key() == "keys") continue;  

        int x = stoi(el.key());                        // Extract the x-value (key)
        int base = stoi(el.value()["base"].get<string>());  // Extract the base
        string y_str = el.value()["value"];            // Extract the y value as a string
        long long y = decodeBase(y_str, base);         // Decode y value using the base

        points[x] = y;  // Store x and decoded y value in the map
    }

    
    vector<pair<int, long long>> selected_points;
    int count = 0;
    for (const auto& point : points) {
        if (count >= k) break;  // We need only k points
        selected_points.push_back(point);
        count++;
    }

    if (selected_points.size() < k) {
        cerr << "Insufficient points provided." << endl;
        return EXIT_FAILURE;
    }

    // Construct the matrix for Gaussian elimination (k equations for (k-1)th degree polynomial)
    vector<vector<double>> matrix(k, vector<double>(k + 1, 0));  

    for (int i = 0; i < k; i++) {
        int x = selected_points[i].first;
        long long y = selected_points[i].second;

        for (int j = 0; j < k; j++) {
            matrix[i][j] = pow(x, k - j - 1);  
        }
        matrix[i][k] = y;  // y value
    }

    // Perform Gaussian elimination to solve for polynomial coefficients
    vector<double> result = gaussianElimination(matrix, k);

    
    cout << "The polynomial is: ";
    for (int i = 0; i < k; i++) {
        if (i > 0) cout << " + ";
        cout << result[i] << "x^" << (k - i - 1);
    }
    cout << " = " << result[k] << endl;

    return 0;
}
