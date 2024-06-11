/*#include <iostream>
#include <string>
#include <stack>
#include <unordered_map>
#include <cctype> // for isdigit
#include <cmath>  // for pow
#include <bitset>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_set>*/
#include <bits/stdc++.h>
#include <complex>

using namespace std;

const int INF = numeric_limits<int>::max();
const double EPS = 1e-9;

class Node{
public:
    int data;
    Node* left, *right;

    Node(int val){
        data = val;
        left = right = NULL;
    }
};

bool isOperator(char ch) {
    return ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^';
}

int precedence(char op) {
    if (op == '^') return 3;       // Highest precedence for power operator
    if (op == '*' || op == '/') return 2;
    if (op == '+' || op == '-') return 1;
    return 0;
}

double applyOp(double a, double b, char op) {
    switch(op) {
        case '+': return a + b;
        case '-': return a - b;
        case '*': return a * b;
        case '/':
            if (b == 0) {
                throw runtime_error("Error: Division by zero");
            }
            return a / b;
        case '^': return pow(a, b); // Power operation
    }
    return 0;
}

double evaluateExpression(const string& expression) {
    stack<double> operands;
    stack<char> operators;
    unordered_map<char, int> precedenceMap;

    precedenceMap['+'] = 1;
    precedenceMap['-'] = 1;
    precedenceMap['*'] = 2;
    precedenceMap['/'] = 2;
    precedenceMap['^'] = 3; // Higher precedence for power operator

    char lastToken = ' '; // Initialize lastToken to a space

    for (int i = 0; i < expression.length(); ++i) {
        if (expression[i] == ' ')
            continue;

        char currentToken = expression[i];

        if (isdigit(currentToken) || currentToken == '.') {
            double operand = 0;
            double decimalMultiplier = 0.1; // To handle digits after decimal point
            bool hasDecimal = false;

            while (i < expression.length() && (isdigit(expression[i]) || expression[i] == '.')) {
                if (expression[i] == '.') {
                    if (hasDecimal) {
                        throw runtime_error("Error: Invalid number format");
                    }
                    hasDecimal = true;
                } else {
                    if (!hasDecimal) {
                        operand = operand * 10 + (expression[i] - '0');
                    } else {
                        operand += (expression[i] - '0') * decimalMultiplier;
                        decimalMultiplier *= 0.1;
                    }
                }
                i++;
            }
            i--; // Adjust index after reading the number
            operands.push(operand);
        }
        else if (currentToken == '(') {
            operators.push(currentToken);
        }
        else if (currentToken == ')') {
            while (!operators.empty() && operators.top() != '(') {
                char op = operators.top();
                operators.pop();
                double operand2 = operands.top();
                operands.pop();
                double operand1 = operands.top();
                operands.pop();
                operands.push(applyOp(operand1, operand2, op));
            }
            if (operators.empty()) {
                throw runtime_error("Error: Unmatched closing parenthesis");
            }
            operators.pop(); // Pop '('
        }
        else if (isOperator(currentToken)) {
            // Check for consecutive operators without operand
            if (lastToken == '+' || lastToken == '-' || lastToken == '*' || lastToken == '/' || lastToken == '^') {
                throw runtime_error("Error: Consecutive operators without operand");
            }

            while (!operators.empty() && precedence(operators.top()) >= precedence(currentToken)) {
                char op = operators.top();
                operators.pop();
                double operand2 = operands.top();
                operands.pop();
                double operand1 = operands.top();
                operands.pop();
                operands.push(applyOp(operand1, operand2, op));
            }
            operators.push(currentToken);
        }

        lastToken = currentToken;
    }

    // Apply remaining operators
    while (!operators.empty()) {
        char op = operators.top();
        operators.pop();
        double operand2 = operands.top();
        operands.pop();
        double operand1 = operands.top();
        operands.pop();
        operands.push(applyOp(operand1, operand2, op));
    }

    if (operands.size() != 1) {
        throw runtime_error("Error: Operand/Operator mismatch");
    }

    return operands.top();
}

// Check if the character is a bitwise operator
bool isBitwiseOperator(char ch) {
    return ch == '&' || ch == '|' || ch == '~' || ch == '<' || ch == '>' || ch == '^';
}

// Define the precedence of each bitwise operator
int bitwisePrecedence(char op) {
    if (op == '~') return 4;  // Highest precedence for bitwise NOT (unary)
    if (op == '^') return 3;  // Higher precedence for bitwise XOR
    if (op == '&') return 2;
    if (op == '|' || op == '<' || op == '>') return 1; // Lowest precedence
    return 0;
}

// Apply the bitwise operation to the operands
int applyBitwiseOp(int a, int b, char op) {
    switch(op) {
        case '&': return a & b;
        case '|': return a | b;
        case '^': return a ^ b;
        case '<': return a << b; // Left shift
        case '>': return a >> b; // Right shift
    }
    return 0;
}

// Evaluate the bitwise expression
int evaluateBitwiseExpression(const string& expression) {
    stack<int> operands;
    stack<char> operators;

    char lastToken = ' '; // Initialize lastToken to a space

    for (int i = 0; i < expression.length(); ++i) {
        if (expression[i] == ' ')
            continue;

        char currentToken = expression[i];

        if (isdigit(currentToken)) {
            int operand = 0;

            while (i < expression.length() && isdigit(expression[i])) {
                operand = operand * 10 + (expression[i] - '0');
                i++;
            }
            i--; // Adjust index after reading the number
            operands.push(operand);
        }
        else if (currentToken == '(') {
            operators.push(currentToken);
        }
        else if (currentToken == ')') {
            while (!operators.empty() && operators.top() != '(') {
                char op = operators.top();
                operators.pop();
                int operand2 = operands.top();
                operands.pop();
                if (op == '~') {
                    operands.push(~operand2);
                } else {
                    int operand1 = operands.top();
                    operands.pop();
                    operands.push(applyBitwiseOp(operand1, operand2, op));
                }
            }
            if (operators.empty()) {
                throw runtime_error("Error: Unmatched closing parenthesis");
            }
            operators.pop(); // Pop '('
        }
        else if (isBitwiseOperator(currentToken)) {
            // Check for consecutive operators without operand
            if (lastToken == '&' || lastToken == '|' || lastToken == '~' || lastToken == '<' || lastToken == '>' || lastToken == '^') {
                throw runtime_error("Error: Consecutive operators without operand");
            }

            while (!operators.empty() && bitwisePrecedence(operators.top()) >= bitwisePrecedence(currentToken)) {
                char op = operators.top();
                operators.pop();
                int operand2 = operands.top();
                operands.pop();
                if (op == '~') {
                    operands.push(~operand2);
                } else {
                    int operand1 = operands.top();
                    operands.pop();
                    operands.push(applyBitwiseOp(operand1, operand2, op));
                }
            }
            operators.push(currentToken);
        }

        lastToken = currentToken;
    }

    // Apply remaining operators
    while (!operators.empty()) {
        char op = operators.top();
        operators.pop();
        int operand2 = operands.top();
        operands.pop();
        if (op == '~') {
            operands.push(~operand2);
        } else {
            int operand1 = operands.top();
            operands.pop();
            operands.push(applyBitwiseOp(operand1, operand2, op));
        }
    }

    if (operands.size() != 1) {
        throw runtime_error("Error: Operand/Operator mismatch");
    }

    return operands.top();
}

// Convert decimal number to binary string
string decimalToBinary(int n) {
    return bitset<32>(n).to_string();
}

// Convert binary string to decimal number
int binaryToDecimal(const string& binary) {
    return stoi(binary, nullptr, 2);
}

// Handle bitwise operations menu
void bitwise() {
    cout << "1. Evaluate expression\n";
    cout << "2. Decimal to binary\n";
    cout << "3. Binary to decimal\n\n";

    int n;
    cin >> n;
    cin.ignore();

    if (n == 1) {
        string expression;
        cout << "Enter a bitwise expression: ";
        getline(cin, expression);

        try {
            int result = evaluateBitwiseExpression(expression);
            cout << "Result: " << result << endl;
        } catch (const exception& e) {
            cerr << e.what() << endl;
        }
    }
    else if (n == 2) {
        int decimal;
        cout << "Enter a decimal number: ";
        cin >> decimal;
        cout << "Binary: " << decimalToBinary(decimal) << endl;
    }
    else if (n == 3) {
        string binary;
        cout << "Enter a binary number: ";
        cin >> binary;
        cout << "Decimal: " << binaryToDecimal(binary) << endl;
    }
}

// Function to perform union of two sets
vector<int> setUnion(const vector<int>& set1, const vector<int>& set2) {
    vector<int> result = set1;
    for (int elem : set2) {
        if (find(set1.begin(), set1.end(), elem) == set1.end()) {
            result.push_back(elem);
        }
    }
    return result;
}

// Function to perform intersection of two sets
vector<int> setIntersection(const vector<int>& set1, const vector<int>& set2) {
    vector<int> result;
    for (int elem : set2) {
        if (find(set1.begin(), set1.end(), elem) != set1.end()) {
            result.push_back(elem);
        }
    }
    return result;
}

// Function to perform difference of two sets (set1 - set2)
vector<int> setDifference(const vector<int>& set1, const vector<int>& set2) {
    vector<int> result;
    for (int elem : set1) {
        if (find(set2.begin(), set2.end(), elem) == set2.end()) {
            result.push_back(elem);
        }
    }
    return result;
}

// Function to perform symmetric difference of two sets
vector<int> setSymmetricDifference(const vector<int>& set1, const vector<int>& set2) {
    vector<int> result = setUnion(set1, set2);
    vector<int> intersect = setIntersection(set1, set2);
    for (int elem : intersect) {
        result.erase(remove(result.begin(), result.end(), elem), result.end());
    }
    return result;
}

// Helper function to read a set from the user
vector<int> readSet() {
    int n;
    cout << "Enter the number of elements in the set: ";
    cin >> n;
    vector<int> set(n);
    cout << "Enter the elements of the set: ";
    for (int i = 0; i < n; ++i) {
        cin >> set[i];
    }
    return set;
}

// Function to generate all subsets of a set
vector<vector<int>> generateSubsets(const vector<int>& set) {
    vector<vector<int>> subsets;
    int n = set.size();
    int subsetCount = 1 << n; // 2^n subsets

    for (int mask = 0; mask < subsetCount; ++mask) {
        vector<int> subset;
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {
                subset.push_back(set[i]);
            }
        }
        subsets.push_back(subset);
    }
    return subsets;
}

// Function to handle set operations menu
void set_operation() {
    cout << "1. Union\n";
    cout << "2. Intersection\n";
    cout << "3. Difference\n";
    cout << "4. Symmetric Difference\n";
    cout << "5. Generate Subsets\n\n";

    int n;
    cin >> n;
    cin.ignore();

    if (n == 5) {
        vector<int> set = readSet();
        vector<vector<int>> subsets = generateSubsets(set);

        cout << "Subsets:\n";
        for (const auto& subset : subsets) {
            cout << "{ ";
            for (int elem : subset) {
                cout << elem << " ";
            }
            cout << "}\n";
        }
        return;
    }

    vector<int> set1 = readSet();
    vector<int> set2 = readSet();

    vector<int> result;
    if (n == 1) {
        result = setUnion(set1, set2);
    }
    else if (n == 2) {
        result = setIntersection(set1, set2);
    }
    else if (n == 3) {
        result = setDifference(set1, set2);
    }
    else if (n == 4) {
        result = setSymmetricDifference(set1, set2);
    }

    cout << "Result: ";
    for (int elem : result) {
        cout << elem << " ";
    }
    cout << endl;
}

// Function to calculate the Greatest Common Divisor (GCD) of two numbers
int gcd(int a, int b) {
    return (b == 0) ? a : gcd(b, a % b);
}

// Function to calculate the Least Common Multiple (LCM) of two numbers
int lcm(int a, int b) {
    return (a * b) / gcd(a, b);
}

// Function to generate prime factors of a number
vector<int> primeFactors(int n) {
    vector<int> factors;
    while (n % 2 == 0) {
        factors.push_back(2);
        n /= 2;
    }
    for (int i = 3; i <= sqrt(n); i += 2) {
        while (n % i == 0) {
            factors.push_back(i);
            n /= i;
        }
    }
    if (n > 2) {
        factors.push_back(n);
    }
    return factors;
}

// Function to perform modular arithmetic for addition, subtraction, and multiplication
int modularArithmetic(int a, int b, int mod, char op) {
    int result;
    switch (op) {
        case '+':
            result = (a % mod + b % mod) % mod;
            break;
        case '-':
            result = (a % mod - b % mod + mod) % mod;
            break;
        case '*':
            result = (a % mod * b % mod) % mod;
            break;
        default:
            cout << "Invalid operation\n";
            return -1;
    }
    return result;
}

// Function to handle number theory operations
void num_theory_operation() {
    cout << "1. Calculate GCD\n";
    cout << "2. Calculate LCM\n";
    cout << "3. Prime Factorization\n";
    cout << "4. Perform Modular Arithmetic\n\n";

    int n;
    cin >> n;
    cin.ignore();

    if (n == 1) {
        int a, b;
        cout << "Enter two numbers: ";
        cin >> a >> b;
        cout << "GCD: " << gcd(a, b) << endl;
    }
    else if (n == 2) {
        int a, b;
        cout << "Enter two numbers: ";
        cin >> a >> b;
        cout << "LCM: " << lcm(a, b) << endl;
    }
    else if (n == 3) {
        int num;
        cout << "Enter a number: ";
        cin >> num;
        cout << "Prime Factors: ";
        vector<int> factors = primeFactors(num);
        for (int factor : factors) {
            cout << factor << " ";
        }
        cout << endl;
    }
    else if (n == 4) {
        int a, b, mod;
        char op;
        cout << "Enter three numbers (a, b, mod) and operation (+, -, *): ";
        cin >> a >> b >> mod >> op;
        cout << "Result of (a " << op << " b) % mod: " << modularArithmetic(a, b, mod, op) << endl;
    }
    else {
        cout << "Invalid choice\n";
    }
}

// Function to calculate factorial
int factorial(int n) {
    if (n == 0 || n == 1) {
        return 1;
    }
    return n * factorial(n - 1);
}

// Function to calculate permutations
int permutations(int n, int r) {
    return factorial(n) / factorial(n - r);
}

// Function to calculate combinations
int combinations(int n, int r) {
    return factorial(n) / (factorial(r) * factorial(n - r));
}

// Function to calculate binomial coefficient
int binomialCoefficient(int n, int r) {
    return factorial(n) / (factorial(r) * factorial(n - r));
}

// Function to calculate basic probability
double basicProbability(int favorable, int total) {
    return (double)favorable / total;
}

// Function to calculate conditional probability P(A|B) = P(A and B) / P(B)
double conditionalProbability(double pAandB, double pB) {
    return pAandB / pB;
}

// Function to calculate probability using Bayes' Theorem
// P(A|B) = (P(B|A) * P(A)) / P(B)
double bayesTheorem(double pBA, double pA, double pB) {
    return (pBA * pA) / pB;
}

// Function to handle combinatorics and probability operations
void combi() {
    cout << "1. Calculate Factorial\n";
    cout << "2. Calculate Permutations (nPr)\n";
    cout << "3. Calculate Combinations (nCr)\n";
    cout << "4. Calculate Binomial Coefficient\n";
    cout << "5. Calculate Basic Probability\n";
    cout << "6. Calculate Conditional Probability\n";
    cout << "7. Apply Bayes' Theorem\n";
    cout << "8. Generate All n-Permutations\n\n";

    int n, r, favorable, total;
    double pAandB, pB, pBA, pA;
    int choice;
    cin >> choice;
    cin.ignore();

    if (choice == 1) {
        cout << "Enter a number: ";
        cin >> n;
        cout << "Factorial: " << factorial(n) << endl;
    }
    else if (choice == 2) {
        cout << "Enter n and r: ";
        cin >> n >> r;
        cout << "Permutations (nPr): " << permutations(n, r) << endl;
    }
    else if (choice == 3) {
        cout << "Enter n and r: ";
        cin >> n >> r;
        cout << "Combinations (nCr): " << combinations(n, r) << endl;
    }
    else if (choice == 4) {
        cout << "Enter n and r: ";
        cin >> n >> r;
        cout << "Binomial Coefficient: " << binomialCoefficient(n, r) << endl;
    }
    else if (choice == 5) {
        cout << "Enter the number of favorable outcomes and total outcomes: ";
        cin >> favorable >> total;
        cout << "Basic Probability: " << basicProbability(favorable, total) << endl;
    }
    else if (choice == 6) {
        cout << "Enter P(A and B) and P(B): ";
        cin >> pAandB >> pB;
        cout << "Conditional Probability P(A|B): " << conditionalProbability(pAandB, pB) << endl;
    }
    else if (choice == 7) {
        cout << "Enter P(B|A), P(A), and P(B): ";
        cin >> pBA >> pA >> pB;
        cout << "Bayes' Theorem P(A|B): " << bayesTheorem(pBA, pA, pB) << endl;
    }
    else if (choice == 8) {
        cout << "Enter n: ";
        cin >>n;
        vector<int> permutation;
        for (int i = 1; i <= n; i++) {
                permutation.push_back(i);
        }
        do {
            for (auto x:permutation) cout<<x<<" ";
            cout<<"\n";
        } while (next_permutation(permutation.begin(),permutation.end()));
    }
    else {
        cout << "Invalid choice\n";
    }
}

// Function to read a matrix from the user
vector<vector<int>> readMatrix(int rows, int cols) {
    vector<vector<int>> matrix(rows, vector<int>(cols));
    cout << "Enter the matrix elements row-wise:\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cin >> matrix[i][j];
        }
    }
    return matrix;
}

// Function to print a matrix
void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }
}

// Function to add two matrices
vector<vector<int>> addMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<int>> result(rows, vector<int>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

// Function to subtract two matrices
vector<vector<int>> subtractMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<int>> result(rows, vector<int>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

// Function to multiply two matrices
vector<vector<int>> multiplyMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int colsB = B[0].size();
    vector<vector<int>> result(rowsA, vector<int>(colsB, 0));
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Function to transpose a matrix
vector<vector<int>> transposeMatrix(const vector<vector<int>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<vector<int>> result(cols, vector<int>(rows));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

// Function to calculate the determinant of a 2x2 matrix
int determinant2x2(const vector<vector<int>>& matrix) {
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

// Function to calculate the determinant of a 3x3 matrix
int determinant3x3(const vector<vector<int>>& matrix) {
    int det = 0;
    det += matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
    det -= matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]);
    det += matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    return det;
}

// Function to calculate matrix power
vector<vector<int>> matrixPower(const vector<vector<int>>& matrix, int power) {
    int n = matrix.size();
    vector<vector<int>> result(n, vector<int>(n, 0));
    vector<vector<int>> base = matrix;

    // Initialize result as identity matrix
    for (int i = 0; i < n; i++) {
        result[i][i] = 1;
    }

    while (power > 0) {
        if (power % 2 == 1) {
            result = multiplyMatrices(result, base);
        }
        base = multiplyMatrices(base, base);
        power /= 2;
    }

    return result;
}

void mat_operation() {
    cout << "Matrix Operations:\n";
    cout << "1. Add Matrices\n";
    cout << "2. Subtract Matrices\n";
    cout << "3. Multiply Matrices\n";
    cout << "4. Transpose Matrix\n";
    cout << "5. Matrix Power\n";
    cout << "6. Determinant\n\n";

    int choice;
    cin >> choice;

    if (choice == 1 || choice == 2 || choice == 3) {
        int rowsA, colsA, rowsB, colsB;
        cout << "Enter number of rows and columns of first matrix: ";
        cin >> rowsA >> colsA;
        vector<vector<int>> A = readMatrix(rowsA, colsA);

        cout << "Enter number of rows and columns of second matrix: ";
        cin >> rowsB >> colsB;
        vector<vector<int>> B = readMatrix(rowsB, colsB);

        if ((choice == 1 || choice == 2) && (rowsA != rowsB || colsA != colsB)) {
            cout << "Error: Matrices must have the same dimensions for addition/subtraction.\n";
            return;
        }

        if (choice == 1) {
            vector<vector<int>> result = addMatrices(A, B);
            printMatrix(result);
        } else if (choice == 2) {
            vector<vector<int>> result = subtractMatrices(A, B);
            printMatrix(result);
        } else if (choice == 3) {
            if (colsA != rowsB) {
                cout << "Error: Number of columns of first matrix must be equal to number of rows of second matrix for multiplication.\n";
                return;
            }
            vector<vector<int>> result = multiplyMatrices(A, B);
            printMatrix(result);
        }
    } else if (choice == 4) {
        int rows, cols;
        cout << "Enter number of rows and columns of the matrix: ";
        cin >> rows >> cols;
        vector<vector<int>> matrix = readMatrix(rows, cols);
        vector<vector<int>> result = transposeMatrix(matrix);
        printMatrix(result);
    } else if (choice == 5) {
        int rows, cols, power;
        cout << "Enter number of rows and columns of the matrix (must be square): ";
        cin >> rows >> cols;
        if (rows != cols) {
            cout << "Error: Matrix must be square for exponentiation.\n";
            return;
        }
        vector<vector<int>> matrix = readMatrix(rows, cols);
        cout << "Enter the power: ";
        cin >> power;
        vector<vector<int>> result = matrixPower(matrix, power);
        printMatrix(result);
    } else if (choice == 6) {
        int rows, cols;
        cout << "Enter number of rows and columns of the matrix (must be 2x2 or 3x3): ";
        cin >> rows >> cols;
        if ((rows != cols) || (rows != 2 && rows != 3)) {
            cout << "Error: Only 2x2 and 3x3 matrices are supported for determinant calculation.\n";
            return;
        }
        vector<vector<int>> matrix = readMatrix(rows, cols);
        int det;
        if (rows == 2) {
            det = determinant2x2(matrix);
        } else {
            det = determinant3x3(matrix);
        }
        cout << "Determinant: " << det << endl;
    } else {
        cout << "Invalid choice.\n";
    }
}

const int N=300005;
vector<pair<int,int>> adj[N];

void printGraph(int n) {
    cout << "\nAdjacency List:\n";
    for (int i = 1; i <= n; ++i) {
        cout << "Node " << i << ":";
        for (const auto& edge : adj[i]) {
            cout << "  (" << edge.first << ", " << edge.second << ") ";
        }
        cout << endl;
    }
}

void BFS(int start) {
    vector<bool> visited(N, false);
    queue<int> q;
    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        cout << node << " ";

        for (auto& neighbor : adj[node]) {
            int nextNode = neighbor.first;
            if (!visited[nextNode]) {
                q.push(nextNode);
                visited[nextNode] = true;
            }
        }
    }
    cout << endl;
}

void DFSprint(int node, vector<bool>& visited) {
    visited[node] = true;
    cout << node << " ";

    for (auto& neighbor : adj[node]) {
        int u = neighbor.first;
        if (!visited[u]) {
            DFSprint(u, visited);
        }
    }
}
void DFS(int start) {
    vector<bool> visited(N, false);
    DFSprint(start, visited);
    cout<<"\n";
}

pair<vector<int>, vector<int>> dijkstra(int source) {
    vector<int> distance(N, INF); // Initialize distances to INF
    vector<int> predecessor(N, -1); // Initialize predecessors to -1
    vector<bool> visited(N, false); // Initialize visited array

    // Priority queue to store {distance, node} pairs
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    distance[source] = 0; // Distance from source to itself is 0
    pq.push({0, source}); // Push the source node into the priority queue

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        // Skip if the node is already visited
        if (visited[u]) continue;

        // Mark the node as visited
        visited[u] = true;

        // Update distances to all the adjacent nodes of u
        for (auto& edge : adj[u]) {
            int v = edge.first;
            int weight = edge.second;
            if (distance[u] + weight < distance[v]) {
                distance[v] = distance[u] + weight;
                predecessor[v] = u; // Update predecessor of v
                pq.push({distance[v], v});
            }
        }
    }

    return {distance, predecessor};
}

void GRAPH() {
    //Storing graph
    cout<<"Store Graph:\n";
    int n,m;
    cout<<"Enter number of nodes: "; cin>>n;
    cout<<"Enter number of edges: "; cin>>m;
    int x,y,w;
    cout << "Enter the edges in the format: node1 node2 weight\n";
    for (int i=0;i<m;i++) {
        cin>>x>>y>>w;
        adj[x].push_back({y,w});
        adj[y].push_back({x,w});
    }

    //Selecting operation on that graph
    cout << "1. Print the graph\n";
    cout << "2. Perform BFS\n";
    cout << "3. Perform DFS\n";
    cout << "4. Find shortest path\n";

    int choice; cin>>choice;
    if (choice==1) {
        printGraph(n);
        cout<<"\n";
    }
    if (choice==2) {
        cout << "Enter the starting node for BFS: ";
        int start; cin >> start;
        BFS(start);
    }
    if (choice==3) {
        cout << "Enter the starting node for DFS: ";
        int start; cin>>start;
        DFS(start);
    }
    if (choice==4) {
        //calculate shortest path from one node to another using dijkstra
        cout << "Enter the source node: ";
        int source; cin >> source;
        cout << "Enter the destination node: ";
        int destination; cin >> destination;

        pair<vector<int>, vector<int>> result = dijkstra(source);
        vector<int> distance = result.first;
        vector<int> predecessor = result.second;

        if (distance[destination] == INF) {
            cout << "There is no path from " << source << " to " << destination << endl;
        } else {
            cout << "Shortest distance from " << source << " to " << destination << " is " << distance[destination] << endl;
            // Print the shortest path
            stack<int> path;
            int current = destination;
            while (current != -1) {
                path.push(current);
                current = predecessor[current];
            }
            cout << "Shortest path: ";
            while (!path.empty()) {
                cout << path.top() << " ";
                path.pop();
            }
            cout << endl;
            }
        }
}



void solveLinearEquations(vector<vector<double>>& mat) {
    int n = mat.size();
    cout << "Solution x:" << endl;
    switch (n) {
        case 2: {
            double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            if (abs(det) < EPS) {
                cout << "No unique solution exists. Determinant is zero." << endl;
                return;
            }
            double x1 = (mat[1][1] * mat[0][2] - mat[0][1] * mat[1][2]) / det;
            double x2 = (mat[0][0] * mat[1][2] - mat[1][0] * mat[0][2]) / det;
            cout << "x1 = " << x1 << ", x2 = " << x2 << endl;
            break;
        }
        case 3: {
            double det = mat[0][0] * mat[1][1] * mat[2][2]
                       + mat[0][1] * mat[1][2] * mat[2][0]
                       + mat[0][2] * mat[1][0] * mat[2][1]
                       - mat[0][2] * mat[1][1] * mat[2][0]
                       - mat[0][1] * mat[1][0] * mat[2][2]
                       - mat[0][0] * mat[1][2] * mat[2][1];
            if (abs(det) < EPS) {
                cout << "No unique solution exists. Determinant is zero." << endl;
                return;
            }
            double x1 = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / det;
            double x2 = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) / det;
            double x3 = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / det;
            cout << "x1 = " << x1 << ", x2 = " << x2 << ", x3 = " << x3 << endl;
            break;
        }
        default:
            cout << "Unsupported number of variables. Please use 2 to 3 variables." << endl;
            break;
    }
}

void solveCubicEquation(double a, double b, double c, double d) {
    // Discriminant
    double delta0 = b * b - 3 * a * c;
    double delta1 = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d;

    // Cubic roots
    double C;
    if (delta0 == 0 && delta1 == 0) {
        C = 0;
    } else {
        double delta = (delta1 * delta1) - (4 * delta0 * delta0 * delta0);
        double u = (cbrt(delta1 + sqrt(delta))) / 2;
        double v = (cbrt(delta1 - sqrt(delta))) / 2;
        C = u + v;
    }

    // Calculate roots
    double x1, x2, x3;
    x1 = (-b - C + sqrt(delta0 + 2 * b * C - 3 * a)) / (3 * a);
    x2 = (-b + (0.5 * (1 + sqrt(-3) * sqrt(-1)) * (delta0 + 2 * b * C - 3 * a))) / (3 * a);
    x3 = (-b + (0.5 * (1 - sqrt(-3) * sqrt(-1)) * (delta0 + 2 * b * C - 3 * a))) / (3 * a);

    cout << "Solutions: x1 = " << x1 << ", x2 = " << x2 << ", x3 = " << x3 << endl;
}

void eqn() {
    //options: Select Vartiable
    cout<<"Select variable: ";
    int n; cin>>n;
    if (n==1) {
        cout<<"Select power: ";
        int power; cin>>power;
        if (power==1) {
            // Equation: ax = c
            double a, c;
            cout << "Enter coefficient 'a': ";
            cin >> a;
            cout << "Enter constant 'c': ";
            cin >> c;
            double x = c / a;
            cout << "Solution: x = " << x << endl;
        }
        if (power==2) {
            // Equation: ax^2 + bx = c
            double a, b, c;
            cout << "Enter coefficient 'a': ";
            cin >> a;
            cout << "Enter coefficient 'b': ";
            cin >> b;
            cout << "Enter constant 'c': ";
            cin >> c;
            double D = b * b - 4 * a * c;
            if (D < 0) {
                cout << "No real solutions." << endl;
            } else {
                double x1 = (-b + sqrt(D)) / (2 * a);
                double x2 = (-b - sqrt(D)) / (2 * a);
                cout << "Solutions: x1 = " << x1 << ", x2 = " << x2 << endl;
            }
        }
        if (power==3) {
            // Equation: ax^3 + bx^2 + cx = d
            double a, b, c, d;
            cout << "Enter coefficient 'a': ";
            cin >> a;
            cout << "Enter coefficient 'b': ";
            cin >> b;
            cout << "Enter coefficient 'c': ";
            cin >> c;
            cout << "Enter constant 'd': ";
            cin >> d;
            solveCubicEquation(a, b, c, d);
        }
        else {
            cout<<"Too much power!\n";
        }
    }
    else if (n>=2 && n<=3) {
        vector<vector<double>> mat(n, vector<double>(n + 1));

        cout << "Enter the coefficients:" << endl;
        for (int i = 0; i < n; ++i) {
            cout << "Equation " << i+1 << ": ";
            for (int j = 0; j <= n; ++j) {
                cin >> mat[i][j];
            }
        }
        solveLinearEquations(mat);
    }
    else cout<<"Too many variables!\n";
}

int main() {
    while (1) {
        cout << "Choose action:\n\n";
        cout << "1. Evaluate arithmetic expression\n";
        cout << "2. Bitwise operation\n";
        cout << "3. Set operation\n";
        cout << "4. Matrix and Determinant\n";
        cout << "5. Number theory operations\n";
        cout << "6. Combinatorics and Probability\n";
        cout << "7. Graphs and Trees\n";
        cout << "8. Solving equation\n\n";

        int n;
        cin >> n;
        cin.ignore();  // Ignore the newline character left in the input buffer

        if (n == 1) {
            string expression;
            cout << "Enter an expression: ";
            getline(cin, expression);

            try {
                double result = evaluateExpression(expression);
                cout << "Result: " << result << endl;
            } catch (const exception& e) {
                cerr << e.what() << endl;
            }
        }
        else if (n==2) {
            bitwise();
            cout<<"\n";
        }
        else if (n==3) {
            set_operation();
            cout<<"\n";
        }
        else if (n==4) {
            mat_operation();
            cout<<"\n";
        }
        else if (n==5) {
            num_theory_operation();
            cout<<"\n";
        }
        else if (n==6) {
            combi();
            cout<<"\n";
        }
        else if (n==7) {
            cout<<"1. Tree operation\n";
            cout<<"2. Graph operation\n";
            int k; cin>>k;
            if(k==1){
                //TREE();
            }
            if (k==2) {
                GRAPH();
            }
            else cout<<"Invalid operation\n";
            cout<<"\n";
        }
        else if (n==8) {
            eqn();
            cout<<"\n";
        }
        else {
            cout<< "Invalid operation\n";
        }
    }

    return 0;
}
