#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include<iomanip>
#include <ctime>
using namespace std;

const int MATRIX_SIZE = 100;
int random_seed = 42;
double eps = 1e-10;

void printVector(const vector<double>& vector) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        cout << setw(10) << vector[i] << " ";
    }
    cout << endl;
}
void printMatrix(const vector<vector<double>>& matrix) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            cout << setw(10) << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void  copyVector(const vector<double>& orig, vector<double>& copy) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        copy[i] = orig[i];
    }
}
void  multiplyVector(vector<double>& F, const vector<vector<double>>& Q) {
    vector<double> ANS(MATRIX_SIZE, 0.0);
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int l = 0; l < MATRIX_SIZE; ++l) {
            ANS[i] += Q[i][l] * F[l];
        }
    }
    copyVector(ANS,F);
}

void  copyMatrix(const vector<vector<double>>& orig, vector<vector<double>>& copy) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            copy[i][j] = orig[i][j];
        }
    }
}
void  multiplyMatrix(vector<vector<double>>& R, const vector<vector<double>>& Q) {
    vector<vector<double>> ANS(MATRIX_SIZE, vector<double>(MATRIX_SIZE, 0.0));
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            for (int l = 0; l < MATRIX_SIZE; ++l) {
                ANS[i][j] += Q[i][l] * R[l][j];
                ANS[i][j] = fabs(ANS[i][j]) < eps ? 0 : ANS[i][j];
            }
        }
    }
    copyMatrix(ANS,R);
}

double maxDist(const vector<double>& vec1, const vector<double>& vec2){
    double max = 0, temp;
    for (int i = 0; i < MATRIX_SIZE; ++i){
        temp = fabs(vec1[i] - vec2[i]);
        max = max > temp ? max : temp;
    }
    return max;
}

void solveSystem(const vector<vector<double>>& A, const vector<double>& F, vector<double>& x) {
    vector<vector<double>> Q(MATRIX_SIZE, vector<double>(MATRIX_SIZE, 0.0));
    vector<vector<double>> R = A;
    vector<double> F1 = F;
    

    //прямой ход
    for (int k = 0; k < MATRIX_SIZE - 1; ++k) {
        vector<double> u(MATRIX_SIZE, 0.0);
        vector<double> v(MATRIX_SIZE, 0.0);
        vector<double> w(MATRIX_SIZE, 0.0);
        
        double normX = 0.0;
        for (int i = k; i < MATRIX_SIZE; ++i) {
            normX += R[i][k] * R[i][k];
        }
        normX = sqrt(normX);
        
        //v - вектор нормали
        //u - вспомогательный х - y или x + y
        u[k] = R[k][k] + (R[k][k] < 0 ? -normX : normX);
        for (int i = k + 1; i < MATRIX_SIZE; ++i) {
            u[i] = R[i][k];
        }
        
        double normU = 0.0;
        for (int i = k; i < MATRIX_SIZE; ++i) {
            normU += u[i] * u[i];
        }
        normU = sqrt(normU);
        
        for (int i = k; i < MATRIX_SIZE; ++i) {
            v[i] = u[i] / normU;
        }
        
        for (int i = 0; i < k; ++i){
            for (int j = i; j < MATRIX_SIZE; ++j){
                Q[i][j] = 0;
                Q[j][i] = 0;
            }
        }
        for (int i = k; i < MATRIX_SIZE; ++i){
            for (int j = k; j < MATRIX_SIZE; ++j){
                Q[i][j] = - 2 * v[i] * v[j];
            }
        }
        for (int i = 0; i < MATRIX_SIZE; ++i){
            Q[i][i] += 1;
        }
        multiplyMatrix(R, Q);
        multiplyVector(F1, Q);       
    }
    
    //обратный ход    
    x[MATRIX_SIZE - 1] = F1[MATRIX_SIZE - 1] / R[MATRIX_SIZE - 1][MATRIX_SIZE - 1];
    for (int i = MATRIX_SIZE - 2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < MATRIX_SIZE; ++j) {
            sum += R[i][j] * x[j];
        }
        x[i] = (F1[i] - sum) / R[i][i];
    }
    
}

double get_rand(){
    srand(random_seed);
    random_seed++;
    return (rand() % 201 - 100.0)/100.0;
}

int main() {
     time_t start, end;
     time(&start);
    vector<vector<double>> A(MATRIX_SIZE, vector<double>(MATRIX_SIZE, 0.0));
    vector<double> F(MATRIX_SIZE, 0.0);
    vector<double> x(MATRIX_SIZE, 0.0);
    vector<double> x_true(MATRIX_SIZE, 0.0);

    ifstream work_file("matrix11.csv");
    string line;
    string tmp;
    for(int i = 0; i < MATRIX_SIZE; ++i){
        getline(work_file, line);
        stringstream stream(line);
        for(int j = 0; j < MATRIX_SIZE; j++){
            getline(stream, tmp, ',');
            A[i][j] = stod(tmp);
        }
    }
    work_file.close();

    for(int i = 0; i < MATRIX_SIZE; i++){
        x_true[i] = get_rand();
    }

    for(int i = 0; i < MATRIX_SIZE; i++){
        double S = 0;
        for(int j = 0; j < MATRIX_SIZE; j++){
            S += A[i][j] * x_true[j];
        }
        F[i] = S;
    }
    
    solveSystem(A, F, x);
    cout << "x_true:" << endl;
    printVector(x_true);
    cout << "x найденные:" << endl;
    printVector(x);

    cout << "Максимум норма погрешности"<<endl;
    cout << maxDist(x, x_true) << endl;
    time(&end);
    double seconds = difftime(end, start);
    cout << "Время работы программы(в секундах):"<< endl;
    cout << seconds << endl;
    return 0;
}
