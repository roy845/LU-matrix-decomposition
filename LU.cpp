#include <iostream>
#include <fstream> // For using files
#include <string>  //for using string operations
#include <vector>  //for vector temp
#include <chrono> // For measuring time of executions

using namespace std;

/**Program to calculate LU decomposition of a square matrix

 Input:Square matrix
 Output:4 files L.txt , L.csv and U.txt , U.csv containing: Lower triangular matrix L , Upper triangular matrix U

**/



//Functions prototypes
void initMatrices(double**& L, double**& U, int n);
void LUdecomposition(double**& A, double**& L, double**& U, int n);
void saveResults(double** matrix, int n, string fileName);
double** readFile(int& n, string fileName);


//Function to initialize matrices: L:lower , U:upper ---> size:nxn
void initMatrices(double**& L, double**& U, int n) {

    L = new double* [n];
    U = new double* [n];
    for (int i = 0; i < n; i++)
    {
        L[i] = new double[n];
        U[i] = new double[n];
    }
}

//LU decomposition fucntion
void LUdecomposition(double**& A, double**& L, double**& U, int n) {

    double sum;
    //first start by initializing matrices L:lower , U:upper --->size:nxn 
    initMatrices(L, U, n);

    for (int i = 0; i < n; i++)
    {
        //Calculation of Upper Triangular matrix
        for (int k = i; k < n; k++)
        {
            // Calculating the Sum of L(i, j) * U(j, k)
            sum = 0.0;
            for (int j = 0; j < i; j++)
                sum = sum + (L[i][j] * U[j][k]);

            //Calculating U(i, k)
            U[i][k] = A[i][k] - sum;

            //Putting zeros in the lower triangular
            if (i < k) {
                U[k][i] = 0;
            }
        }

        //Calculation of lower Triangular matrix
        for (int k = i; k < n; k++)
        {
            if (i == k)
                L[i][i] = 1; //Putting ones in the main diagonal
            else
            {
                //Calculating the Sum of L(k, j) * U(j, i)
                sum = 0.0;
                for (int j = 0; j < i; j++)
                    sum = sum + (L[k][j] * U[j][i]);

                //Calculating L(k, i)
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
            //Putting zeros in the upper triangular
            if (k > i) {
                L[i][k] = 0;
            }
        }
    }
}

//Function to read matrix from txt file
double** readFile(int& n, string fileName) {


    string line;
    string delimiter = ",";
    size_t pos = 0;
    string num_str;

    vector<double> temp;
    // Open the file to read
    ifstream dataFile(fileName);
    //if fail to read
    if (dataFile.fail()) {
        cerr << "Error opening file: " << fileName << endl;
        exit(1);
    }

    cout << "Reading the data from file: " << fileName << endl;

    //get data from file line by line and count the number of rows.
    while (getline(dataFile, line))
    {
        //number of rows in file
        n++;


        while ((pos = line.find(delimiter)) != string::npos)// --> get numbers from line
        {

            num_str = line.substr(0, pos);//-->number in line

            temp.push_back(stod(num_str));//-->push back the number in line to vector
            line.erase(0, pos + delimiter.length());

        }

        temp.push_back(stod(line));
    }
    cout << "Finished reading the data from file: " << fileName << endl;

    cout << "Number of lines read:" << n << endl;

    //Create matrix and fill it with data from file
    double** matrix = new double* [n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            matrix[i][j] = temp[j + i * n];


        }
    }

    return matrix;
}

//Function to save results in file
void saveResults(double** matrix, int n, string fileName) {

    cout << "Saving result in: " << fileName << endl;
    ofstream output(fileName);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            //if not the end of the line write the number to file with -->,
            if (j != n - 1) {
                output << matrix[i][j] << ",";
            }
            //make sure not to write to file -->, at the end of the line
            if (j == n - 1) {
                output << matrix[i][j];
            }
        }
        //write new line at the end of the current line
        output << "\n";
    }
}



int main() {

    //initializations
    int n = 0; //Matrix size
    double** A; //Matrix to decompose
    double** L; //Lower Triangular matrix
    double** U; //Upper Triangular matrix
    string fileName = "data_10.txt";

    //Measuring execution time of reading matrix from txt file
    auto start0 = chrono::high_resolution_clock::now();
    A = readFile(n, fileName); //Reading matrix from txt file
    auto elapsed0 = chrono::high_resolution_clock::now() - start0;//how much time elapsed from start0
    auto secondsTimeRead = std::chrono::duration_cast<std::chrono::seconds>(elapsed0).count();
    cout << "time readFile: " << secondsTimeRead << "s" << endl;



    //Measuring execution time of LU decomposition 
    auto start1 = chrono::high_resolution_clock::now();
    cout << "Starting LU Decomposition..." << endl;
    LUdecomposition(A, L, U, n);
    auto elapsed1 = chrono::high_resolution_clock::now() - start1;
    auto secondsLU = std::chrono::duration_cast<std::chrono::seconds>(elapsed1).count();
    cout << "time LU decomposition: " << secondsLU << "s" << endl;



    //Measuring execution time of saving results to file
    auto start2 = chrono::high_resolution_clock::now();
    saveResults(L, n, "L.txt");
    saveResults(U, n, "U.txt");
    saveResults(L, n, "L.csv");
    saveResults(U, n, "U.csv");
    auto elapsed2 = chrono::high_resolution_clock::now() - start2;
    auto secondsSaveResults = std::chrono::duration_cast<std::chrono::seconds>(elapsed2).count();
    cout << "time saveResults: " << secondsSaveResults << "s" << endl;


    cout << "Finished" << endl;



    char ch = getchar();
    return 0;
}


