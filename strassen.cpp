//------------------------strassen.cpp-----------------------------------------
// Author: Noah Nguyen
// Created: May 27, 2023
//-----------------------------------------------------------------------------
// Implementing the Strassen's Algorithm
//-----------------------------------------------------------------------------

#include <iostream>
#include <vector>
using namespace std;

/** 2 n x n matrices are presented by
 *      vector<vector<int>> A
 *      vector<vector<int>> B
 *      A.size() == B.size() == n
 */

//--------------------------subtractMatrix()-----------------------------------
// perform matrix subtraction
// C = A - B
//-----------------------------------------------------------------------------
vector<vector<int>> subtractMatrix(vector<vector<int>> A,
                                   vector<vector<int>> B) {
  int n = A.size();
  vector<vector<int>> C(n, vector<int>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = A[i][j] - B[i][j];
  return C;
} // end of subtractMatrix

//--------------------------subtractMatrix()-----------------------------------
// perform matrix addition
// C = A + B
//-----------------------------------------------------------------------------
vector<vector<int>> addMatrix(vector<vector<int>> A, vector<vector<int>> B) {
  int n = A.size();
  vector<vector<int>> C(n, vector<int>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = A[i][j] + B[i][j];
  return C;
} // end of addMatrix

//-------------------------strassen()------------------------------------------
// perform multiplication using strassen algorithm
//-----------------------------------------------------------------------------
vector<vector<int>> strassen(vector<vector<int>> A, vector<vector<int>> B) {
  int n = A.size();
  // make sure both are n size
  // return null vector
  if (n != B.size()) {
    cerr << "Invalid matrix" << endl;
    return vector<vector<int>>();
  }
  // matrix size 1 and 2 -- base case (2x2 matrix)
  // compute sing standard method because the overhead of further recursion or
  // the Strassen's method wouldn't be beneficial for such small-sized matrices.
  if (n <= 2) {
    vector<vector<int>> C(n, vector<int>(n));
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
    return C;
  } else {
    // splitting matrix
    int new_size = n / 2;
    vector<vector<int>> a11(new_size, vector<int>(new_size)),
        a12(new_size, vector<int>(new_size)),
        a21(new_size, vector<int>(new_size)),
        a22(new_size, vector<int>(new_size)),
        b11(new_size, vector<int>(new_size)),
        b12(new_size, vector<int>(new_size)),
        b21(new_size, vector<int>(new_size)),
        b22(new_size, vector<int>(new_size));

    // dividing matrices into 4 sub-matrices
    for (int i = 0; i < new_size; i++) {
      for (int j = 0; j < new_size; j++) {
        a11[i][j] = A[i][j];
        a12[i][j] = A[i][j + new_size];
        a21[i][j] = A[i + new_size][j];
        a22[i][j] = A[i + new_size][j + new_size];

        b11[i][j] = B[i][j];
        b12[i][j] = B[i][j + new_size];
        b21[i][j] = B[i + new_size][j];
        b22[i][j] = B[i + new_size][j + new_size];
      }
    }

    // compute 7 products
    vector<vector<int>> p1 = strassen(a11, subtractMatrix(b12, b22)),
                        p2 = strassen(addMatrix(a11, a12), b22),
                        p3 = strassen(addMatrix(a21, a22), b11),
                        p4 = strassen(a22, subtractMatrix(b21, b11)),
                        p5 = strassen(addMatrix(a11, a22), addMatrix(b11, b22)),
                        p6 = strassen(subtractMatrix(a12, a22),
                                      addMatrix(b21, b22)),
                        p7 = strassen(subtractMatrix(a11, a21),
                                      addMatrix(b11, b12));

    // calculate 4 quadrants
    vector<vector<int>> c11 = addMatrix(subtractMatrix(addMatrix(p5, p4), p2),
                                        p6),
                        c12 = addMatrix(p1, p2), c21 = addMatrix(p3, p4),
                        c22 = subtractMatrix(
                            subtractMatrix(addMatrix(p5, p1), p3), p7);

    // calculate the final result
    vector<vector<int>> C(n, vector<int>(n));
    for (int i = 0; i < new_size; i++) {
      for (int j = 0; j < new_size; j++) {
        C[i][j] = c11[i][j];
        C[i][j + new_size] = c12[i][j];
        C[i + new_size][j] = c21[i][j];
        C[i + new_size][j + new_size] = c22[i][j];
      }
    }
    return C;
  }

} // end of strassen

int main() {
  //-----------------------test 1----------------------------------------------
  vector<vector<int>> A1 = {{1, 2}, {3, 4}};
  vector<vector<int>> B1 = {{5, 6}, {7, 8}};

  cout << "Matrix A1: \n";
  for (int i = 0; i < A1.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < A1[i].size(); j++) {
      cout << A1[i][j] << " ";
    }
    cout << "]\n";
  }

  cout << "\nMatrix B1: \n";
  for (int i = 0; i < B1.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < B1[i].size(); j++) {
      cout << B1[i][j] << " ";
    }
    cout << "]\n";
  }

  vector<vector<int>> C1 = strassen(A1, B1);

  cout << "\nMatrix C1 = A1 * B1: \n";
  for (int i = 0; i < C1.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < C1[i].size(); j++) {
      cout << C1[i][j] << " ";
    }
    cout << "]\n";
  }
  //---------------------------------------------------------------------------

  //----------------------test 2-----------------------------------------------
  vector<vector<int>> A2 = {
      {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
  vector<vector<int>> B2 = {
      {17, 18, 19, 20}, {21, 22, 23, 24}, {25, 26, 27, 28}, {29, 30, 31, 32}};

  cout << "Matrix A2: \n";
  for (int i = 0; i < A2.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < A2[i].size(); j++) {
      cout << A2[i][j] << " ";
    }
    cout << "]\n";
  }

  cout << "\nMatrix B2: \n";
  for (int i = 0; i < B2.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < B2[i].size(); j++) {
      cout << B2[i][j] << " ";
    }
    cout << "]\n";
  }

  vector<vector<int>> C2 = strassen(A2, B2);

  cout << "\nMatrix C2 = A2 * B2: \n";
  for (int i = 0; i < C2.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < C2[i].size(); j++) {
      cout << C2[i][j] << " ";
    }
    cout << "]\n";
  }
  //---------------------------------------------------------------------------

  //-------------------------------test 3--------------------------------------
  vector<vector<int>> A3(10, vector<int>(10));
  vector<vector<int>> B3(10, vector<int>(10));
  int count = 1;
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      A3[i][j] = count;
      B3[i][j] = count + 100; // just to create different matrices
      count++;
    }
  }
  vector<vector<int>> C3 = strassen(A3, B3);

  cout << "\nMatrix C3 = A3 * B3: \n";
  for (int i = 0; i < C3.size(); i++) {
    cout << "[ ";
    for (int j = 0; j < C3[i].size(); j++) {
      cout << C3[i][j] << " ";
    }
    cout << "]\n";
  }
  //---------------------------------------------------------------------------
  return 0;
}