#include <iostream>
#include <stdio.h>
using namespace std;



class QRDescomposition 
{ 
    // Access specifier 
    public: 
  
    // Data Members 
    double data[3][3]; //original matrix
    double A[3][3]; //copy of the original matrix
    double Q[3][3]; // Q matrix (QR descomposition)
    double Qt[3][3]; // Q matrix transposed
    double R[3][3]; // R matrix (QR descomposition)
    double QtxA[3][3]; // Qt*A matrix 
    double AxQ[3][3]; // A*Q matrix 

    // Member Functions() 

    //print the diagonal (when QR iteration has sufficient iterations the diagonal of the matrix aproximates to the eigenvalues)
    void print_Diagonal() 
    { 
        for(int j = 0; j<3; j++){
            cout << data[j][j] << ",";
        }
        cout << endl;
    } 

    //this function prints a 3x3 matrix
    void print_Matrix(double matrix[3][3]) 
    { 
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                cout << matrix[i][j] << ",";
            }   
            cout << "\t" <<endl;
        }
        cout << endl;
    } 


    //this function prints the original matrix 
    void print_Data() 
    { 
        cout << "Printing Data "<< endl;
        print_Matrix(data);
    } 

    //this function prints a copy of the original matrix 
    void print_A() 
    { 
        cout << "Printing A "<< endl;
        print_Matrix(A);
    } 

    void print_Q() 
    { 
        cout << "Printing Q "<< endl;
        print_Matrix(Q);
    } 

    void print_R() 
    { 
        cout << "Printing R "<< endl;
        print_Matrix(R);
    } 

    void print_Qt() 
    { 
        cout << "Printing Qt "<< endl;
        print_Matrix(Qt);
    } 

    void print_QtxA() 
    { 
        cout << "Printing QtxA "<< endl;
        print_Matrix(QtxA);
    } 

    void print_AxQ() 
    { 
        cout << "Printing AxQ "<< endl;
        print_Matrix(AxQ);
    }  

    //this function set the data of the object
    void set_Data(double input[3][3]) 
    { 
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                A[i][j] = 0;
                Q[i][j] = 0;
                Qt[i][j] = 0;
                R[i][j] = 0;
                QtxA[i][j] = 0;
                AxQ[i][j] = 0;
            }   
            
        }

       for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                data[i][j] = input[i][j];
                A[i][j] = input[i][j];
            } 
        }

        gram_Schmidt();
    }

    //this function computes the square root of a number
    double sqrt(double num){
        if( num < 0 ) return -1;
        if( num == 0 ) return 0;

        double xhi = num;
        double xlo = 0;
        double guess = num/2; 

        while(guess * guess != num){
            if(guess * guess > num){
                xhi = guess;
            }else{
                xlo = guess;
            }
            double new_guess = (xhi + xlo)/2;
            if(new_guess == guess){
                break;
            }
            guess = new_guess;
        }
        return guess;
    }

    // this function descompose the matrix in the Q and R matrixes
    void gram_Schmidt(){
        int k, i, j;

        for (k=0; k<3; k++){
            R[k][k]=0; // equivalent to sum = 0
            
            for (i=0; i<3; i++)
                R[k][k] = R[k][k] + A[i][k] * A[i][k]; // rkk = sqr(a0k) + sqr(a1k) + sqr(a2k) 
            R[k][k] = sqrt(R[k][k]);  // ||a||
            //cout << endl << "R"<<k<<k<<": " << R[k][k];
            
            for (i=0; i<3; i++) {
                Q[i][k] = A[i][k]/R[k][k];
                //cout << " q"<<i<<k<<": "<<Q[i][k] << " ";
            }
            
            for(j=k+1; j<3; j++) {
                R[k][j]=0;
                for(i=0; i<3; i++) R[k][j] += Q[i][k] * A[i][j];
                //cout << endl << "r"<<k<<j<<": " <<R[k][j] <<endl;
                for (i=0; i<3; i++) A[i][j] = A[i][j] - R[k][j]*Q[i][k];
                //for (i=0; i<3; i++) cout << "a"<<j<<": " << A[i][j]<< " ";
            }
        }

        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                Qt[j][i] = Q[i][j];
            }
        }
    }
}; 

//this function multiplies two matrixes
void mat_Mult(double A[3][3], double B[3][3], double C[3][3]){
    for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                C[i][j] = 0;
            }   
            
        }

    for ( int i = 0; i < 3; i++ )
    {
        for ( int j = 0; j < 3; j++ )
        {
            for ( int k = 0; k < 3; k++ ) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

//QR iteration function
void QR_Algorithm(QRDescomposition& qrDescomp, int count){

    //qrDescomp.print_Data();
    //qrDescomp.print_A();
    //qrDescomp.print_Q();
    //qrDescomp.print_Qt();
    //qrDescomp.print_R();
    //qrDescomp.print_QtxA();
    //qrDescomp.print_AxQ();

    qrDescomp.print_Diagonal();

    mat_Mult(qrDescomp.Qt, qrDescomp.data, qrDescomp.QtxA);
    mat_Mult(qrDescomp.QtxA, qrDescomp.Q, qrDescomp.AxQ);

    if(count > 1){
        QRDescomposition qrDescomp_1;
        qrDescomp_1.set_Data(qrDescomp.AxQ);
        QR_Algorithm(qrDescomp_1, --count);
    }

    
}

int main (){
    //Identity matrix
    double I[3][3] = {  
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };

    //Matrix
    double M[3][3] = {  
            {14, -7, -2},
            {-18, 0, 12},
            {-1, -2, 16}
        };

    
    QRDescomposition qrDes;

    qrDes.set_Data(M);

    //funct(qrDes);
    QR_Algorithm(qrDes, 200);

    return 0;

}