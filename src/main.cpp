#include <iostream>
#include <fstream>
#include <iomanip>
#define _matrix_class
using namespace std;

template <class T> class matrix
{
private:
    int s;  // Rozmiar macierzy
    T* A;      // Elementy
public:
    int n, m;
    matrix(int r, int c);                    // Konstruktor tworzacy macierz nxm
    matrix(int r);                           // Konstruktor tworzacy macierz nxn
    matrix(string filename);                 // Konstruktor tworzacy macierz z pliku po podaniu nazwy utworzy sie plik, do ktorego chcemy wrzucic macierz
    ~matrix();                               // Destruktor
    T getv(int r, int c);                    // Pobiera element na podanej pozycji
    void setv(int r, int c, T v);            // Ustawia element na podanej pozycji
    int rows();                              // Zwraca liczbę wierszy macierzy
    int cols();                              // Zwraca liczbę kolumn macierzy
    void print();                            // Wyswietla macierz
    void multiply(matrix<T>& a, matrix<T>& b);           //Mnozy dwie macierze
    void add(matrix<T>& a, matrix<T>& b);               //Dodaje dwie macierze
    void subtract(matrix<T>& a, matrix<T>& b);          //Odejmuje dwie macierze
    void store(string filename);                       //Zachowuje w pliku o podanej nazwie, jesli takiego nie ma to go tworzy
};
template <class T> matrix<T>::matrix(int r)
{
    n = r;
    m = r;
    s = n * m;
    A = new T[s];
    for (int i = 0; i < s; i++)
    {
        A[i] = 0;
    }
}

template <class T> matrix<T>::matrix(int r, int c)
{
    n = r;
    m = c;
    s = n * m;
    A = new T[s];
    for (int i = 0; i < s; i++)
    {
        A[i] = 0;
    }
}

template <class T> matrix<T>::matrix(string filename)
{
    ifstream fin(filename);  // tworzenie pliku o podanej nazwie

    //fin.open("dane.txt", ifstream::in);
    fin >> n >> m;

    s = n * m;

    A = new T[s];

    for (int i = 0; i < s; i++) fin >> A[i];

    fin.close();  //wydaje mi sie, ze nie trzeba close, ale nie jestem pewna
}

template <class T> matrix<T>::~matrix()
{
    delete[] A;
}

template <class T> T matrix<T>::getv(int r, int c)
{
    return A[r * m + c];
}

template <class T> void matrix<T>::setv(int r, int c, T v)
{
    A[r * m + c] = v;
}

template <class T> int matrix<T>::rows()
{
    return n;
}


template <class T> int matrix<T>::cols()
{
    return m;
}

template <class T> void matrix<T>::print()
{
    cout << "Macierz:" << n << "x" << m << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            cout << setw(5) << A[i * m + j];
        cout << endl;
    }
}

template <class T> void matrix<T>::multiply(matrix<T>& a, matrix<T>& b)
{
    if (a.m == b.n)
    {
        T sum;
        int i, j, k;

        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
            {
                sum = 0;
                for (k = 0; k < a.m; k++) sum += a.getv(i, k) * b.getv(k, j);
                A[i * m + j] = sum;
            }
    }
    else
    {
        cout << "Can't multiply two matrix if number of columns in first matrix not equals to rows of second matrix, change the size" << endl;
        cout << "The matrix you see on screen is first matrix initiated by constructor" << endl;
    }
}

template <class T> void matrix<T>::add(matrix<T>& a, matrix<T>& b)
{
    T sum;
    int i, j;
    if ((a.n == b.n) && (a.m == b.m))
    {
        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
            {
                sum = 0;
                sum += a.getv(i, j) + b.getv(i, j);
                A[i * m + j] = sum;
            }
    }
    else
    {
        cout << "Can't add two matrix in another size, change the size" << endl;
        cout << "The matrix you see on screen is first matrix initiated by constructor" << endl;
    }
}

template <class T> void matrix<T>::subtract(matrix<T>& a, matrix<T>& b)
{
    T sum;
    int i, j;
    if ((a.n == b.n) && (a.m == b.m))
    {
        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
            {
                sum = 0;
                sum += a.getv(i, j) - b.getv(i, j);
                A[i * m + j] = sum;
            }
    }
    else
    {
        cout << "Can't subtract two matrix in another size, change the size" << endl;
        cout << "The matrix you see on screen is first matrix initiated by constructor" << endl;
    }
}
template <class T> void matrix<T>::store(string filename)
{
    ofstream fin(filename);  //tworzy plik o podanej nazwie, jezeli nie istnieje

   // fin.open(filename, ifstream::out); // Otwieramy strumień do odczytu, nie potrzebujemy jezeli chcemy tworzyc plik
    fin << "Macierz:" << n << "x" << m << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            fin << setw(5) << A[i * m + j];
        fin << endl;
    }

    fin.close();   //teoretycznie niepotrzebne fin.close(), ale nie jestem pewna
}
// Program główny
//---------------

int main()
{
    cout << "----------------------------------------------------------------" << endl;
    cout << "          Gabriela Wielgus, Projekt 1 - JiPP, Klasy, Macierze          " << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << "-------------------------------------------------------------------------------------" << endl;
    cout << "test 1. zainicjowanie dwoch macierzy zerami o podanych wymiarach n x m gdzie n!=m \n podalam A.m = B.n i A.n = B.m dla celow wykonania mnozenia, dodawanie i odejmowanie w tym przypadku jest bezcelowe,\n nie zgodne z zalozeniami\n";
    cout << "-------------------------------------------------------------------------------------" << endl;
    matrix<double> A(8,4);
    matrix<double> B(4, 8);// A.m == B.n, dla mozliwosc mnozenia tych dwoch macierzy
    matrix<double> C(A.n, B.m); //wynik mnozenia macierzy A i B
    matrix<double> D(A.n, A.m); //wynik dodawania macierzy a i b
    matrix<double> E(A.n, A.m); // wynik odejmowania macierzy a i b
    cout << "Macierz A: " << endl;
    cout << "Metoda rows(); dla macierzy A\t" << A.rows() << endl << endl;
    cout << "Metoda cols(); dla macierzy A\t" << A.cols() << endl << endl;
    A.print(); 
    cout << "Macierz B: " << endl;
    B.print();
    cout << "Ustawienie elementu macierzy A: r=1, c=1, v=30 " << endl << endl;
    A.setv(1, 1, 30.);
    A.print();
    cout << "Metoda get(); dla macierzy A\t" << A.getv(1, 1) << endl << endl;;
    cout << "Wynik dodawania macierzy: " << endl;
    D.add(A, B);
    D.print();
    cout << "Wynik odejmowania macierzy: " << endl;
    E.subtract(A, B);
    E.print();
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierz A:" << endl;
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
            A.setv(i, j, 30.);
    }
    A.print();
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierz B:" << endl;
    for (int i = 0; i < B.rows(); i++)
    {
        for (int j = 0; j < B.cols(); j++)
            B.setv(i, j, 2.);
    }
    B.print();
    C.multiply(A, B);
    cout << "Wynik mnozenia macierzy A i B" << endl;
    C.print();
    cout << endl << endl;
    A.store("Macierz_A.txt"); //tworzy plik i do niego zapisuje
    B.store("Macierz_B.txt");
    cout << "-------------------------------------------------------------------------------------" << endl;
    cout << "test 2. zainicjowanie dwoch macierzy kwadratowych n x m gdzie n=m=8 \n";
    cout << "-------------------------------------------------------------------------------------" << endl;
    matrix<double> F(8);
    matrix<double> G(8);
    matrix<double> H(F.n, G.m); 
    matrix<double> I(F.n, G.m); 
    matrix<double> J(F.n, G.m); 
    cout << "Macierz F: " << endl;
    cout << "Metoda rows(); dla macierzy A\t" << F.rows() << endl << endl;
    cout << "Metoda cols(); dla macierzy A\t" << F.cols() << endl << endl;
    F.print();
    cout << "Macierz G: " << endl;
    G.print();
    cout << "Ustawienie elementu macierzy F: r=1, c=1, v=30 " << endl << endl;
    F.setv(1, 1, 30.);
    F.print();
    cout << "Metoda get(); dla macierzy F\t" << F.getv(1, 1) << endl << endl;
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierzy F:" << endl;
    for (int i = 0; i < F.rows(); i++)
    {
        for (int j = 0; j < F.cols(); j++)
            F.setv(i, j, 5.);
    }
    F.print();
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierz G:" << endl;
    for (int i = 0; i < G.rows(); i++)
    {
        for (int j = 0; j < G.cols(); j++)
            G.setv(i, j, 7.);
    }
    G.print();
    cout << "Wynik mnozenia macierzy: " << endl;
    H.multiply(F, G);
    H.print();
    cout << endl << endl;
    cout << "Wynik dodawania macierzy: " << endl;
    I.add(F, G);
    I.print();
    cout << "Wynik odejmowania macierzy: " << endl;
    J.subtract(F, G);
    J.print();

    F.store("Macierz_F.txt");
    G.store("Macierz_G.txt");
    cout << "-------------------------------------------------------------------------------------" << endl;
    cout << "test 3. zainicjowanie dwoch macierzy n x m gdzie n!=m z wybranych plikow \n";
    cout << "-------------------------------------------------------------------------------------" << endl;
    matrix<double> K("dane.txt");
    matrix<double> L("dane2.txt");
    matrix<double> M("dane3.txt"); 
    matrix<double> N(K.n, L.m); 
    matrix<double> O(K.n, L.m); 
    matrix<double> P(K.n, L.m); 
    cout << "Macierz z pliku dane.txt : " << endl;
    K.print();
    cout << "Macierz z pliku dane2.txt : " << endl;
    L.print();
    cout << "Macierz z pliku dane3.txt (plik pusty): " << endl;
    M.print();
    cout << "Metoda rows(); dla macierzy dane.txt (K)\t" << K.rows() << endl << endl;
    cout << "Metoda cols(); dla macierzy dane.txt (K)\t" << K.cols() << endl << endl;
    cout << "Wynik mnozenia macierzy: " << endl;
    N.multiply(K, L);
    N.print();
    cout << endl << endl;
    cout << "Wynik dodawania macierzy: " << endl;
    O.add(K, L);
    O.print();
    cout << "Wynik odejmowania macierzy: " << endl;
    P.subtract(K, L);
    P.print();

    N.store("Wynik_mnozenia_macierze_plik.txt");
    O.store("Wynik_dodawania_macierze_plik.txt");
    P.store("Wynik_odejmowania_macierze_plik.txt");

  
    return 0;
}
