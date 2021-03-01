#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

template <class T> class matrix
{
private:
    T** A;                                    // Elementy
    int n, m;
public:
    matrix(const matrix<T>& other);          // Konstruktor kopiujacy
    matrix(int r, int c);                    // Konstruktor tworzacy macierz nxm
    matrix(int r);                           // Konstruktor tworzacy macierz nxn
    matrix(string filename);                 // Konstruktor tworzacy macierz z pliku po podaniu nazwy utworzy sie plik, do ktorego chcemy wrzucic macierz
    ~matrix();                               // Destruktor
    T getv(int r, int c);                    // Pobiera element na podanej pozycji
    void setv(int r, int c, T v);            // Ustawia element na podanej pozycji
    int rows();                              // Zwraca liczbę wierszy macierzy
    int cols();                              // Zwraca liczbę kolumn macierzy
    void print();                            // Wyswietla macierz
    matrix<T> multiply(const matrix<T>& a);       // Mnozy dwie macierze
    matrix<T> add(matrix<T>& a);            // Dodaje dwie macierze
    matrix<T> subtract(matrix<T>& a);       // Odejmuje dwie macierze
    void store(string filename);             // Zachowuje w pliku o podanej nazwie, jesli takiego nie ma to go tworzy
};
template<class T> matrix<T>::matrix(const matrix<T>& other) 
{
    n = other.n;
    m = other.m;

    A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[m];
    if (other.A)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                A[i][j] = other.A[i][j];
    }
}

template <class T> matrix<T>::matrix(int r)
{
    n = r;
    m = r;
    A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[m];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A[i][j] = 0.0;
}

template <class T> matrix<T>::matrix(int r, int c)
{
    n = r;
    m = c;
    A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[m];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A[i][j] = 0.0;
}

template <class T> matrix<T>::matrix(string filename)
{
    ifstream fin(filename);
   
	    fin >> n >> m;

	    A = new double* [n];
	    for (int i = 0; i < n; i++)
		A[i] = new double[m];

	    for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
		    fin >> A[i][j];

	    fin.close();
    
}

template <class T> matrix<T>::~matrix()
{
    for (int i = 0; i < n; i++)
        delete[] A[i];
    delete[] A;
}

template <class T> T matrix<T>::getv(int r, int c)
{
    if ((r >= 0 && r <= n) && (c >= 0 && c <= m))
    {
         return A[r][c];
    }
    else
    {
        cout << "Can't get number of matrix, wrong number of rows or columns size\n" << endl;
        exit(1);
    }
}

template <class T> void matrix<T>::setv(int r, int c, T v)
{
    if ((r >= 0 && r <= n) && (c >= 0 && c <= m))
    {
        A[r][c] = v;
    }
    else
    {
        cout << "Can't set number of matrix, wrong number of rows or columns size\n" << endl;
    }
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


    cout << "Matrix:" << n << "x" << m << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << setw(5) << A[i][j];
        }
        cout << endl;
    }
}

template <class T> matrix<T> matrix<T>::multiply(const matrix<T>& a)
{
    T sum = 0;
    int r = n;
    int c = a.m;
    matrix result(r, c);
    if (m != a.n)
    {
        cout << "\nCan't multiply two matrix if number of columns in first matrix not equals to rows of second matrix, change the size\n" << endl;
        exit(1);
    }
    else
    {
        int i, j, k;

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < a.m; j++)
            {
                sum = 0;
                for (k = 0; k < m; k++)

                    sum += A[i][k] * (a.A[k][j]);

                result.A[i][j] = sum;
            }
        }
    }
    return result;
}

template <class T> matrix<T> matrix<T>::add(matrix<T>& a)
{
    int r = n;
    int c = m;
    matrix result(r, c);
    if ((n != a.n) || (m != a.m))
    {
        cout << "\nCan't add two matrix in another size, change the size\n" << endl;
        exit(1);
    }
    else
    {
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                result.A[i][j] = (*this).A[i][j] + a.A[i][j];
    }
    return result;
}

template <class T> matrix<T> matrix<T>::subtract(matrix<T>& a)
{
    int r = n;
    int c = m;
    matrix result(r, c);
    if ((n != a.n) || (m != a.m))
    {
        cout << "\nCan't add two matrix in another size, change the size\n" << endl;
        exit(1);
    }
    else
    {
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                result.A[i][j] = (*this).A[i][j] - a.A[i][j];
    }
    return result;
}

template <class T> void matrix<T>::store(string filename)
{
    ofstream fin(filename);  //tworzy plik o podanej nazwie, jezeli nie istnieje

   // fin.open(filename, ifstream::out); // Otwieramy strumień do odczytu, nie potrzebujemy jezeli chcemy tworzyc plik
    fin << n << " " << m << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            fin << setw(5) << A[i][j];
        fin << endl;
    }

    fin.close();   
}


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
    cout << "Macierz A: \n" << endl;
    cout << "Metoda rows(); dla macierzy A\t" << A.rows() << endl << endl;
    cout << "Metoda cols(); dla macierzy A\t" << A.cols() << endl << endl;
    A.print(); 
    cout << "Macierz B: \n" << endl;
    B.print();
    cout << "Ustawienie elementu macierzy A: r=1, c=1, v=30 \n" << endl << endl;
    A.setv(1, 1, 30.);
    A.print();
    //cout << "Metoda get(); dla macierzy A\t" << endl << endl;
    //A.getv(1, 100000000); //zly rozmiar
    //cout << "Wynik dodawania macierzy A i B: \n" << endl; //zly rozmiar macierzy, wyswietla error i exit()
    //matrix<double> adding(A.add(B));
    //adding.print();
    //cout << "Wynik odejmowania macierzy A i B: \n" << endl; //zly rozmiar macierzy, wyswietla error i exit()
    //matrix<double> subtracting(A.subtract(B));
    //subtracting.print();
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierz A:\n" << endl;
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
            A.setv(i, j, 30.);
    }
    A.print();
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierz B:\n" << endl;
    for (int i = 0; i < B.rows(); i++)
    {
        for (int j = 0; j < B.cols(); j++)
            B.setv(i, j, 2.);
    }
    B.print();
   
    cout << "Wynik mnozenia macierzy A i B\n" << endl;
    matrix<double> multiplying(A.multiply(B)); //przekazanie macierzy D, ktora jest zwracana podczas operacji do konstruktora kopiujacego
    multiplying.print();
    cout << endl << endl;
    A.store("Macierz_A.txt"); //tworzy plik i do niego zapisuje
    B.store("Macierz_B.txt");

    cout << "-------------------------------------------------------------------------------------" << endl;
    cout << "test 2. zainicjowanie dwoch macierzy kwadratowych n x m gdzie n=m=8 \n";
    cout << "-------------------------------------------------------------------------------------" << endl;
    matrix<double> F(8);
    matrix<double> G(8);
   
    cout << "Macierz F: \n" << endl;
    cout << "Metoda rows(); dla macierzy F\t" << F.rows() << endl << endl;
    cout << "Metoda cols(); dla macierzy F\t" << F.cols() << endl << endl;
    F.print();
    cout << "Macierz G: \n" << endl;
    G.print();
    cout << "Ustawienie elementu macierzy F: r=1, c=1, v=30 \n" << endl << endl;
    F.setv(1, 1, 30.);
    F.print();
    cout << "Metoda get(); dla macierzy F\t" << F.getv(1, 1) << endl << endl;
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierzy F:\n" << endl;
    for (int i = 0; i < F.rows(); i++)
    {
        for (int j = 0; j < F.cols(); j++)
            F.setv(i, j, 5.);
    }
    F.print();
    cout << "Zmiana elementow macierzy, dla lepszego zobrazowania dzialania, nowa macierz G:\n" << endl;
    for (int i = 0; i < G.rows(); i++)
    {
        for (int j = 0; j < G.cols(); j++)
            G.setv(i, j, 7.);
    }
    G.print();
    cout << "Wynik mnozenia macierzy F i G: \n" << endl;
    matrix<double> fileM(F.multiply(G));
    fileM.print();
    cout << endl << endl;
    cout << "Wynik dodawania macierzy F i G: \n" << endl;
    matrix<double> fileA(F.add(G));
    fileA.print();
    cout << "Wynik odejmowania macierzy F i G: \n" << endl;
    matrix<double> fileS(F.subtract(G));
    fileS.print();
    F.store("Macierz_F.txt");
    G.store("Macierz_G.txt");
    cout << "-------------------------------------------------------------------------------------" << endl;
    cout << "test 3. zainicjowanie dwoch macierzy n x m gdzie n!=m z wybranych plikow \n";
    cout << "-------------------------------------------------------------------------------------" << endl;
    matrix<double> K("dane.txt");
    matrix<double> L("dane2.txt");
    matrix<double> M("dane3.txt"); 
    cout << "Macierz z pliku dane.txt : \n" << endl;
    K.print();
    cout << "Macierz z pliku dane2.txt : \n" << endl;
    L.print();
    cout << "Macierz z pliku dane3.txt (plik pusty): \n" << endl;
    M.print();
    cout << "Metoda rows(); dla macierzy dane.txt (K)\t" << K.rows() << endl << endl;
    cout << "Metoda cols(); dla macierzy dane.txt (K)\t" << K.cols() << endl << endl;
    cout << "Wynik mnozenia macierzy: \n" << endl;
    matrix<double> multiplyingM(K.multiply(L)); //zly rozmiar - komunikat 
    multiplyingM.print();
    cout << endl << endl;
    cout << "Wynik dodawania macierzy: \n" << endl;
    matrix<double> addingM(K.add(L));
    addingM.print();
    cout << "Wynik odejmowania macierzy: \n" << endl;
    matrix<double> subtractingM(K.subtract(L));
    subtractingM.print();

    K.store("Macierz_K.txt");
    L.store("Macierz_L.txt");
   
    return 0;
}

