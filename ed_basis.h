#ifndef BASIS_H
#define BASIS_H

#define l_int long long
#define PI 3.14159265358979323846

// random number
double ran_num();
// bit operation, default l_int, 'int' for length of int 
int numOfBit1(const l_int& a);
void findBit1(const l_int& a, const int& n, int* b);
void print_binary(const l_int& a, const int& n);

l_int Power(const int& m, const int& n);
l_int Factorial(const int& m);
l_int Combination(const int& m, const int& n);

class Basis
{
private:
    int L;
    int nop;
    l_int Dim;
    l_int* State;

public:
    Basis(int _L, int _nop);
    ~Basis();

    int get_L(void) { return L; }
    int get_nop(void) { return nop; }
    l_int get_Dim(void) { return Dim; }

    l_int get_state(const l_int& index);
    l_int get_index(const l_int& s);

    double Cal_chargesign(const l_int& s, const int& i);
    double Cal_chargesign(const l_int& s, const int& i, const int& j);
    // debug
    void PrintStates();
};
#endif
