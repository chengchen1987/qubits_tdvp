#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>
#include <sstream>

#include "ed_basis.h"

using namespace std;

// bitwise operations
int numOfBit1(const l_int& b)
{
    int a = b;
    int cnt = 0;
    while (a != 0)
    {
        ++cnt;
        a &= (a - 1);
    }
    return cnt;
}

// find first n '1's of a 
void findBit1(const l_int& a, const int& n, int* b)
{
    int x = 0;
    int i = 0;
    while (x < n) {
        if (((a >> i) & 1) == 1) {
            b[x] = i;
            x++;
        }
        i++;
    }
}

void print_binary(const l_int& a, const int& n)
{
    cout << "  ";
    for (int ix = 0; ix < n; ix++) cout << ((a >> (n - ix - 1)) & 1);
    cout << "  ";
}

l_int Power(int& m, int& n)
{
    l_int ans = 1;
    for (int i = 0; i < n; i++) {
        ans = ans * (long long)m;
    }
    return ans;
}

l_int Factorial(const int& m)
{
    l_int ans = 1;
    for (int i = 1; i < m + 1; i++) {
        ans = ans * i;
    }
    return ans;
}

l_int Combination(const int& m, const int& n)
{
    l_int  ans = 1;
    int minmn = min(n, m - n);
    int maxmn = max(n, m - n);
    for (int i = 0; i < minmn; i++) {
        ans = ans * (maxmn + 1 + i);
    }
    ans = ans / Factorial(minmn);
    return ans;
}

// basis 
Basis::Basis(int _L, int _nop) :
    L(_L),
    nop(_nop),
    Dim(Combination(L,nop))
{
    State = new l_int[Dim];
    
    l_int j = 0;
    int totalcharge;
    for (l_int s = 0; s < (1<<L); s++) {
        totalcharge = numOfBit1(s);
        if (totalcharge == nop) {
            State[j] = s;
            j++;
        }
    }
    // debug
    //  // PrintStates();
}

Basis::~Basis() 
{
    delete[] State;
}

l_int Basis::get_state(const l_int& s)
{
    return State[s];
}

l_int Basis::get_index(const l_int& s) 
{
    if (s < 0) return -1;
    l_int bmin = 0, bmax = Dim;
    l_int b;
    while (bmin <= bmax) {
        b = bmin + (bmax - bmin) / 2;
        l_int aux = State[b];
        if (s == aux) {
            return b;
        }
        else if (s > aux) {
            bmin = b + 1;
        }
        else if (s < aux) {
            bmax = b - 1;
        }
    }
    return -1;
}

double Basis::Cal_chargesign(const l_int& numk, const int& i) 
{
    l_int numk1 = (numk >> (i + 1));
    int cnt = 0;
    while (numk1 != 0)
    {
        ++cnt;
        numk1 &= (numk1 - 1);
    }
    double sign = cnt % 2 == 0 ? 1.0 : -1.0;
    return sign;
}

double Basis::Cal_chargesign(const l_int& numk, const int& i, const int& j) 
{
    double sign;
    if (j < i) {
        cout << "Error! In function Cal_chargesign(const int &numk, const int &i, const int &j) there must be i <= j";
        exit(-100);
    }
    else if (j < i + 1) {
        sign = 1;
    }
    else {
        int cnt = 0;
        for (int site = i + 1; site < j; site++) {
            cnt += (numk >> site) & 1;
        }
        sign = cnt % 2 == 0 ? 1.0 : -1.0;
    }
    return sign;
}

void Basis::PrintStates() 
{
    ofstream ofd("states");
    for (l_int i = 0; i < Dim; i++) {
        ofd << setw(4) << i << setw(8) << State[i];
        ofd << "  ";
        for (int j = 0; j < L; j++) {
            ofd << ((State[i] >> (L - 1 - j)) & 1);
        }
        ofd << endl;
    }
    ofd.close();
}
