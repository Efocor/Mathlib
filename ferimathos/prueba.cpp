// main.cpp
#include <iostream>
#include "ferimathos.h"

int main() {
    int a = 12;
    int b = 18;

    std::cout << "MCD(" << a << ", " << b << ") = " << BibliotecaMatematica::mcd(a, b) << std::endl;
    std::cout << "MCM(" << a << ", " << b << ") = " << BibliotecaMatematica::mcm(a, b) << std::endl;

    int n = 17;
    std::cout << n << (BibliotecaMatematica::esPrimo(n) ? " es" : " no es") << " un número primo." << std::endl;

    int valorFactorial = 5;
    std::cout << "Factorial de " << valorFactorial << " es " << BibliotecaMatematica::factorial(valorFactorial) << std::endl;

    int base = 2;
    int exponente = 10;
    int modulo = 1000000007;
    std::cout << base << "^" << exponente << " mod " << modulo << " = " << BibliotecaMatematica::potenciaModular(base, exponente, modulo) << std::endl;

    int numero = 7;
    int mod = 11;
    std::cout << "Inverso de " << numero << " mod " << mod << " = " << BibliotecaMatematica::inversoModular(numero, mod) << std::endl;

    // ¡Agrega más pruebas para otras funciones aquí!

    return 0;
}

//g++ ferimathos.cpp prueba.cpp -o programa
//programa.exe
