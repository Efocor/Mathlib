// ferimathos.h
#ifndef FERIMATHOS_H
#define FERIMATHOS_H

namespace BibliotecaMatematica {

    // Función para calcular el máximo común divisor (MCD)
    int mcd(int a, int b);

    // Función para calcular el mínimo común múltiplo (MCM)
    int mcm(int a, int b);

    // Función para calcular el máximo común divisor extendido (Algoritmo de Euclides extendido)
    int euclidesExtendido(int a, int b, int &x, int &y);

    // Función para verificar si un número es primo
    bool esPrimo(int n);

    // Función para verificar si un número es primo
    bool esPrimoMer(int n);

    // Función para calcular el factorial de un número entero no negativo (implementación iterativa)
    unsigned long long factorialIterativo(int n);

    // Función para calcular el factorial de un número entero no negativo (implementación recursiva)
    int factorial(int n);

    // Función para calcular la potencia de un número entero
    int potencia(int base, int exponente);

    // Función para calcular la potencia modular de un número entero
    int potenciaModular(int base, int exponente, int modulo);

    // Función para calcular el inverso modular de un número entero
    int inversoModular(int a, int modulo);

    // Función para verificar si un número es un número perfecto
    bool esNumeroPerfecto(int n);

    // Función para calcular el número de combinaciones (número de formas de elegir k elementos de un conjunto de n elementos)
    unsigned long long combinaciones(int n, int k);

    // Función para calcular el número de permutaciones (número de formas de organizar k elementos de un conjunto de n elementos)
    unsigned long long permutaciones(int n, int k);

    // Función para calcular el residuo de una congruencia lineal (residuo de ax ≡ b (mod m))
    int residuoCongruenciaLineal(int a, int b, int m);

    // Función para verificar si un número es un número de Fibonacci
    bool esNumeroFibonacci(int n);

    // Función para calcular el número de Euler (phi de Euler)
    int eulerPhi(int n);

    // Función para calcular el número de Euler (phi de Euler)
    int eulerPhi2(int n);

    // Función para verificar si un número es un número de Lucas
    bool esNumeroLucas(int n);

    // Función para calcular el coeficiente binomial (n sobre k)
    unsigned long long coeficienteBinomial(int n, int k);

    // Verifica si un número es un número de Harshad (número que es divisible por la suma de sus dígitos)
    bool esNumeroHarshad(int n);

    // Calcula el número de Stirling de segunda clase (S(n, k))
    int numeroStirlingSegundaClase(int n, int k);

    // Calcula el número de Bell (número de particiones de un conjunto de n elementos)
    int numeroBell(int n);

    // Calcula la función phi de Möbius (μ(n))
    int funcionPhiMobius(int n);

    // Exponencial de un número
    double exponencial(double x);

    // Función auxiliar para calcular el logaritmo natural (ln) utilizando la serie de Taylor
    double logaritmoNatural(double x);

    // Función auxiliar para calcular el producto de dos números de punto flotante sin errores de redondeo
    double productoSinError(double a, double b);

    // Función auxiliar para calcular el valor absoluto de un número de punto flotante
    double valorAbsoluto(double num);

    // Función auxiliar para calcular la fracción decimal de un número de punto flotante
    double fraccionDecimal(double num);

    // Función auxiliar para calcular la parte entera de un número de punto flotante
    double parteEntera(double num);

    // Función para calcular la potencia de un número (base) elevado a un exponente (exponente) para enteros
    int potenciaEntera(int base, int exponente);

    // Función para calcular la potencia de un número (base) elevado a un exponente (exponente)
    double potencia(double base, double exponente);

    // Función para calcular el polinomio de Chebyshev de segunda clase U_n(x)
    int polinomioChebyshevSegundaClase(int n, int x);

    // Función para calcular el polinomio de Chebyshev de primera clase T_n(x)
    int polinomioChebyshevPrimeraClase(int n, int x);

    // Números de Fibonacci usando iteración
    int fibonacci(int n);

    // Calcula el número de Carmichael (número compuesto que cumple a^(n-1) ≡ 1 (mod n) para todo a coprimo con n)
    bool esNumeroCarmichael(int n);

    // Verifica si un número es poderoso
    bool esNumeroPoderoso(int n);

    // Función para calcular la potencia de un número (base) elevado a un exponente (exponente)
    double potencia(double base, double exponente);

    // Función para calcular la potencia de un número (base) elevado a un exponente (exponente) para enteros
    int potenciaEntera(int base, int exponente);

    // Función auxiliar para calcular la parte entera de un número de punto flotante
    double parteEntera(double num);

    // Función auxiliar para calcular la fracción decimal de un número de punto flotante
    double fraccionDecimal(double num);

    // Función auxiliar para calcular el producto de dos números de punto flotante sin errores de redondeo
    double productoSinError(double a, double b);

    // Función auxiliar para calcular el logaritmo natural (ln) utilizando la serie de Taylor
    double logaritmoNatural(double x);

    // Función auxiliar para calcular el exponencial (exp) utilizando la serie de Taylor
    double exponencial(double x);

    // Función de Ackermann
    int ackermann(int m, int n);

    // Función para calcular el módulo de un número entero (a % b) de forma segura.
    int modulo(int a, int b);

    // Función para calcular el inverso multiplicativo módulo 'm'
    int inversoMultiplicativo(int a, int modulo);

    // Función para calcular la exponenciación rápida
    int potenciaRapida(int base, int exponente);

    // Test de Miller-Rabin para verificar si un número es probablemente primo
    bool esProbablementePrimo(int n, int iteraciones);

    // Generador lineal congruente (LCG) para generar números aleatorios
    int lcg(int seed, int a, int c, int m);

    // Calcula el número de Smith (número cuya suma de dígitos es igual a la suma de dígitos de sus factores primos)
    bool esNumeroSmith(int n);
    
    // Algoritmo extendido de Euclides
    int mcdExtendido(int a, int b, int& x, int& y);

    // Función para calcular el número de Wagstaff W(n) = 3^(n) - 1
    int numeroWagstaff(int n);

    // Función para calcular el número de Mersenne M(p) = 2^p - 1
    int numeroMersenne(int p);

    // Función para verificar si un número es un número de Mersenne (M(p) = 2^p - 1)
    bool esNumeroMersenne(int p);

    // Teorema chino del resto, encuentra el número x tal que: x ≡ a[i] (mod m[i])
    int teoremaChinoResto(const int residuos[], const int modulos[], int n);

    // Devuelve un puntero a un arreglo dinámico de enteros que contiene todos los números primos.
    int* cribaEratostenes(int n, int& numPrimos);

    // Función phi de Euler: Calcula la cantidad de enteros positivos menores o iguales a n que son coprimos con n.
    int phiEuler(int n);

    // Calcula el número de Bernoulli Bn utilizando la fórmula recursiva de Bernoulli
    double numeroBernoulli(int n);

    // Calcula el número de Leyland (número de la forma x^y + y^x)
    bool esNumeroLeyland(int n);

    // Función para verificar si un número es una potencia de dos
    bool esPotenciaDeDos(int n);

    // Función para calcular el número de Lucas-Lehmer L(n) = 2^(n) - 1
    int numeroLucasLehmer(int n);

    // Función para calcular el número de Fermat F(n) = 2^(2^n) + 1
    int numeroFermat(int n);

    // Función para verificar si un número es una raíz primitiva módulo n.
    bool esRaizPrimitiva(int g, int n);

    // Función para calcular la función φ de Euler para un número n.
    int eulerPhi(int n);

    // Función para calcular el símbolo de Jacobi.
    int simboloJacobi(int a, int n);

    // Función Zeta de Riemann (ζ(s))
    double funcionZetaRiemann(double s, int iteraciones = 1000);

    // Función auxiliar para calcular el factorial de un número de punto flotante
    double factorial(double n);

    // Función para calcular la función gamma de Euler (Γ) utilizando la fórmula de Euler
    double funcionGammaEuler(double x);

    // Función para calcular la función gamma de Stirling (Γ) utilizando aproximaciones de π y e
    double funcionGammaStirling(double x, int iteracionesPi = 1000000, int iteracionesE = 100);

    // Función para calcular la función gamma de Lanczos (Γ)
    double funcionGammaLanczos(double x, int iteracionesPi = 1000000);

    // Aproximación de π utilizando la fórmula de Leibniz para π/4
    double aproximacionPi(int iteraciones = 1000000);

    // Aproximación de e utilizando la fórmula de Euler para e
    double aproximacionE(int iteraciones = 100);

    // Función para calcular el valor de 'sin(x)' utilizando la expansión en serie de Taylor
    double sin(double x, int iteraciones = 10);

    // Función para calcular el valor de 'cos(x)' utilizando la expansión en serie de Taylor
    double cos(double x, int iteraciones = 10);

    // Función para calcular el polinomio de Chebyshev de segunda clase U_n(x)
    int polinomioChebyshevSegundaClase(int n, int x);

    // Función para verificar el teorema de Fermat para exponentes pequeños (n <= 20)
    bool verificarTeoremaFermat(int n);

    // Función para calcular la función gamma de Euler (Γ) utilizando la fórmula de Euler
    double funcionGammaEuler(double x);

    // Función para calcular la función gamma de Legendre (Γ) utilizando la fórmula de Legendre
    double funcionGammaLegendre(double x);

    // Función para calcular la función eta de Riemann (η) utilizando la serie de Dirichlet
    Complejo funcionEtaRiemann(const Complejo& s, int iteraciones = 1000);

    // Función para calcular la función zeta de Riemann (ζ) utilizando la serie de Dirichlet
    Complejo funcionZetaRiemann(const Complejo& s, int iteraciones = 1000);

} // namespace BibliotecaMatematica

#endif // FERIMATHOS_H