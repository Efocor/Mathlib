// ferimathos.cpp
// #include "ferimathos.h"
#include "ferimathos.h"

namespace BibliotecaMatematica {

    // Implementación de las funciones

    // Función para calcular π utilizando la serie de Leibniz
    double calcularPi(int iteraciones = 100000) {
        double pi = 0.0;
        for (int i = 0; i < iteraciones; ++i) {
            int signo = (i % 2 == 0) ? 1 : -1;
            double termino = static_cast<double>(signo) / (2 * i + 1);
            pi += termino;
        }
        pi *= 4.0;
        return pi;
    }

    // Función para calcular el valor absoluto de un número
    double abs(double x) {
        return (x < 0) ? -x : x;
    }

    // Función para calcular la raíz cuadrada de un número utilizando el método de aproximaciones sucesivas (método de Herón)
    double sqrt(double x, int iteraciones = 100) {
        if (x < 0) {
            // La raíz cuadrada no está definida para números negativos
            return 0.0;
        }

        double aproximacion = x; // Tomamos la aproximación inicial como el propio número
        for (int i = 0; i < iteraciones; ++i) {
            aproximacion = 0.5 * (aproximacion + x / aproximacion);
        }

        return aproximacion;
    }

    // Función para calcular el valor absoluto de un número entero
    int abs(int x) {
        return (x < 0) ? -x : x;
    }

    int mcd(int a, int b) {
        if (b == 0) {
            return a;
        }
        return mcd(b, a % b);
    }

    int mcm(int a, int b) {
        return a / mcd(a, b) * b;
    }

    // Función para calcular el máximo común divisor extendido (Algoritmo de Euclides extendido)
    int euclidesExtendido(int a, int b, int &x, int &y) {
        if (a == 0) {
            x = 0;
            y = 1;
            return b;
        }

        int x1, y1;
        int gcd = euclidesExtendido(b % a, a, x1, y1);

        x = y1 - (b / a) * x1;
        y = x1;

        return gcd;
    }


    bool esPrimo(int n) {
        if (n <= 1) {
            return false;
        }
        for (int i = 2; i * i <= n; ++i) {
            if (n % i == 0) {
                return false;
            }
        }
        return true;
    }

    unsigned long long factorialIterativo(int n) {
        unsigned long long resultado = 1;
        for (int i = 1; i <= n; ++i) {
            resultado *= i;
        }
        return resultado;
    }

    int factorial(int n) {
        if (n == 0 || n == 1)
            return 1;
        else
            return n * factorial(n - 1);
    }

    int potencia(int base, int exponente) {
        if (exponente == 0) {
            return 1;
        }
        int resultado = 1;
        while (exponente > 0) {
            if (exponente % 2 == 1) {
                resultado *= base;
            }
            base *= base;
            exponente /= 2;
        }
        return resultado;
    }

    int potenciaModular(int base, int exponente, int modulo) {
        int resultado = 1;
        base %= modulo;
        while (exponente > 0) {
            if (exponente % 2 == 1) {
                resultado = (resultado * base) % modulo;
            }
            base = (base * base) % modulo;
            exponente /= 2;
        }
        return resultado;
    }

    int inversoModular(int a, int modulo) {
        int m0 = modulo;
        int y = 0, x = 1;

        if (modulo == 1) {
            return 0;
        }

        while (a > 1) {
            int q = a / modulo;
            int t = modulo;
            modulo = a % modulo, a = t;
            t = y;
            y = x - q * y;
            x = t;
        }

        if (x < 0) {
            x += m0;
        }

        return x;
    }

    bool esNumeroPerfecto(int n) {
        int sumaDivisores = 0;
        for (int i = 1; i <= n / 2; ++i) {
            if (n % i == 0) {
                sumaDivisores += i;
            }
        }
        return sumaDivisores == n;
    }

    unsigned long long combinaciones(int n, int k) {
        if (k > n - k) {
            k = n - k;
        }
        unsigned long long resultado = 1;
        for (int i = 0; i < k; ++i) {
            resultado *= n - i;
            resultado /= i + 1;
        }
        return resultado;
    }

    unsigned long long permutaciones(int n, int k) {
        unsigned long long resultado = 1;
        for (int i = 0; i < k; ++i) {
            resultado *= n - i;
        }
        return resultado;
    }

    int residuoCongruenciaLineal(int a, int b, int m) {
        int x = 0, y = 0;
        int g = euclidesExtendido(a, m, x, y);
        if (b % g != 0) {
            return -1;
        }
        x = (x * (b / g)) % m;
        if (x < 0) {
            x += m;
        }
        return x;
    }

    bool esNumeroFibonacci(int n) {
        int a = 0, b = 1;
        while (b <= n) {
            if (b == n) {
                return true;
            }
            int c = a + b;
            a = b;
            b = c;
        }
        return false;
    }

    // Calcula el número euler
    int eulerPhi(int n) {
        int resultado = n;
        for (int i = 2; i * i <= n; ++i) {
            if (n % i == 0) {
                while (n % i == 0) {
                    n /= i;
                }
                resultado -= resultado / i;
            }
        }
        if (n > 1) {
            resultado -= resultado / n;
        }
        return resultado;
    }

    bool esNumeroLucas(int n) {
        int a = 2, b = 1;
        while (b <= n) {
            if (b == n) {
                return true;
            }
            int c = a + b;
            a = b;
            b = c;
        }
        return false;
    }

    // Calcula el coeficiente binomial (número de combinaciones de n elementos tomados de k en k)
    unsigned long long coeficienteBinomial(int n, int k) {
        if (k > n - k) {
            k = n - k;
        }
        unsigned long long resultado = 1;
        for (int i = 0; i < k; ++i) {
           resultado *= n - i;
           resultado /= i + 1;
        }
        return resultado;
    }

    // Verifica si un número es un número de Harshad (número que es divisible por la suma de sus dígitos)
    bool esNumeroHarshad(int n) {
        int sumaDigitos = 0;
        int numero = n;
        while (numero > 0) {
            sumaDigitos += numero % 10;
            numero /= 10;
        }
        return n % sumaDigitos == 0;
    }

    // Calcula el número de Stirling de segunda clase (S(n, k))
    int numeroStirlingSegundaClase(int n, int k) {
        if (n == k || k == 1) {
            return 1;
        }
        return k * numeroStirlingSegundaClase(n - 1, k) + numeroStirlingSegundaClase(n - 1, k - 1);
    }

    // Calcula el número de Bell (número de particiones de un conjunto de n elementos)
    int numeroBell(int n) {
        int bell[n + 1][n + 1];
        bell[0][0] = 1;
        for (int i = 1; i <= n; ++i) {
            bell[i][0] = bell[i - 1][i - 1];
            for (int j = 1; j <= i; ++j) {
                bell[i][j] = bell[i - 1][j - 1] + bell[i][j - 1];
            }
        }
        return bell[n][0];
    }

    // Calcula la función phi de Möbius (μ(n))
    int funcionPhiMobius(int n) {
        if (n == 1) {
            return 1;
        }
        int resultado = 0;
        for (int i = 2; i * i <= n; ++i) {
            if (n % i == 0) {
                n /= i;
                if (n % i == 0) {
                    return 0;
                }
                resultado = -resultado + funcionPhiMobius(n);
                break;
            }
        }
        if (n > 1) {
            resultado = -resultado + 1;
        }
        return resultado;
    }

    // Función auxiliar para calcular el exponencial (exp) utilizando la serie de Taylor
    double exponencial(double x) {
        const int iteracionesMaximas = 100; // Número máximo de iteraciones para evitar bucle infinito

        double resultado = 1.0;
        double termino = 1.0;
        int n = 1;

        while (n < iteracionesMaximas) {
            termino *= x / n;
            resultado += termino;
            n++;
        }

        return resultado;
    }

    // Función auxiliar para calcular el logaritmo natural (ln) utilizando la serie de Taylor
    double logaritmoNatural(double x) {
        if (x <= 0.0) {
            // El logaritmo no está definido para valores no positivos
            return 0.0;
        }

        const int iteracionesMaximas = 100; // Número máximo de iteraciones para evitar bucle infinito

        double resultado = 0.0;
        double termino = (x - 1) / x;
        int n = 1;

        while (n < iteracionesMaximas) {
            if (n % 2 == 1) {
                resultado += termino;
            } else {
                resultado -= termino;
            }

            termino *= (x - 1) / x;
            n++;
        }

        return resultado;
    }

    // Función auxiliar para calcular el producto de dos números de punto flotante sin errores de redondeo
    double productoSinError(double a, double b) {
        double resultado = 0.0;
        while (b > 0.0) {
            if (b >= 1.0) {
                resultado += a;
                b -= 1.0;
            }
            a *= 2.0;
            b *= 2.0;
        }
        return resultado;
    }

    // Función auxiliar para calcular el valor absoluto de un número de punto flotante
    double valorAbsoluto(double num) {
        return (num < 0) ? -num : num;
    }

    // Función auxiliar para calcular la fracción decimal de un número de punto flotante
    double fraccionDecimal(double num) {
        return num - static_cast<long long>(num);
    }

    // Función auxiliar para calcular la parte entera de un número de punto flotante
    double parteEntera(double num) {
        return num - fraccionDecimal(num);
    }

    // Función para calcular la potencia de un número (base) elevado a un exponente (exponente) para enteros
    int potenciaEntera(int base, int exponente) {
        if (base == 0) {
            return 0;
        }
        if (exponente == 0) {
            return 1;
        }

        int resultado = 1;
        while (exponente > 0) {
            if (exponente % 2 == 1) {
                resultado *= base;
            }
            base *= base;
            exponente /= 2;
        }
        return resultado;
    }

    // Función para calcular la potencia de un número (base) elevado a un exponente (exponente)
    double potencia(double base, double exponente) {
        if (base == 0.0) {
            return 0.0;
        }
        if (exponente == 0.0) {
            return 1.0;
        }

        // Si el exponente es entero, utiliza el algoritmo de exponenciación binaria
        double resultado = 1.0;
        if (fraccionDecimal(exponente) == 0.0) {
            int expEntero = static_cast<int>(valorAbsoluto(exponente));
            if (exponente < 0) {
                base = 1.0 / base;
            }
            while (expEntero > 0) {
                if (expEntero % 2 == 1) {
                    resultado = productoSinError(resultado, base);
                }
                base = productoSinError(base, base);
                expEntero /= 2;
            }
        } else {
            // Si el exponente es fraccional, utiliza la fórmula: base^exponente = exp(exponente * log(base))
            resultado = exponencial(exponente * logaritmoNatural(base));
        }

        return resultado;
    }

    // Función para calcular la función zeta de Riemann (ζ(s))
    double funcionZetaRiemann(double s, int iteraciones = 1000) {
        double resultado = 0.0;
        for (int n = 1; n <= iteraciones; ++n) {
            resultado += 1.0 / (potenciaEntera(n, s));
        }
        return resultado;
    }

    // Verifica si un número es poderoso
    bool esNumeroPoderoso(int n) {
        for (int p = 2; p * p <= n; ++p) {
            if (n % (p * p) == 0) {
                return true;
            }
        }
        return false;
    }

    // Calcula el símbolo de Legendre (a/p)
    int simboloLegendre(int a, int p) {
        a %= p;
        if (a == 0) {
            return 0;
        }
        if (a == 1) {
            return 1;
        }
        if (a == 2) {
            if (p % 8 == 1 || p % 8 == 7) {
                return 1;
            } else {
                return -1;
            }
        }
        if (a == p - 1) {
            if (p % 4 == 1) {
                return 1;
            } else {
                return -1;
            }
        }
        if (mcd(a, p) != 1) {
            return 0;
        }
        if (a % 2 == 0) {
            int signo = (p * p - 1) / 8;
            return signo * simboloLegendre(a / 2, p);
        }
        if ((a - 1) * (p - 1) / 4 % 2 == 0) {
            return simboloLegendre(p % a, a);
        } else {
            return -simboloLegendre(p % a, a);
        }
    }

    // Calcula el símbolo de Jacobi (a/n)
    int simboloJacobi(int a, int n) {
        if (n <= 0 || n % 2 == 0) {
            return 0;
        }
        if (a == 0 || a == 1) {
            return a;
        }
        if (a < 0) {
            if ((n - 1) / 2 % 2 == 0) {
                return simboloJacobi(-a, n);
            } else {
                return -simboloJacobi(-a, n);
            }
        }
        if (a % 2 == 0) {
            if ((n * n - 1) / 8 % 2 == 0) {
                return simboloJacobi(a / 2, n);
            } else {
                return -simboloJacobi(a / 2, n);
            }
        }
        if (a % 4 == 3 && n % 4 == 3) {
            return -simboloJacobi(n, a);
        } else {
            return simboloJacobi(n, a);
        }
    }

    // Calcula el número de Carmichael (número compuesto que cumple a^(n-1) ≡ 1 (mod n) para todo a coprimo con n)
    bool esNumeroCarmichael(int n) {
        if (esPrimo(n)) {
            return false;
        }
        for (int a = 2; a < n; ++a) {
            if (mcd(a, n) == 1 && potenciaModular(a, n - 1, n) != 1) {
                return false;
            }
        }
        return true;
    }

    // Calcula el número de Smith (número cuya suma de dígitos es igual a la suma de dígitos de sus factores primos)
    bool esNumeroSmith(int n) {
        int sumaDigitos = 0;
        int numero = n;
        while (numero > 0) {
            sumaDigitos += numero % 10;
            numero /= 10;
        }
        int sumaDigitosFactores = 0;
        for (int i = 2; i * i <= n; ++i) {
            if (n % i == 0) {
                int sumaDigitosFactor = 0;
                while (n % i == 0) {
                    n /= i;
                    int numero = i;
                    while (numero > 0) {
                        sumaDigitosFactor += numero % 10;
                        numero /= 10;
                    }
                }
                sumaDigitosFactores += sumaDigitosFactor;
            }
        }
        if (n > 1) {
            int sumaDigitosFactor = 0;
            while (n > 0) {
                sumaDigitosFactor += n % 10;
                n /= 10;
            }
            sumaDigitosFactores += sumaDigitosFactor;
        }
        return sumaDigitos == sumaDigitosFactores;
    }

    // Calcula el número de Leyland (número de la forma x^y + y^x)
    bool esNumeroLeyland(int n) {
        for (int x = 2; x * x <= n; ++x) {
            for (int y = 2; y < x; ++y) {
                if (potencia(x, y) + potencia(y, x) == n) {
                    return true;
                }
            }
        }
        return false;
    }

    // Calcula el número de Bernoulli Bn utilizando la fórmula recursiva de Bernoulli
    double numeroBernoulli(int n) {
        if (n == 0) {
            return 1.0;
        }
        if (n == 1) {
            return -0.5;
        }
        double resultado = 0.0;
        for (int k = 0; k < n; ++k) {
            resultado += combinaciones(n + 1, k) * numeroBernoulli(k) / (n + 1 - k);
        }
        return 1.0 - resultado;
    }

    // Calcula el número de Euler utilizando la fórmula de Euler
    double numeroEuler(int n) {
        double resultado = 0.0;
        for (int k = 1; k <= n; ++k) {
            resultado += funcionPhiMobius(k) / potenciaEntera(2, k);
        }
        return 1.0 - resultado;
    }

    // Función phi de Euler: Calcula la cantidad de enteros positivos menores o iguales a n que son coprimos con n.
    int phiEuler(int n) {
        int result = n;
        for (int i = 2; i * i <= n; ++i) {
            if (n % i == 0) {
                while (n % i == 0) {
                    n /= i;
                }
                result -= result / i;
            }
        }
        if (n > 1) {
            result -= result / n;
        }
        return result;
    }

     // Función para calcular la criba de Eratóstenes y encontrar todos los números primos menores o iguales a n.
    // Devuelve un puntero a un arreglo dinámico de enteros que contiene todos los números primos.
    int* cribaEratostenes(int n, int& numPrimos) {
        bool* esPrimo = new bool[n + 1];

        for (int i = 0; i <= n; ++i) {
            esPrimo[i] = true;
        }
        esPrimo[0] = esPrimo[1] = false;

        for (int p = 2; p * p <= n; ++p) {
            if (esPrimo[p]) {
                // Marcar los múltiplos de p como no primos.
                for (int i = p * p; i <= n; i += p) {
                    esPrimo[i] = false;
                }
            }
        }

        // Contar los números primos.
        numPrimos = 0;
        for (int i = 2; i <= n; ++i) {
            if (esPrimo[i]) {
                numPrimos++;
            }
        }

        // Almacenar los números primos en un arreglo dinámico.
        int* primos = new int[numPrimos];
        int idx = 0;
        for (int i = 2; i <= n; ++i) {
            if (esPrimo[i]) {
                primos[idx++] = i;
            }
        }

        delete[] esPrimo;
        return primos;
    }

        // Función para calcular el módulo de un número entero (a % b) de forma segura.
    int modulo(int a, int b) {
        int result = a % b;
        return result >= 0 ? result : result + b;
    }

    // Función para calcular el inverso modular de un número entero 'a' módulo 'm'.
    int inversoModular(int a, int m) {
        int m0 = m, t, q;
        int x0 = 0, x1 = 1;

        if (m == 1) {
            return 0;
        }

        while (a > 1) {
            q = a / m;
            t = m;
            m = a % m, a = t;
            t = x0;
            x0 = x1 - q * x0;
            x1 = t;
        }

        if (x1 < 0) {
            x1 += m0;
        }

        return x1;
    }

    // Función para resolver un sistema de congruencias simultáneas utilizando el teorema chino del resto.
    // Se asume que los módulos son coprimos entre sí.
    int teoremaChinoResto(const int residuos[], const int modulos[], int n) {
        int resultado = 0;
        int M = 1;
        for (int i = 0; i < n; ++i) {
            M *= modulos[i];
        }
        for (int i = 0; i < n; ++i) {
            int Mi = M / modulos[i];
            resultado += residuos[i] * inversoModular(Mi, modulos[i]) * Mi;
            resultado = modulo(resultado, M);
        }
        return resultado;
    }

    // es primo de mersenne(?)
    bool esPrimoMer(int n) {
        if (n <= 1) {
            return false;
        }
        for (int i = 2; i * i <= n; ++i) {
            if (n % i == 0) {
                return false;
            }
        }
        return true;
    }

    bool esPotenciaDeDos(int n) {
        return n > 0 && (n & (n - 1)) == 0;
    }

    bool esNumeroMersenne(int p) {
        return esPrimo(p) && esPotenciaDeDos(p + 1);
    }
  
    // Función para calcular el número de Mersenne M(p) = 2^p - 1
    int numeroMersenne(int p) {
        return potenciaEntera(2, p) - 1;
    }

    // Función para calcular el número de Fermat F(n) = 2^(2^n) + 1
    int numeroFermat(int n) {
        return potenciaEntera(2, potenciaEntera(2, n)) + 1;
    }

    // Función para calcular el número de Lucas-Lehmer L(n) = 2^(n) - 1
    int numeroLucasLehmer(int n) {
        return potenciaEntera(2, n) - 1;
    }

    // Función para calcular el número de Wagstaff W(n) = 3^(n) - 1
    int numeroWagstaff(int n) {
        return potenciaEntera(3, n) - 1;
    }

    // Función para calcular la función φ de Euler para un número n.
    int eulerPhi2(int n) {
        int result = n;
        for (int p = 2; p * p <= n; ++p) {
            if (n % p == 0) {
                while (n % p == 0) {
                    n /= p;
                }
                result -= result / p;
            }
        }
        if (n > 1) {
            result -= result / n;
        }
        return result;
    }

    // Función para verificar si un número es una raíz primitiva módulo n.
    bool esRaizPrimitiva(int g, int n) {
        if (g <= 0 || n <= 1) {
            return false;
        }

        int phi = eulerPhi(n);

        // Utilizamos un arreglo para verificar si los residuos son distintos.
        bool* residuos = new bool[phi];
        for (int i = 0; i < phi; ++i) {
            residuos[i] = false;
        }

        int resultado = 1;
        for (int i = 0; i < phi; ++i) {
            resultado = (resultado * g) % n;
            if (residuos[resultado]) {
                delete[] residuos;
                return false;
            }
            residuos[resultado] = true;
        }

        delete[] residuos;
        return true;
    }

    // Algoritmo extendido de Euclides
    int mcdExtendido(int a, int b, int& x, int& y) {
        if (a == 0) {
            x = 0;
            y = 1;
            return b;
        }
        int x1, y1;
        int gcd = mcdExtendido(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return gcd;
    }

    // Exponenciación rápida módulo 'm'
    int potenciaModular(int base, int exponente, int modulo) {
        int resultado = 1;
        base %= modulo;
        while (exponente > 0) {
            if (exponente % 2 == 1) {
                resultado = (resultado * base) % modulo;
            }
            base = (base * base) % modulo;
            exponente /= 2;
        }
        return resultado;
    }

    // Generador lineal congruente (LCG) para generar números aleatorios
    int lcg(int seed, int a, int c, int m) {
        return (a * seed + c) % m;
    }

    // Test de Miller-Rabin para verificar si un número es probablemente primo
    bool esProbablementePrimo(int n, int iteraciones) {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0) return false;

        int d = n - 1;
        while (d % 2 == 0) d /= 2;

        // Semilla para el generador LCG
        int seed = 42; // Puedes cambiar el número 42 por cualquier otra semilla

        for (int i = 0; i < iteraciones; ++i) {
            int a = 2 + lcg(seed, 1, 1, n - 4); // Generar 'a' aleatoriamente en el rango [2, n - 2]
            seed = lcg(seed, 3, 1, n); // Actualizar la semilla para la siguiente iteración
            int x = potenciaModular(a, d, n);

            if (x == 1 || x == n - 1) continue;

            while (d != n - 1) {
                x = (x * x) % n;
                d *= 2;

                if (x == 1) return false;
                if (x == n - 1) break;
            }

            if (x != n - 1) return false;
        }

        return true;
    }

    // Gauss
    // Función para calcular la suma de los primeros n enteros
    int sumaEnteros(int n) {
        return (n * (n + 1)) / 2;
    }

    // Función para calcular la suma de los cuadrados de los primeros n enteros
    int sumaCuadrados(int n) {
        return (n * (n + 1) * (2 * n + 1)) / 6;
    }

    // Luego:
    // Función para calcular la exponenciación rápida
    int potenciaRapida(int base, int exponente) {
        int resultado = 1;
        while (exponente > 0) {
            if (exponente % 2 == 1) {
                resultado *= base;
            }
            base *= base;
            exponente /= 2;
        }
        return resultado;
    }

    // Función para calcular la exponenciación rápida módulo 'm'
    int potenciaModular(int base, int exponente, int modulo) {
        int resultado = 1;
        base %= modulo;
        while (exponente > 0) {
            if (exponente % 2 == 1) {
                resultado = (resultado * base) % modulo;
            }
            base = (base * base) % modulo;
            exponente /= 2;
        }
        return resultado;
    }

    // Función para calcular el inverso multiplicativo módulo 'm'
    int inversoMultiplicativo(int a, int modulo) {
        int m0 = modulo;
        int y = 0, x = 1;
        if (modulo == 1) {
            return 0;
        }
        while (a > 1) {
            int q = a / modulo;
            int t = modulo;
            modulo = a % modulo, a = t;
            t = y;
            y = x - q * y;
            x = t;
        }
        if (x < 0) {
            x += m0;
        }
        return x;
    }

    // Función de Ackermann
    int ackermann(int m, int n) {
        if (m == 0) return n + 1;
        if (n == 0) return ackermann(m - 1, 1);
        return ackermann(m - 1, ackermann(m, n - 1));
    }

    // Números de Fibonacci usando iteración
    int fibonacci(int n) {
        if (n <= 0) return 0;
        if (n == 1) return 1;

        int prev = 0;
        int current = 1;
        for (int i = 2; i <= n; ++i) {
            int temp = current;
            current = prev + current;
            prev = temp;
        }

        return current;
    }

    // Función para calcular el polinomio de Chebyshev de primera clase T_n(x)
    int polinomioChebyshevPrimeraClase(int n, int x) {
        if (n == 0) return 1;
        if (n == 1) return x;
        
        int Tn_minus_1 = x;
        int Tn_minus_2 = 1;
        int Tn = 0;

        for (int i = 2; i <= n; ++i) {
            Tn = 2 * x * Tn_minus_1 - Tn_minus_2;
            Tn_minus_2 = Tn_minus_1;
            Tn_minus_1 = Tn;
        }

        return Tn;
    }

    // Función para calcular el polinomio de Chebyshev de segunda clase U_n(x)
    int polinomioChebyshevSegundaClase(int n, int x) {
        if (n == 0) return 1;
        if (n == 1) return 2 * x;

        int Un_minus_1 = 2 * x;
        int Un_minus_2 = 1;
        int Un = 0;

        for (int i = 2; i <= n; ++i) {
            Un = 2 * x * Un_minus_1 - Un_minus_2;
            Un_minus_2 = Un_minus_1;
            Un_minus_1 = Un;
        }

        return Un;
    }

    // Función para verificar el teorema de Fermat para exponentes pequeños (n <= 20)
    bool verificarTeoremaFermat(int n) {
        if (n <= 2) {
            // No hay soluciones para n <= 2
            return false;
        }

        for (int x = 1; x <= 20; ++x) {
            for (int y = 1; y <= 20; ++y) {
                for (int z = 1; z <= 20; ++z) {
                    int lhs = potenciaModular(x, n, z) + potenciaModular(y, n, z);
                    int rhs = potenciaModular(z, n, z);
                    if (lhs == rhs) {
                        // Se encontró una solución para la ecuación x^n + y^n = z^n
                        return true;
                    }
                }
            }
        }

        // No se encontraron soluciones para la ecuación x^n + y^n = z^n
        return false;
    }

    // Función para calcular el valor de 'sin(x)' utilizando la expansión en serie de Taylor
    double sin(double x, int iteraciones = 10) {
        double resultado = 0.0;
        for (int n = 0; n < iteraciones; ++n) {
            double coeficiente = (n % 2 == 0) ? 1.0 : -1.0;
            resultado += coeficiente * potenciaEntera(x, 2 * n + 1) / factorial(2 * n + 1);
        }
        return resultado;
    }

    // Función para calcular el valor de 'cos(x)' utilizando la expansión en serie de Taylor
    double cos(double x, int iteraciones = 10) {
        double resultado = 0.0;
        for (int n = 0; n < iteraciones; ++n) {
            double coeficiente = (n % 2 == 0) ? 1.0 : -1.0;
            resultado += coeficiente * potenciaEntera(x, 2 * n) / factorial(2 * n);
        }
        return resultado;
    }

    // Estructura para representar números complejos
    struct Complejo {
        double real;
        double imag;
    };

    // Función para calcular la función zeta de Riemann (ζ) utilizando la serie de Dirichlet
    Complejo funcionZetaRiemann(const Complejo& s, int iteraciones = 1000) {
        Complejo resultado = {0.0, 0.0};
        for (int n = 1; n <= iteraciones; ++n) {
            double realPart = 1.0 / potenciaEntera(n, s.real) * cos(s.imag * logaritmoNatural(n));
            double imagPart = 1.0 / potenciaEntera(n, s.real) * sin(s.imag * logaritmoNatural(n));
            resultado.real += realPart;
            resultado.imag += imagPart;
        }
        return resultado;
    }

    // Función para calcular la función eta de Riemann (η) utilizando la serie de Dirichlet
    Complejo funcionEtaRiemann(const Complejo& s, int iteraciones = 1000) {
        Complejo resultado = {0.0, 0.0};
        Complejo signo = {1.0, 0.0};
        for (int n = 1; n <= iteraciones; ++n) {
            double realPart = signo.real / potenciaEntera(n, s.real) * cos(s.imag * logaritmoNatural(n));
            double imagPart = signo.imag / potenciaEntera(n, s.real) * sin(s.imag * logaritmoNatural(n));
            resultado.real += realPart;
            resultado.imag += imagPart;
            signo = {-signo.real, -signo.imag};
        }
        return resultado;
    }

    // Función para calcular la función gamma de Legendre (Γ) utilizando la fórmula de Legendre
    double funcionGammaLegendre(double x) {
        if (x <= 0.0) {
            // La función gamma de Legendre no está definida para valores no positivos
            return 0.0;
        }

        const int iteracionesMaximas = 100; // Número máximo de iteraciones para evitar bucle infinito

        double resultado = 1.0;
        double termino = 1.0;
        int n = 1;

        while (n < iteracionesMaximas) {
            termino *= (x + n) / exponencial(n);
            resultado += termino;
            n++;
        }

        return resultado / x;
    }

    // Función para calcular la función gamma de Euler (Γ) utilizando la fórmula de Euler
    double funcionGammaEuler(double x) {
        if (x <= 0.0) {
            // La función gamma de Euler no está definida para valores no positivos
            return 0.0;
        }

        const int iteracionesMaximas = 100; // Número máximo de iteraciones para evitar bucle infinito

        double resultado = 1.0;
        double termino = 1.0;
        int n = 1;

        while (n < iteracionesMaximas) {
            termino *= (x + n) / n;
            resultado *= termino;
            n++;
        }

        return resultado;
    }

    // Aproximación de π utilizando la fórmula de Leibniz para π/4
    double aproximacionPi(int iteraciones = 1000000) {
        double pi = 0.0;
        for (int i = 0; i < iteraciones; ++i) {
            pi += (i % 2 == 0) ? 1.0 / (2 * i + 1) : -1.0 / (2 * i + 1);
        }
        return 4.0 * pi;
    }

    // Aproximación de e utilizando la fórmula de Euler para e
    double aproximacionE(int iteraciones = 100) {
        double e = 1.0;
        double factorial = 1.0;
        for (int i = 1; i <= iteraciones; ++i) {
            factorial *= i;
            e += 1.0 / factorial;
        }
        return e;
    }

    // Función para calcular la función gamma de Stirling (Γ) utilizando aproximaciones de π y e
    double funcionGammaStirling(double x, int iteracionesPi = 1000000, int iteracionesE = 100) {
        if (x <= 0.0) {
            // La función gamma de Stirling no está definida para valores no positivos
            return 0.0;
        }

        double pi = aproximacionPi(iteracionesPi);
        double e = aproximacionE(iteracionesE);

        return sqrt(2.0 * pi / x) * potenciaEntera(x / e, x);
    }

    // Función auxiliar para calcular el factorial de un número de punto flotante
    double factorial(double n)
    {
        if (n <= 1.0) {
            return 1.0;
        }
        return n * factorial(n - 1.0);
    }

    // Función para calcular la función gamma de Lanczos (Γ)
    double funcionGammaLanczos(double x, int iteracionesPi = 1000000) {
        if (x <= 0.0) {
            // La función gamma de Lanczos no está definida para valores no positivos
            return 0.0;
        }

        double pi = aproximacionPi(iteracionesPi);
        const double raizDosPi = sqrt(2.0 * pi);

        // Los coeficientes de la expansión en serie de Stirling
        const double coeficientes[] = {
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };

        double suma = coeficientes[0];
        for (int i = 1; i < 9; ++i) {
            suma += coeficientes[i] / (x + i);
        }

        double potencia = x + 7.5;
        double exponente = x + 0.5;
        return raizDosPi * potenciaEntera(exponente, potencia) * exponencial(-exponente) * suma;
    }

} // namespace BibliotecaMatematica
