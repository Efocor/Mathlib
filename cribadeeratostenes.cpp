// Por Felipe Alexander Correa Rodríguez

// Criba de Eratóstenes:
// Se trata de un código que permite encontrar todos los números primos hasta un límite dado.
// Complejidad: O(n log log n)
// Encontrar esos primos sirve para resolver problemas como:
// - Encontrar los divisores de un número.
// - Encontrar la factorización de un número.
// - Encontrar el número de divisores de un número.

// La idea de este código es probar todos los números desde 2 hasta el límite, y para cada número primo encontrado,
// marcar todos sus múltiplos como no primos. Al final, los números que no hayan sido marcados como no primos serán primos.

// Librerías
#include <iostream> 
#include <cstring> // Para usar memset

using namespace std;

const int MAX_LIMITE = 1000000; // Límite máximo para la Criba de Eratóstenes.

// Función auxiliar para mostrar un mensaje de error y salir del programa.
void mostrarError(const string& mensaje) {
    cout << "Error: " << mensaje << endl;
    exit(1);
}

// Función para calcular la raíz cuadrada (New-Raph Meth).
double raizCuadrada(int x) {
    if (x < 0) {
        mostrarError("No se puede calcular la raíz cuadrada de un número negativo.");
    }

    double aprox = x;
    double diff = 1.0;

    while (diff > 1e-9) { // 1e-9 = 10^-9
        double siguiente = 0.5 * (aprox + x / aprox);
        diff = abs(aprox - siguiente);
        aprox = siguiente;
    }

    return aprox;
}

// Función para aplicar la Criba de Eratóstenes con optimización de la Criba de Euler.
int cribaDeEratostenesOptimizada(int limite, int primos[]) {
    if (limite < 0 || limite > MAX_LIMITE) {
        mostrarError("El límite debe estar en el rango de 0 a " + to_string(MAX_LIMITE) + ".");
    }

    const int TAM_BITSET = (MAX_LIMITE + 1) / 64 + 1;
    unsigned long long esPrimo[TAM_BITSET];
    memset(esPrimo, 0xFF, sizeof(esPrimo)); // Inicializar todos los bits en 1.

    int menorFactor[MAX_LIMITE + 1];
    for (int i = 0; i <= limite; i++) {
        menorFactor[i] = i;
    }

    int cantPrimos = 0;
    esPrimo[0] &= ~(1ULL);
    esPrimo[1 / 64] &= ~(1ULL << (1 % 64)); // 1 no es primo.

    int limiteRaiz = 1 + static_cast<int>(raizCuadrada(limite));
    for (int p = 2; p <= limiteRaiz; p++) {
        if (esPrimo[p / 64] & (1ULL << (p % 64))) { // Si p es primo.
            primos[cantPrimos++] = p; // Agregar p a la lista de primos.
            for (int i = p * p; i <= limite; i += p) {
                esPrimo[i / 64] &= ~(1ULL << (i % 64));
                menorFactor[i] = min(menorFactor[i], p);
            }
        }
    }

    for (int i = limiteRaiz + 1; i <= limite; i++) {
        if (esPrimo[i / 64] & (1ULL << (i % 64))) { // Si i es primo.
            primos[cantPrimos++] = i; // Agregar i a la lista de primos.
            menorFactor[i] = i;
        }
    }

    return cantPrimos;
}

int main() {
    cout << "Criba de Eratóstenes con Optimización de la Criba de Euler" << endl;

    int limite;
    cout << "Ingrese el límite para encontrar los números primos (0 <= limite <= " << MAX_LIMITE << "): ";
    cin >> limite;

    int primos[MAX_LIMITE];
    int cantPrimos = cribaDeEratostenesOptimizada(limite, primos);

    cout << "Los números primos hasta " << limite << " son: ";
    for (int i = 0; i < cantPrimos; i++) {
        cout << primos[i] << " "; // Debería imprimir ej: 2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 ...
    }
    cout << endl;

    return 0;
}
