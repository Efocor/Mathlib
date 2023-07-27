#include <iostream>
#include <vector>

using namespace std;

// Función para calcular el máximo común divisor (MCD) de dos números utilizando el algoritmo de Euclides.
int mcd(int a, int b) {
    if (b == 0)
        return a;
    return mcd(b, a % b);
}

// Función para calcular el inverso multiplicativo de 'a' módulo 'b'.
int inversoModular(int a, int b) {
    a = a % b;
    for (int x = 1; x < b; x++) {
        if ((a * x) % b == 1)
            return x;
    }
    return 1;
}

// Función para aplicar el Teorema Chino del Resto y obtener la solución única del sistema de congruencias.
int teoremaChinoResto(const vector<int>& residuos, const vector<int>& modulos) {
    int n = residuos.size();
    
    // Validar que el número de congruencias sea mayor o igual a 2.
    if (n < 2) {
        cout << "Error: Debe ingresar al menos dos congruencias para aplicar el Teorema Chino del Resto." << endl;
        return -1; // Valor inválido para indicar un error.
    }

    int productoModulos = 1; // Producto de todos los módulos.

    // Paso 1: Calcular el producto de todos los módulos.
    for (int i = 0; i < n; i++) {
        // Validar que los módulos sean positivos.
        if (modulos[i] <= 0) {
            cout << "Error: Los módulos deben ser números positivos." << endl;
            return -1;
        }

        productoModulos *= modulos[i];
    }

    // Paso 2: Calcular las constantes Mi y xi para cada congruencia.
    vector<int> constantesMi(n), constantesXi(n);
    for (int i = 0; i < n; i++) {
        constantesMi[i] = productoModulos / modulos[i];

        // Validar que los módulos sean primos entre sí dos a dos.
        for (int j = i + 1; j < n; j++) {
            if (mcd(modulos[i], modulos[j]) != 1) {
                cout << "Error: Los módulos deben ser primos entre sí para aplicar el Teorema Chino del Resto." << endl;
                return -1;
            }
        }

        // Calcular el inverso multiplicativo módulo modulos[i].
        constantesXi[i] = inversoModular(constantesMi[i], modulos[i]);
    }

    // Paso 3: Calcular el resultado utilizando la fórmula del Teorema Chino del Resto.
    int resultado = 0;
    for (int i = 0; i < n; i++) {
        // Aplicar la fórmula para obtener el resultado parcial.
        resultado = (resultado + (residuos[i] * constantesMi[i] * constantesXi[i]) % productoModulos) % productoModulos;
    }

    // El resultado debe estar en el rango [0, productoModulos-1].
    if (resultado < 0)
        resultado += productoModulos;

    return resultado;
}

int main() {
    cout << "Teorema Chino del Resto:" << endl;
    cout << "Resultado matemático que proporciona una solución única para un sistema de congruencias lineales," << endl;
    cout << "útil cuando se trabaja con grandes números o cuando se necesita resolver problemas que involucran múltiples módulos de forma más eficiente." << endl << endl;

    int n;
    cout << "Ingrese la cantidad de congruencias: ";
    cin >> n;

    vector<int> residuos(n), modulos(n);

    cout << "Ingrese las congruencias en el formato 'a i mod m i':" << endl;
    for (int i = 0; i < n; i++) {
        cout << "Congruencia " << i + 1 << ": ";
        cin >> residuos[i] >> modulos[i];
    }

    // Calcular la solución utilizando el Teorema Chino del Resto.
    int resultado = teoremaChinoResto(residuos, modulos);

    if (resultado != -1) {
        cout << "La solución para el sistema de congruencias es: " << resultado << endl;
    }

    return 0;
}
