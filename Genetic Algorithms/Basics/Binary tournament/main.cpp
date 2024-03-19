#include <iostream>
#include <vector>

// Función de torneo binario, devuelve los índices de aquellos individuos aptos
std::vector<int> seleccion_torneo_binario(const std::vector<int>& baraja, const std::vector<int>& aptitudes) {
    std::vector<int> sobrevivientes;
    for (size_t i = 0; i < baraja.size(); i += 2) {
        if (i + 1 < baraja.size()) {
            int indice1 = baraja[i] - 1;
            int indice2 = baraja[i + 1] - 1;
            if (aptitudes[indice1] > aptitudes[indice2]) {
                sobrevivientes.push_back(baraja[i]);
            } else {
                sobrevivientes.push_back(baraja[i + 1]);
            }
        }
    }
    return sobrevivientes;
}

int main() {
    std::vector<int> aptitudes = {178, 350, 245, 368, 122, 100};
    std::vector<int> baraja = {2, 5, 1, 3, 4, 6};
    std::vector<int> sobrevivientes = seleccion_torneo_binario(baraja, aptitudes);

    std::cout << "Individuos seleccionados como sobrevivientes: ";
    for (int indice : sobrevivientes) {
        std::cout << indice << " ";
    }

    return 0;
}
