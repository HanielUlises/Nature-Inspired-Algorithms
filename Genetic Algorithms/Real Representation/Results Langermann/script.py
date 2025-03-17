import numpy as np

# Función para procesar un solo archivo
def procesar_archivo(archivo):
    fitness_values = []
    with open(archivo, 'r') as f:
        for line in f:
            # Buscamos la línea que contiene la palabra "Fitness"
            if 'Fitness' in line:
                # Extraemos el valor de Fitness
                parts = line.split('|')
                fitness = float(parts[2].split(':')[1].strip())
                fitness_values.append(fitness)
    return fitness_values

# Listas para almacenar todos los resultados de los 10 archivos
todos_fitness = []

# Procesamos los archivos results0.txt, results1.txt, ..., results9.txt
for i in range(10):
    archivo = f'Genetic Algorithms/Real Representation/Results/results{i}.txt'
    fitness_values = procesar_archivo(archivo)
    todos_fitness.append(fitness_values)

# Convertir todos los resultados a un solo array para análisis
todos_fitness = np.array(todos_fitness)

# Calcular las estadísticas
mejor = np.min(todos_fitness)  # Mejor fitness (el más cercano a 0 si es un problema de minimización)
peor = np.max(todos_fitness)   # Peor fitness (el más alejado de 0 si es un problema de minimización)
media = np.mean(todos_fitness) # Media de los fitness
desviacion_estandar = np.std(todos_fitness)  # Desviación estándar

# Mostrar los resultados
print(f"Mejor fitness: {mejor}")
print(f"Peor fitness: {peor}")
print(f"Media: {media}")
print(f"Desviacion estandar: {desviacion_estandar}")
