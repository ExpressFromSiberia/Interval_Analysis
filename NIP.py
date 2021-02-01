from Algorithm import ProccesCycle

import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np


# Вывод ответа - вектора решений
def PrintAnswer(answer):
    # Блок вывода решения
    Dim = int(len(answer) // 2)
    for i in range(0, Dim):
        u_0 = answer[2*i]
        u_1 = answer[2*i + 1]
        print(i + 1, ": [", u_0, ",", u_1, "]", end="")
        if u_0 is not None and u_1 is not None:
            if u_0 <= u_1:
                print("    ->")
            else:
                print("     <-")
        else:
            print()


# Блок считывания данных из файла
def GetDataFromFile():
    with open("interval_data_256x16.txt") as fi:
        filelist = [line.strip('\n') for line in fi]

        # Считываем первую строку - строку параметров
        # lines_count - количество строк матрицы
        # columns_count - количество столбцов матрицы
        # Eps - точность искомого решения
        # Tau -  релаксационный параметр
        # IterLim - максимальное количество итераций
        # MatrixCount - количество матриц, необходимых для нахождения решения
        # DetLim - порог для определителя матриц, необходимых для нахождения решения
        lines_count, columns_count, Eps, Tau, IterLim, MatrixCount, DetLim = filelist[0].split()
        lines_count = int(lines_count)
        columns_count = int(columns_count)
        Eps = float(Eps)
        Tau = float(Tau)
        IterLim = int(IterLim)
        MatrixCount = int(MatrixCount)
        DetLim = float(DetLim)

        # Инициализируем память под данные
        C = [888 for _ in range(lines_count * 2 * columns_count)]
        d = [None for _ in range(2 * lines_count)]
        # Обработка последней строки - строки правой части
        last_str = list(map(float, filelist[-1].split()))
        for i, j in zip(range(0, lines_count, 1), range(0, 2 * lines_count, 2)):
            d[2 * i] = last_str[j]
            d[2 * i + 1] = last_str[j + 1]

        # Обработка для интервальной матрицы
        for i in range(2, len(filelist) - 2):
            T = list(map(float, filelist[i].split()))
            i -= 2
            k = 0
            for j in range(0, len(T), 2):
                C[i * 2 * columns_count + 2 * k] = T[j]
                C[i * 2 * columns_count + 2 * k + 1] = T[j + 1]
                k += 1

    return lines_count, columns_count, Eps, Tau, IterLim, MatrixCount, DetLim, C, d


# Разделяем исходную матрицу на четверти
def SeparationMatrix(Matrix, d, columns_count, lines_count):
    quarter_matrix_array = [[None for _ in range(0, len(Matrix) // 4)], [None for _ in range(0, len(Matrix) // 4)],
                            [None for _ in range(0, len(Matrix) // 4)], [None for _ in range(0, len(Matrix) // 4)]]
    half_d_array = [[None for _ in range(0, len(d) // 2)], [None for _ in range(0, len(d) // 2)]]

    for i in range(0, lines_count // 2):
        for j in range(0, columns_count):
            quarter_matrix_array[0][i*columns_count + j] = Matrix[i*2*columns_count + j]
            quarter_matrix_array[1][i*columns_count + j] = Matrix[i*2*columns_count + j + columns_count]
            quarter_matrix_array[2][i*columns_count + j] = Matrix[i*2*columns_count + lines_count + j]
            quarter_matrix_array[3][i * columns_count + j] = Matrix[i * 2 * columns_count + lines_count
                                                                    + j + columns_count]
        for k in range(2):
            half_d_array[0][i*2 + k] = d[i*2 + k]
            half_d_array[1][i*2 + k] = d[i*2 + lines_count + k]

    return quarter_matrix_array, half_d_array


def main():
    lines_count, columns_count, Eps, Tau, IterLim, MatrixCount, DetLim, C, d = GetDataFromFile()
    plotMatrix = np.array(C).reshape((lines_count, 2*columns_count))

    plt.figure()
    sbn.heatmap(plotMatrix, )
    plt.title('Matrix')
    plt.xlabel('Columns')
    plt.ylabel('Rows')
    plt.savefig("Matrix.png")

    quarter_matrix_array, half_d_array = SeparationMatrix(C, d, columns_count, lines_count)

    answer = [[None for _ in range(0, columns_count)], [None for _ in range(0, columns_count)],
              [None for _ in range(0, columns_count)], [None for _ in range(0, columns_count)]]

    for i in range(0, 4):
        answer[i] = ProccesCycle(quarter_matrix_array[i], half_d_array[i // 2], columns_count // 2,
                              lines_count // 2, Tau, Eps, IterLim, DetLim)
        print("\n\nВектор решений полученный для", i+1,"четверти исходной матрицы \n")
        PrintAnswer(answer[i])

    # answer = ProccesCycle(C, d, columns_count, lines_count, Tau, Eps, IterLim, DetLim)
    # PrintAnswer(answer)

    return 0


if __name__ == "__main__":
    main()
