from math import fabs
from Verification import VerificationAnswer
import numpy as np


# Метод Гаусса
def solver(xx, DIM, F):
    ip = [999 for _ in range(DIM)]
    ier = 0
    ip[DIM - 1] = 1

    for k in range(0, DIM - 1):
        m = k
        p = fabs(F[k * DIM + k])
        # Ищем максимальный элемент в строке k
        for i in range(k + 1, DIM):
            q = fabs(F[i * DIM + k])
            if p < q:
                m = i
                p = q

        ip[k] = m
        p = F[m * DIM + k]

        if m != k:
            ip[DIM - 1] = -ip[DIM - 1]
            F[m * DIM + k] = F[k * DIM + k]
            F[k * DIM + k] = p

        if p == 0:
            ier = k
            ip[DIM - 1] = 0
            break
        p = 1. / p

        for i in range(k + 1, DIM):
            F[i * DIM + k] *= -p

        for j in range(k + 1, DIM):
            p = F[m * DIM + j]
            F[m * DIM + j] = F[k * DIM + j]
            F[k * DIM + j] = p
            if p != 0:
                for i in range(k + 1, DIM):
                    F[i * DIM + j] += F[i * DIM + k] * p

    if ier != 0 or F[DIM * DIM - 1] == 0:
        print("\n\r   Sorry, the interval matrix of the equation")
        print("\n is not absolutely regular.")
        print("\n\r   Try to change it a little. \n")

    for k in range(0, DIM - 1):
        m = ip[k]
        q = xx[m]
        xx[m] = xx[k]
        xx[k] = q
        for i in range(k + 1, DIM):
            xx[i] += F[i * DIM + k] * q

    for j in range(0, DIM - 1):
        k = DIM - j - 1
        xx[k] /= F[k * DIM + k]
        q = -xx[k]
        for i in range(0, k):
            xx[i] += F[i * DIM + k] * q

    xx[0] /= F[0]
    return 0


# Субдифференциальный метод Ньютона
def Algorithm(Dim, DIM, C, F, d, xx, Tau, Eps, IterLim):
    x = [999 for _ in range(DIM)]
    check = solver(xx, DIM, F)
    ni = 0

    while True:
        ni += 1
        r = 0

        for i in range(0, DIM):
            x[i] = xx[i]
            for j in range(0, DIM):
                F[i * DIM + j] = 0

        for i in range(0, Dim):
            s_0 = 0
            s_1 = 0
            for j in range(0, Dim):
                g_0 = C[i * DIM + 2 * j]
                g_1 = C[i * DIM + 2 * j + 1]
                h_0 = x[j]
                h_1 = x[j + Dim]

                if g_0 * g_1 > 0:
                    l = 0 if g_0 > 0 else 2
                else:
                    l = 1 if g_0 <= g_1 else 3
                if h_0 * h_1 > 0:
                    m = 1 if h_0 > 0 else 3
                else:
                    m = 2 if h_0 <= h_1 else 4

                tmp = 4 * l + m
                if tmp == 1:
                    t_0 = g_0 * h_0
                    t_1 = g_1 * h_1
                    F[i * DIM + j] = g_0
                    F[(i + Dim) * DIM + j + Dim] = g_1
                if tmp == 2:
                    t_0 = g_1 * h_0
                    t_1 = g_1 * h_1
                    F[i * DIM + j] = g_1
                    F[(i + Dim) * DIM + j + Dim] = g_1
                if tmp == 3:
                    t_0 = g_1 * h_0
                    t_1 = g_0 * h_1
                    F[i * DIM + j] = g_1
                    F[(i + Dim) * DIM + j + Dim] = g_0
                if tmp == 4:
                    t_0 = g_0 * h_0
                    t_1 = g_0 * h_1
                    F[i * DIM + j] = g_0
                    F[(i + Dim) * DIM + j + Dim] = g_0
                if tmp == 5:
                    t_0 = g_0 * h_1
                    t_1 = g_1 * h_1
                    F[i * DIM + j + Dim] = g_0
                    F[(i + Dim) * DIM + j + Dim] = g_1
                if tmp == 6:
                    u_0 = g_0 * h_1
                    v_0 = g_1 * h_0
                    u_1 = g_0 * h_0
                    v_1 = g_1 * h_1
                    if u_0 < v_0:
                        t_0 = u_0
                        F[i * DIM + j + Dim] = g_0
                    else:
                        t_0 = v_0
                        F[i * DIM + j] = g_1
                    if u_1 > v_1:
                        t_1 = u_1
                        F[(i + Dim) * DIM + j] = g_0
                    else:
                        t_1 = v_1
                        F[(i + Dim) * DIM + j + Dim] = g_1
                if tmp == 7:
                    t_0 = g_1 * h_0
                    t_1 = g_0 * h_0
                    F[i * DIM + j] = g_1
                    F[(i + Dim) * DIM + j] = g_0
                if tmp == 8:
                    t_0 = 0
                    t_1 = 0
                if tmp == 9:
                    t_0 = g_0 * h_1
                    t_1 = g_1 * h_0
                    F[i * DIM + j + Dim] = g_0
                    F[(i + Dim) * DIM + j] = g_1
                if tmp == 10:
                    t_0 = g_0 * h_1
                    t_1 = g_0 * h_0
                    F[i * DIM + j + Dim] = g_0
                    F[(i + Dim) * DIM + j] = g_0
                if tmp == 11:
                    t_0 = g_1 * h_1
                    t_1 = g_0 * h_0
                    F[i * DIM + j + Dim] = g_1
                    F[(i + Dim) * DIM + j] = g_0
                if tmp == 12:
                    t_0 = g_1 * h_1
                    t_1 = g_1 * h_0
                    F[i * DIM + j + Dim] = g_1
                    F[(i + Dim) * DIM + j] = g_1
                if tmp == 13:
                    t_0 = g_0 * h_0
                    t_1 = g_1 * h_0
                    F[i * DIM + j] = g_0
                    F[(i + Dim) * DIM + j] = g_1
                if tmp == 14:
                    t_0 = 0
                    t_1 = 0
                if tmp == 15:
                    t_0 = g_1 * h_1
                    t_1 = g_0 * h_1
                    F[i * DIM + j + Dim] = g_1
                    F[(i + Dim) * DIM + j + Dim] = g_0
                if tmp == 16:
                    u_0 = g_0 * h_0
                    v_0 = g_1 * h_1
                    u_1 = g_0 * h_1
                    v_1 = g_1 * h_0
                    if u_0 > v_0:
                        t_0 = u_0
                        F[i * DIM + j] = g_0
                    else:
                        t_0 = v_0
                        F[i * DIM + j + Dim] = g_1
                    if u_1 < v_1:
                        t_1 = u_1
                        F[(i + Dim) * DIM + j + Dim] = g_0
                    else:
                        t_1 = v_1
                        F[(i + Dim) * DIM + j] = g_1
                s_0 += t_0
                s_1 += t_1
            xx[i] = t_0 = s_0 - d[2 * i]
            xx[i + Dim] = t_1 = s_1 - d[2 * i + 1]
            t_0 = fabs(t_0)
            t_1 = fabs(t_1)
            r += t_0 if t_0 > t_1 else t_1

        check = solver(xx, DIM, F)
        q = 0
        for i in range(0, DIM):
            xx[i] = x[i] - xx[i] * Tau
            q += fabs(xx[i])
        if q == 0:
            q = 1

        if r / q < Eps or ni > IterLim:
            # print("Number of iterations = ", ni)
            # print("1-norm of defect of approximate solution = ", r)
            break


# Удаление столбцов из матрицы (точечной), которые малоинформативны
def RemoveColumn(A, lines_count, columns_count):
    resColVal = 0
    colSum = [0 for _ in range(columns_count)]
    remaining_columns = []

    check = np.array(A).reshape((lines_count, columns_count))

    for i in range(0, columns_count):
        for j in range(0, lines_count):
            colSum[i] += A[i + j*columns_count]

    for i in range(0, len(colSum)):
        if abs(colSum[i]) > 0:
            # print("The sum of the values in the", i, "column = ", abs(colSum[i]), "\n")
            resColVal += 1
            remaining_columns.append(i)

    remA = [None for _ in range(lines_count * resColVal)]
    k = 0

    check_A = np.array(A).reshape((lines_count, columns_count))

    for i in range(0, lines_count):
        for j in remaining_columns:
                remA[k] = A[j + i*columns_count]
                k += 1

    return remA, resColVal, remaining_columns


# "Отрезаем" квадратную матрицу по заданной комбинации строк
def GetCuttingMatrix(A, lines_count, columns_count, combination):
    cuttingA = [None for _ in range(0, len(combination) * columns_count)]
    k = 0

    for i in combination:
        for j in range(0, columns_count):
            cuttingA[k] = A[i * columns_count + j]
            k += 1

    return cuttingA


# Получение интервальной квадратной матрицы, полученной из изначальной
def GetIntervalSquareMatrix(C, columns_count, lines_combination, remaining_columns):
    squareMatrix_lines = len(lines_combination)
    squareMatrix_columns = 2 * len(remaining_columns)
    squareA = [None for _ in range(squareMatrix_lines * squareMatrix_columns)]

    for i in range(0, squareMatrix_lines):
        k = 0
        for j in range(0, squareMatrix_columns, 2):
            squareA[i * squareMatrix_columns + j] = C[lines_combination[i] * 2 * columns_count
                                                      + 2 * remaining_columns[k]]
            squareA[i * squareMatrix_columns + j + 1] = C[
                lines_combination[i] * 2 * columns_count + 2 * remaining_columns[k] + 1]
            k += 1

    return squareA, squareMatrix_lines, squareMatrix_columns


# Получение соответствующих конкретной матрице значений вектора правой части
def GetRightHandVector(d, lines_combination):
    specialRightHandVector = [None for _ in range(2 * len(lines_combination))]
    k = 0

    for i in range(0, 2 * len(lines_combination), 2):
        specialRightHandVector[i] = d[2 * lines_combination[k]]
        specialRightHandVector[i + 1] = d[2 * lines_combination[k] + 1]
        k += 1
    return specialRightHandVector


# Тут я вынес создание матрицы F
def GetFMatrix(C, Dim, DIM):
    F = [None for _ in range(DIM * DIM)]
    for i in range(0, Dim):
        for j in range(0, Dim):
            p = 0.5 * (C[i * DIM + 2 * j] + C[i * DIM + 2 * j + 1])
            if p >= 0:
                F[i * DIM + j] = p
                F[(i + Dim) * DIM + j + Dim] = p
                F[(i + Dim) * DIM + j] = 0
                F[i * DIM + j + Dim] = 0
            else:
                F[i * DIM + j] = 0
                F[(i + Dim) * DIM + j + Dim] = 0
                F[(i + Dim) * DIM + j] = p
                F[i * DIM + j + Dim] = p
    return F


def GetxxVector(d, Dim, DIM):
    xx = [None for _ in range(DIM)]
    for i, j in zip(range(0, Dim, 1), range(0, DIM, 2)):
        xx[i] = d[j]
        xx[i + Dim] = d[j + 1]
    return xx


# Получаем массив индексов, которые определяют квадратные матрицы с определителем, выше заданного порога
def GetArrayMatrix(remC, lines_count, resColVal, DetLim, first_line, last_line):
    # Оставляем только те комбинации строк, которые образуют матрицы с "хорошим" определителем
    k = 0
    check_remC = np.array(remC).reshape((lines_count, resColVal))
    for i in range(last_line, lines_count + 1):
        lines_combination = [i for i in range(first_line + k, last_line + k)]
        cuttingMatrix = GetCuttingMatrix(remC, lines_count, resColVal, lines_combination)
        check_cuttingMatrix = np.array(cuttingMatrix).reshape((resColVal, resColVal))
        tmpDet = np.linalg.det(np.array(cuttingMatrix).reshape((resColVal, resColVal)))
        # print("Взяты строки с ", first_line + k, " по ", last_line + k, ". Определитель = ", "%.e" % tmpDet)
        k += 1

        if abs(tmpDet) >= DetLim:
            yield lines_combination, list(range(resColVal))
        else:
            rem_cuttingC, rem_column_cuttingC, remaining_columns = RemoveColumn(cuttingMatrix, resColVal, resColVal)
            if 0 < rem_column_cuttingC < resColVal:
                for l, c in GetArrayMatrix(rem_cuttingC, last_line - first_line, rem_column_cuttingC, DetLim,
                                           0, rem_column_cuttingC):
                    yield np.array(lines_combination)[l].tolist(), np.array(remaining_columns)[c].tolist()


# Вылонение цикла для получение решений по квадратным матрицам
def ProccesCycle(C, d, columns_count, lines_count, Tau, Eps, IterLim, DetLim):

    # Получаем для основой матрицы точечный вариант
    midC = GetMidMatrix(C)
    check_midC = np.array(midC).reshape((lines_count, columns_count))

    # Удаляем "вредные" столбцы
    remC, resColVal, remaining_columns = RemoveColumn(midC, lines_count, columns_count)
    check_remC = np.array(remC).reshape((lines_count, resColVal))

    squarexx = [999 for _ in range(2 * len(remaining_columns))]
    answer = [None for _ in range(2 * columns_count)]
    if len(remC) == 0:
        return answer

    # Запускаем цикл для массива подматриц и получаем решения
    first_line = 0
    last_line = len(remaining_columns)
    for i, (result_lines_combination, result_columns_combination) in enumerate(
            GetArrayMatrix(remC, lines_count, resColVal, DetLim, first_line, last_line)):

        # Получаем матрицу и всё для неё по комбинации строк
        columns = np.array(remaining_columns)[result_columns_combination].tolist()
        squareA, squareMatrix_lines, squareMatrix_columns = GetIntervalSquareMatrix(C, columns_count,
                                                                                    result_lines_combination,
                                                                                    columns)
        squareF = GetFMatrix(squareA, squareMatrix_lines, squareMatrix_columns)
        squared = GetRightHandVector(d, result_lines_combination)
        squarexx = GetxxVector(squared, squareMatrix_lines, 2*squareMatrix_lines)

        check = solver(squarexx, squareMatrix_columns, squareF)
        # Решаем для этой матрицы систему
        Algorithm(squareMatrix_lines, squareMatrix_columns, squareA, squareF, squared, squarexx, Tau, Eps, IterLim)
        MergeAnswers(answer, squarexx, columns)
        ver_answers = []
        for c in columns:
            ver_answers += [answer[2*c], answer[2*c + 1]]

        verification_result = VerificationAnswer(squareA, ver_answers, squared)
        # if verification_result < Eps:
        #     print("\n Verification success \n")
        # else:
        #     print("\n Verification failed \n")
    return answer


# Получение матрицы midA
def GetMidMatrix(A):
    midA = [None for _ in range(len(A) // 2)]
    k = 0
    for i in range(0, len(A), 2):
        midA[k] = 0.5 * (A[i] + A[i + 1])
        k += 1
    return midA


# Объединяем старые решения и полученное новое
def MergeAnswers(allAnswers, newAnswers, columns):

    # Делаем правильные проекции векторов
    for i in range(0, len(newAnswers), 2):
        if newAnswers[i] > newAnswers[i + 1]:
            newAnswers[i + 1], newAnswers[i] = newAnswers[i], newAnswers[i + 1]

    # Пересекаем ответы
    for i, c in enumerate(columns):
        j, k = 2*i, 2*c
        allAnswers[k] = max(allAnswers[k], newAnswers[j]) if allAnswers[k] is not None else newAnswers[j]
        allAnswers[k + 1] = min(allAnswers[k + 1], newAnswers[j + 1]) if allAnswers[k + 1] is not None\
            else newAnswers[j + 1]
