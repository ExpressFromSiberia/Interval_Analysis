# Получение положительной части вещественного числа
def PositivePart(a):
    return max(a, 0)


# Получение отрицательной части вещественного числа
def NegativePart(a):
    return max(-a, 0)


# Сложение интервалов внутри одного интервального вектора
def IntervalAddition(Vector):
    resInt = [0, 0]
    for i in range(0, len(Vector), 2):
        resInt[0] += Vector[i]
        resInt[1] += Vector[i + 1]

    return resInt


# Умножение интервальных векторов с помощью формул Лакеева
def MultiplicationVecVec(firstVec, secondVec):
    if len(firstVec) != len(secondVec):
        print("In VerificationInclusion: The dimensions of the vectors do not match! \n")
        return 0

    Dim = len(firstVec)
    resultVec = []

    for i in range(0, Dim, 2):
        firstVec_ielem_left = firstVec[i]
        firstVec_ielem_right = firstVec[i + 1]
        secondVec_ielem_left = secondVec[i]
        secondVec_ielem_right = secondVec[i + 1]

        resultVec.append(max(PositivePart(firstVec_ielem_left)*PositivePart(secondVec_ielem_left),
                             NegativePart(firstVec_ielem_right)*NegativePart(secondVec_ielem_right)) -
                         max(PositivePart(firstVec_ielem_right)*NegativePart(secondVec_ielem_left),
                             NegativePart(firstVec_ielem_left)*PositivePart(secondVec_ielem_right)))

        resultVec.append(max(PositivePart(firstVec_ielem_right)*PositivePart(secondVec_ielem_right),
                             NegativePart(firstVec_ielem_left)*NegativePart(secondVec_ielem_left)) -
                         max(PositivePart(firstVec_ielem_left)*NegativePart(secondVec_ielem_right),
                             NegativePart(firstVec_ielem_right)*PositivePart(secondVec_ielem_left)))

    return resultVec


# Умножение интервальной матрицы на интервальный вектор
def MultiplicationMatVec(Matrix, Vector):
    if len(Matrix) != len(Vector)*len(Vector) // 2:
        print("In MultiplicationMatVec: The dimensions of the vector and the matrix did not match! \n")
        return 0

    Dim = len(Vector)
    resultVec = []

    for i in range(0, int(Dim/2)):
        tmp_vec = MultiplicationVecVec(Matrix[i*Dim : i*Dim + Dim], Vector)
        resultVec += IntervalAddition(tmp_vec)
    return resultVec
