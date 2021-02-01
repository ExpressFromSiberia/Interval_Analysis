from scipy import linalg
from Interval_Operations import MultiplicationMatVec


# Проверка включения интервального вектора Vec_a в интервальный вектор Vec_b
def VerificationInclusion(Vec_a, Vec_b):
    if len(Vec_a) != len(Vec_b):
        print("In VerificationInclusion: The dimensions of the vectors do not match! \n")
        return 0
    verification_result = 1

    for i in range(0, len(Vec_a), 2):
        if (Vec_a[i] < Vec_b[i]) or (Vec_a[i + 1] > Vec_b[i + 1]):
            verification_result = 0
            return verification_result

    return verification_result


# Проверка полученного решения через сумму норм разностей
def VerificationNorm(Ax, b):

    Dim = len(Ax) // 2
    left_parts_vec = [None for _ in range(Dim)]
    right_parts_vec = [None for _ in range(Dim)]

    for i, j in zip(range(0, Dim, 1), range(0, 2*Dim, 2)):
        left_parts_vec[i] = Ax[j] - b[j]
        right_parts_vec[i] = Ax[j + 1] - b[j + 1]

    verification_result = linalg.norm(left_parts_vec) + linalg.norm(right_parts_vec)

    return verification_result


# Проверка полученного ответа через систему Ax = b
def VerificationAnswer(A, x, b):
    # Меняю представление вектора x
    Dim = len(x) // 2
    tmp_x = x.copy()
    for i, j in zip(range(0, Dim, 1), range(0, len(x), 2)):
        tmp_x[j] = x[i]
        tmp_x[j + 1] = x[i + Dim]

    left_sideVec = MultiplicationMatVec(A, tmp_x)
    verification_result = VerificationNorm(left_sideVec, b)

    return verification_result