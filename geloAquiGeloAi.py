import numpy as np
# import matplotlib.pyplot as plt

class coca(object):
    def __init__(self, dx, dy, dimensions, alpha, isolate, dt, K):
        self.dx = dx
        self.dy = dy
        self.dimensions = dimensions

        self.dx2 = dx**2
        self.dy2 = dy**2
        self.alpha = alpha
        self.isolate = isolate
        self.dt = dt
        self.f0x = alpha * dt / self.dx2
        self.f0y = alpha * dt / self.dy2
        self.K = K

    def create_temp_matrix(self, temp_int, temp_ext):
        """
        dimenions is a tuple with the width and the height in that order
        """
        nx = int(self.dimensions[0]/self.dx) + 1
        ny = int(self.dimensions[1]/self.dy) + 1
        matrix_temp = np.zeros([nx, ny])
        for i in range(matrix_temp.shape[0]):
            for j in range(matrix_temp.shape[1]):
                if i == 0:
                    matrix_temp[i][j] = temp_ext[0]
                elif i == nx - 1:
                    matrix_temp[i][j] = temp_ext[2]
                elif j == 0:
                    matrix_temp[i][j] = temp_ext[1]
                elif j == ny - 1:
                    matrix_temp[i][j] = temp_ext[3]
                else:
                    matrix_temp[i][j] = temp_int
        return matrix_temp

    def Eq1(self, T0, T1, Tmenos1):
        return (-2*T0 + Tmenos1 + T1) / self.dx2

    def Eq2(self, T0, T1, Tmenos1):
        return (-2*T0 + Tmenos1 + T1) / self.dy2
        
    def Eq3(self, T0, T1, Tmenos1):
        return (-2*T0 + Tmenos1 + T1) / self.dz2

    def newT(self, dTx, dTy, dTz, T_atual):
        return self.dt * (self.alpha*(dTx + dTy + dTz)) + T_atual

    def fluxo_esquerda(self, T0, T_1x, dT, T_1y, T_menos1y):
        return self.f0x * (2 * T_1x - 2 * self.dx * dT + T_1y + T_menos1y) + (1 - 4 * self.f0x) * T0

    def fluxo_direita(self, T0, dT, T_menos1x, T_1y, T_menos1y):
        return self.f0x * (2 * T_menos1x - 2 * self.dx * dT  + T_1y + T_menos1y) + (1 - 4 * self.f0x) * T0

    def fluxo_topo(self, T0, T_1x, T_menos1x, T_1y, dT):
        return self.f0y * (2 * T_1y - 2 * self.dy * dT + T_1x + T_menos1x) + (1 - 4 * self.f0y) * T0

    def fluxo_base(self, T0, T_1x, T_menos1x, dT, T_menos1y):
        return self.f0y * (2 * T_menos1y - 2 * self.dy * dT + T_1x + T_menos1x) + (1 - 4 * self.f0y) * T0

    def q2LinhaX(self, T_1x, T_menos1x):
        return (self.K * ((T_1x - T_menos1x) / (2 * self.dx)))

    def q2LinhaY(self, T_1y, T_menos1y):
        return (self.K * ((T_1y - T_menos1y) / (2 * self.dy)))

    def passo(self, matrix):
        copy = np.copy(matrix)
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if (i == 0 and j == matrix.shape[1] - 1) or (i == 0 and j == 0) or (i == matrix.shape[0] - 1 \
                    and j == matrix.shape[1] - 1) or (i == matrix.shape[0] - 1 and j == 0):
                    # print(i, j, "canto")
                    pass
                elif i == 0:
                    if self.isolate[0]:
                        dT = self.q2LinhaY(matrix[i + 1][j], matrix[i - 1][j])
                    else:
                        dT = 0
                    # print(i, j, "topo")
                    copy[i, j] = self.fluxo_topo(matrix[i][j], matrix[i][j + 1], matrix[i][j - 1], matrix[i + 1][j], dT)

                elif j == 0:
                    if self.isolate[1]:
                        dT = self.q2LinhaX(matrix[i][j + 1], matrix[i][j - 1])
                    else:
                        dT = 0
                    # print(i, j, "esquerda")
                    copy[i, j] = self.fluxo_esquerda(matrix[i][j], matrix[i][j + 1], dT, matrix[i + 1][j], matrix[i - 1][j])

                elif i == matrix.shape[0] - 1:
                    if self.isolate[2]:
                        dT = self.q2LinhaY(matrix[i + 1][j], matrix[i - 1][j])
                    else:
                        dT = 0
                    # print(i, j, "base")
                    copy[i, j] = self.fluxo_base(matrix[i][j], matrix[i][j + 1], matrix[i][j - 1], dT, matrix[i - 1][j])

                elif j == matrix.shape[1] - 1:
                    if self.isolate[3]:
                        dT = self.q2LinhaX(matrix[i][j + 1], matrix[i][j - 1])
                    else:
                        dT = 0
                    # print(i, j, "direita")
                    copy[i, j] = self.fluxo_direita(matrix[i][j], dT, matrix[i][j - 1], matrix[i + 1][j], matrix[i - 1][j])

                else:
                    dTx = self.Eq1(matrix[i][j], matrix[i + 1][j], matrix[i - 1][j])
                    dTy = self.Eq2(matrix[i][j], matrix[i][j + 1], matrix[i][j - 1])

                    dTz = 0
                    # print(i, j, "meio")
                    copy[i, j] = self.newT(dTx, dTy, dTz, matrix[i][j])

        return copy

    def solve(self, matrix, tMax):
        list_t = np.arange(0, tMax, self.dt)
        for t in list_t:
            matrix = self.passo(matrix)
        return matrix

    def checkErro(self, matrizNossa, matrizGabarito, erro):
        for i in range(matrizGabarito.shape[0]):
            for j in range(matrizGabarito.shape[1]):
                if (matrizGabarito[i][j] - matrizNossa[i][j] > erro):
                    # Return 0 - Tá errado
                    return 0
        # Se chegou até aqui, tá tudo dentro do erro
        return 1

    def rodaErro(self, matrizNossa, matrizGabarito, erro):
        if self.checkErro(self, matrizNossa, matrizGabarito, erro):
            print("Tá certo, gg wp")
            return
        print("TÁ ERRADO")

    def parserText(self, matriz):
        stringF = ""
        for i in range(matriz.shape[0]):
            stringF += "[   "
            for j in range(matriz.shape[1]):
                stringF += str(round(matriz[i][j], 5)) + "    "
            stringF += "]\n"
        with open("saida.txt", "w") as txt:
            txt.write(str(stringF))

dx = 0.10
dy = 0.10
dimensions = [0.40, 0.40]
k = 0.23 # condutividade termica
rho = 2.7e-6 # densidade
Cp = 897 # calor especifico
alpha = k / (rho * Cp)
isolate = [1, 0, 1, 1] # Se a parede é isolada recebe 0, caso contrario 1
temp_ext = [150, 0, 0, 50]
dt = 1e-2

# dados tarefa 2 - Aula 20
gabaritoTarefa2 = [[150,     150,     150,     150,     0], 
                   [12.9929, 12.9981, 13.1632, 16.5732, 50], 
                   [0.5977,  0.6035,  0.7945,  4.8822,  50], 
                   [0.0189,  0.0244,  0.2061,  4.1606,  50], 
                   [0,       0,       0,       0,        0]]

erroTarefa2 = 2.8039e-04



coke = coca(dx, dy, dimensions, alpha, isolate, dt, k)

object_temp = 0

matrix_temp = coke.create_temp_matrix(object_temp, temp_ext)

tMax = 10
solved = coke.solve(matrix_temp, tMax)
coke.rodaErro(solved, gabaritoTarefa2, erroTarefa2)

coke.parserText(solved)