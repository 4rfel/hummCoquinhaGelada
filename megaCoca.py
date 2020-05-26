import numpy as np
# import matplotlib.pyplot as plt

class coca(object):
    def __init__(self, dx, dy, dimensions, alpha, dt, K):
        self.dx = dx
        self.dy = dy
        self.dimensions = dimensions

        self.dx2 = dx**2
        self.dy2 = dy**2
        self.alpha = alpha
        self.dt = dt
        self.f0x = alpha * dt / self.dx2
        self.f0y = alpha * dt / self.dy2
        self.K = K

    def create_temp_matrix(self, temp_int, temp_ext):
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

    def passo(self, matrix):
        copy = np.copy(matrix)
        for i in range(1, matrix.shape[0]-1):
            for j in range(1,matrix.shape[1]-1):
                dTx = self.Eq1(matrix[i][j], matrix[i + 1][j], matrix[i - 1][j])
                dTy = self.Eq2(matrix[i][j], matrix[i][j + 1], matrix[i][j - 1])
                dTz = 0
                copy[i, j] = self.newT(dTx, dTy, dTz, matrix[i][j])

        return copy

    def solve(self, matrix, tMax):
        list_t = np.arange(0, tMax, self.dt)
        for t in list_t:
            matrix = self.passo(matrix)
        return matrix

    def parserText(self, matriz):
        stringF = ""
        for i in range(len(matriz)):
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
alpha = 0.25
temp_ext = [150, 0, 0, 0]
dt = 1e-2

coke = coca(dx, dy, dimensions, alpha, dt, k)

object_temp = 0

matrix_temp = coke.create_temp_matrix(object_temp, temp_ext)

tMax = 10
solved = coke.solve(matrix_temp, tMax)

print(solved)

coke.parserText(solved)