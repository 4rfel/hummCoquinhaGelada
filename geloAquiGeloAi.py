import numpy as np

class coca(object):
    def __init__(self, dx, dy, dimenions, K, ro, Cp, q):
        self.dx = dx
        self.dy = dy
        self.dimensions = dimenions

        self.dx2 = dx**2
        self.dy2 = dy**2
        rocp = ro * Cp
        self.k1 = K/rocp
        self.k2 = q / rocp

    def create_temp_matrix(self, temp_int, temp_ext):
        """
        dimenions is a tuple with the width and the height in that order
        """
        nx = int(self.dimensions[0]/self.dx) + 2
        ny = int(self.dimensions[1]/self.dy) + 2
        matrix_temp = np.zeros([nx, ny])
        for i in range(matrix_temp.shape[0]):
            for j in range(matrix_temp.shape[1]):
                if i == 0:
                    matrix_temp[i][j] = temp_ext[0]
                elif j == 0:
                    matrix_temp[i][j] = temp_ext[1]
                elif i == nx - 1:
                    matrix_temp[i][j] = temp_ext[2]
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

    def newT(self, dTx, dTy, dTz):
        return self.k1*(dTx + dTy + dTz) + self.k2

    def passo(self, matrix):
        copy = np.copy(matrix)
        for i in range(1, matrix.shape[0]-1):
            for j in range(1, matrix.shape[1]-1):
                dTx = self.Eq1(matrix[i][j], matrix[i + 1][j], matrix[i - 1][j])
                dTy = self.Eq2(matrix[i][j], matrix[i][j + 1], matrix[i][j - 1])
                dTz = 0
                copy[i, j] = self.newT(dTx, dTy, dTz)
        return copy

    def solve(self, matrix, tMax, dt):
        list_t = np.arange(0, tMax, dt)
        for t in list_t:
            matrix = self.passo(matrix)
        return matrix




dx = 0.1
dy = 0.1
dimensions = [0.5, 0.5]

K = 1
ro = 1
Cp = 1
q = 1

object_temp = 25

temp_ext = [100, 0, 0, 0]

tMax = 10
dt = 1e-2

coke = coca(dx, dy, dimensions, K, ro, Cp, q)

matrix_temp = coke.create_temp_matrix(object_temp, temp_ext)
# print(matrix_temp)

solved = coke.solve(matrix_temp, tMax, dt)
# print(solved)
