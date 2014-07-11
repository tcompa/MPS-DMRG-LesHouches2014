import numpy

class Ham_Heis:
    def __init__(N, d, J, Jz, h):
        self.d = d
        self.J = J
        self.Jz = Jz
        self.h = h
        self.N = N

    def __getitem__(self, i):
        if i == 0:
            return None
        elif i == 1:
            return self.row
        elif i == N:
            return self.column
        else:
            return self.matrix

    def read_ham_Heis_mpo(self):

        identity = numpy.identity(2)

        S_plus = numpy.zeros((2, 2))
        S_plus[0][1] = 1
        
        S_minus = numpy.zeros((2, 2))
        S_minus[1][0] = 1
        
        S_zed = numpy.identity(2)
        S_zed[1][1] = -1

        mat_W = [None]

        row = numpy.zeros((1, 5, 2, 2))
        column = numpy.zeros((5, 1, 2, 2))
        matrix = numpy.zeros((5, 5, 2, 2))
        # fill row and column
        row[0][0][:, :] = - self.h * S_zed[:, :]
        row[0][1][:, :] = self.J / 2. * S_minus[:, :]
        row[0][2][:, :] = self.J / 2. * S_plus[:, :]
        row[0][3][:, :] = self.Jz * S_zed[:, :]
        row[0][4][:, :] = identity[:, :]

        column[0][0] = identity[:, :]
        column[1][0] = S_plus[:, :]
        column[2][0] = S_minus[:, :]
        column[3][0] = S_zed[:, :]
        column[4][0] = - self.h * S_zed[:, :]
        
        matrix[:, 0, :, :] = column[:, 0]
        matrix[4, :, :, :] = row[0, :]

        mat_W.append(row[:, :])
        for site in range(1, self.N-1):
            mat_W.append(matrix[:, :, :, :])

        self.column = column
        self.row = row
        self.matrix = matrix
#endclass


def read_ham_Heis_mpo(N, d, J, Jz, h):
    print '[read_ham_Heis_mpo] begin'
    identity = numpy.identity(2)

    S_plus = numpy.zeros((2, 2))
    S_plus[0][1] = 1
    
    S_minus = numpy.zeros((2, 2))
    S_minus[1][0] = 1
    
    S_zed = numpy.identity(2)
    S_zed[1][1] = -1

    #mat_W = numpy.zeros((N+1, 5, 5, 2, 2))
    mat_W = [None]

    row = numpy.zeros((1, 5, 2, 2))
    column = numpy.zeros((5, 1, 2, 2))
    matrix = numpy.zeros((5, 5, 2, 2))
    # fill row and column
    row[0][0][:, :] = - h * S_zed[:, :]
    row[0][1][:, :] = J / 2. * S_minus[:, :]
    row[0][2][:, :] = J / 2. * S_plus[:, :]
    row[0][3][:, :] = Jz * S_zed[:, :]
    row[0][4][:, :] = identity[:, :]

    column[0][0] = identity[:, :]
    column[1][0] = S_plus[:, :]
    column[2][0] = S_minus[:, :]
    column[3][0] = S_zed[:, :]
    column[4][0] = - h * S_zed[:, :]
    
    matrix[:, 0, :, :] = column[:, 0]
    matrix[4, :, :, :] = row[0, :]

    mat_W.append(row[:, :])
    for site in range(2, N):
        mat_W.append(matrix[:, :, :, :])
    mat_W.append(column[:, :])

    print '[read_ham_Heis_mpo] mat_W[0]:', type(mat_W[0])
    for w in range(1, len(mat_W)):
        print '[read_ham_Heis_mpo] mat_W[%i].shape:' % w, mat_W[w].shape
    print '[read_ham_Heis_mpo] end'
    return mat_W

if __name__ == "__main__":
    print("running")
    mat_w = read_ham_Heis_mpo(3, 2, 1, 10, 5)

    print("0 -----------------------------------------------------------")
    print (mat_w[0])
    print("1 -----------------------------------------------------------")
    print (mat_w[1])
    print("2 -----------------------------------------------------------")
    print (mat_w[2])
    print("3 -----------------------------------------------------------")
    print (mat_w[3])
    print("4 -----------------------------------------------------------")
