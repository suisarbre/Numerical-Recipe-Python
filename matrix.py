class Matrix:
    def __init__(self, data):
        self.data = data
        self.rows = len(data)
        self.cols = len(data[0])
        
        if not all(len(row) == self.cols for row in data):
            raise ValueError("Inconsistent row lengths")
        
    def __str__(self):
        return '\n'.join(' '.join(str(cell) for cell in row) for row in self.data)
    
    def __repr__(self):
        return f'Matrix\n({self})'
    
    def __getitem__(self, idx):
        return self.data[idx]
    
    def __setitem__(self, idx, new_row):
        self.data[idx] = new_row
    
    def __add__(self,other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrices must have the same dimensions")
        return Matrix([[self.data[i][j] + other.data[i][j] for j in range(self.cols)] for i in range(self.rows)])
    
    def __sub__(self,other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrices must have the same dimensions")
        return Matrix([[self.data[i][j] - other.data[i][j] for j in range(self.cols)] for i in range(self.rows)])
    
    def __mul__(self,other):
        if self.cols != other.rows:
            raise ValueError("Number of columns in the first matrix must be equal to the number of rows in the second matrix")
        return Matrix([[sum(self.data[i][k]*other.data[k][j] for k in range(self.cols)) for j in range(other.cols)] for i in range(self.rows)])
    
    def __eq__(self,other):
        if self.rows != other.rows or self.cols != other.cols:
            return False
        return all(self.data[i][j] == other.data[i][j] for j in range(self.cols) for i in range(self.rows))
    
    def __len__(self):
        return self.rows
    
    def transpose(self):
        return Matrix([[self.data[j][i] for j in range(self.rows)] for i in range(self.cols)])
    
    def row(self, idx):
        return self.data[idx]
    
    def col(self, idx):
        return [row[idx] for row in self.data]
    
from Chapter2.gaussj import gaussj    
if __name__ == "__main__":
    m1 = Matrix([
        [2, 1, -1],
        [-3, -1, 2],
        [-2, 1, 2]
    ])
    m2 = Matrix([
        [8, -11, -3],
        [1, 1, 1],
        [2, 2, 2]
    ])
    gaussj(m1, m2)
    print(m1)
    print(m2)
    