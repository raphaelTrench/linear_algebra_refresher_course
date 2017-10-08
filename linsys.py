from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from hyperplane import Hyperplane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Hyperplane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = system.dimension

        col = 0
        for row in range(num_equations):
            while col < num_variables:
                c = MyDecimal(system[row].normal_vector[col])
                if c.is_near_zero():
                    swap_succeeded = system.swapWithNonZeroBelow(row, col)
                    if not swap_succeeded:
                        col += 1
                        continue

                system.clearCoefficientsBelow(row, col)
                col += 1
                break

        return system

    def compute_rref(self):
        tf = self.compute_triangular_form()
        index = tf.indices_of_first_nonzero_terms_in_each_row()
        num_equations = len(tf)        

        for i in range(num_equations)[::-1]:
            j = index[i]
            if (j<0):
                continue
            tf.scaleRowToMakeCoefficientOne(i, j)
            tf.clearCoefficientsAbove(i, j)
        return tf

    def compute_solution(self):
        try:
            return self.GE_Solution()
        except Exception as e:
            if(str(e) == self.NO_SOLUTIONS_MSG):
                return str(e)
            else:
                raise e  
    
    def raise_exception_if_contradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not(constant_term.is_near_zero()):
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise e

    def raise_exception_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if(num_pivots < num_variables):
           raise Exception(self.INF_SOLUTIONS_MSG)

    def GE_Solution(self):
        rref = self.compute_rref()
        rref.raise_exception_if_contradictory_equation()

        direction_vectors = rref.extract_vectors_for_parametrization()
        basepoint = rref.extract_basepoint_for_parametrization()
        return Parametrization(basepoint,direction_vectors);

    def extract_vectors_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_var_indices = set(range(num_variables)) - set(pivot_indices)

        direction_vectors = []

        for free_var in free_var_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for i,p in enumerate(self.planes):
                pivot_var = pivot_indices[i]
                if(pivot_var < 0):
                    break
                vector_coords[pivot_var] = -p.normal_vector[free_var]
            direction_vectors.append(Vector(vector_coords))

        return direction_vectors


    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()

        basepoint_coords = [0] * num_variables

        for i,p in enumerate(self.planes):
            pivot_var = pivot_indices[i]
            if(pivot_var<0):
                break
            basepoint_coords[pivot_var] = p.constant_term

        return Vector(basepoint_coords)

    def scaleRowToMakeCoefficientOne(self,row,col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)

    def clearCoefficientsAbove(self,row,col):
        for row_to_be_added_to in range(row)[::-1]:
            n = self[row_to_be_added_to].normal_vector
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, row_to_be_added_to)

    def swapWithNonZeroBelow(self,row,col):
        num_equations = len(self)

        for k in range(row + 1, num_equations):
            coefficient = MyDecimal(self[k].normal_vector[col])
            if not coefficient.is_near_zero():
                self.swap_rows(row, k)
                return True
        return False    

    def clearCoefficientsBelow(self,row,col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector[col])

        for row_to_be_added_to in range(row + 1, num_equations):
            n = self[row_to_be_added_to].normal_vector
            gamma = n[col]
            alpha = -gamma / beta
            self.add_multiple_times_row_to_row(alpha, row, row_to_be_added_to)

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        aux = self[row].normal_vector.times_scalar(coefficient)
        aux2 = self[row].constant_term * coefficient
        self[row] = Hyperplane(normal_vector=aux, constant_term=aux2)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        row1 = self[row_to_add].normal_vector.times_scalar(coefficient)
        k1 = self[row_to_add].constant_term * coefficient
        row2 = self[row_to_be_added_to].normal_vector.plus(row1)
        k2 = self[row_to_be_added_to].constant_term
        k = k1 + k2
        self[row_to_be_added_to] = Hyperplane(normal_vector=row2,constant_term=k)        

   

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

class Parametrization(object):

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = (
        'The basepoint and direction vectors should all live in the same '
        'dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output

p1 = Hyperplane(normal_vector=Vector(['0.786', '0.786', '8.123', '1.111', '-8.363']),
                constant_term='-9.955')
p2 = Hyperplane(normal_vector=Vector(['0.131', '-0.131', '7.05', '-2.813', '1.19']),
                constant_term='-1.991')
p3 = Hyperplane(normal_vector=Vector(['9.015', '-5.873', '-1.105', '2.013', '-2.802']),
                constant_term='-3.982')

system = LinearSystem([p1, p2, p3])
print (system.compute_solution())